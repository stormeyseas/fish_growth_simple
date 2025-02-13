### Functions modified from Baldan et al 2018 R package for aquaculture. 
### https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195732
### https://github.com/cran/RAC/tree/master/R

library(matrixStats)
library(here)


# CREATE REPO STURCTURE

create_repo <- \(this_path){
  
  for(dir in c("data", "figures", "explore")){
    
    this_dir <- file.path(this_path, dir)
    
    if(!dir.exists(this_dir))
      
      dir.create(this_dir)
    
  }
}






#creates the same structure for all species in data folder

create_species_folders  <-  \(this_species){
  
  this_dir <- sprintf(here("data/%s"), this_species)
  
  if(!dir.exists(this_dir))
    
    dir.create(this_dir)
  
  for(this_sub_dir in c("forcings",  "management", "params", "figures", "data_products")){
    
    dir.create(file.path(this_dir, this_sub_dir), recursive = TRUE, showWarnings = FALSE)
    
  } #end for loop
  
  for(this_figures_dir in c("inputs", "outputs")){
    
    if(dir.exists(file.path(this_dir, "figures"))){
      
      dir.create(file.path(this_dir, "figures", this_figures_dir), recursive = TRUE, showWarnings = FALSE)
      
    } #end if loop
    
  } #end for loop
  
}







# DATA LOADER FUNCTION
# This function tidies, gapfills and standardises input data.

data_loader <- \(Path, Farm_id) {
 

  # Reads forcing files
  Ttem = read.csv(sprintf(file.path(Path, "forcings/Water_temperature_%s.csv"), Farm_id), header = FALSE)       # Reading the temperature time series (daily series) data


  #Extracts vectors from the forcing files
  timeT=as.matrix(Ttem[,1])                     # Vector of the times of Temperature measurements
  Temperature=as.double(as.matrix(Ttem[,2]))    # Vector of  water temperature time series (daily series)

  forcings=list(timeT,Temperature)

  return(forcings)
}


# Function that solves the population dynamics equations including discontinuities
Pop_fun <- \(Nseed, mort,times) {
  
  # Integration times
  ti=times[1]
  tf=times[2]
  timestep=times[3]
  
  # Initial condition and vectors initialization
  N=as.vector(matrix(0,nrow=ti))  # Initialize vector N
  N[ti]=Nseed                     # Impose initial condition
  dN=as.vector(matrix(0,nrow=ti)) # Initialize vector dN
  
  for (t in ti:tf){ # for cycle that solves population ODE with Euler method
    
    dN[t]=-mort*N[t]             # individuals increment
    N[t+1]=N[t]+dN[t]*timestep   # Individuals at time t+1

    
   #Taking out the management alterations for now
        
    # for (i in 1:length(manag[,1])) {  # For cycle that adjusts N according with management strategies
    #   if (t==manag[i,1]) {              # if statement to check if it is the time to adjust N
    #     N[t+1]=N[t]+manag[i,2]
    #     
    #   } # close if
    # } # close for
  } # close for
  
  output=N
  return(output)
} # close function




#PREPROCESSOR FUNCTION - pulls in input data and gap fills where necessary to ensure consistent and complete time series

preprocess <- \(Path, Forcings, Feed_type, Stocking, Farm_id){

  cat("Data preprocessing")
  
  # Extracts forcings values from the list
  timeT <- Forcings[[1]]
  Temperature <- Forcings[[2]]
  
  
  # Read forcings and parameters from .csv files
  Param_matrix <- read.csv(sprintf(file.path(Path,"params/Parameters_%s.csv"), Feed_type), sep = ",")           # Reading the matrix containing parameters and their description
  Food <- read.csv(sprintf(file.path(Path, "forcings/Food_characterization_%s.csv"), Feed_type), sep = ",", header = FALSE)    # Reading the food composition (Proteins, Lipids, Carbohydrates) data
  
  # Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
  Spp_param <-  as.matrix(Param_matrix[c(1:21, 26),3])           # Vector containing all parameters
  Spp_param <-  suppressWarnings(as.numeric(Spp_param))
  Dates <- Param_matrix[22:23,3]                      # Vector containing the starting and ending date of the simulation
  #IC=as.double(as.matrix(Param_matrix[24,3]))        # Initial weight condition - hashed out for randomization in population model
  CS <- as.double(as.matrix(Param_matrix[25,3]))      # Commercial size
  Food <- as.double(as.matrix(Food[,1]))              # Food composition (Proteins, Lipids, Carbohydrates) data
  
  # Prepare data for ODE solution
  #t0 <- min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeG[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) #  Minimum starting date for forcings and observations
  timestep <- 1                                        # Time step for integration [day]
  ti <- as.numeric(gsub("day_", "", Dates[1]))   # Start of integration [day]
  tf <- as.numeric(gsub("day_", "", Dates[2]))    # End of integration [day]
  #weight=as.vector(matrix(0,nrow=ti))               # Initialize vector weight
  #weight[ti]=IC                                     # Weight initial value [g]
  times <- cbind(ti, tf, timestep)                    # Vector with integration data
  
  # Food composition vector
  Pcont <- Food[1]       # [-] Percentage of proteins in the food
  Lcont <- Food[2]       # [-] Percentage of lipids in the food
  Ccont <- Food[3]       # [-] Percentage of carbohydrates in the food
  
  # Read files with population parameters and management strategies (not used here)

Pop_matrix <- read.csv(file.path(Path,"params/Population.csv"), sep = ",") |> 
  mutate(Value = case_when(Quantity == "Nseed" ~ Stocking,
                           TRUE ~ Value)) # Reading the matrix containing population parameters and their description
  
  #Management <- read.csv(file.path(Path,"management/Management.csv"), sep = ",")   # Reading the matrix containing seeding and harvesting management
  
  # Extract population parameters
  meanW <- as.double(as.matrix(Pop_matrix[1,3]))      # [g] Dry weight average
  deltaW <- as.double(as.matrix(Pop_matrix[2,3]))     # [g] Dry weight standard deviation
  Wlb <- as.double(as.matrix(Pop_matrix[3,3]))        # [g] Dry weight lower bound
  meanImax <- as.double(as.matrix(Pop_matrix[4,3]))   # [l/d gDW] Clearance rate average
  deltaImax <- as.double(as.matrix(Pop_matrix[5,3]))  # [l/d gDW] Clearance rate standard deviation
  Nseed <-as.double(as.matrix(Pop_matrix[6,3]))      # [-] number of seeded individuals
  mortmyt <- as.double(as.matrix(Pop_matrix[7,3]))    # [1/d] natural mortality rate
  nruns <- as.double(as.matrix(Pop_matrix[8,3]))      # [-] number of runs for population simulation
  
  Pop_param <- c(meanW, deltaW, Wlb, meanImax, deltaImax, Nseed, mortmyt, nruns)

  
  # Population differential equation solution - uses first order Runge Kutta integration (Eulers method) to estimate population number at subsequent time step given an initial slope (N*mortality)
  
  N <- Pop_fun(Nseed, mortmyt, times)
  
  
  # Print to screen inserted parameters
  cat(" \n")
  cat('The model will be executed with the following parameters:\n');
  cat(" \n")
  for (i in 1:21){
    cat(paste0(toString(Param_matrix[i,2]), ": ", toString(Param_matrix[i,3]), " " ,toString(Param_matrix[i,4])),"\n")
  }
  
  cat(" \n")
  cat("Integration is performed between ", toString(Dates[1]), " and ", toString(Dates[2]),"\n")
  cat(" \n")
  cat("The food has the following composition: \n")
  cat(toString(Pcont*100),"% proteins\n")
  cat(toString(Lcont*100),"% lipids\n")
  cat(toString(Ccont*100),"% carbohydrates\n")
  cat(" \n")
  cat('Commercial size is ', toString(CS)," g")
  cat(" \n")
  
  
  # Print to screen population characteristics
  
  cat(" \n")
  cat('The population is simulated by assuming that initial weight and maximum ingestion rate are normally distributed:\n');
  cat(" \n")
  for (i in 1:5){
    cat(paste0(toString(Pop_matrix[i,2]), ": ", toString(Pop_matrix[i,3]), " " ,toString(Pop_matrix[i,4])),"\n")
  }
  
  cat(" \n")
  cat("The population is initially composed by ", unique(Stocking), " Individuals\n")
  cat(" \n")
  cat("The mortality rate is:", toString(Pop_matrix[7,3]),'1/d\n' )
  
  
  
  # Plot to file inserted forcing functions
  
  cat(" \n")
  cat("Forcings are represented in graphs available at the following folder:\n")
  cat(file.path(Path,"/figures/inputs\n"))
  
  
  
  # Plot Temperature forcing
  
  ggplot(data = data.frame(Day = seq(ti:tf), Temperature = Temperature))+
    aes(x = Day, y = Temperature)+
    geom_line()+
    labs(y = bquote(Temperature~degree~C), title = paste(as_label(this_species), as_label(Farm_id)))+
    theme_bw()+
    theme(title = element_text(face="bold", size=8))
  
  ggsave(filename = file.path(Path, sprintf("figures/inputs/Water_temperature_%s.jpeg", Farm_id)), dpi=150, width = 12, height = 8, units="cm")
  

#Plot population change
ggplot(data = data.frame(Day = seq(ti:tf), Stocked_fish = N[ti:tf]))+
    aes(x = Day, y = Stocked_fish)+
    geom_line()+
    labs(y = "Number of individuals", title = paste(as_label(this_species), as_label(Farm_id)))+
    theme_bw()+
    theme(title = element_text(face="bold", size=8))
  
  
ggsave(filename = file.path(Path, sprintf("figures/inputs/Population_%s.jpeg", Farm_id)), dpi=150, width = 12, height = 8, units="cm")
  

  
  output <- list(Spp_param, Pop_param, Temperature, Food,times, Dates, N,CS)
  return(output)
}





ind_equations <- function(Path, Pop_param, Spp_param, Temp, Food, times, N){
  
  
  # Parameters definition
  ingmax=rnorm(1, mean = Pop_param[4], sd = Pop_param[5])     # [g/d] Maximum ingestion rate
  alpha=Spp_param[2]         # [-] Feeding catabolism coefficient
  betaprot=Spp_param[3]      # [-] Assimilation coefficient for protein
  betalip=Spp_param[4]       # [-] Assimilation coefficient for lipid
  betacarb=Spp_param[5]      # [-] Assimilation coefficient for carbohydrates
  epsprot=Spp_param[6]       # [J/gprot] Energy content of protein
  epslip=Spp_param[7]        # [J/glip] Energy content of lipid
  epscarb=Spp_param[8]       # [J/gcarb] Energy content of carbohydrate
  epsO2=Spp_param[9]         # [J/gO2] Energy consumed by the respiration of 1g of oxygen
  pk=Spp_param[10]           # [1/day] Temperature coefficient for the fasting catabolism
  k0=Spp_param[11]           # [1/Celsius degree]  Fasting catabolism at 0 Celsius degree
  m=Spp_param[12]            # [-] Weight exponent for the anabolism
  n=Spp_param[13]            # [-] Weight exponent for the catabolism
  betac=Spp_param[14]        # [-]  Shape coefficient for the H(Tw) function
  Tma=Spp_param[15]          # [Celsius degree] Maximum lethal temperature  
  Toa=Spp_param[16]          # [Celsius degree] Optimal temperature
  Taa=Spp_param[17]          # [Celsius degree] Lowest feeding temperature
  omega=Spp_param[18]        # [gO2/g] Oxygen consumption - weight loss ratio
  a=Spp_param[19]            # [J/gtissue] Energy content of fish tissue
  k=Spp_param[20]            # [-] Weight exponent for energy content
  eff=Spp_param[21]          # [-] Food ingestion efficiency
  fcr = Spp_param[22]
  
  
  # Food composition definition
  Pcont=Food[1]       # [-] Percentage of proteins in the food
  Lcont=Food[2]       # [-] Percentage of lipids in the food
  Ccont=Food[3]       # [-] Percentage of carbohydrates in the food
  
  # Weight and resource
  ti = times[1]
  tf = times[2]                                 # number of days 
  dt = times[3]                               # time step
  # isave = 1
  # nsave  <- floor(tmax/(dt*isave)) # no. of time slots to save
  days<-(1:tf)*dt
  
  weight <- rep(0,tf)
  weight[ti]<- rnorm(1, mean = Pop_param[1], sd = Pop_param[2]) # start size in grams
  # 
  # # vectors of the other important growth metrics: intake, functional response, assimilated energy and growth increment
  fgT <- rep(0, tf)
  frT <- rep(0, tf)
  Tfun <- rep(0, tf)
  ing <- rep(0, tf)
  resource <- rep(0, tf)
  G <- rep(0, tf)
  ingvero <- rep(0, tf)
  epstiss <- rep(0, tf)
  assE <- rep(0, tf)
  Pexc <- rep(0, tf)
  Lexc <- rep(0, tf)
  Cexc <- rep(0, tf)
  exc <- rep(0, tf)
  Pwst <- rep(0, tf)
  Lwst <- rep(0, tf)
  Cwst <- rep(0, tf)
  wst <- rep(0, tf)
  anab <- rep(0, tf)
  catab <- rep(0, tf)
  metab <- rep(0, tf)
  O2 <- rep(0, tf)
  NH4 <- rep(0, tf)
  dw <- rep(0, tf)
  
  # EQUATIONS
  for (i in 1:(tf-1)){
    
    # Forcing temperature
    fgT[i]= exp(betac*(Temp[i]-Toa))*((Tma-Temp[i])/(Tma-Toa))^(betac*(Tma-Toa))   # Optimum Temperature dependence for ingestion
    frT[i]= exp(pk*Temp[i])                                                       # Exponential Temperature dependence for catabolism
    # Tfun[i]=cbind(fgT[i], frT[i])                                               # Output with temperature limitation functions
    
    # Ingested mass
    ing[i]=ingmax*(weight[i]^m)*fgT[i]   # [g/d] Potential ingestion rate
    resource[i] = 0.066*weight[i]^0.75 #from Alex's code
    G[i] = eff*resource[i]
    
    # # Lowest feeding temperature threshold
    if (Temp[i]<Taa) {
      ing[i]=0
    }
    
    # reduced feeding at temps higher than optimal temp
    if (Temp[i]>Toa) {
      ing[i] = 0
    }
    
    # Available food limitation
    if (ing[i] > G[i]) {
      ingvero[i] = G[i]         # [g/d] Actual ingestion rate
    }  else {
      ingvero[i] = ing[i]     # [g/d] Actual ingestion rate
    }
    
    # Energy content of somatic tissue [J/g] Source: Lupatsch et al. (2003)
    epstiss[i] = a*weight[i]^k
    
    # Ingested energy
    diet = Pcont*epsprot*betaprot+Lcont*epslip*betalip+Ccont*epscarb*betacarb # [J/g] Energy content of the ingested food
    assE[i] = ingvero[i]*diet # [J/d] Ingested energy
    
    # Compute excretion (faeces)
    Pexc[i] = (1-betaprot)*Pcont*ingvero[i]  # Excreted proteins [g/d]
    Lexc[i] = (1-betalip)*Lcont*ingvero[i]   # Excreted lipids [g/d]
    Cexc[i] = (1-betacarb)*Ccont*ingvero[i]  # Excreted carbohydrates [g/d]
    # exc[i]=cbind(Pexc[i],Lexc[i],Cexc[i])        # Output with excretion values
    
    # Compute waste (this is uneaten feed only - does not include solid faecal matter)
    Pwst[i]=(resource[i]-ingvero[i])*Pcont     # Proteins to waste [g/d]
    Lwst[i]=((G[i]/eff)-ingvero[i])*Lcont     # Lipids to waste [g/d]
    Cwst[i]=((G[i]/eff)-ingvero[i])*Ccont     # Carbohydrates to waste [g/d]
    # wst[i]=cbind(Pwst[i],Cwst[i],Lwst[i])        # Output with waste values
    
    # Metabolism terms
    anab[i]=assE[i]*(1-alpha)                    # Net anabolism [J/d]
    catab[i]=epsO2*k0*frT[i]*(weight[i]^n)*omega    # Fasting catabolism [J/d]
    # metab[i]=cbind(anab[i],catab[i])                # Output with metabolic rates
    
    # O2 and NH4 produced
    O2[i]=catab[i]/epsO2          # O2 consumed [g02/d]
    NH4[i]=O2[i]*0.06             # NH4 produced [gN/d]
    
    # Mass balance
    dw[i] = (anab[i]-catab[i])/epstiss[i] # weight increment [g/d]
    
    weight[i+1] = weight[i] + dw[i]*dt
  }
  # Function outputs
  output=cbind(days, weight, biomass = weight*N[ti:tf], dw, epstiss, Pexc_total = Pexc*N[ti:tf], Lexc_total = Lexc*N[ti:tf], Cexc_total = Cexc*N[ti:tf], Pwst_total = Pwst*N[ti:tf], Lwst_total = Lwst*N[ti:tf], Cwst_total = Cwst*N[ti:tf], ing_total = ing*N[ti:tf], ingvero_total = ingvero*N[ti:tf], anab, catab, O2, NH4_total = NH4*N[ti:tf], resource_total = resource*N[ti:tf], fgT, frT, Temp)
  return(output) 
}



loop <- function(Path, Pop_param, Spp_param, Temp, Food, times, N, Farm_id){
  
  cat("\n")
  cat("Running population simulation for", Farm_id)
  cat("\n")
  
  
  ti = times[1]
  tf = times[2]
  nruns = Pop_param[8]
  
  
  #initiate matrices to fill for each population iteration
  
  weight_mat <- matrix(data = 0, nrow = nruns, ncol = tf)
  biomass_mat <- matrix(data = 0, nrow = nruns, ncol = tf)
  dw_mat <- matrix(data = 0, nrow = nruns, ncol = tf)
  epistiss_mat <- matrix(data = 0, nrow = nruns, ncol = tf)
  Pexc_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  Lexc_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  Cexc_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  Pwst_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  Lwst_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  Cwst_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  ingvero_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  anab_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  catab_mat <- matrix(data = 0, nrow = nruns, ncol = tf)
  O2_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  NH4_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  resource_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  fgT_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  frT_mat <- matrix(data = 0, nrow = nruns, ncol = tf) 
  
  
  for(n in 1:nruns){
    
    ind_output <- ind_equations(Path = Path, Pop_param = Pop_param, Spp_param = Spp_param, Temp = Temp, Food = Food, times = times,N = N)
    
    #append to matrix
    weight_mat[n,] <- ind_output[,2]
    biomass_mat[n,] <- ind_output[,3]
    dw_mat[n,] <- ind_output[,4]
    epistiss_mat[n,] <- ind_output[,5]
    Pexc_mat[n,] <-  ind_output[,6]
    Lexc_mat[n,] <-  ind_output[,7]
    Cexc_mat[n,] <-  ind_output[,8]
    Pwst_mat[n,] <-  ind_output[,9]
    Lwst_mat[n,] <-  ind_output[,10]
    Cwst_mat[n,] <-  ind_output[,11]
    ingvero_mat[n,] <-  ind_output[,13]
    anab_mat[n,] <-  ind_output[,14]
    catab_mat[n,] <-  ind_output[,15]
    O2_mat[n,] <-  ind_output[,16]
    NH4_mat[n,] <-  ind_output[,17]
    resource_mat[n,] <- ind_output[,18]
    fgT_mat[n,] <-  ind_output[,19]
    frT_mat[n,] <-  ind_output[,20]
    
  }
  
  weight_stat = cbind(colMeans(weight_mat), colSds(weight_mat))
  biomass_stat = cbind(colMeans(biomass_mat), colSds(biomass_mat))
  dw_stat = cbind(colMeans(dw_mat), colSds(dw_mat))
  epistiss_stat = cbind(colMeans(epistiss_mat), colSds(epistiss_mat))
  Pexc_stat = cbind(colMeans(Pexc_mat), colSds(Pexc_mat))
  Lexc_stat = cbind(colMeans(Lexc_mat), colSds(Lexc_mat))
  Cexc_stat = cbind(colMeans(Cexc_mat), colSds(Cexc_mat))
  Pwst_stat = cbind(colMeans(Pwst_mat), colSds(Pwst_mat))
  Lwst_stat = cbind(colMeans(Lwst_mat), colSds(Lwst_mat))
  Cwst_stat = cbind(colMeans(Cwst_mat), colSds(Cwst_mat))
  ingvero_stat = cbind(colMeans(ingvero_mat), colSds(ingvero_mat))
  anab_stat = cbind(colMeans(anab_mat), colSds(anab_mat))
  catab_stat = cbind(colMeans(catab_mat), colSds(catab_mat))
  O2_stat = cbind(colMeans(O2_mat), colSds(O2_mat))
  NH4_stat = cbind(colMeans(NH4_mat), colSds(NH4_mat))
  resource_stat = cbind(colMeans(resource_mat), colSds(resource_mat))
  fgT_stat = cbind(colMeans(fgT_mat), colSds(fgT_mat))
  frT_stat = cbind(colMeans(frT_mat), colSds(frT_mat))
  
  
  
  out_loop <- list(weight_stat, biomass_stat, dw_stat, 
                   epistiss_stat, 
                   Pexc_stat, Lexc_stat, Cexc_stat,
                   Pwst_stat, Lwst_stat, Cwst_stat, 
                   ingvero_stat, anab_stat, catab_stat,
                   O2_stat, NH4_stat, resource_stat,
                   fgT_stat, frT_stat)
  
  return(out_loop)
}




post_process <- function(Path, Farm_id, Feed_type, out_loop, times, N, CS) {
  
  cat('Data post-processing\n')
  cat('\n')
  
  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end
  
  # Extracts outputs from the output list
  W_stat=out_loop[[1]]
  Bio_stat = out_loop[[2]]
  Dw_stat = out_loop[[3]]
  Epistiss_stat = out_loop[[4]]
  Pexc_stat=out_loop[[5]]
  Lexc_stat=out_loop[[6]]
  Cexc_stat=out_loop[[7]]
  Pwst_stat=out_loop[[8]]
  Lwst_stat=out_loop[[9]]
  Cwst_stat=out_loop[[10]]
  Ingvero_stat=out_loop[[11]]
  A_stat=out_loop[[12]]
  C_stat=out_loop[[13]]
  NH4_stat=out_loop[[14]]
  O2_stat=out_loop[[15]]
  Resource_stat = out_loop[[16]]
  fgT_stat=out_loop[[17]]
  frT_stat=out_loop[[18]]
  
  
  # Days to commercial size
  
  # Lower bound
  foo <- function(w,S){which(w>S)[1]}
  arg=as.data.frame(W_stat[,1]-W_stat[,2])
  days <- apply(arg,1,foo,S=CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex)==0) {
    Lb_daysToSize="Not reaching the commercial size"
  }else{  Lb_daysToSize <- min(NonNAindex)
  }
  
  # Mean
  foo <- function(w,S){which(w>S)[1]}
  arg=as.data.frame(W_stat[,1])
  days <- apply(arg,1,foo,S=CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex)==0) {
    Mean_daysToSize="Not reaching the commercial size"
  }else{  Mean_daysToSize <- min(NonNAindex)
  }
  
  # Upper bound
  foo <- function(w,S){which(w>S)[1]}
  arg=as.data.frame(W_stat[,1]+W_stat[,2])
  days <- apply(arg,1,foo,S=CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex)==0) {
    Ub_daysToSize="Not reaching the commercial size"
  }else{  Ub_daysToSize <- min(NonNAindex)
  }
  
  # List containing days to size
  daysToSize<-as.list(cbind(Ub_daysToSize,Mean_daysToSize,Lb_daysToSize))
  
  
  
  out_post = list(W_stat, Bio_stat, Dw_stat, Epistiss_stat, Pexc_stat, Lexc_stat, Cexc_stat, Pwst_stat, Lwst_stat, Cwst_stat, A_stat, C_stat, NH4_stat,O2_stat,  Resource_stat, fgT_stat, frT_stat, daysToSize, N)
  
 
  
  
  # PLOT RESULTS
  #shorter and longer version of days due to lag being used in some functions
  
  days <- seq(from = ti, to = tf-1, by = 1) # create a dates vector to plot results
  days2 <- seq(ti, by = 1, length = tf-ti+1) # create a dates vector to plot results
  
  
  
  
  
  # Plot weight
  
  weight_df = data.frame(days = days2, weight = W_stat[,1], lower_bound = W_stat[,1]-W_stat[,2], upper_bound = W_stat[,1]+W_stat[,2])
  
  #plot weight    
  ggplot(data = weight_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = weight))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Weight (g)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/weight_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  
  # plot biomass
  
  
  biomass_df = data.frame(days = days2, biomass = Bio_stat[,1], lower_bound = Bio_stat[,1]-Bio_stat[,2], upper_bound = Bio_stat[,1]+Bio_stat[,2])
  
  
  ggplot(data = biomass_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = biomass))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Total farm biomass (tonnes)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = \(label) label/1e+6)
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/biomass_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  # plot dW
  
  
  dw_df = data.frame(days = days, dw = Dw_stat[,1][-730], lower_bound = Dw_stat[,1][-730]-Dw_stat[,2][-730], upper_bound = Dw_stat[,1][-730]+Dw_stat[,2][-730])
  
  
  ggplot(data = dw_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = dw))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = bquote(Delta~ Body~weight~g~day^1))+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/dw_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  
  #plot energy content of somatic tissue
  
  epistiss_df = data.frame(days = days, epistiss = Epistiss_stat[,1][-730], lower_bound = Epistiss_stat[,1][-730]-Epistiss_stat[,2][-730], upper_bound = Epistiss_stat[,1][-730]+Epistiss_stat[,2][-730])
  
  
  ggplot(data = epistiss_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = epistiss))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Energy content of somatic tissue (J/g)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))
  
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/epistiss_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  # plot excretion
  
  excretion_df <- 
    rbind(
      data.frame(days = days2, 
                 mean_excretion = Pexc_stat[,1], 
                 lower_bound = Pexc_stat[,1]-Pexc_stat[,2], 
                 upper_bound = Pexc_stat[,1]+Pexc_stat[,2],
                 nutrient = "Protein"),
      data.frame(days = days2, 
                 mean_excretion = Lexc_stat[,1], 
                 lower_bound = Lexc_stat[,1]-Lexc_stat[,2], 
                 upper_bound = Lexc_stat[,1]+Lexc_stat[,2],
                 nutrient = "Lipid"),
      data.frame(days = days2, 
                 mean_excretion = Cexc_stat[,1], 
                 lower_bound = Cexc_stat[,1]-Cexc_stat[,2], 
                 upper_bound = Cexc_stat[,1]+Cexc_stat[,2],
                 nutrient = "Carbohydrates")
    )
  
  
  ggplot(data = excretion_df |> filter(days != 730))+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = nutrient), alpha = 0.2)+
    geom_line(aes(x = days, y = mean_excretion, colour = nutrient))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Excretion (kg/day)")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.title = element_blank(),
          legend.position.inside= c(0.2, 0.8),
          legend.background = element_rect(fill="transparent"))+
    scale_y_continuous(labels = \(label) label/1e+3)
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/excretion_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  # plot wasted feed
  
  feed_waste_df <- 
    rbind(
      data.frame(days = days2, 
                 mean_waste = Pwst_stat[,1], 
                 lower_bound = Pwst_stat[,1]-Pwst_stat[,2], 
                 upper_bound = Pwst_stat[,1]+Pwst_stat[,2],
                 nutrient = "Protein"),
      data.frame(days = days2, 
                 mean_waste = Lwst_stat[,1], 
                 lower_bound = Lwst_stat[,1]-Lwst_stat[,2], 
                 upper_bound = Lwst_stat[,1]+Lwst_stat[,2],
                 nutrient = "Lipid"),
      data.frame(days = days2, 
                 mean_waste = Cwst_stat[,1], 
                 lower_bound = Cwst_stat[,1]-Cwst_stat[,2], 
                 upper_bound = Cwst_stat[,1]+Cwst_stat[,2],
                 nutrient = "Carbohydrates")
    )
  
  ggplot(data = feed_waste_df |> filter(days != 730))+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = nutrient), alpha = 0.2)+
    geom_line(aes(x = days, y = mean_waste, colour = nutrient))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Feed waste (kg/day)")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.25, 0.8),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_blank())+
    scale_y_continuous(labels = \(label) label/1e+3)
  
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/feed_waste_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  #plot ingestion
  
  
  ingvero_df = data.frame(days = days, ingvero = Ingvero_stat[,1][-730], lower_bound = Ingvero_stat[,1][-730]-Ingvero_stat[,2][-730], upper_bound = Ingvero_stat[,1][-730]+Ingvero_stat[,2][-730])
  
  
  
  ggplot(data = ingvero_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = ingvero))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Ingestion (kg/day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = \(label) label/1e+3)
  
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/actual_ingestion_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  metab_df =  rbind(
    data.frame(days = days, value = A_stat[,1][-730], lower_bound = A_stat[,1][-730]-A_stat[,2][-730], upper_bound = A_stat[,1][-730]+A_stat[,2][-730], metabolic_f = "Anabolism"),
    data.frame(days = days, value = C_stat[,1][-730], lower_bound = C_stat[,1][-730]-C_stat[,2][-730], upper_bound = C_stat[,1][-730]+C_stat[,2][-730], metabolic_f = "Catabolism")
  )
  
  ggplot(data = metab_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = metabolic_f), alpha = 0.2)+
    geom_line(aes(x = days, y = value, colour = metabolic_f))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Net matbolic rate (J/day)")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.25, 0.8),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_blank())+
    scale_y_continuous(labels = \(label) label/1e+3)
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/metabolic_rate_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  #Plot O2 consumption
  
  O2_df = data.frame(days = days, o2_consumption = O2_stat[,1][-730], lower_bound = O2_stat[,1][-730]-O2_stat[,2][-730], upper_bound = O2_stat[,1][-730]+O2_stat[,2][-730])
  
  ggplot(data = O2_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = o2_consumption))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Oxygen consumption (kg/day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = \(labels) labels/1e+3)
  
  
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/o2_consumption_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  #plot ammonia
  
  NH4_df = data.frame(days = days, nh4_production = NH4_stat[,1][-730], lower_bound = NH4_stat[,1][-730]-NH4_stat[,2][-730], upper_bound = NH4_stat[,1][-730]+NH4_stat[,2][-730])
  
  
  ggplot(data = NH4_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = nh4_production))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Ammonia production (kg N/ day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = \(labels) labels/1e+3)
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/nh4_production_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  resource_df = data.frame(days = days, feed_resource = Resource_stat[,1][-730], lower_bound = Resource_stat[,1][-730]-Resource_stat[,2][-730], upper_bound = Resource_stat[,1][-730]+Resource_stat[,2][-730])
  
  
  ggplot(data = resource_df)+
    geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = feed_resource))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Available feed (kg/day)")+
    guides(alpha = "none", fill = "none", colour = "none")+
    theme(text=element_text(size=8))+
    scale_y_continuous(labels = \(labels) labels/1e+3, limits = c(0,9000000))
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/available_feed_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  
  
  
  #plot temperature functions
  
  temp_f_df =  rbind(
    data.frame(days = days, value = fgT_stat[,1][-730], lower_bound = fgT_stat[,1][-730]-fgT_stat[,2][-730], upper_bound = fgT_stat[,1][-730]+fgT_stat[,2][-730], metabolic_f = "Anabolism"),
    data.frame(days = days, value = frT_stat[,1][-730], lower_bound = frT_stat[,1][-730]-frT_stat[,2][-730], upper_bound = frT_stat[,1][-730]+frT_stat[,2][-730], metabolic_f = "Catabolism")
  )
  
  
  
  ggplot(data = temp_f_df)+
    #geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = metabolic_f), alpha = 0.2)+
    geom_line(aes(x = days, y = value, colour = metabolic_f))+
    theme_bw()+
    labs(x = "Production cycle (days)", y = "Temperature dependence function")+
    theme(text=element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.25, 0.5),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_blank())+
    scale_y_continuous(labels = \(label) label/1e+3)
  
  ggsave(filename = file.path(Path, sprintf("figures/outputs/%s/temp_function_%s.jpeg", Feed_type, Farm_id)), dpi = 150, width = 12, height = 8, units="cm")
  
  
  #save data
  
  qsave(x= weight_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/weight_output_%s.qs", Feed_type, Farm_id)))
  qsave(x = biomass_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/biomass_output_%s.qs", Feed_type, Farm_id)))
  qsave(x = dw_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/dw_output_%s.qs", Feed_type, Farm_id)))
  qsave(epistiss_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/epistiss_output_%s.qs", Feed_type, Farm_id)))
  qsave(excretion_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/excretion_output_%s.qs", Feed_type, Farm_id)))
  qsave(feed_waste_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/feed_waste_output_%s.qs", Feed_type, Farm_id)))
  qsave(ingvero_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/ingvero_output_%s.qs", Feed_type, Farm_id)))
  qsave(metab_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/metabolic_rate_output_%s.qs", Feed_type, Farm_id)))
  qsave(O2_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/O2_consumption_output_%s.qs", Feed_type, Farm_id)))
  qsave(NH4_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/NH4_production_output_%s.qs", Feed_type, Farm_id)))
  qsave(resource_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/feed_available_output_%s.qs", Feed_type, Farm_id)))
  qsave(temp_f_df, file = file.path(Path, sprintf("data_products/model_outputs/%s/temp_function_output_%s.qs", Feed_type, Farm_id)))
  qsave(daysToSize, file = file.path(Path, sprintf("data_products/model_outputs/%s/days_to_size_output_%s.qs", Feed_type, Farm_id)))

}



model_run <- function(Path, Forcings, Feed_type, Stocking, Farm_id){
  
  cat(" \n")
  cat('Population bioenergetic model for ', str_to_sentence(basename(Path)), Farm_id, "\n")
  cat(" \n")
  
  
  out_pre <- preprocess(Path = Path, Forcings = Forcings, Feed_type = Feed_type, Stocking = Stocking, Farm_id = Farm_id)
  
  Spp_param=out_pre[[1]]
  Pop_param=out_pre[[2]]
  Temp=out_pre[[3]]
  Food=out_pre[[4]]
  #IC=out_pre[[5]]
  times=out_pre[[5]]
  Dates=out_pre[[6]]
  N=out_pre[[7]]
  CS=out_pre[[8]]
  
  # loop through multiple iterations of farm level dynamics
  out_loop <- loop(Path = Path, Spp_param = Spp_param, Pop_param = Pop_param, Temp = Temp, Food = Food, times = times, N = N, Farm_id = Farm_id)
  
  #plot and save model outputs
  out_post <- post_process(Path = Path, Farm_id = Farm_id, Feed_type = Feed_type, out_loop = out_loop, times = times, N = N, CS = CS)
  
  return(out_post)
  
}



