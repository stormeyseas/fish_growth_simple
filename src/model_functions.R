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

 data_loader <- \(this_path, this_farm_id) {
 

  # Reads forcing files
  Ttem = read.csv(sprintf(file.path(this_path, "forcings/Water_temperature_%s.csv"), this_farm_id), header = FALSE)       # Reading the temperature time series (daily series) data
  #DaF =  read.csv(file.path(this_path, "forcings/Feeding.csv"), sep = ",", header = FALSE)  # Reading the individual feeding dose time series (daily series) data


  #Extracts vectors from the forcing files
  timeT=as.matrix(Ttem[,1])                     # Vector of the times of Temperature measurements
  Temperature=as.double(as.matrix(Ttem[,2]))    # Vector of  water temperature time series (daily series)
  # timeF=as.matrix(DaF[,1])                      # Vector of the times of feeding dose
  # Feed=as.double(as.matrix(DaF[,2]))               # Vector of the individual feeding dose time series (daily series)
  #Dates=Param_matrix[22:23,3]                   # Vector containing the starting and ending day of the simulation


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

preprocess <- \(this_path, these_forcings, this_feed, this_stocking_N){

  cat("Data preprocessing")
  
  # Extracts forcings values from the list
  timeT <- these_forcings[[1]]
  Temperature <- these_forcings[[2]]
  # timeF <- these_forcings[[3]]
  # Feed <- these_forcings[[4]]
  
  
  # Read forcings and parameters from .csv files
  Param_matrix <- read.csv(file.path(this_path,"params/Parameters.csv"), sep = ",")           # Reading the matrix containing parameters and their description
  Food <- read.csv(sprintf(file.path(this_path, "forcings/Food_characterization_%s.csv"), this_feed), sep = ",", header = FALSE)    # Reading the food composition (Proteins, Lipids, Carbohydrates) data
  
  # Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
  Param <-  as.matrix(Param_matrix[1:21,3])           # Vector containing all parameters
  Param <-  suppressWarnings(as.numeric(Param))
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

Pop_matrix <- read.csv(file.path(this_path,"params/Population.csv"), sep = ",") |> 
  mutate(Value = case_when(Quantity == "Nseed" ~ this_stocking_N,
                           TRUE ~ Value)) # Reading the matrix containing population parameters and their description
  
  #Management <- read.csv(file.path(this_path,"management/Management.csv"), sep = ",")   # Reading the matrix containing seeding and harvesting management
  
  # Extract population parameters
  meanW <- as.double(as.matrix(Pop_matrix[1,3]))      # [g] Dry weight average
  deltaW <- as.double(as.matrix(Pop_matrix[2,3]))     # [g] Dry weight standard deviation
  IC<- deltaW
  Wlb <- as.double(as.matrix(Pop_matrix[3,3]))        # [g] Dry weight lower bound
  meanImax <- as.double(as.matrix(Pop_matrix[4,3]))   # [l/d gDW] Clearance rate average
  deltaImax <- as.double(as.matrix(Pop_matrix[5,3]))  # [l/d gDW] Clearance rate standard deviation
  Nseed <-as.double(as.matrix(Pop_matrix[6,3]))      # [-] number of seeded individuals
  mortmyt <- as.double(as.matrix(Pop_matrix[7,3]))    # [1/d] natural mortality rate
  nruns <- as.double(as.matrix(Pop_matrix[8,3]))      # [-] number of runs for population simulation
  

  #Not using realistic management values at present
  
    # # Prepare management values
  # manag <- as.matrix(matrix(0,nrow=length(Management[,1]),ncol=2))
  # 
  # for (i in 1:length(Management[,1])) {
  #   manag[i,1]=as.numeric(as.Date(Management[i,1], "%d/%m/%Y"))-t0
  #   if ((Management[i,2])=="h") {
  #     manag[i,2]=-as.numeric(Management[i,3])
  #   } else {
  #     manag[i,2]=as.numeric(Management[i,3])
  #   }
  # }
  # 
  # 
  
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
  cat("The population is initially composed by ", unique(this_stocking_N), " Individuals\n")
  cat(" \n")
  cat("The mortality rate is:", toString(Pop_matrix[7,3]),'1/d\n' )
  
  
  #TAKING OUT MANAGEMENT FOR NOW
  # Print to screen management actions
  # cat(" \n")
  # cat('The population is managed according with following list (h:harvesting s:seeding):\n');
  # cat(" \n")
  # for (i in 1:length(Management[,1])){
  #   cat(paste0(toString(Management[i,1])," ", toString(Management[i,2]), " " ,toString(Management[i,3])),"individuals\n")
  # }
  # 
  # cat(" \n")
  # cat("The individual model will be executed ", toString(nruns), " times in order to simulate a population\n")
  # cat(" \n")
  
  
  
  # Plot to file inserted forcing functions
  
  cat(" \n")
  cat("Forcings are represented in graphs available at the following folder:\n")
  cat(file.path(this_path,"/figures/inputs\n"))
  
  
  
  # Plot Temperature forcing
  
  ggplot(data = data.frame(Day = seq(ti:tf), Temperature = Temperature))+
    aes(x = Day, y = Temperature)+
    geom_line()+
    labs(y = bquote(Temperature~degree~C), title = paste(as_label(this_species), as_label(this_farm_id)))+
    theme_bw()+
    theme(title = element_text(face="bold", size=8))
  
  ggsave(filename = file.path(this_path, sprintf("figures/inputs/Water_temperature_%s.jpeg", this_farm_id)), dpi=150, width = 12, height = 8, units="cm")
  
  
  # Masking feed supply for now
  # # Plot feeding rate forcing
  # Gintsave=Gint[(ti+1):tf]
  # filepath <- file.path(this_path,"/figures/inputs/Feeding.jpeg")
  # jpeg(filepath,800,600)
  # days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)
  # plot(days, Gintsave, ylab="Feed (g/d)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
  # labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  # axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  # dev.off()
  

#Plot population change
ggplot(data = data.frame(Day = seq(ti:tf), Stocked_fish = N[ti:tf]))+
    aes(x = Day, y = Stocked_fish)+
    geom_line()+
    labs(y = "Number of individuals", title = paste(as_label(this_species), as_label(this_farm_id)))+
    theme_bw()+
    theme(title = element_text(face="bold", size=8))
  
  
ggsave(filename = file.path(this_path, sprintf("figures/inputs/Population_%s.jpeg", this_farm_id)), dpi=150, width = 12, height = 8, units="cm")
  

  
  output <- list(Param, Temperature, Food, IC, times, Dates, N,CS)
  return(output)
}





# MODEL EQUATIONS  - calculates the growth and outputs for the bioenergetic model


equations <- \(Param, N, Temp, Food, weight){
  
  # Parameters definition
  # Parameters definition
  ingmax=Param[1]        # [g/d] Maximum ingestion rate
  alpha=Param[2]         # [-] Feeding catabolism coefficient
  betaprot=Param[3]      # [-] Assimilation coefficient for protein
  betalip=Param[4]       # [-] Assimilation coefficient for lipid
  betacarb=Param[5]      # [-] Assimilation coefficient for carbohydrates
  epsprot=Param[6]       # [J/gprot] Energy content of protein
  epslip=Param[7]        # [J/glip] Energy content of lipid
  epscarb=Param[8]       # [J/gcarb] Energy content of carbohydrate
  epsO2=Param[9]         # [J/gO2] Energy consumed by the respiration of 1g of oxygen
  pk=Param[10]           # [1/day] Temperature coefficient for the fasting catabolism
  k0=Param[11]           # [1/Celsius degree]  Fasting catabolism at 0 Celsius degree
  m=Param[12]            # [-] Weight exponent for the anabolism
  n=Param[13]            # [-] Weight exponent for the catabolism
  betac=Param[14]        # [-]  Shape coefficient for the H(Tw) function
  Tma=Param[15]          # [Celsius degree] Maximum lethal temperature 
  Toa=Param[16]          # [Celsius degree] Optimal temperature 
  Taa=Param[17]          # [Celsius degree] Lowest feeding temperature
  omega=Param[18]        # [gO2/g] Oxygen consumption - weight loss ratio
  a=Param[19]            # [J/gtissue] Energy content of fish tissue
  k=Param[20]            # [-] Weight exponent for energy content
  eff=Param[21]          # [-] Food ingestion efficiency
  
  
  
  # Food composition definition
  Pcont=Food[1]       # [-] Percentage of proteins in the food
  Lcont=Food[2]       # [-] Percentage of lipids in the food
  Ccont=Food[3]       # [-] Percentage of carbohydrates in the food
  
  
  # EQUATIONS
  
  #G = weight*ingmax #total feed provided
  
  # Forcing temperature
  fgT=((Tma-Temp)/(Tma-Toa))^(betac*(Tma-Toa))*exp(betac*(Temp-Toa)) # Optimum Temperature dependance for ingestion [-]
  frT= exp(pk*Temp)                                                  # Exponential Temperature dependance for catabolism [-]
  Tfun=cbind(fgT, frT)                                               # 
  
  # Ingested mass
  ing=ingmax*(weight^m)*fgT  # Potential ingestion rate  [g/d]
  
  G=ing/eff   # Food ingestion efficiency is used to assume that fish are fed to satiation but more is provided because they do not feed with absolute efficiency
  
  # Lowest feeding temperature threshold
  if (is.na(Temp)) print("Temp")
  
  if (Temp<Taa) {
    ing=0
  }
  
  # Available food limitation
  if (ing>G) {
    ingvero=G      # [g/d] Actual ingestion rate
  } else {
    ingvero=ing  # [g/d] Actual ingestion rate
  }
  
  # Energy content of somatic tissue [J/g]. Lupatsch et al. (2003)
  epstiss=a*weight^k
  
  # Ingested energy
  diet=Pcont*epsprot*betaprot+Lcont*epslip*betalip+Ccont*epscarb*betacarb  # [J/g]
  assE=ingvero*diet   # [J/d]
  
  # Compute excretion
  Pexc=(1-betaprot)*Pcont*ingvero*N/1e3  # Excreted proteins [kg/d]
  Lexc=(1-betalip)*Lcont*ingvero*N/1e3   # Excreted lipids [kg/d]
  Cexc=(1-betacarb)*Ccont*ingvero*N/1e3  # Excreted carbohydrates [kg/d]
  exc=cbind(Pexc,Lexc,Cexc)
  
  # Compute waste
  Pwst=((G)-ingvero)*Pcont*N/1e3  # Proteins to waste [kg/d]
  Lwst=((G/eff)-ingvero)*Lcont*N/1e3  # Lipids to waste [kg/d]
  Cwst=((G/eff)-ingvero)*Ccont*N/1e3  # Carbohydrates to waste [kg/d]
  wst=cbind(Pwst,Cwst,Lwst)
  
  # Metabolism terms
  anab=assE*(1-alpha)                    # Net anabolic rate [J/d]
  catab=epsO2*k0*frT*(weight^n)*omega    # Fasting catabolic rate [J/d]
  metab=cbind(anab,catab)
  
  # O2 consumed and NH4 produced
  O2=catab/epsO2*N/1e3          # O2 consumed [kg02/d]
  NH4=O2*0.06*N/1e3             # NH4 produced [kgN/d]
  
  # Mass balance
  dw = (anab-catab)/epstiss   # Weight increment [g/d]
  
  # Function outputs
  output=list(dw,exc,wst,ing,ingvero,Tfun,metab, O2, NH4, G)
  return(output)
}


#Runge-Kutte integration

RKsolver <- \(Param, Temperature, Food, IC, times, N){ #taken G out as an input and included it as a more dynamic function of ingestion (adjusted by efficiency)
  
  # Integration extremes definition
  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end
  timestep=times[3]     # Timestep for integration
  
  # Initial condition definition
  weight=as.vector(matrix(0,nrow=ti))      # Initialize vector w
  weight[ti]=IC                            # Define initial condition
  
  # initialize outputs
  wst=as.matrix(matrix(0,nrow=ti,ncol=3))      # Initialize food to waste vector
  exc=as.matrix(matrix(0,nrow=ti,ncol=3))      # Initialize excretion vector
  ing=as.matrix(matrix(0,nrow=ti,ncol=1))             # Initialize potential ingestion vector
  ingvero=as.matrix(matrix(0,nrow=ti,ncol=1))         # Initialize actual ingestion vector
  tfun=as.matrix(matrix(0,nrow=ti,ncol=2))     # Initialize temperature limitations vector
  metab=as.matrix(matrix(0,nrow=ti,ncol=2))   # Initialize metabolic rates vector
  O2=as.matrix(matrix(0,nrow=ti))       # Initialize oxygen consumption rates vector
  NH4=as.matrix(matrix(0,nrow=ti))      # Initialize ammonia release rates vector
  feed_available = as.matrix(matrix(0,nrow=ti))  # Initialize feed availability rates vector
  
  for (t in ti:(tf-1)) {
    
    # Compute Runge-Kutta increments
    
    # 1
    Tapp=Temperature[t]
    #Gapp=G[t]
    Napp=N[t]
    output<-equations(Param = Param, N = Napp, Temp = Tapp,  Food = Food, weight = weight[t]) # Taken Gapp out because G is resolved in equations
    dw=unlist(output[1])
    k1=timestep*dw
    
    # 2
    Tapp=approx(seq(from=1,to=tf,by=timestep),Temperature,xout=(t+timestep/2))
    #Gapp=approx(seq(from=1,to=tf,by=timestep),G,xout=(t+timestep/2))
    output<-equations(Param, Napp, Tapp$y, Food, weight = weight[t]+k1/2)
    dw=unlist(output[1])
    k2=timestep*dw;
    
    # 3
    Tapp=approx(seq(from=1,to=tf,by=timestep),Temperature,xout=(t+timestep/2))
    #Gapp=approx(seq(from=1,to=tf,by=timestep),G,xout=(t+timestep/2))
    output<-equations(Param, Napp, Tapp$y,Food, weight[t]+k2/2)
    dw=unlist(output[1])
    k3=timestep*dw;
    
    # 4
    Tapp=Temperature[t+timestep]
    #Gapp=G[t+timestep]
    output<-equations(Param = Param, N = Napp, Temp = Tapp,  Food = Food, weight = weight[t]+k3)
    dw=unlist(output[1])
    k4=timestep*dw;
    
    # Compute weight at t+1 using Runge-Kutta increments
    weight[t+timestep]=weight[t]+(k1+2*k2+2*k3+k4)/6;
    
    # Compute the other outputs of the model
    output<-equations(Param = Param, N = N[t+timestep], Temp = Temperature[t+timestep], Food = Food, weight = weight[t+timestep])
    
    # Extracts outputs from the output list
    excretion=output[[2]]
    waste=output[[3]]
    ingestion=unlist(output[4])
    ingestionvero=unlist(output[5])
    temperaturefun=output[[6]]
    metabolism=output[[7]]
    oxygen=output[[8]]
    ammonia=output[[9]]
    feed_input = output[[10]]
    
    # Outputs creation
    wst=rbind(wst, waste)
    exc=rbind(exc, excretion)
    ing=rbind(ing, ingestion)
    ingvero=rbind(ingvero, ingestionvero)
    tfun=rbind(tfun, temperaturefun)
    metab=rbind(metab, metabolism)
    O2=rbind(O2, oxygen)
    NH4=rbind(NH4, ammonia)
    feed_input = rbind(feed_available, feed_input)
    
  }  # Close cycle
  
  output=list(weight,exc,wst,ing,ingvero,tfun,metab, O2, NH4, feed_input)
  return(output) # Bream_pop_RKsolver output
  
} # Close function






# MONTE-CARLO SIMULATION

loop <- \(Param, Tint, Food, IC, times, N, this_path) {
  
  cat("Population processing\n")
  
  ti=times[1]
  tf=times[2]
  t0=times[4]
  
  # Read files with population parameters and management strategies
  Pop_matrix <- read.csv(file.path(this_path,"params/Population.csv"), sep = ",")   # Reading the matrix containing population parameters and their description
  #Management <- read.csv(file.path(this_path,"management/Management.csv"), sep = ",")   # Reading the matrix containing seeding and harvesting management
  
  # Extract population parameters
  meanW<-as.double(as.matrix(Pop_matrix[1,3]))      # [g] Dry weight average
  deltaW<-as.double(as.matrix(Pop_matrix[2,3]))     # [g] Dry weight standard deviation
  Wlb<-as.double(as.matrix(Pop_matrix[3,3]))        # [g] Dry weight lower bound
  meanImax<-as.double(as.matrix(Pop_matrix[4,3]))   # [l/d gDW] Clearence rate average
  deltaImax<-as.double(as.matrix(Pop_matrix[5,3]))  # [l/d gDW] Clearance rate standard deviation
  Nseed<-as.double(as.matrix(Pop_matrix[6,3]))      # [-] number of seeded individuals
  mortmyt<-as.double(as.matrix(Pop_matrix[7,3]))    # [1/d] natural mortality rate
  nruns<-as.double(as.matrix(Pop_matrix[8,3]))      # [-] number of runs for population simulation
  
  # Prepare management values
  # manag <- as.matrix(matrix(0,nrow=length(Management[,1]),ncol=2))
  # for (i in 1:length(Management[,1])) {
  #   manag[i,1]=as.numeric(as.Date(Management[i,1], "%d/%m/%Y"))-t0
  #   if ((Management[i,2])=="h") {
  #     manag[i,2]=-as.numeric(Management[i,3])
  #   } else {
  #     manag[i,2]=as.numeric(Management[i,3])
  #   }
  # }
  
  # Vectors initialization
  saveIC=as.vector(matrix(0,nrow=nruns))            # Initialize initial conditions records: saves perturbated initial conditions for each run
  saveImax=as.vector(matrix(0,nrow=nruns))          # Initialize maximum ingestion rate records: saves perturbated maximum ingestion rate for each run
  W=as.matrix(matrix(0,nrow=nruns,ncol=tf))         # initialize weight vector
  Pexc=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize excreted proteines vector
  Lexc=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize excreted lipids vector
  Cexc=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize excreted carbohydrates vector
  Pwst=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize proteines to waste vector
  Lwst=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize lipids to waste vector
  Cwst=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize carbohydrates to waste vector
  ingestion=as.matrix(matrix(0,nrow=nruns,ncol=tf)) # Initialize actual ingestion vector
  A=as.matrix(matrix(0,nrow=nruns,ncol=tf))         # Initialize anabolic rate vector
  C=as.matrix(matrix(0,nrow=nruns,ncol=tf))         # Initialize catabolic rate vector
  O2=as.matrix(matrix(0,nrow=nruns,ncol=tf))        # Initialize oxygen consumption rate vector
  NH4=as.matrix(matrix(0,nrow=nruns,ncol=tf))       # Initialize ammonium release rate vector
  feed_available = as.matrix(matrix(0,nrow=nruns,ncol=tf)) # Initialise feeding vector
  
  # Population LOOP
  
  pb <- txtProgressBar(min = 0, max = nruns, style = 3)
  
  for (ii in 1:nruns){
    
    
    # Weight initialization
    IC=rnorm(1,meanW,deltaW)     # [g] initial weight extracted from a normal distribution
    IC=max(IC, Wlb)              # Lower bound for weight distribution
    saveIC[ii]=IC                # Saves initial condition values on W for each run
    
    # Maximum clearance rate initialization
    Imax=rnorm(1,meanImax,deltaImax)  # [l/d gDW] Maximum ingestion rate extracted from a normal distribution
    Imax=max(Imax,0)                  # Forces maximum ingestion rate to be positive
    saveImax[ii]=Imax                 # Saves initial condition values on Imax for each run
    
    # Perturb the parameters vector
    Param[1]=Imax
    
    # Solves ODE with perturbed parameters
    output<-RKsolver(Param, Tint,  Food, IC, times, N)
    
    # Unlist outputs
    weight=unlist(output[1])
    exc=output[[2]]
    wst=output[[3]]
    ing=unlist(output[4])
    ingvero=unlist(output[5])
    Tfun=output[[6]]
    metab=output[[7]]
    oxygen=output[[8]]
    ammonium=output[[9]]
    feed_input = output[[10]]
    
    
    
    # Saves results of each run to compute statistics
    W[ii,1:length(weight)]=weight           # Tissue dry weight [g]
    
    Pexc[ii,1:length(exc[,1])]=exc[,1]      # Excreted proteins [kg/d]
    Lexc[ii,1:length(exc[,1])]=exc[,2]      # Excreted lipids [kg/d]
    Cexc[ii,1:length(exc[,1])]=exc[,3]      # Excreted carbohydrates [kg/d]
    
    Pwst[ii,1:length(wst[,1])]=wst[,1]      # Proteins to waste [kg/d]
    Lwst[ii,1:length(wst[,1])]=wst[,2]      # Lipids to waste [kg/d]
    Cwst[ii,1:length(wst[,1])]=wst[,3]      # Carbohydrates to waste [kg/d]
    
    ingestion[ii,1:length(ingvero)]=t(ingvero)  # Actual ingested food [g/d]
    
    
    A[ii,1:length(metab[,1])]=metab[,1]     # Net anabolism [J/d]
    C[ii,1:length(metab[,1])]=metab[,2]     # Fasting catabolism [J/d]
    
    O2[ii,1:length(oxygen)]=oxygen           # Tissue dry weight [kgO2/d]
    NH4[ii,1:length(ammonium)]=ammonium      # Tissue dry weight [kgN/d]
    feed_available[ii,1:length(feed_input)]=feed_input      # Tissue dry weight [kgN/d]
    
    setTxtProgressBar(pb, ii)
    
  } # Close population loop
  
 close(pb)
  
  # Temperaure limitation functions
  
  fgT=Tfun[,1] # Optimum  dependance from temperature for ingestion
  frT=Tfun[,2] # Exponential dependance from temperature for catabolism
  
  # Statistics computation:
  # Statistics vectors contain the mean and the standard deviation of a variable
  
  W_stat=t(rbind(colMeans(W), colSds(W)))
  
  Pexc_stat=t(rbind(colMeans(Pexc), colSds(Pexc)))
  Lexc_stat=t(rbind(colMeans(Lexc), colSds(Lexc)))
  Cexc_stat=t(rbind(colMeans(Cexc), colSds(Cexc)))
  
  ingestion_stat=t(rbind(colMeans(ingestion), colSds(ingestion)))
  
  Pwst_stat=t(rbind(colMeans(Pwst), colSds(Pwst)))
  Lwst_stat=t(rbind(colMeans(Lwst), colSds(Lwst)))
  Cwst_stat=t(rbind(colMeans(Cwst), colSds(Cwst)))
  
  A_stat=t(rbind(colMeans(A), colSds(A)))
  C_stat=t(rbind(colMeans(C), colSds(C)))
  
  O2_stat=t(rbind(colMeans(O2), colSds(O2)))
  NH4_stat=t(rbind(colMeans(NH4), colSds(NH4)))
  feed_stat = t(rbind(colMeans(feed_available), colSds(feed_available)))
  
  output=list(W_stat,Pexc_stat,Lexc_stat,Cexc_stat,ingestion_stat,Pwst_stat,Lwst_stat,Cwst_stat,A_stat,C_stat,fgT,frT, O2_stat, NH4_stat, feed_stat)
  return(output)
  
}




# POST PROCESSING - plots the model outputs

post_process <- function(this_path, this_farm_id, output,times,Dates,N, CS) {
  
  cat('Data post-processing\n')
  cat('\n')
  
  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end
  
  # Extracts outputs from the output list
  W_stat=output[[1]]
  Pexc_stat=output[[2]]
  Lexc_stat=output[[3]]
  Cexc_stat=output[[4]]
  ingestion_stat=output[[5]]
  Pwst_stat=output[[6]]
  Lwst_stat=output[[7]]
  Cwst_stat=output[[8]]
  A_stat=output[[9]]
  C_stat=output[[10]]
  fgT=output[[11]]
  frT=output[[12]]
  O2_stat=output[[13]]
  NH4_stat=output[[14]]
  feed_stat = output[[15]]
  
  # Adjusts results acoording with integration extremes
  # now day 1 coincides with ti
  weightSave=W_stat[ti:tf,]
  
  PexcSave=Pexc_stat[ti:tf,]
  LexcSave=Lexc_stat[ti:tf,]
  CexcSave=Cexc_stat[ti:tf,]
  
  ingestionSave=ingestion_stat[(ti+1):tf,]
  
  PwstSave=Pwst_stat[ti:tf:tf,]
  LwstSave=Lwst_stat[ti:tf,]
  CwstSave=Cwst_stat[ti:tf,]
  
  ASave=A_stat[ti:tf,]
  CSave=C_stat[ti:tf,]
  
  fgT=fgT[(ti+1):tf]
  frT=frT[(ti+1):tf]
  
  O2Save=O2_stat[(ti+1):tf,]
  NH4Save=NH4_stat[(ti+1):tf,]
  
  feedSave = feed_stat[(ti+1):tf,]
  
  N=N[ti:tf]
  
  # Days to commercial size
  
  # Lower bound
  foo <- function(w,S){which(w>S)[1]}
  arg=as.data.frame(weightSave[,1]-weightSave[,2])
  days <- apply(arg,1,foo,S=CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex)==0) {
    Lb_daysToSize="Not reaching the commercial size"
  }else{  Lb_daysToSize <- min(NonNAindex)
  }
  
  # Mean
  foo <- function(w,S){which(w>S)[1]}
  arg=as.data.frame(weightSave[,1])
  days <- apply(arg,1,foo,S=CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex)==0) {
    Mean_daysToSize="Not reaching the commercial size"
  }else{  Mean_daysToSize <- min(NonNAindex)
  }
  
  # Upper bound
  foo <- function(w,S){which(w>S)[1]}
  arg=as.data.frame(weightSave[,1]+weightSave[,2])
  days <- apply(arg,1,foo,S=CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex)==0) {
    Ub_daysToSize="Not reaching the commercial size"
  }else{  Ub_daysToSize <- min(NonNAindex)
  }
  
  # List containing days to size
  daysToSize<-as.list(cbind(Ub_daysToSize,Mean_daysToSize,Lb_daysToSize))
  
  output=list(weightSave,PexcSave,LexcSave,CexcSave,ingestionSave,PwstSave,LwstSave,CwstSave,ASave,CSave,fgT,frT,O2Save, NH4Save, N,daysToSize)
  
  # Plot results
  days <- seq(from = ti, to = tf-1, by = 1) # create a dates vector to plot results
  days2 <- seq(ti, by = 1, length = tf-ti+1) # create a dates vector to plot results
  
  # Plot weight

  weight = data.frame(days = days2, weight = weightSave[,1], lower_bound = weightSave[,1]-weightSave[,2], upper_bound = weightSave[,1]+weightSave[,2])

#plot weight    
ggplot(data = weight)+
  geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "salmon", colour = "grey50", linetype = "dashed")+
    geom_line(aes(x = days, y = weight))+
  theme_bw()+
  labs(x = "Production cycle (days)", y = "Weight (g)")+
  guides(alpha = "none", fill = "none", colour = "none")+
  theme(text=element_text(size=8))

ggsave(filename = file.path(this_path, sprintf("figures/outputs/weight_%s.jpeg", this_farm_id)), dpi = 150, width = 12, height = 8, units="cm")

  
  # plot excretion

excretion <- 
  rbind(
  data.frame(days = days, 
             mean_excretion = PexcSave[,1], 
             lower_bound = PexcSave[,1]-PexcSave[,2], 
             upper_bound = PexcSave[,1]+PexcSave[,2],
             nutrient = "Protein"),
  data.frame(days = days, 
             mean_excretion = LexcSave[,1], 
             lower_bound = LexcSave[,1]-LexcSave[,2], 
             upper_bound = LexcSave[,1]+LexcSave[,2],
             nutrient = "Lipid"),
  data.frame(days = days, 
             mean_excretion = CexcSave[,1], 
             lower_bound = CexcSave[,1]-CexcSave[,2], 
             upper_bound = CexcSave[,1]+CexcSave[,2],
             nutrient = "Carbohydrates")
  )

ggplot(data = excretion)+
  geom_ribbon(aes(x = days, ymin = lower_bound, ymax = upper_bound, fill = nutrient), alpha = 0.2)+
  geom_line(aes(x = days, y = mean_excretion, colour = nutrient))+
  theme_bw()+
  labs(x = "Production cycle (days)", y = "Excretion (kg/day)")+
  guides(alpha = "none", fill = "none", colour = "none")+
  theme(text=element_text(size=8))

ggsave(filename = file.path(this_path, sprintf("figures/outputs/excretion_%s.jpeg", this_farm_id)), dpi = 150, width = 12, height = 8, units="cm")



  # plot wasted feed
  filepath=file.path(this_path,"/figures/outputs/wasted_feed.jpeg")
  jpeg(filepath,800,600)
  Lub=LwstSave[,1]+LwstSave[,2]
  Pub=PwstSave[,1]+PwstSave[,2]
  Cub=CwstSave[,1]+CwstSave[,2]
  Llb=as.matrix(matrix(0,nrow=length(Lub),ncol=1))
  Plb=as.matrix(matrix(0,nrow=length(Pub),ncol=1))
  Clb=as.matrix(matrix(0,nrow=length(Cub),ncol=1))
  for (i in 1:length(LwstSave[,1]-LwstSave[,2])){
    Llb[i]=max(LwstSave[i,1]-LwstSave[i,2],0)
    Plb[i]=max(PwstSave[i,1]-PwstSave[i,2],0)
    Clb[i]=max(CwstSave[i,1]-CwstSave[i,2],0)
  }
  maxub=max(Lub,Pub,Cub)
  plot(days,LwstSave[,1],ylab="Wasted feed (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
  polygon(c(days,rev(days)),c(Llb,rev(Lub)),col="grey75",border=FALSE)
  lines(days,LwstSave[,1],lwd=2,col="red")
  polygon(c(days,rev(days)),c(Plb,rev(Pub)),col="grey75",border=FALSE)
  lines(days,PwstSave[,1],lwd=2,col="green")
  polygon(c(days,rev(days)),c(Clb,rev(Cub)),col="grey75",border=FALSE)
  lines(days,CwstSave[,1],lwd=2,col="blue")
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  legend("topleft",c("Proteins","Lipids","Carbohydrates"),fill=c("red","green","blue"))
  dev.off()
  
  # plot ingested food
  filepath=file.path(this_path,"figures/outputs/actual_ingestion.jpeg")
  jpeg(filepath,800,600)
  ub=ingestionSave[,1]+ingestionSave[,2]
  lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
  for (i in 1:length(ingestionSave[,1]-ingestionSave[,2])){
    lb[i]=max(ingestionSave[i,1]-ingestionSave[i,2],0)
  }
  maxub=max(ingestionSave[,1]+ingestionSave[,2])
  plot(days,ingestionSave[,1],ylab="Ingested food (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
  polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
  lines(days,ingestionSave[,1],lwd=2,col="red")
  lines(days,lb,col="blue")
  lines(days,ub,col="blue")
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()
  
  # plot limitation functions
  filepath=file.path(this_path,"figures/outputs/temperature_response.jpeg")
  jpeg(filepath,800,600)
  ub=max(max(fgT),max(frT))
  plot(days,fgT,ylab="Temperature response function",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
  lines(days,frT,col="blue")
  legend("topright",c("Anabolism","Catabolism"),fill=c("red","blue"))
  labDates <- seq(ti:tf)
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()
  
  # plot metabolic rates
  filepath=file.path(this_path,"/figures/outputs/metabolism.jpeg")
  jpeg(filepath,800,600)
  Aub=ASave[,1]+ASave[,2]
  Cub=CSave[,1]+CSave[,2]
  Alb=as.matrix(matrix(0,nrow=length(Aub),ncol=1))
  Clb=as.matrix(matrix(0,nrow=length(Cub),ncol=1))
  for (i in 1:length(LwstSave[,1]-LwstSave[,2])){
    Alb[i]=max(ASave[i,1]-ASave[i,2],0)
    Clb[i]=max(CSave[i,1]-CSave[i,2],0)
  }
  maxub=max(Aub,Cub)
  plot(days,ASave[,1],ylab="Metabolic rates (J/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
  polygon(c(days,rev(days)),c(Alb,rev(Aub)),col="grey75",border=FALSE)
  lines(days,ASave[,1],lwd=2,col="red")
  polygon(c(days,rev(days)),c(Clb,rev(Cub)),col="grey75",border=FALSE)
  lines(days,CSave[,1],lwd=2,col="blue")
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  legend("topleft",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
  dev.off()
  
  # plot population dynamics
  filepath=file.path(this_path,"/figures/outputs/Population.jpeg")
  jpeg(filepath,800,600)
  plot(days2, N, ylab="Number of individuals", xlab="", xaxt = "n",type="l",cex.lab=1.4)
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()
  
  # plot O2 consumption
  filepath=file.path(this_path,"/figures/outputs/O2_consumption.jpeg")
  jpeg(filepath,800,600)
  ub=O2Save[,1]+O2Save[,2]
  lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
  for (i in 1:length(O2Save[,1]-O2Save[,2])){
    lb[i]=max(O2Save[i,1]-O2Save[i,2],0)
  }
  maxub=max(O2Save[,1]+O2Save[,2])
  plot(days,O2Save[,1],ylab="O2 consumption (kgO2/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
  polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
  lines(days,O2Save[,1],lwd=2,col="red")
  lines(days,lb,col="blue")
  lines(days,ub,col="blue")
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()
  
  # plot NH4 production
  filepath=file.path(this_path,"/figures/outputs/NH4_release.jpeg")
  jpeg(filepath,800,600)
  ub=NH4Save[,1]+NH4Save[,2]
  lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
  for (i in 1:length(NH4Save[,1]-NH4Save[,2])){
    lb[i]=max(NH4Save[i,1]-NH4Save[i,2],0)
  }
  maxub=max(NH4Save[,1]+NH4Save[,2])
  plot(days,NH4Save[,1],ylab="NH4 release (kgN/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
  polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
  lines(days,NH4Save[,1],lwd=2,col="red")
  lines(days,lb,col="blue")
  lines(days,ub,col="blue")
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()
  
  # Results save
  filepath=file.path(this_path,"/data_products/weight.csv")
  write.csv(weightSave,filepath)
  
  filepath=file.path(this_path,"data_products/faeces_production_Proteins.csv")
  write.csv(PexcSave,filepath)
  
  filepath=file.path(this_path,"/data_products/faeces_production_Lipids.csv")
  write.csv(LexcSave,filepath)
  
  filepath=file.path(this_path,"/data_products/faeces_production_Carbohydrates.csv")
  write.csv(CexcSave,filepath)
  
  filepath=file.path(this_path,"data_products/wasted_feed_Proteins.csv")
  write.csv(PwstSave,filepath)
  
  filepath=file.path(this_path,"/data_products/wasted_feed_Lipids.csv")
  write.csv(LwstSave,filepath)
  
  filepath=file.path(this_path,"/data_products/wasted_feed_Carbohydrates.csv")
  write.csv(CwstSave,filepath)
  
  filepath=file.path(this_path,"/data_products/actual_ingestion.csv")
  write.csv(ingestionSave,filepath)
  
  filepath=file.path(this_path,"/data_products/anabolic_rate.csv")
  write.csv(ASave,filepath)
  
  filepath=file.path(this_path,"/data_products/catabolic_rate.csv")
  write.csv(CSave,filepath)
  
  tfun=cbind(fgT,frT)
  
  filepath=file.path(this_path,"/data_products/temperature_response.csv")
  write.csv(tfun,filepath)
  
  filepath=file.path(this_path,"/data_products/O2_consumption.csv")
  write.csv(O2Save,filepath)
  
  filepath=file.path(this_path,"/data_products/NH4_release.csv")
  write.csv(NH4Save,filepath)
  
  filepath=file.path(this_path,"/data_products/population.csv")
  write.csv(N,filepath)
  
  filepath=file.path(this_path,"/data_products/Days_to_commercial_size.csv")
  write.csv(daysToSize,filepath)
  
  return(output)
  
}













# Bioenergetic population model

model <- \(this_path, forcings, this_species){
  
  rm(list=ls())           # Clean workspace
  
  cat('Population bioenergetic model for ', str_to_sentence(this_species), "\n")
  cat(" \n")
  

  # Run the preprocessor for the first time to print to screen parameters and forcing selected
  out_pre<-preprocess(this_path = this_path, these_forcings = forcings)
  
  # While cycle to repeat the pre-processing until correct inputs are inserted
  selector="y"
  
  while (identical(selector,"y")=="TRUE") {
    cat(" \n")
    selector=readline("Do you want to change the inputs? [y/n]")
    
    if (identical(selector,"n")=="TRUE") {break}
    
    cat(" \n")
    cat("Insert forcings and parameters in the following folder\n")
    cat(file.path(this_path))
    cat(" \n")
    cat("Type y if you entered the correct inputs\n")
    cat("The data will be preprocessed again")
    selector=readline(" ")
    
    out_pre<-preprocess(this_path,forcings)
    selector="y"
  }
  
  # Extract preprocessor outputs
  Param=out_pre[[1]]
  Temp=out_pre[[2]]
  Food=out_pre[[3]]
  IC=out_pre[[4]]
  times=out_pre[[5]]
  Dates=out_pre[[6]]
  N=out_pre[[7]]
  CS=out_pre[[8]]
  
  
  
  
  # Manages population
  out_RKsolver <- loop(Param, Tint, Food, IC, times, N, this_path)
  
  # Post-process data
  out_post<-post_process(this_path, out_RKsolver, times, Dates,N, CS)
  
  cat(" ")
  cat("End")
  
  return(out_post)
  
  
}



