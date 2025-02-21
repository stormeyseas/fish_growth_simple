#Template rasters


base_rast <- rast(res = 0.5)
values(base_rast) <- 1:ncell(base_rast)


gall_rast <- rast(ext(-20037508, 20037508, -10000000, 10000000), 
                  resolution = 10000, 
                  crs = "ESRI:54016")

values(gall_rast) <- 1:ncell(gall_rast)


