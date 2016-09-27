#############################
#
#
#
#
#
############################
install.packages('pblappy')
install.packages('sp')
install.packages('rgdal')
install.packages('gstat')

require(pblappy)
require(sp)
require(rgdal)
require(gstat)
require(parallel)
# read in the pcp data from the mod  

setwd('~/Desktop/pcp')
all_files = list.files('.')
all_frames = lapply(all_files,read.csv)
mod_pcp_all = do.call(rbind, all_frames)

mod_pcp = mod_pcp_all[,-c(6,7,8,9,10,11,12)]
dates = paste(mod_pcp$year, mod_pcp$month, mod_pcp$day, sep = '-')
dates = as.Date(dates)
mod_pcp$dates = dates
mod_pcp$year = NULL
mod_pcp$month = NULL
mod_pcp$day = NULL

mod_order = mod_pcp[order(mod_pcp$dates),]
mod_order$sp = as.factor(mod_order$dates)

dates = mod_order$dates

pcp_days = split(mod_order, mod_order$sp)

directory = "/Users/grad/Desktop/lseg"
layer_name = "P6Beta_v3_LSegs_081516_Albers"
lseg = readOGR(dsn = directory, layer = layer_name)

lseg_day_pcp = function(pcp_days, lseg){
  ids = lseg@data$FIPS_NHL
  lseg@data = merge(lseg@data, pcp_days, by.x = 'FIPS_NHL', by.y = 'seg', all = TRUE)
  lseg@data$sp = NULL
  lseg@data$dates[is.na(lseg@data$dates)] = lseg@data$dates[1]
  lseg@data$precip.in.[is.na(lseg@data$precip.in.)] = 0
  lseg@data$FIPS_NHL = as.character(ids)
  return(lseg)
  
}

lseg_pcp = lapply(pcp_days, lseg_day_pcp, lseg)

################################
# plots
lseg_plot = lseg
lseg_plot@data$dates = NULL
spplot(lseg_plot)
###############################

directory = "/Users/grad/Desktop/umces/cso/css"
layer_name = "CSS_final_jan14_clean"
cso = readOGR(dsn = directory, layer = layer_name)

### read in the land segment shapefile.
extract_pcp_mod = function(lseg, cso, pcp_unit = 'in', threshold = 0.01, write_csv = FALSE){
  ## prject vor to same as the input shape.. 
  lseg <- spTransform(lseg, proj4string(cso))
  ### get the precip data ... 
  ## extracts the precip value for each cso area, id? 
  p <- over(cso, lseg[2], layer = 3, fun = mean, na.rm = FALSE)
  
  # Make the output data frame, all ids, one day ... 
  output1 = data.frame(id = cso$NPDES_ID, pcp = p$precip.in., ac = cso$ACRES )
  
  if(pcp_unit == 'mm'){
    print('currently not supported')
  } 
  if(pcp_unit == 'in') {
    output1$pcp[which(output1$pcp <= 0.01)] <- 0
  }
  
  if(write_csv == TRUE){
    # here is one-days output......
    write.csv(out, paste('station_pcp', 'cso_mod_pcp.csv', sep = '_'), row.names = FALSE)
  }
  
  return(output1)
}

#  extraction the pcp ... 
mod_pcp = pblapply(lseg_pcp, extract_pcp_mod, cso)

mod_all_pcp = do.call('rbind', mod_pcp)

write.csv(mod_pcp, '~/Desktop/pcp/modeled_pcp.csv')

# outside the other functions, because this makes it CBP specific ....
cso_est = function(x){
  mm_out = x / 10
  inches = mm_out*0.03937
  cso_est = 1567*inches^2 - 46.955*inches + 1309.8
  cso_est = cso_est / 1e6
  return(cso_est)
}

#calculate cso output on outp
output$cso = cso_est(output$pcp) * output$ac 

out1 = aggregate.data.frame(output$pcp, by = list(output$id), FUN = sum)
out2 = aggregate.data.frame(output$ac, by = list(output$id), FUN = sum)
out = merge(out1, out2, by = 'Group.1')
colnames(out) = c('id', 'pcp', 'ac')

outp = out[order(out$id),]

# write it separte,
write.csv(outp, paste(getwd(), 'output_pcp.csv', sep = '/'), row.names = FALSE)
# end call 

## then apply the cso equation and aggregation ... 


##########################
#observing the results 
cso_orig = read.csv('~/Desktop/cso.csv')
tetra_tech = cso_orig[,5]

dat = data.frame(tetra_tech, mod)

plot_ecdf(dat)

NSE(mod, tetra_tech)

library(ggplot2)
ggplot(dat, aes(x=mod, y=tetra_tech)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  xlab('mod')


plot_ecdf <- function(dat, legend_names = colnames(dat), title = 'ECDFs', xlab = 'value', ylab = 'probability', out_type = 'display', file_name = 'ecdf-plot', out_dir = getwd()){
  require(ggplot2)
  require(reshape2)
  
  ecdfs <- list()
  for(i in 1:ncol(dat)){
    ecdfs[[i]] = ecdf(dat[,i])
  }
  
  ecdf_frames <- list()
  for(i in 1:length(ecdfs)){
    now_dat <- data.frame(dat[,i], ecdfs[[i]](dat[,i]))
    ecdf_frames[[i]] <- now_dat[order(now_dat[,1]),]
  }
  
  for(i in 1:length(ecdf_frames)){
    colnames(ecdf_frames[[i]]) <- c('x', 'y')
  }
  
  melter <- data.frame(matrix(ncol = length(ecdf_frames) + 1, nrow = nrow(ecdf_frames[[1]])))
  melter[,1] <- ecdf_frames[[1]][,2]
  
  for(i in 1:length(ecdf_frames)){
    melter[,i+1] <- ecdf_frames[[i]][,1]
  }
  
  colnames(melter) <- c('X1', legend_names)
  
  data_long <- melt(melter, id='X1')
  
  plot_ecdf <- ggplot(data = data_long,
                      aes(x = value, y = X1, colour = variable)) +
    labs(title = title, x = xlab, y = ylab) +
    geom_line()
  
  if(out_type == 'display'){
    plot(plot_ecdf)
    print('Printed plot to display device')
  }
  
  if(out_type == 'file'){
    invisible(plot_ecdf, file=paste0(paste(out_dir, file_name, sep = '/'), ".png"), plot=plot_qq)
    print(paste('Saved plot to ', out_dir, sep = ''))
  }
}