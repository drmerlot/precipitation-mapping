######################################################################
#                                                                    #
#           Written By Andrew R Sommerlot                            #
#           Version 0.1.0                                            #
#           September 16 2016                                        #
#           <andy.sommerlot@gmail.com>                               #
#           In Brief: A Group of functions gathering precip data     #
#              from the GHCN network, making thiesen polygons, and   #
#              extracting aggregated precip to input polygons.       #
#                                                                    #
######################################################################

## first is the ghcn gathering function .....
#Requires annoying install due to archieved packages ....
install.packages(c('R.utils', 'R.oo', 'R.methodsS3', 'abind', 'ncdf'))
path_to_file = '~/Desktop/umces/cso/ncdf_1.6.7.tar.gz'
install.packages(path_to_file, repos = NULL, type="source")
path_to_file = '~/Desktop/umces/cso/GhcnDaily_1.5.tar.gz'
install.packages(path_to_file, repos = NULL, type="source")



getGHCN <- function(lat, long, NSdeg=0.5, EWdeg=0.5, StartYear= 2000, EndYear=(as.POSIXlt(Sys.time())$year+1900), StnOpt = "Any" ){
  require(SWATmodel)
  require(GhcnDaily)
  require(plyr)
  
  read_ghcn_raw=function(GHCNfilename,startyear=StartYear,elements=c("TMIN","TMAX","PRCP")){
    require(Hmisc)
    suppressWarnings(rm(list=c("alldata")))
    suppressWarnings(rm(list=objects(pattern="allghcndata")))
    
    for (fileline in readLines(GHCNfilename)){
      stn=substr(fileline,1,11)
      year=substr(fileline,12,15)
      mo=substr(fileline,16,17)
      yearmo=paste(year,mo,sep="")
      element=as.character(substr(fileline,18,21))
      
      if(year<startyear | (length(grep(element,elements))<1)){next()}
      
      dy=1:monthDays(paste(year,mo,"01",sep="/"))
      date=as.Date(format="%Y%m/%d",paste(yearmo,"/",dy,sep=""))
      tempdata=data.frame(date,var=as.numeric(as.character(substring(fileline,c(seq(22,((length(dy)-1)*8+22),8)),c(seq(26,((length(dy)-1)*8+26),8))))))
      names(tempdata)=c("date",element)
      
      if(!(exists(paste("allghcndata",element,sep="")))){
        assign(paste("allghcndata",element,sep=""),tempdata)
      } else {
        assign(paste("allghcndata",element,sep=""),rbind(get(paste("allghcndata",element,sep="")),tempdata))
      }
    }
    
    for (dfnames in objects(pattern="allghcndata")){
      print(dfnames)
      
      if(exists("alldata")){alldata=merge(alldata,get(dfnames),by="date",all=T)}
      else {alldata=get(dfnames)}
    }
    
    return(alldata)
  }
  
  if (!("ghcndailytmp.dly" %in% list.files()) ) downloadDailyInventory()  # only download if needed
  junk=readDailyInventory()
  SubJunk <- subset(junk,Lat>(lat-NSdeg)&Lat<(lat+NSdeg)&Lon>(long-EWdeg)&Lon<(long+EWdeg)&LastYear>=EndYear & FirstYear < StartYear & (Element =="TMAX" | Element=="TMIN" | Element=="PRCP"))

  if (StnOpt == "Full"){   #  Only returns stations that have max/min temp and precip
    GHCNstns <- subset(count(subset(SubJunk,Lat>(lat-NSdeg)&Lat<(lat+NSdeg)&Lon>(long-EWdeg)&Lon<(long+EWdeg)& LastYear>=EndYear & FirstYear < StartYear & (Element =="TMAX" | Element=="TMIN" | Element=="PRCP"))[1]) ,freq==3)
  } else GHCNstns <- SubJunk  #  Otherwise returns all stations with any information
  
  if (nrow(GHCNstns) < 1){ 
    print("no GHCN stations within your search window, you could try making NSdeg or EWdeg larger to expand it")
  } else {
    GHCN <-list()
  
    if (!("ghcnd-stations.txt" %in% list.files()) ) download.file(GHCN.DAILY.METADATA.URL, destfile="ghcnd-stations.txt")
    ghcnSta = readLines("ghcnd-stations.txt") 
    
    for (b in unique(GHCNstns$Id)){
      download.file(paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/", b, ".dly", sep=""), "ghcndailytmp.dly")
      GHCN[[b]] <- read_ghcn_raw("ghcndailytmp.dly")
    }	
  
    Stns <- names(GHCN)    
    Lats <- SubJunk$Lat[sapply(Stns, function(x) min(grep(x, SubJunk$Id)))]
    Lons <- SubJunk$Lon[sapply(Stns, function(x) min(grep(x, SubJunk$Id)))]
    subst =   ghcnSta[sapply(Stns, function(n) grep(n, ghcnSta))] 
    ghcnElev = 	as.numeric(sapply(strsplit(subst, " +"), function(s) s[4]))
    StationInfo = data.frame(ID=Stns, Lat=Lats, Lon=Lons, Elev=ghcnElev)
    
    #GHCN[[b]]$Lat <- SubJunk$Lat[min(which(SubJunk$Id == b))]    ##  Add a column reporting decimal latitude
    #GHCN[[b]]$Lon <- SubJunk$Lon[min(which(SubJunk$Id == b))]	 ##  Add a column reporting decimal longitude
    #GHCN[[b]]$Elev <- as.numeric(ghcnElev)          ##  Add a column reporting decimal longitude
  
    GHCN[['StationInfo']] = StationInfo
    
    return(GHCN)	
  }
}


### read output of GHCN functio
read_ghcn = function(test_gage, cut_date = "2015-12-31", write_csv = TRUE){

  all_dat = test_gage
  all_dat[length(all_dat)] <- NULL 

  datefilter = function(dataf, cut_date){
    if(length(which(dataf$date == as.Date(as.Date(cut_date)))) != 0){
      return(dataf)
    }else{
      return(NA)
    }
  }

  data = lapply(all_dat, datefilter, cut_date = cut_date)
  locs = test_gage[[length(test_gage)]][!is.na(data),]
  data = data[!is.na(data)]

  cutdate = function(dataf, cut_date){
    cut <- which(dataf$date == as.Date(cut_date))
    cut_frame <- dataf[1:cut,]
    return(cut_frame)
  }


  make.na = function(dataf, val){
    dataf[dataf == val] <- NA
    return(dataf)
  }

  data = lapply(data, FUN = cutdate, cut_date = cut_date)
  data = lapply(data, make.na, val = -9999)

  # remove any data with less than the max rows 
  max_r = max(unlist(lapply(data, nrow)))

  maxonly = function(dataf, max_r){
    if(nrow(dataf) != max_r){
      return(NA)
    } else{
      return(dataf)
    }
  }

  data = lapply(data, maxonly, max_r = max_r)
  locs = locs[!is.na(data),]
  data = data[!is.na(data)]

  getpcp = function(dataf){
    pcp = dataf$PRCP
    return(pcp)
  }

  pcp = lapply(data, getpcp)

  ## loop through the rows will get all stations values per day.. 
  pcp_df = as.data.frame(pcp)
  #dates = all_dat[[1]]$date[1:which(all_dat[[1]]$date == cut_date)]
  dates = data[[1]]$date

  if(write_csv == TRUE) {
    write.csv(pcp_df, 'pcp_df_84_16_st_ghcn.csv', row.names = FALSE)
    write.csv(dates, 'dates_84_16_st_ghcn.csv', row.names = FALSE)
    write.csv(locs, 'locs_84_16_st_ghcn.csv', row.names = FALSE)
  }
  
  pcp <- split(pcp_df, seq(nrow(pcp_df)))
  
  all_out = list(pcp, dates, locs)
  names(all_out) = c('pcp', 'dates', 'locs')
  
  return(all_out)
}


##################################
# process

make_df = function(pcp_l, test) {
  id = test$locs$ID
  x = test$locs$Lon
  y = test$locs$Lat
  df <- data.frame(id = test$locs$ID, precip = as.numeric(pcp_l), x = test$locs$Lon, y = test$locs$Lat)
  return(df)
}


# remove NA values.. 
remove_nas = function(pcp1){
  pcp1 = na.omit(pcp1)
}

## don't need for thiesen method, or any shapefile defined sp distribution. 
drop_zeros = function(pcp1, drop_frac = 0.8, write_csv = FALSE){
  zeros = which(pcp1$precip == 0)
  drops = sample(zeros, round(length(zeros)* drop_frac), replace=FALSE)
  pcp1 = pcp1[-drops,]
  
  if(write_csv = TRUE){
    write.csv(pcp1, 'pcp1.csv', row.names =  FALSE)
  }
  return(pcp1)
}

#################################
## spatial functions....
#install.packages('gstat')
# Load precipitation data
#gauges$precip[which(is.na(gauges$precip))] = 0
# Put data into a Spatial Data Frame
gauge_spatial = function(pcp1){
  require(sp)
  require(rgdal)
  require(gstat)
  gauges = pcp1
  colnames(gauges) = c("gauge_id","precip","x","y")
  gauges.xy <- gauges[c("x","y")]
  coordinates(gauges) <- gauges.xy
  # Set the correct CRS for the gauges data
  proj4string(gauges) <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
  
  return(gauges)
}

plot_gauges = function(gauges, item = 1){
  plot(all_gauges[[item]], pch=16, col='blue', cex=all_gauges[[item]]$precip/30)
}


## now the thiessen

thiessen_poly <- function(pcp1, plot_out = FALSE) {
  require(dismo)
  points = pcp1[,c(3,4)]
  vor <- voronoi(points)
  proj4string(vor) <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
  vor$precip <- pcp1$precip
  vor$id <- pcp1$id
  
  if(plot_out == TRUE){
    spplot(vor, zcol="precip", col.regions=colorRampPalette(brewer.pal(9,"Blues"))(100))
  }
  return(vor)
}

## now the functions for the spatial extraction.
# !! call !!

# define the TT cso equation .... 
# outside the other functions, because this makes it CBP specific ....
cso_est = function(x){
  mm_out = x / 10
  inches = mm_out*0.03937
  cso_est = 1567*inches^2 - 46.955*inches + 1309.8
  #cso_est[which(x <= threshold)] <- 0
  cso_est = cso_est / 1e6
  return(cso_est)
}

ghcn_default_mm = function(x){ 
  mm_out = x / 10
  return(mm_out)
}

mm2in = function(x) {
  inches = x*0.03937
}


extract_pcp = function(vor, cso, pcp_unit = 'mm', threshold = 0.01, write_csv = FALSE){
  ## prject vor to same as the input shape.. 
  vor <- spTransform(vor, proj4string(cso))
  ### get the precip data ... 
  ## extracts the precip value for each cso area, id? 
  p <- over(cso, vor, fun = mean, na.rm = FALSE)

  # Make the output data frame, all ids, one day ... 
  output1 = data.frame(id = cso$NPDES_ID, pcp = p$precip, ac = cso$ACRES)
  output1$pcp = ghcn_default_mm(output1$pcp)
  
  if(pcp_unit == 'in'){
    output1$pcp = mm2in(output1$pcp)
    output1$pcp[which(output1$pcp <= 0.01)] <- 0
  } 
  if(pcp_unit == 'mm') {
    output1$pcp[which(output1$pcp <= 0.01)] <- 0
  }
  
  if(write_csv == TRUE){
    # here is one-days output......
    write.csv(out, paste('station_pcp', 'cso_overflow.csv', sep = '_'), row.names = FALSE)
  }

  return(output1)
}

finish_pcp <- function(out, dates, write_csv = TRUE){ 
  for( i in 1:length(dates)){
    out[[i]]$date = dates[i]
  }
  output = do.call('rbind', out)
  if( write_csv == TRUE) {
    write.csv(output, paste(getwd(), 'output_pcp.csv', sep = '/'), row.names = FALSE)
  }
  return(output)
}

plot_extract_pcp = function(vor, dates, cso, date_lab = dates[[1]], pad.x = .5, pad.y = .5) {
  # reproject the data onto a "longlat" projection
  vor2 = vor[[which(dates == date_lab)]]
  subsetTransform <- spTransform(vor2, CRS("+proj=longlat"))
  cso_trans <- spTransform(cso, CRS("+proj=longlat"))
  subsetTransform$id = rep(1:nrow(subsetTransform))
  
  # determine the bounding box of the spatial object
  b <- bbox(subsetTransform)
  
  # get and plot a map
  cb_watershed <- ggmap(get_map(location = b, maptype = "satellite", zoom = 6))
  
  #
  subsetTransformFortified <- fortify(subsetTransform, region = "id")
  subsetTransformFortified <- merge(subsetTransformFortified,
                                    subsetTransform@data, by = "id")
  # convert to inches
  subsetTransformFortified$precip = subsetTransformFortified$precip / 10 * 0.03937
  
  cso_overlay = fortify(cso_trans, region = "NPDES_ID")
  cso_overlay  <- merge(cso_overlay,
                        cso_trans@data, by.x = "id", by.y = "NPDES_ID")

  big_plot = cb_watershed + geom_polygon(data = subsetTransformFortified,
                              aes(x = long, y = lat, group = group,
                                  fill = precip), alpha = 0.5 ) +
    scale_x_continuous(limits = c(b[1, 1] - pad.x, b[1,2] + pad.x)) +
    scale_y_continuous(limits = c(b[2,1] - pad.y, b[2,2] + pad.y)) +
    #labs(title= paste("Thiesen Poly GHCN Precipitation", date_lab, sep = ' ')) +
    scale_fill_gradient(name="Precip in", limits=c(0,max(subsetTransformFortified$precip)), low="white", high="#0000ff") + 
    geom_polygon(data = cso_overlay,
                 aes(x = long, y = lat, group = group, fill = shape), 
                 fill = 'darkred', alpha = 0.9, show.legend = TRUE) 
  print(big_plot)
}