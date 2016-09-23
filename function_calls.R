# calls for cso precip.. 


setwd('/Users/grad/Desktop/umces/cso/precip')

test_gage = getGHCN(lat = 39.9, lon = -77.18, NSdeg = 3, EWdeg = 3, StartYear = 1984, EndYear = 2015, StnOpt = "Full")

test = read_ghcn(test_gage, cut_date = "2015-12-31", write_csv = FALSE)

pcp_l = test$pcp

tl = lapply(pcp_l, make_df, test)

require(pbapply)
tl = pblapply(tl, remove_nas)

# Have a look
require(pbapply)
all_gauges = pblapply(tl, gauge_spatial)

# look at a gauge dist
## !! call
plot_gauges(all_gauges, item = 2000)
# !! end call

## !! call
require(pbapply)
vor_l = pblapply(all_gauges, thiessen_poly)

# !! call, define outside to save time ...
direcotry = "/Users/grad/Desktop/umces/cso/css"
layer_name = "CSS_final_jan14_clean"
cso = readOGR(dsn = directory, layer = layer_name)
## !! end call

# !! call 
cso = readOGR(dsn = directory, layer = layer_name)
require(pbapply)

out = mclapply(vor_l, extract_pcp, cso = cso,  write_csv = FALSE)
# !! end call 

# !! call 
output = finish_pcp(out, test$dates, write_csv = FALSE)

#order it if you want 
outp = output[order(output$id),]

# cso conv function. 
# define the TT cso equation .... 
# outside the other functions, because this makes it CBP specific ....
cso_est = function(x, threshold = 0.01){
  mm_out = x / 10
  inches = mm_out*0.03937
  cso_est = 1567*inches^2 - 46.955*inches + 1309.8
  cso_est[which(x <= threshold)] <- 0
  cso_est = cso_est / 1e6
  return(cso_est)
}


#calculate cso output on outp
outp$cso = cso_est(outp$pcp) * outp$ac 

# write it separte,
write.csv(outp, paste(getwd(), 'output_pcp.csv', sep = '/'), row.names = FALSE)
# end call 

##plot the min and max number of stations, and min and max rainfall
# number of stations 
lens = lapply(tl, nrow)

mins = which(lens == min(unlist(lens)))
plot.min = mins[[length(mins)]]
date.min = dates[[plot.min]]

maxs = which(lens == max(unlist(lens)))
plot.max = maxs[[length(maxs)]]
date.max = dates[[plot.max]]


theme_set(theme_gray(base_size = 22))
plot_extract_pcp(vor = vor_l, dates = dates, cso = cso, date_lab = dates[[which(dates == date.min)]])




