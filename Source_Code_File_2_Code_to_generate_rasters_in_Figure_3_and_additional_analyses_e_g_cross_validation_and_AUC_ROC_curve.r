######################################################
### Code to generate Figure 3 and related analysis
### #################################################

# N.b. the raster files were further edited in ArcGIS
# for improved aesthetic appearance


############
##libraries
###########
require(raster)
require(rgdal)
library(INLA)
library(RColorBrewer)
library(cvTools)
library(zoo)
library(boot)
library(pROC)
library(ROCR)



###########
#data load
##########
#define/load dataset containing location and time info

dd <- "INSERT PATH/OBJECT CONTAINING LOCATION AND TIME INFORMATION "

#define/load Rt (response)
Rt <- "INSERT PATH/OBJECT CONTAINING RESPONSE (REPRODUCTION NUMBERS)"
dd$Rt <- Rt


#Rt<-load("C:/Users/ir515/Dropbox/El_Sal_transmission_chain_reconstruction/spatial_netrate_V4_2908/out/sd_mid_vivax_yunnan.RData")
dd$x <- "COLUMN/VECTOR OF LONGITUDES"
dd$y <- "COLUMN/VECTOR OF LATITUDES"

#load covariate rasters
lsf <- c(acc, elev, lst_day, lst_night, lst_diff, precip, urban)




###################
#functions required
###################
#auc

fun.auc <-
  function(pred, obs) {
    # Run the ROCR functions for AUC calculation
    ROC_perf <- performance(prediction(pred, obs), "tpr", "fpr")
    ROC_sens <- performance(prediction(pred, obs), "sens", "spec")
    ROC_err <- performance(prediction(pred, labels = obs), "err")
    ROC_auc <- performance(prediction(pred, obs), "auc")
    # AUC value
    AUC <- ROC_auc@y.values[[1]] # AUC
    # Mean sensitivity across all cutoffs
    x.Sens <- mean(as.data.frame(ROC_sens@y.values)[, 1])
    # Mean specificity across all cutoffs
    x.Spec <- mean(as.data.frame(ROC_sens@x.values)[, 1])
    # Sens-Spec table to estimate threshold cutoffs
    SS <-
      data.frame(SENS = as.data.frame(ROC_sens@y.values)[, 1],
                 SPEC = as.data.frame(ROC_sens@x.values)[, 1])
    # Threshold cutoff with min difference between Sens and Spec
    SS_min_dif <-
      ROC_perf@alpha.values[[1]][which.min(abs(SS$SENS - SS$SPEC))]
    # Threshold cutoff with max sum of Sens and Spec
    SS_max_sum <-
      ROC_perf@alpha.values[[1]][which.max(rowSums(SS[c("SENS", "SPEC")]))]
    # Min error rate
    Min_Err <- min(ROC_err@y.values[[1]])
    # Threshold cutoff resulting in min error rate
    Min_Err_Cut <-
      ROC_err@x.values[[1]][which(ROC_err@y.values[[1]] == Min_Err)][1]
    # Kick out the values
    round(cbind(
      AUC,
      x.Sens,
      x.Spec,
      SS_min_dif,
      SS_max_sum,
      Min_Err,
      Min_Err_Cut
    ),
    3)
  }

# AUC plot function
fun.aucplot <- function(pred, obs, title) {
  # Run the AUC calculations
  ROC_perf <- performance(prediction(pred, obs), "tpr", "fpr")
  ROC_sens <- performance(prediction(pred, obs), "sens", "spec")
  ROC_err <- performance(prediction(pred, labels = obs), "err")
  ROC_auc <-
    performance(prediction(pred, obs), "auc")  # Spawn a new plot window (Windows OS)
  # Plot the curve
  plot(
    ROC_perf,
    colorize = T,
    print.cutoffs.at = seq(0, 1, by = 0.1),
    lwd = 3,
    las = 1,
    main = title
  )
  # Add some statistics to the plot
  text(1,
       0.25,
       labels = paste("Npres = ", sum(obs == 1), sep = ""),
       adj = 1)
  text(1,
       0.20,
       labels = paste("Nabs = ", sum(obs == 0), sep = ""),
       adj = 1)
  text(1,
       0.15,
       labels = paste("AUC = ", round(ROC_auc@y.values[[1]], digits = 2), sep =
                        ""),
       adj = 1)
  text(1,
       0.10,
       labels = paste("Sens = ", round(mean(
         as.data.frame(ROC_sens@y.values)[, 1]
       ), digits = 2), sep = ""),
       adj = 1)
  text(1,
       0.05,
       labels = paste("Spec = ", round(mean(
         as.data.frame(ROC_sens@x.values)[, 1]
       ), digits = 2), sep = ""),
       adj = 1)
}

#data shaping functions

#latlong to xyz
ll.to.xyz <- function(ll) {
  if (is.null(colnames(ll))) {
    colnames(ll) <- c('longitude', 'latitude')
  }
  if (colnames(ll)[1] == 'x' & colnames(ll)[2] == 'y') {
    colnames(ll) <- c('longitude', 'latitude')
  }
  if (colnames(ll)[1] == 'lon' & colnames(ll)[2] == 'lat') {
    colnames(ll) <- c('longitude', 'latitude')
  }
  ll[, 'longitude'] <- ll[, 'longitude'] * (pi / 180)
  ll[, 'latitude'] <- ll[, 'latitude'] * (pi / 180)
  x = cos(ll[, 'latitude']) * cos(ll[, 'longitude'])
  y = cos(ll[, 'latitude']) * sin(ll[, 'longitude'])
  z = sin(ll[, 'latitude'])
  return(cbind(x, y, z))
}

#temporal data clean
temporalInfo <- function(f) {
  library(zoo)
  start_dates <-
    as.yearmon(paste(f[, "year_start"], "-", f[, "month_start"], sep = ""))
  end_dates <-
    as.yearmon(paste(f[, "year_end"], "-", f[, "month_end"], sep = ""))
  mid_dates <- (start_dates + end_dates) / 2
  return(as.numeric(mid_dates))
}

#colour matching
match.cols <- function(val, n) {
  colfunc <-
    colorRampPalette(
      c(
        "snow1",
        "snow2",
        "snow3",
        "seagreen",
        "orange",
        "firebrick",
        "darkred"
      ),
      space = "rgb",
      bias = 5
    )#colors
  #	colfunc <- colorRampPalette(c("blue","cyan","yellow","orange","red"), space = "rgb",bias=1)
  
  col <-
    data.frame(val = seq(min(val), max(val), length.out = n), col = colfunc(n))
  out <- rep(NA, length(col))
  for (i in 1:length(val)) {
    out[i] <- as.character(col[which.min(abs(col$val - val[i])), 'col'])
  }
  return(out)
}

#alternative latlong to xyz
ll2xyz <- function (longlat, radius = 6371) {
  # check inputs
  stopifnot(is.finite(radius) & radius > 0)
  stopifnot(inherits(longlat, 'matrix') |
              inherits(longlat, 'data.frame'))
  stopifnot(ncol(longlat) == 2)
  
  # extract required columns
  longitude <- longlat[, 1]
  latitude <- longlat[, 2]
  
  # convert
  ans <- data.frame(
    x = radius * cos(latitude) * cos(longitude),
    y = radius * cos(latitude) * sin(longitude),
    z = radius * sin(latitude)
  )
  
  return (ans)
  
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#xyz to latlong
xyz2ll <- function (xyz, radius = 6371) {
  # check inputs
  stopifnot(is.finite(radius) & radius > 0)
  stopifnot(inherits(longlat, 'matrix') |
              inherits(longlat, 'data.frame'))
  stopifnot(ncol(xyz) == 3)
  
  # extract columns
  x <- xyz[, 1]
  y <- xyz[, 2]
  z <- xyz[, 3]
  
  # convert
  ans <- data.frame(longitude = atan2(y, x),
                    latitude = asin(z / radius))
  
  return (ans)
  
}
#trimming
trim2 <- function(x, out = "raster") {
  if (!any(out == c("matrix", "raster")))
    stop("output must be a matrix or raster")
  if (class(x) == "matrix" &
      out == "raster")
    stop("if you supply a matrix, you must use out='matrix'")
  if (class(x) == "RasterLayer") {
    if (out == "raster") {
      cres <- 0.5 * res(x)
      crs <- projection(x)
      y <- x
    }
    x <- matrix(as.array(x), nrow = nrow(x), ncol = ncol(x))
  }
  if (class(x) != "matrix") {
    stop("x must be a matrix or raster")
  } else {
    r.na <- c.na <- c()
    for (i in 1:nrow(x))
      r.na <- c(r.na, all(is.na(x[i,])))
    for (i in 1:ncol(x))
      c.na <- c(c.na, all(is.na(x[, i])))
    r1 <-
      1 + which(diff(which(r.na)) > 1)[1]
    r2 <- nrow(x) -  which(diff(which(rev(r.na))) > 1)[1]
    c1 <-
      1 + which(diff(which(c.na)) > 1)[1]
    c2 <- ncol(x) - which(diff(which(rev(c.na))) > 1)[1]
    x <- x[r1:r2, c1:c2]
    if (out == "raster") {
      xs <- xFromCol(y, col = c(c1, c2)) + c(-1, 1) * cres[1]
      ys <- yFromRow(y, row = c(r2, r1)) + c(-1, 1) * cres[2]
      x <- raster(
        x,
        xmn = xs[1],
        xmx = xs[2],
        ymn = ys[1],
        ymx = ys[2],
        crs = crs
      )
    }
  }
  return(x)
}


#colour scale function
colfunc <-
  colorRampPalette(
    c(
      "snow1",
      "snow2",
      "snow3",
      "seagreen",
      "orange",
      "firebrick",
      "darkred"
    ),
    space = "rgb",
    bias = 1
  )

###########
#data shape
###########

#Stack and extract raster
st <- stack(lsf)
cn <- calc(st, sum)
cn[!is.na(cn)] = 1
#trim2(x=cn)
NAvalue(cn) = -9999
wh <- extract(cn, cbind(dd$x, dd$y))

#cropping and formatting stack
NAvalue(st) = -9999
st <- crop(st, cn)
#st[[1]]=asinh(st[[1]])

#cells from xy coords - match to covariate rasters
dd$cell <- cellFromXY(cn, cbind(dd$x, dd$y))
covs <- st[dd$cell]
covs <- as.data.frame(covs)

#make matrix of response and covs
mat <- cbind(dd$Rt, covs)

#exclude missingness
Rt <- Rt[complete.cases(mat)]
covs <- covs[complete.cases(mat),]

#df <- data.frame(matrix(unlist(dd), nrow=nrow(dd), byrow=T),stringsAsFactors=FALSE)
dd <- dd[complete.cases(mat),]

###########################
##exploratory vizualisation
###########################


par(mfrow = c(1, 2))
cols <- match.cols(Rt, 1000)
image(cn, col = 'grey')
points(dd$x, dd$y, pch = 16, col = cols)

image(cn, col = 'grey')
points(dd$x[dd$I == 1], dd$y[dd$I == 1], pch = 16, col = 'red')
points(dd$x[dd$I == 0], dd$y[dd$I == 0], pch = 16, col = 'blue')

##############
#INLA model
##############

# make mesh

mesh1d = inla.mesh.1d(
  seq(2010, 2017, by = 1),
  interval = c(2010, 2017),
  degree = 2,
  boundary = c('free')
)

xyz <- as.data.frame(ll.to.xyz(cbind(dd$x, dd$y)))
mesh = inla.mesh.2d(
  loc = cbind(xyz['x'], xyz['y'], xyz['z']),
  cutoff = 0.005,
  #should be 0.003
  min.angle = c(25, 25),
  max.edge = c(0.005, 1)
)

# spde function
spde = inla.spde2.matern(mesh, alpha = 2)

A.est =
  inla.spde.make.A(
    mesh,
    loc = cbind(xyz[, 'x'], xyz[, 'y'], xyz[, 'z']),
    group = dd$Date.symptoms.year,
    group.mesh = mesh1d
  )

#-- Create index matrix --#
field.indices =
  inla.spde.make.index("field", n.spde = mesh$n, n.group = mesh1d$m)

est.cov <- as.list(as.data.frame(covs))
est.cov$Intercept = 1
covariate.names <- names(est.cov)
Rtb = Rt
Rtb[Rt > 0] = 1

#Rtb = log(Rt+1)

stack.est =
  inla.stack(
    data = list(response = Rtb),
    A = list(A.est, 1),
    effects =
      list(c(field.indices),
           c(est.cov)),
    tag = "est",
    remove.unused = FALSE,
    compress = FALSE
  )

# Define formula

formula <- as.formula(paste(
  paste("response ~ -1 + "),
  # paste("f(field, model=spde) + ",sep=""),
  paste(
    "f(field, model=spde,group=field.group, control.group=list(model='ar1')) + ",
    sep = ""
  ),
  paste(covariate.names, collapse = '+'),
  sep = ""
))
stack.est <- stack.est

# Call INLA and get results #
mod.pred =   inla(
  formula,
  data = inla.stack.data(stack.est),
  family = "binomial",
  control.predictor = list(
    A = inla.stack.A(stack.est),
    compute = TRUE,
    quantiles = NULL
  ),
  control.compute = list(cpo = TRUE, dic = TRUE, config =
                           TRUE),
  keep = FALSE,
  verbose = TRUE,
  #num.threads=5,
  control.inla = list(
    strategy = 'gaussian',
    int.strategy = 'eb',
    verbose = TRUE,
    fast = TRUE,
    dz = 1,
    step.factor = 1,
    stupid.search = FALSE
  )
)

index_lp = inla.stack.index(stack.est, "est")$data
lp = mod.pred$summary.linear.predictor$mean[index_lp]

##############
## Diagnostics
##############
#AUC
fun.auc(plogis(lp), Rtb)
fun.aucplot(plogis(lp), Rtb, "AUC")

#######################
##visualisation of maps
#######################

samp = inla.posterior.sample(100, mod.pred)

#by year
st_pred <- stack()
st_sd <- stack()
for (year in 2011:2016) {
  ind <- cn
  pred_val <- getValues(ind)#get values again
  w <- is.na(pred_val) #find NAs again
  index <- 1:length(w)
  index <- index[!w]
  pred_locs <-
    xyFromCell(ind, 1:ncell(ind))  #get prediction locations
  pred_locs <- pred_locs[!w,] #remove NA cells
  
  colnames(pred_locs) <- c('longitude', 'latitude')
  pred_locs <- ll.to.xyz(pred_locs) #get locations
  pred_mat <- matrix(nrow = nrow(pred_locs), ncol = length(names))
  Xp <- st[index]
  Xp = as.data.frame(Xp)
  Xp$intercept = 1
  A = inla.spde.make.A(
    mesh,
    loc = cbind(pred_locs[, 'x'], pred_locs[, 'y'], pred_locs[, 'z']),
    group = rep(year, nrow(pred_locs)),
    group.mesh = mesh1d
  )
  #A= inla.spde.make.A(mesh, loc=cbind(pred_locs[,'x'],pred_locs[,'y'],pred_locs[,'z']))
  
  lp <-
    as.matrix(Xp) %*% (mod.pred$summary.fixed$mean) + drop(A %*% mod.pred$summary.random$field$mean)
  
  P <- cn
  P[!w] = plogis(lp)
  print(year)
  st_pred <- addLayer(st_pred, P)
  
  #standard dev
  S <- stack()
  for (xx in 1:100) {
    coeff <- c()
    for (xxx in 1:length(covariate.names)) {
      coeff[xxx] <-
        samp[[xx]]$latent[grep(covariate.names[xxx], rownames(samp[[xx]]$latent)),]
    }
    field <-
      samp[[xx]]$latent[grep('field', rownames(samp[[xx]]$latent)),]
    
    Z <- cn
    Z[!w] = (drop(coeff %*% t(Xp)) + drop(A %*% field))
    S <- addLayer(S, Z)
  }
  lt <- function(x)
    exp(x) / (1 + exp(x))
  Qi <-
    calc(
      S,
      fun = function(x) {
        quantile(x, probs = c(.025, .975), na.rm = TRUE)
      }
    )
  
  M <- calc(S, mean)
  Sd <- calc(lt(S), sd)
  st_sd <- addLayer(st_sd, Sd)
  
  
}

# keeps scales standard
for (i in 1:nlayers(st_pred)) {
  st_pred[[i]][1:5] = 1
  st_pred[[i]][6:10] = 0
  
}

for (i in 1:nlayers(st_sd)) {
  st_sd[[i]][1:5] = 0.35
  st_sd[[i]][6:10] = 0
  
}
dev.off()
year <- c(2011:2016)
par(mfrow = c(2, 3))
for (i in 1:nlayers(st_pred)) {
  plot(st_pred[[i]], col = colfunc(1000),  axes = FALSE)
  writeRaster(st_pred[[i]], paste0(year[i], '.tif'))
}

par(mfrow = c(2, 3))
for (i in 1:nlayers(st_sd)) {
  plot(st_sd[[i]],  col = colfunc(1000), axes = FALSE)
  writeRaster(st_sd[[i]], paste0('sd', year[i], '.tif'))
}


par(mfrow = c(1, 2))

P[2] = 0
P[1] = 1
#colors

plot(P, col = colfunc(1000))