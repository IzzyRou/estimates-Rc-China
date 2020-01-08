#################################################################################
## Code generating Figure 5  and expected serial interval quoted in manuscript #
################################################################################

####################
# Define functions #
####################

#truncated normal
rtnorm <- function(n, mean, sd, a = 0, b = 1){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk
t_col <- function(color, percent = 10, name = NULL) {
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  ## Save the color
  invisible(t.col)
  
}
## END
l<-function(t,As){
  t[t<0]<-0
  As*t*exp(-0.5*As*t*t)}

line<-function(t,shift=15, A){lines(t,l(t-shift,A),type ="l", lwd=0.5, col= mycol)}

mycol<-t_col("red",percent =50)

# illustration of SI prior with parameters used in main text- this is from a truncated normal distribution, normal uniform bounded distribution version commented out

#incubation period
shift<-15
#Alpha parameter
As<-0.003

t<- c(1:200)

plot(t,l(t-shift,As),type ="l", lwd=0.5, col= "red", xlab = "Time", ylab="Probability density", ylim=c(0,0.12), cex.lab=1.5, cex.axis =1.5 )

#we set a prior for Alpha with mean 0.003 and sd 0.01,using a truncated normal distribution
#As<-runif(300,0.00001,0.02)
As<-rtnorm(1000,0.003,0.01)
#As<-rtnorm(1000,0.002,0.1)


quant<-quantile(As,c(0.025,0.975,0.5))
As_2.5<-quant[1]
As_97.5<-quant[2]
As_0.5<-quant[3]

#2.5% quantile 
t[which.max(l(t-shift,As_2.5))]

#97.5% quantile
t[which.max(l(t-shift,As_97.5))]

#median
t[which.max(l(t-shift,As_0.5))]

for(i in 1:length(As)){
  line(t=t,A=As[i])
}
As<-0.003
lines(t,l(t-shift,As), type ="l", lwd=5, col= "black")
lines(t,l(t-shift,As_2.5), type ="l", lwd=5,lty="dashed", col= "black")
lines(t,l(t-shift,As_97.5), type ="l", lwd=5, lty= "dashed", col= "black")
title("Serial Interval Distribution", cex.main=2)