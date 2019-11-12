##############################################################################
# bootstrap percentile method
##############################################################################
bootstrap <- function(data, nboot, alpha){
  n=length(data)   # size of original sample
  xbar = mean(data)   # mean of original sample
  tmpdata = sample(data,n*nboot,replace = TRUE)   # re-sample from original data
  bootstrap_sample = matrix(tmpdata,nrow = n,ncol = nboot)  # arrange re-sampled data into (n x nboot) matrix
  bsmeans = colMeans(bootstrap_sample)
  deltastar = bsmeans - xbar
  sorteddeltastar = sort(deltastar)
  c_value = (1-alpha)*0.5
  d2=sorteddeltastar[as.integer(nboot*c_value)]
  d1=sorteddeltastar[as.integer(nboot*(1-c_value))]
  CI = xbar-c(d1,d2)
  return(CI)
}


library(latex2exp)

##############################################################################
# processing input file
##############################################################################
file_processor <- function(infile, nfrom, nto, nby, num_runs){
  input <- read.table(file = infile, header = FALSE)
  output <- data.frame(n= integer(0), id= integer(0), time= integer(0))
  nLines = (nto/nby)*num_runs
  for (n in seq(nfrom, nto, nby)){
    for (r in 1:nLines){
      if (input[r,1] == n){
        output <- rbind(output, input[r,])
      }
    }
  }
  return(output)  
}

##############################################################################
# compare two dataset (2)
##############################################################################
two_plotter <- function(file1, file2, 
                               nfrom, nto, nby, num_runs, 
                               legend1, legend2){
  
  data_set1 <- file_processor(file1, nfrom, nto, nby)
  num_n <- nrow(data_set1)/num_runs   # number of different values of n
  sample1 <- matrix(data_set1[,3], nrow = num_runs, ncol = num_n)
  y1 <-  colMeans(sample1)
  
  data_set2 <- file_processor(file2, nfrom, nto, nby)
  sample2 <- matrix(data_set2[,3], nrow = num_runs, ncol = num_n)
  y2 <-  colMeans(sample2)
  
  matplot(seq(nfrom, nto, nby), cbind (y1, y2), pch = c(2,17), 
           xlab = "Length of bitstring", ylab = "Empirical runtime", 
          col="black", cex=1.5)
  require(latex2exp)
  
  legend("topleft",c(TeX(legend1), TeX(legend2)),
         col="black", pch=c(2,17), 
         lwd=c(2), cex=1.3)
  
  #title("Testing UMDA on LeadingOnes with/without Elitism\
  #      param settings: lambda=n and mu=sqrt(n)")
}

##############################################################################
# compare two dataset (3)
##############################################################################
three_plotter <- function(file1, file2, file3,
                        nfrom, nto, nby, num_runs, 
                        legend1, legend2, legend3){
  
  data_set1 <- file_processor(file1, nfrom, nto, nby)
  num_n <- nrow(data_set1)/num_runs   # number of different values of n
  sample1 <- matrix(data_set1[,3], nrow = num_runs, ncol = num_n)
  y1 <-  colMeans(sample1)
  
  data_set2 <- file_processor(file2, nfrom, nto, nby)
  sample2 <- matrix(data_set2[,3], nrow = num_runs, ncol = num_n)
  y2 <-  colMeans(sample2)
  
  data_set3 <- file_processor(file3, nfrom, nto, nby)
  sample3 <- matrix(data_set3[,3], nrow = num_runs, ncol = num_n)
  y3 <-  colMeans(sample3)
  
  matplot(seq(nfrom, nto, nby), cbind (y1, y2, y3), pch = c(2,17,16), 
          xlab = "Length of bitstring", ylab = "Empirical runtime", 
          col="black", cex=1.5)
  
  legend("topleft",c(latex2exp(legend1),latex2exp(legend2), latex2exp(legend3)),
         col="black", pch=c(2,17,16), 
         lwd=c(2), cex=1.3)
  
  #title("Testing UMDA on LeadingOnes with/without Elitism\
   #     param settings: lambda=n and mu=sqrt(n)")
}


##############################################################################
# plot data with confident intervals
##############################################################################
plotter <- function(file_name, nfrom, nto, nby, alpha, nboot, num_runs){
  
  data <- file_processor(file_name, nfrom, nto, nby, num_runs)
  num_n <- nrow(data)/num_runs  # number of different n
  sample_data <- matrix(data[,3], nrow = num_runs, ncol = num_n)
  mean_times = colMeans(sample_data)
  
  # upper half of confidence interval
  uci <- seq(from=1, to=num_n, by=1)
  
  # lower half of confidence interval
  lci <- seq(from=1, to=num_n, by=1)
  
  for(i in 1:num_n){
    ci <- bootstrap(sample_data[,i], nboot, alpha)
    uci[i] <- ci[2]
    lci[i] <- ci[1]
  }
  
  x <- seq(from=nfrom, to=nto, by=nby) #x data
  y <- mean_times   #y data
  
  # plot the graph
  require(plotrix)
  plotCI(x, y, ui=uci, li=lci, 
         xlab = "Length of bitstring", ylab = "Empirical runtime",
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col="red")
  
  #title("Empirical Runtime with 95% confidence interval using 100 samples")
}

##############################################################################
# Regression
##############################################################################
regression <- function(file_name, nfrom, nto, nby, alpha, nboot, num_runs){
  
  data <- file_processor(file_name, nfrom, nto, nby, num_runs)
  num_n <- nrow(data)/num_runs  # number of different n
  sample_data <- matrix(
    #sapply(data[,3], function(x) log2(x)), 
    data[,3],
    nrow = num_runs, 
    ncol = num_n
  )
  mean_times = colMeans(sample_data)
  
  # upper half of confidence interval
  uci <- seq(from=1, to=num_n, by=1)
  
  # lower half of confidence interval
  lci <- seq(from=1, to=num_n, by=1)
  
  for(i in 1:num_n){
    ci <- bootstrap(sample_data[,i], nboot, alpha)
    uci[i] <- ci[2]
    lci[i] <- ci[1]
  }
  
  x <- seq(from=nfrom, to=nto, by=nby) #x data
  y <- mean_times   #y data
  
  # n^2
  m1 <- nls(y~a*x^2)
  r1 <- cor(y, predict(m1))
  summary(m1)
  cor(y, predict(m1))
  
  # n^{3/2}
  m2 <- nls(y~a*x^(3/2))
  r2 <- cor(y, predict(m2))
  summary(m2)
  cor(y, predict(m2))
  
  # n^2*log n
  m3 <- nls(y~a*(x^2)*log2(x))
  r3 <- cor(y, predict(m3))
  summary(m3)
  cor(y, predict(m3))
  
  # plot the graph
  require(plotrix)
  plotCI(x, y, ui=uci, li=lci, 
         xlab = "Length of bitstring", ylab = "Number of fitness evaluations",
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col="red")
  
  # add non-regression lines
  lines(x, predict(m1), lty="solid",lwd=2)
  lines(x, predict(m2), lty="dotted", lwd=3)
  lines(x, predict(m3), lty="dashed", lwd=3)
  
  
  # add legend
  legend('topleft',c(TeX("$4.2905\\cdot n^2$"), TeX("$122.83\\cdot n^{3/2}$"), TeX("$0.4407\\cdot n^2\\log n$")),
         lty=c("solid", "dotted","dashed"),
         lwd=c(2,3,3), cex=1.5, bty='n')
  
  #title("Empirical Runtime with 95% confidence interval using 100 samples")
}

##############################################################################
# plot probabilistic model
##############################################################################
model_plotter <- function(file_name, par_x, par_y, n ){
  data <- read.table(file = file_name, header = FALSE)
  gfrom <- 1 
  gby <- nrow(data)/(par_x*par_y)
  rows <- nrow(data)  #rows
  cols <- ncol(data)
  
  par(mfrow=c(par_y, par_x)) # all plots on one page
  
  # add lines to the graph
  for(i in seq.int(from=gfrom, to=rows, by=gby)){
    heading = paste("Generation ",i) 
    plot(seq(from=1, to=n, by=1), type='n', ylim=c(0,1), main=heading, xlab = 'positions', ylab = 'probability')
    lines(seq(from=1, to=n, by=1), data[i, 3:cols], type='l', col='red')
  }
}

##############################################################################
# Plot the variance of the bit-string
##############################################################################
var_plotter <- function(file_name_with_model){
  data <- read.table(file = file_name_with_model, header = FALSE)
  var_matrix <- data.frame(variance= double(0))
  cols = ncol(data)  # number of columns
  n <- cols-2
  rows = nrow(data)  # number of rows
  for (r in 1:rows){
    var <- 0
    for (c in 3:cols){
      var <- var + (1-data[r,c])*data[r,c]
    }
    var_matrix <- rbind(var_matrix, var)
  }
  x <- 1:rows
  y <- var_matrix[,1]
  max_var <- max(var_matrix[,1])
  min_var <- round(min(var_matrix[,1]), digits = 4)
  plot(x, type='n', ylim=c(0,n/4), xlab = 'generations', ylab = 'variance')
  lines(x, y, type='l', col='red')
  legend('topright',c(paste("max-value: ", max_var), paste("min-value: ", min_var)),
                 lwd=c(3), cex=1.5, bty='n')
  
  #title(sprintf("Variance of the bitstring over generations for n=%d \
  #        parameter settings: lambda=n and mu=sqrt(n)",n))
}

##############################################################################
# var_k plotter
##############################################################################
var_k_plotter <- function(var_k_file_name, n){
  var_k <- read.table(file = var_k_file_name, header = FALSE)
  x <- 1:nrow(var_k)
  y <- var_k[,1]
  
  max_var <- max(var_k[,1])
  min_var <- round(min(var_k[,1]), digits = 4)
  
  plot(x, type = "n",ylim = c(0,n/4), xlab = "generations", ylab = "variance of k-region")
  lines(x, y, type = "l", col="red")
  legend('topright',c(paste("max-value: ", max_var), paste("min-value: ", min_var)),
         lwd=c(3), cex=1.5, bty='n')
  
  #title(sprintf("Variance of k-region over generations for n=%d \
  #        parameter settings: lambda=n and mu=sqrt(n)",n))
} 

##############################################################################
# plotting all marginal probabilities over generations 
##############################################################################
model_generation_plotter <- function(file_name, n){
  data <- read.table(file = file_name, header = FALSE)
  x <- seq(from=1, to=nrow(data), by=1)
  plot(x, ylim = c(0,1), type = "n", xlab = "generations", ylab = "marginal probability")
  for (i in 2:ncol(data)){
    lines(x, data[,i], type = "l")
  }
  #title(sprintf("Marginal probabilities over generations for n=%d \
  #        parameter settings: lambda=n and mu=sqrt(n)log(n)",n))
}


##############################################################################
# print k-j relationship
##############################################################################

kj_plotter <- function(file_name, n){
  data <- read.table(file = file_name, header = FALSE)
  x <- data[,1]
  y <- data[,2]
  plot(0:n, ylim = c(0,n+1), type = "n", xlab = "width of k-region (k)", ylab = "current level (j)")
  lines(x, y, type = "l", col="red", lwd=3)
}


##############################################################################
# Working Space
##############################################################################

#kj_plotter("kj_small_mu.log", 10000)
#kj_plotter("kj_large_mu.log", 10000)


#multilines_plotter("om_non_elitism.log", "om_elitism.log",
#                   100, 2000, 100,
#                   100,
#                   "Non-Elitism", "Elitism (keep the best)")

#plotter("om_non_elitism.log", 100, 2000, 100, 0.95, 100, 100)
#regression("om_non_elitism.log", 100, 2000, 100, 0.95, 100, 100)

#model_plotter("lo_non_elitism_model.log", 2, 5, 100)
#var_plotter("om_non_elitism_model_2.log")
#var_k_plotter("var_k.log", 10000)
#model_generation_plotter("om_large_mu_n_100.log", 100)
regression("pbil_los_large_2_1000", 50, 400, 50, 0.95, 100, 100)

title("Empirical runtime of PBIL on OneMax \
with 95% confidence interval using 100 samples \
          parameter settings: lambda=n, mu=sqrt(n), eta=0.1")

three_plotter("jump_small_r_large_mu.log", 
              "jump_medium_r_large_mu.log",
              "jump_large_r_large_mu.log",
              100, 2000, 100, 100, "r=n/4", "r=n/2", "r=3n/4")

two_plotter("jump_medium_r_small_mu.log", "jump_medium_r_large_mu.log",
                        100, 2000, 100, 100, 
            '$\\mu = \\sqrt{n}\\log(n)$', 
            '$\\mu = \\sqrt{n}$')

title(TeX("Empirical runtime of UMDA on JUMP with $r = n/2$"))



model_generation_plotter("pbil_om_100_model", 100)

title("PBIL on OneMax \
      parameters: n=100, mu=sqrt(n), lambda=10*sqrt(n) and eta=0.1")


#######################################################################################
setEPS()
postscript(file="exp2.eps", width = 14.0, height = 7.0)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

setEPS()
postscript(file="exp1.eps", width = 14.0, height = 7.0)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
file_name = "pbil_los_small_2_1000"
nfrom = 50
nto = 1000
nby = 50
num_runs <- 100
nboot = 100
alpha = .95

data <- file_processor(file_name, nfrom, nto, nby, num_runs)
num_n <- nrow(data)/num_runs  # number of different n
sample_data <- matrix(
  #sapply(data[,3], function(x) log2(x)), 
  data[,3],
  nrow = num_runs, 
  ncol = num_n
)
mean_times = colMeans(sample_data)

# upper half of confidence interval
uci <- seq(from=1, to=num_n, by=1)

# lower half of confidence interval
lci <- seq(from=1, to=num_n, by=1)

for(i in 1:num_n){
  ci <- bootstrap(sample_data[,i], nboot, alpha)
  uci[i] <- ci[2]
  lci[i] <- ci[1]
}

x <- seq(from=nfrom, to=nto, by=nby) #x data
y <- mean_times   #y data

# n^2
#m1 <- nls(y~a*x^2)
#r1 <- cor(y, predict(m1))
#summary(m1)
#cor(y, predict(m1))

# n^{3/2}
#m2 <- nls(y~a*x^(3/2))
#r2 <- cor(y, predict(m2))
#summary(m2)
#cor(y, predict(m2))

# n^2*log n
#m3 <- nls(y~a*(x^2)*log2(x))
#r3 <- cor(y, predict(m3))
#summary(m3)
#cor(y, predict(m3))

# plot the graph
require(plotrix)
plotCI(x, y, ui=uci, li=lci, 
       xlab = "Problem size", ylab = "Actual Number of Fitness Evaluations",
       cex.lab=2.0, cex.axis=2.0, cex.main=1.5, cex.sub=1.5, col="red")

# add non-regression lines
#lines(x, predict(m1), lty="solid",lwd=4)
#lines(x, predict(m2), lty="dotted", lwd=5)
#lines(x, predict(m3), lty="dashed", lwd=4)

lines(x, 4.2905*x^2, lty="solid",lwd=4)
lines(x, 122.83*x^{3/2}, lty="dotted", lwd=5)
lines(x, 0.4407*x^2*log2(x), lty="dashed", lwd=4)

grid(nx=NULL, ny = NULL, col = "gray48", lty = "dotted",
     lwd = 2, equilogs = TRUE)

# add legend
legend('topleft',c(TeX("$4.2905\\cdot n^2$"), TeX("$122.83\\cdot n^{1.5}$"), TeX("$0.4407\\cdot n^2\\log n$")),
       lty=c("solid", "dotted","dashed"),
       lwd=c(4,4,4),
       bty='n',
       cex = 2.0,
       box.lwd = par("lwd"), box.lty = par("lty"), box.col = par("fg"))
dev.off()


dev.off()
###########################################################################################

file_name = "pbil_los_small_2_1000"
nfrom = 50
nto = 1000
nby = 50
num_runs <- 100

data <- file_processor(file_name, nfrom, nto, nby, num_runs)
num_n <- nrow(data)/num_runs  # number of different n
sample_data <- matrix(
  #sapply(data[,3], function(x) log2(x)), 
  data[,3],
  nrow = num_runs, 
  ncol = num_n
)
mean_times = colMeans(sample_data)

x <- seq(50,1000,50)
y <- mean_times

high_y <- apply(sample_data,2,max)
low_y <- apply(sample_data,2,min)

plot(x,y, type="n",ylim = c(0,1.6*max(high_y)), xlab = "Problem Size", ylab = "Average Runtime")

lines(x, high_y,col = 'pink')
lines(x, low_y, col = 'pink')

polygon(c(x, rev(x)), c(high_y, rev(low_y)),
        col = "pink", border = NA)

lines(x,y, type = 'l', col = "red")

# models
m1 <- function(n, a){a*n*n*log2(n)}
#for (a in seq(0.4,1.0,0.05)){
  lines(x, m1(x,.7), col='blue')
#}
