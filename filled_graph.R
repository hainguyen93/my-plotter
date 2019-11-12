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

#############################################################################################
#parameters for runtime
file_name = "pbil_los_small_2_1000"
nfrom = 50
nto = 1000
nby = 50
num_runs <- 100
y_scale <- 1.0


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

mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

plot(x, y, type="n",ylim = c(0,y_scale*max(high_y)), 
     xlab = "Problem Size", ylab = "Actual Running Time",
     cex.lab = 1.5, cex.axis = 1.5)

grid(nx=NULL, ny = NULL, col = "gray", lty = "dotted",
     lwd = 2, equilogs = TRUE)

#lines(x, high_y,col = 'pink')
#lines(x, low_y, col = 'pink')

polygon(c(x, rev(x)), c(high_y, rev(low_y)),
        col = "grey", border = NA)

#average runtime
lines(x,y, type = 'l', lwd = 2, col = "red")

# model n^2
m1 <- function(n, a){a*n^2}
lines(x, m1(x,4.2905), col='blue', lwd = 2, pch=6)

# model n^1.5
m2 <- function(n,a){a*n^1.5}
lines(x, m2(x,122.83), col='orange', lwd=2)

#model n^2log n
m3 <- function(n,a){a*n^2*log2(n)}
lines(x, m3(x,0.4407), col = 'yellow', lwd=2)

axis(side=2, at=seq(0,1000,100))

# add legend
legend('topleft',c(TeX("$6\\~n^2$"), TeX("PBIL Act. Runtime")),
       lty=c("solid"),
       col = c("blue","red"),
       lwd=c(2),
       cex = 1.5,
       box.lwd = par("lwd"), box.lty = par("lty"), box.col = par("fg"))



setEPS()
postscript(file="pbil_los_2_1000.eps", width = 7.0, height = 8.0)
#add plot here
dev.off()

