pre_processor <- function(file_name, eta_from, eta_to, eta_by){
  data <- read.table(file = file_name, header = FALSE)
  output <- data.frame(eta=numeric(0), id= integer(0), time= integer(0))
  tol = 1e-5
  for (eta in seq(eta_from, eta_to, eta_by)){
    for (row in 1:nrow(data)){
      if (abs(data[row,1]- eta) <= tol){
        output <- rbind(output, data[row,])
      }
    }
  }
  return(output)  
}

file_name = "eta_small_100"
eta_from = 0.05
eta_to = 1.00
eta_by = 0.05

data <- read.table(file = file_name, header = FALSE)
nrow(data)

processed_data = pre_processor(file_name,eta_from,eta_to,eta_by)

row <- 1
eta <- 0.05

#lambda <- function(x){2*exp(1.0)*10/x^2}

d  <- data.frame("0.05" = sapply(processed_data[seq(row,row+99),3], function(x) log2(x)))
#d  <- data.frame("0.05" = sapply(processed_data[seq(row,row+99),3], function(x) log2(x)))
#d["0.1] <- sapply(processed_data[seq(101,200,1),3], function(x) log2(x))
#d[toString(0.15)] <- sapply(processed_data[seq(201,300,1),3], function(x) log2(x))

#row <- row + 100

for (eta in seq(0.1, eta_to, eta_by)){
  row <- row + 100
  print(eta)
  d[toString(eta)] <- sapply(processed_data[seq(row,row+99),3], function(x) log2(x))
}


mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
boxplot(d, 
        names=c(seq(eta_from,eta_to,eta_by)),
        ylab = "Log Number of Fitness Evaluations",
        xlab = "Smoothing Parameter",
        col = "red",
        cex.lab=2.0, cex.axis=2.0)
grid(nx=NULL, ny = NULL, col = "gray48", lty = "dotted",
     lwd = 2, equilogs = TRUE)
#title("Experiement with smoothing parameter\\
#      settings: n=100, mu = sqrt(n), lambda = (2e/eta^2)mu")

setEPS()
postscript(file="exp2-runtime.eps", width = 14.0, height = 7.0)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))

dev.off()


pdf(file = "eta_2_1000.pdf")

dev.off()




####################################################################
#file_name = "eta_small_100"
file_name = "pbil_1_los_fix_1000_eta_2"
eta_from = 0.05
eta_to = 1.00
eta_by = 0.05

data <- read.table(file = file_name, header = FALSE)
nrow(data)

processed_data = pre_processor(file_name,eta_from,eta_to,eta_by)

row <- 1
eta <- 0.05

lmbda <- 10*1000^0.5

d  <- data.frame("0.05" = sapply(processed_data[seq(row,row+99),3], function(x) x/lmbda))
#d  <- data.frame("0.05" = sapply(processed_data[seq(row,row+99),3], function(x) log2(x)))
#d["0.1] <- sapply(processed_data[seq(101,200,1),3], function(x) log2(x))
#d[toString(0.15)] <- sapply(processed_data[seq(201,300,1),3], function(x) log2(x))

#row <- row + 100

for (eta in seq(0.1, eta_to, eta_by)){
  row <- row + 100
  print(eta)
  #d[toString(eta)] <- processed_data[seq(row,row+99),3]
  d[toString(eta)] <- sapply(processed_data[seq(row,row+99),3], function(x) x/lmbda)
}


mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
boxplot(d,
        names=c(seq(eta_from,eta_to,eta_by)),
        ylab = "Actual Number of Generations",
        xlab = "Smoothing Parameter",
        col = "red",
        cex.lab=2.0, cex.axis=2.0)
grid(nx=NULL, ny = NULL, col = "gray48", lty = "dotted",
     lwd = 2, equilogs = TRUE)

title("Parameters: Fixed lambda/mu = 10, mu=sqrt(n) and n=1000")

