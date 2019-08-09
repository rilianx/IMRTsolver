library(ggplot2)
options(width=120)


getData <- function (data, evals, nrep) {
  res <- as.data.frame(matrix(NA, nrow=length(evals), 
                              ncol=nrep))
  for (k in 1:nrep) {
    for (v in 1:length(evals)) {
      aa <- (data[[k]][,1] <= evals[v])
      s <- 1
      if (sum(aa) > 0)
        s <- max(which(data[[k]][,1] <= evals[v]))     
      res[v,k] <- data[[k]][s,"F"] 
    }
  }
  return(res)
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)< 2){
  stop("You must provide files")
}


# First argument is the test name 
outname     <- args[1]
file.prefix <- args[2]
save.folder <- args[3]
ntests     <- as.numeric(args[4])
cat ("Processing ", ntests, " tests\n")

# Next the name of the tests
test.names <- c()
a <- 4
for (i in 1:ntests){
  a <- a + 1
  test.names <- c(test.names, args[a])
}

# Finally the files sufix
file.sufix <- c()
for (i in (a+1):length(args)) {
  file.sufix <- c(file.sufix, args[i])
}

nrep <-30
instance.names <- read.table("instances.txt", stringsAsFactors=FALSE)[,1]
ninstances     <- length(instance.names)

data  <- list() # data by instance/tests
mean.data <- list()
data.files <- list()

# instances loop
for (i in 1:ninstances) {
  data[[i]]  <- list()
  mean.data[[i]] <- as.data.frame(matrix(NA, ncol=5, nrow=0))
  # Get data of the different tests
  for (j in 1:ntests){
    aux.data <- list() # data from the convergence file
    evals    <- c() # all number of evaluations recorded
    # Get data from different repetitions
    for (k in 1:nrep) {
      # get filename
      file <- paste("../output/",instance.names[i],
                    "_", file.prefix, file.sufix[j],
                     "_", k, ".conv.traj", sep="")
      print(file)
      aux.data[[k]] <- read.table(file, sep=";", colClasses="numeric",
                                  skip=1, stringsAsFactors=FALSE)
      colnames(aux.data[[k]]) <- c("eval", "time", "F")
      evals <- unique(c(evals, aux.data[[k]][,"eval"]))
    }

    evals <- sort(evals)
    data[[i]][[j]] <- getData (aux.data, evals, nrep)

    mm <- as.numeric(rowMeans(data[[i]][[j]]))
    #ss <- as.numeric(apply(data[[i]][[j]], 1, sd))
    #m1 <- mm-ss
    #m2 <- mm+ss
    m1 <- as.numeric(apply(data[[i]][[j]], 1, min))
    m2 <- as.numeric(apply(data[[i]][[j]], 1, max))
    tt <- rep(test.names[j], length(evals))

    aux <- as.data.frame(cbind(evals,tt,mm, m1, m2), 
                         stringsAsFactors=FALSE)
    mean.data[[i]] <- rbind(mean.data[[i]], aux)
  }
   

  colnames(mean.data[[i]]) <- c("eval","test", "F", "min", "max")
  mean.data[[i]]$eval <- as.numeric(mean.data[[i]]$eval)
  mean.data[[i]]$F <- as.numeric(mean.data[[i]]$F)
  mean.data[[i]]$min <- as.numeric(mean.data[[i]]$min)
  mean.data[[i]]$max <- as.numeric(mean.data[[i]]$max)

  #for (j in 1:length(test.names)) {
  #  dd <- mean.data[[i]][mean.data[[i]]$test == test.names[j],]
  #  p <- ggplot(data=dd, aes(x=eval, y=mean))
  #  p <- p + geom_line()
  #  p <- p + geom_ribbon(aes(ymin=min, ymax=max), color="grey70",alpha=0.4)
  #  p <- p + ggtitle(paste("ILS", outname, " - ", test.names[j]))
  #  p <- p + theme(plot.title = element_text(hjust = 0.5))
  #  ggsave(paste(save.folder,"/",outname, "-", test.names[j], "-", instance.names[i], "-conv.pdf", sep=""), plot = p, width=10, height=5)    
  #}

  p <- ggplot(data=mean.data[[i]], aes(x=eval, y=F, color=test))
  p <- p + geom_line(aes(x=eval, y=F, color=test))
  p <- p + geom_ribbon(aes(ymin=min, ymax=max, fill=test), color="grey70",alpha=0.4)
  p <- p + coord_cartesian(xlim = c(0,10000), ylim = c(30,200))
  p <- p + ggtitle(paste("ILS", outname, " - ", instance.names[j])) 

  ggsave(paste(save.folder,"/",outname, "-", instance.names[i], 
               "-all-conv.pdf", sep=""), 
         plot = p, width=10, height=5)
}



