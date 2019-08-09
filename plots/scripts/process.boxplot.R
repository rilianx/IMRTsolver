# data.matrix: groups of data should be on columns, that is each column will be a boxplot
# data.labels: names of the data to be printed
# plot.mar: exterior margin for the boxplot c(bottom, left, top, right)
# plot.sides.mar: margins for  (axis title, axis labels, axis line) default c(3, 1, 0)
# plot.title: title of the plot
# plot.yaxis.label: Label of y axis
# file.name: filename to save the plot WITHOUT EXTENSION
# file.type: type of time (png, eps, pdf)
# pdf.dim: dimensions of the pdf c(width, height)
# plot.axis.style: option las of par. Axis text sty 0: always parallel to the axis [_default_], 1: always horizontal, 2: always perpendicular to the axis, 3: always vertical.
# outliers: Remove or not outliers
# colors: Colors for the boxplot
# add.points: Print all points with jitter
# add.text: add text in the location specified by at (check the code)
# axis.style  0: always parallel to the axis, 1: always horizontal, 2: always perpendicular to the axis, 3: always vertical.
do.boxplot <- function (data.matrix, data.labels=c(), plot.mar=c(6,14,4,4), 
                        plot.sides.mar=c(3, 1, 0), plot.title="", plot.yaxis.label="", 
                        file.name="output", file.type="pdf", plot.axis.style=1, 
                        plot.size.axis=4, plot.size.title=3, plot.text.yaxis=plot.size.axis, 
                        vertical.x.axis=1, outliers=TRUE, colors=c(), add.points=FALSE, 
                        add.text=NULL, line.add.text=9, at.add.text=1.1, pdf.dim=c(20,8), 
                        small.x.axis=FALSE){
  
  #starting output
  if(file.type == "eps")
     postscript(file = paste(file.name,"-bxp.eps", sep=""), width=500)
  if(file.type == "png")
     #png(file = paste(file.name,"-bxp.png", sep=""), width = 560, height = 320, units = "px", pointsize = 8, bg = "white")
     png(file = paste(file.name,"-bxp.png", sep=""), width = 1060, height = 700, units = "px", pointsize = 16, bg = "white")
  #long 
  #png(file = paste(file.name,"-bxp.png", sep=""), width = 800, height = 320, units = "px", pointsize = 8, bg = "white")

  if(file.type == "pdf")
    cairo_pdf(file = paste(file.name,"-bxp.pdf", sep=""), width=pdf.dim[1], height=pdf.dim[2])

  #getting min and max
  plot.min <- min(apply(data.matrix, 2, min))
  plot.max <- max(apply(data.matrix, 2, max))
  
  #plot characteristics
  par(las=plot.axis.style, mar=plot.mar, cex.axis=plot.size.axis+1, cex.main=plot.size.title, lwd=5, mgp=plot.sides.mar)
  
  #plot
  boxplot(data.matrix, main=plot.title, xaxt="n", outline=outliers, col=colors, lwd=5)

  if(add.points){
    for(i in 1:ncol(data.matrix)){
       mj <- jitter(rep(i,nrow(data.matrix)), factor=10/i)
       points(mj, data.matrix[,i], pch=20, col=rgb(0,0,0,.2) , cex=5)
    }
  }
  
  dist.y<- max(nchar(as.character(plot.min)), nchar(as.character(plot.max)))-1

  #y axis
  if(length(data.labels)==ncol(data.matrix)){
    #For two line label:
    #mtext(plot.yaxis.label, side=2, line=(plot.mar[2]-8), cex=plot.text.yaxis+0.5, las=0)
    #For onle line label:
    #mtext(plot.yaxis.label, side=2, line=(plot.mar[2]-3.5), cex=plot.text.yaxis+0.5, las=0)
    mtext(plot.yaxis.label, side=2, line=(plot.mar[2]-4), cex=plot.text.yaxis+1, las=0)
  }

  ## HERE add.text option
  if(!is.null(add.text)){
    #side=1 for down
    #side=3 for up
    #at=1 for left
    #at=2 for right
    #mtext(add.text,side=3,line=-4, cex=plot.size.axis+0.5, at=2)
    #mtext(add.text,side=3,line=-3, cex=3, at=1)
    mtext(add.text,side=3,line=-line.add.text, cex=2, at=at.add.text)
  
  }

  rr<-3
  if(vertical.x.axis == 2) rr<- 4
  #x axis
  
  if(small.x.axis)
    axis(1, at=c(1:length(data.labels)), labels=data.labels, 
         line=plot.mar[1] - rr, tick=FALSE, 
         las=vertical.x.axis, cex.axis=plot.size.axis)
  else
    axis(1, at=c(1:length(data.labels)), labels=data.labels, 
         line=plot.mar[1] - rr+0.5, tick=FALSE, 
         las=vertical.x.axis, cex.axis=plot.size.axis+1)
  
  dev.off()
}

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

if (length(args)< 2){
  stop("You must provide files")
}


nrep <-30
instance.names <- read.table("instances.txt", 
                             stringsAsFactors=FALSE)[,1]
ninstances     <- length(instance.names)

# First argument is the test name 
outname <- args[1]
# Second argument is the folder in which saving
folder <- args[2]
# Third argument the number of tests
ntests <- as.numeric(args[3])
cat ("Reading ", ntests, " tests\n")

# Next the name of the tests
test.names <- c()
a = 3
for (i in 1:ntests){
  a = a + 1
  test.names <- c(test.names, args[a])
}

# Finally the get and read files
files <- c()
data  <- list() #data directly from convergence file
for (i in 1:ntests) {
  a     <- a + 1
  file  <- args[a]
  files <- c(files, file)
  cat ("Reading ", test.names[i], ": ", file, " \n")
  data[[i]] <- read.table(file = file, 
                          header = TRUE, 
                          sep = ",", 
                          comment.char = "$", 
                          stringsAsFactors = FALSE)
}


i.data  <- list() #data by instance
mi.data <- matrix(NA, nrow = length(instance.names), 
                  ncol = ntests) # data mean per instance and test
#results in long format
results <- as.data.frame(matrix(NA, nrow = 0, ncol = 4,
                         dimnames = list(c(), c("instance", "test", "sig", "F"))))

# Instance loop
for (t in 1:ninstances) {
  i.data[[t]] <- matrix(NA, nrow=nrep, ncol=ntests)
  last        <- t*nrep
  initial     <- (last-nrep+1)
  # Put data in per instance format
  for (i in 1:ntests) {
    i.data[[t]][1:nrep,i] <- as.numeric(data[[i]][initial:last,"F"])
  }
  # Mean per instance
  mi.data[t,] <- colMeans(i.data[[t]])

  # Making statistical tests
  st   <- "" 
  pval <- c()
  for (i in 1:(ntests-1)) {
    for (j in (i+1):ntests) {
      if (sum(i.data[[t]][,i] - i.data[[t]][,j]) == 0) {
        pval <- c(pval, 1)
        next
      }
      res <- wilcox.test(x = i.data[[t]][,i], 
                         y = i.data[[t]][,j],  
                         paired = TRUE)
      pval <- c(pval, res$p.value)
      st   <- paste(st, test.names[i], " vs ", 
                    test.names[j], "= ", 
                    res$p.value, "\n", sep="")
    }
  }
  pval <- p.adjust(pval, method="bonferroni")

  ## Making large plot results
  cc <- colMeans(i.data[[t]]) 
  wb <- which.min(cc)

  # Check if there is a test that is better!
  k  <-0
  all.b <- TRUE
  for (i in 1:(ntests-1)) {
    for (j in (i+1):ntests) {
      k <- k+1
      if (i!=wb && j!=wb) next
      if (pval[k] >= 0.05) all.b <- FALSE
    }
  }
   
  # Put results in long format
  for (i in 1:ntests) {
    sig <- "-"
    if (wb == i && all.b) sig <- "+"
    for(j in 1:nrep) {
      k = nrow(results) + 1
      results[k,] <- c(instance.names[t], test.names[i], sig, i.data[[t]][j,i]) 
    }
  }
  
  #Make individial plots
  do.boxplot (i.data[[t]], data.labels=test.names, 
              plot.title=paste("ILS", outname, "-", instance.names[t]),
              plot.yaxis.label="F", 
              file.name=paste(folder,"/",outname,"-",instance.names[t],sep=""), 
              file.type="pdf", outliers=TRUE, add.points=TRUE, add.text=st)
}

## Make overall plots
pval <- c()
st   <- ""
# Make statistical tests
for (i in 1:(ntests-1)) {
  for (j in (i+1):ntests) {
    if (sum(mi.data[,i] - mi.data[,j]) == 0) {
      pval <- c(pval, 1)
      next
    }
    res <- wilcox.test(x = mi.data[,i], 
                       y = mi.data[,j], 
                       paired=TRUE)
    pval <- c(pval, res$p.value)
  }
}
pval <- p.adjust(pval, method="bonferroni")
# Print the tests
k <- 1
for (i in 1:(ntests-1)) {
  for (j in (i+1):ntests) {
    st <- paste(st, test.names[i], " vs ", test.names[j], "= ", pval[k], "\n", sep="")
    k <- k+1
  }
}

do.boxplot (mi.data, data.labels=test.names, 
            plot.title=paste("ILS", outname), 
            plot.yaxis.label="mean F", 
            file.name=paste(folder,"/",outname,"-summary",sep=""), 
            file.type="pdf", outliers=TRUE, add.points=TRUE, 
            add.text=st)

## Make all instance plot
results$F <- as.numeric(results$F)

#save(file="file.Rdata", results)

#p <- ggplot(data = results, aes(x=instance, y=F)) 
#p <- p + geom_boxplot(aes(fill=test)) + coord_flip() 
#p <- p + facet_wrap( ~ instance, scales="free")

p <- ggplot(data = results, aes(x=test, y=F)) 
p <- p + geom_boxplot(aes(fill=sig)) 
p <- p + facet_wrap( ~ instance, scales="free", ncol = 5) 
p <- p + ggtitle(paste("ILS", outname)) 
p <- p + guides(fill=guide_legend(title="Wilcoxon test")) 
p <- p + xlab("Neighborhood") 
p <- p + theme(plot.title = element_text(hjust = 0.5))

ggsave(paste(folder,"/",outname, "-all-bxp.pdf", sep=""), 
       plot = p, width=10, height=5)

