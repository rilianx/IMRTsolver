
library(stringr)

sets <- c("0_70_140_210_280", "10_80_150_220_290", "15_85_155_225_295", "20_90_160_230_300", "25_95_165_235_305", "30_100_170_240_310", "35_105_175_245_315", "40_110_180_250_320", "45_115_185_255_325", "50_120_190_260_330", "55_125_195_265_335", "5_75_145_215_285", "60_130_200_270_340", "65_135_205_275_345")

prefix <- "~/instances/dao/Equidistantes/"
datfiles <- c("DDM_BLADDER.dat", "DDM_PTVHD.dat", "DDM_RECTUM.dat")

location <- "~/instances/dao/instances/"

allinstances <- c()

for (s in sets) {
  datlines <- paste(prefix,s, "_",datfiles,sep="")
  aline <- gsub("_", " ", s)
  angles <- unlist(strsplit(s, "_"))
  coordfiles <- paste(prefix,"Coordinates/CoordinatesBeam_",angles,".txt", sep="")

  for (a in 1:length(angles)) 
    coordfiles[a] <- paste(angles[a],";",coordfiles[a],sep="")

  filename <- paste(location, s, "_data.txt", sep="")
  fileConn<-file(filename)
  writeLines(c(aline,datlines), fileConn)
  close(fileConn)

  aux <- paste("--file-dep ", filename, sep="")

  filename <- paste(location, s, "_coordinates.txt", sep="")
  fileConn<-file(filename)
  writeLines(coordfiles, fileConn)
  close(fileConn)

  aux <- paste(aux, " --file-coord ", filename, sep="")

  allinstances <- c(allinstances,aux)
}

filename <- "instances.txt"
fileConn<-file(filename)
writeLines(allinstances, fileConn)
close(fileConn)

