##########################
#   XCMS-processing      #
##########################
#
#by Francisco Jose de Novais
#This script to perform processing data (.netCDF, .mzMXL or mzData) to MetaboAnalyst.
#
#If you need the package:
# source("http://bioconductor.org/biocLite.R")
# library("BiocInstaller")
# biocLite("mzR")
#biocLite("xcms")
#biocLite("faahKO")
#biocLite("multtest")
#biocLite ("CAMERA")
#
#
library("xcms")
#
#
mzfiles <- list.files("/PASTEPATH/" ,recursive = TRUE, full=T) # input files mzXML or mXML (step 1)
mzfiles #check if all files were included
xset <- xcmsSet(mzfiles,polarity = "Positive//Negative",method="centWave//profile", peakwidth=c(5,20),ppm=15, snthresh=5,lockMassFreq = 554.2615)  # peak picking/detection/ UPLC (GUO et al, 2015 for details)(step 2)
  xset 
xsg <- group(xset, bw=3)    # peak alignment with plots (add sleep class = number) (step 3.1)  
  xsg
xsg1 <- retcor(xsg, method="loess",span=0.67,family = "gaussian", plottype="deviation") # retention time correction (step 3.2) also plot the result
  xsg1
xsg2 <- group(xsg1, bw=3)   # re-align (step 3.3)
  xsg2
xsg3 <- fillPeaks(xsg2)  # filling in missing peak data (step 4) ##Need a good memory
  xsg3
#
#reporttab <- diffreport (xsg2, "HRFI", "LRFI", "STATSNEG", 10, metlin=0.15)# for statistical tab ##Need a good memory
#  reporttab
#write.csv(file="../../../../../../reportttab.csv", reporttab) # save statistical with compounds
#check<-peakTable(xsg3, filebase="peakList")
#
dat <- groupval(xsg2, "medret", "into")  # get peak intensity matrix (step 5)
  dat <- rbind(group = as.character(phenoData(xsg2)$class), dat)  # add group label
write.csv (file="../NEGATIVE/teste.csv",dat)  # save the data to CSV file
#
dat1 <- groupnames(xsg3, mzdec = 4, rtdec = 4)#head para o arquivo final
write.csv (file="../../../../../../head.csv",dat1)
#
#############################
#     CAMERA-ANNOTATION     #
#############################
#
library("CAMERA")#Annotation
library(faahKO)
#Create an xsAnnotate object
xsa <- xsAnnotate(xsg2,polarity = "negative//positive")
#Group after RT value of the xcms grouped peak
xsaF <- groupFWHM(xsa, perfwhm=0.6)
#Annotate isotopes
xsaFI <- findIsotopes(xsaF, mzabs = 0.01)
#Verify grouping
xsaC <- groupCorr(xsaFI,cor_eic_th = 0.75 )
#Annotate adducts
xsaFA <- findAdducts(xsaFI, polarity="negative//positive")
#Get final peaktable and store on harddrive
peaklist <- getPeaklist(xsaFA)
write.csv(peaklist,file="../mzXML-Centwave/result_CAMERA.csv")
plotEICs(xsaFA, pspec=2, maxlabel=5)
#
#
#GUO, L. et al. Three plasma metabolite signatures for diagnosing high altitude pulmonary edema.Sci Rep. 2015 Oct 13;5:15126. doi: 10.1038/srep15126.
#
# CREATING PLOTS
#TIC
getTICs <- function(xcmsSet=NULL,files=NULL, pdfname="TICs.pdf",rt=c("raw","corrected")) {
  if (is.null(xcmsSet)) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("./", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
      files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)
  } else {
    files <- filepaths(xcmsSet)
  }
  
  N <- length(files)
  TIC <- vector("list",N)
  
  for (i in 1:N) {
    cat(files[i],"n")
    if (!is.null(xcmsSet) && rt == "corrected")
      rtcor <- xcmsSet@rt$corrected[[i]] else 
        rtcor <- NULL
    TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
  }
  
  pdf(pdfname,w=16,h=10)
  cols <- rainbow(N)
  lty = 1:N
  pch = 1:N
  xlim = range(sapply(TIC, function(x) range(x[,1])))
  ylim = range(sapply(TIC, function(x) range(x[,2])))
  plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time", ylab = "TIC")
  for (i in 1:N) {
    tic <- TIC[[i]]
    points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
  }
  legend("topright",paste(basename(files)), col = cols, lty = lty, pch = pch)
  dev.off()
  
  invisible(TIC)
}
#NEGATIVE MODE
getTICs(xcmsSet=xsg1, pdfname="NEGATIVETICs.pdf",rt="corrected")
#POSITIVE MODE
getTICs(xcmsSet=pxsg3, pdfname="POSITIVETICs.pdf",rt="corrected")
#
mzfiles <- list.files('C:/Users/' ,recursive = TRUE, full=T)
#
#That's all Folks!!
q()
