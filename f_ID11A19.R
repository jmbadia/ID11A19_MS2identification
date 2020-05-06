#+++ ABOUT ------------------------------
#(based on structuredCode_Template.R)

#OBJECTIVE: contain all miscelanious functions  

#TO IMPROVE:

#+++ COMMON VARIABLES-----------------------------
options(stringsAsFactors=FALSE)

##+++ PACKGS & SOURCES ------------------------------
if(!require(lsa)){install.packages("lsa"); require(lsa)}
if(!require(ggplot2)){install.packages("ggplot2"); require(ggplot2)}
if(!require(pbapply)){install.packages("pbapply"); require(pbapply)}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require(mzR)) BiocManager::install("mzR")

#_----------------
#+++ FUNCTIONS----------

#print all the info refered to REFspectras or to spectras of REFidmetabolites
infoSpectra <- function(REFidspectras=NULL, REFidmetabolites=NULL){
  
  if(all(!is.null(REFidspectras),!is.null(REFidmetabolites))){
    print("ERROR: one & only one of arguments are allowed ")
    return()
  }
  
  if((is.null(REFidspectras) & !is.null(REFidmetabolites))){
    REFidspectras <- df_EspectreMetabolit[df_EspectreMetabolit$idmetabolit%in%REFidmetabolites,"idespectre"]
    }
  if(!is.null(REFidspectras)){
    cbind(df_metaespectres[match(REFidspectras,df_metaespectres$idespectre),],df_metametabolits[match(df_EspectreMetabolit[df_EspectreMetabolit$idespectre%in%REFidspectras,"idmetabolit"],df_metametabolits$idmetabolit),])
  }
}

#plot a REF spectra and a UNK spectra
plotSpectra <- function(unk, ref, resultIdent){
  rmatch <- list_fragments$espectre[[which(list_fragments$idespectre==ref)]]
  
  spectredesconegut <- UNK$Spectra$spectra[[which(UNK$Spectra$idUNKSpectra == unk)]]
  #relativitzem
  rmatch["intensity",] <- 100*rmatch["intensity",]/max(rmatch["intensity",])
  spectredesconegut["intensity",] <- 100*spectredesconegut["intensity",]/max(spectredesconegut["intensity",])
  #print(i)
  #print(spectredesconegut);print(rmatch)
  metadata <- resultIdent[resultIdent$UNKidspectra==unk & resultIdent$REFidspectra==ref,]
  df1<-data.frame(x = rmatch[1,],y = rmatch[2,])
  p<-ggplot() +
    geom_linerange(data=df1,aes(x=x, ymax=y, ymin=0),colour="red")+
    coord_cartesian(ylim = c(max(df1$y), -max(df1$y)))
  df2<-data.frame(x = spectredesconegut[1,], y = spectredesconegut[2,])
  p <- p + geom_linerange(data=df2,aes(x=x, ymax=-y, ymin=0))
  text_x <- max(c(df1$x,df2$x))
  text_y <- 95
  p <- p + 
    annotate(geom="text", x=min(c(df1$x,df2$x)), y=-text_y, label=paste0("SC: ", metadata$score,"\nM: ", metadata$MATCHmassNum,"/",metadata$UNKmassNum,"\ncoSim: ", metadata$cossim), hjust = 0, vjust=1)+
    annotate(geom="text", x=min(c(df1$x,df2$x)), y=text_y, label=paste0("REFidmetab.: ", metadata$REFidmetabolite,"\nname: ", metadata$REFmetaboliteName,"\nMmi: ", round(metadata$REFmonoisotMW,3)),color="red",hjust = 0, vjust=0)+
    annotate(geom="text", x=text_x, y=text_y, label=paste0("REFidspectra: ",ref,"\nadduct: ", metadata$REFprecAdduct,"\nprecMZ: ", round(metadata$REFprecMZ,3)),color="red",hjust = 1, vjust=0) + 
    annotate(geom="text", x=text_x, y=-text_y, label=paste0("UNKidspectra: ",unk,"\nadduct?: ", metadata$assmdUNKAdduct,"\nprecMZ: ", round(metadata$UNKprecMZ,3)),hjust = 1, vjust=1)
  p
  #print(ggplotly(p))
}  

binSpectra <- function(spectra, decimals2bin){
  keepRownames<-rownames(spectra[[1]])
  # round spectral masses and sum their intensities up when the rounded mass match
  spectra<-pblapply(spectra, function(x){
    x["mass-charge",]<-round(x["mass-charge",], decimal_bin)
    a<-t(x)
    x<-t(aggregate(a[,"intensity"],by=list(a[,"mass-charge"]),sum))
    rownames(x)<-keepRownames
    return(x)
  })
  return(spectra)
}

#summarize rows with same file into one, grouping its scans and removing the duplicated ones  
optimizeScansList <- function(scans){
  uniqueFilenames <- unique(scans$filename)
  optimizedScans <- scans[0, ]#copy only structure of scans
  for(id_file in seq_along(uniqueFilenames)){
    optimizedScans[id_file,"filename"] <- uniqueFilenames[id_file] 
    posSameFile <- scans$filename==optimizedScans$filename[id_file]
    scansSameFile <- unlist(scans$acquisitionNum[posSameFile])
    #in case appears among the scans a NA, it results in just a NA
    if(any(is.na(scansSameFile))) optimizedScans$acquisitionNum[id_file] <- NA
    else optimizedScans$acquisitionNum[id_file] <- I(list(unique(scansSameFile)))
  }
  return(optimizedScans)
}
