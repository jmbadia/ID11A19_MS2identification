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
if (!require(mzR)) BiocManager::install("mzR")

#_----------------
#+++ FUNCTIONS----------

#print all the info refered to REFspectras or to spectras of REFidmetabolites
infoSpectra <- function(REFidspectras=NULL, REFidmetabolites=NULL){
  
  if(all(!is.null(REFidspectras), !is.null(REFidmetabolites))){
    print("ERROR: one & only one of arguments are allowed ")
    return()
  }
  
  if((is.null(REFidspectras) & !is.null(REFidmetabolites))){
    REFidspectras <- DB$df_spectraMetabolite[DB$df_spectraMetabolite$ID_metabolite$ID_metabolite %in% REFidmetabolites,"ID_spectra"]
    }
  if(!is.null(REFidspectras)){
    cbind(DB$df_spectra[match(REFidspectras, DB$df_spectra$ID_spectra), ], DB$df_metabolite[match(DB$df_spectraMetabolite[DB$df_spectraMetabolite$ID_spectra %in% REFidspectras,"ID_metabolite"], DB$df_metabolite$ID_metabolite),])
  }
}

#plot a REF spectra and a UNK spectra
plotSpectra <- function(unkIDspectra, refIDspectra, UNKlist_fragments, REFlist_fragments, resultIdent){
  rmatch <- REFlist_fragments$spectra[[which(REFlist_fragments$ID_spectra==refIDspectra)]]
  
  spectredesconegut <- UNKlist_fragments$spectra[[which(UNKlist_fragments$idUNKSpectra == unkIDspectra)]]
  #relativitzem
  rmatch["intensity",] <- 100*rmatch["intensity",]/max(rmatch["intensity",])
  spectredesconegut["intensity",] <- 100*spectredesconegut["intensity",]/max(spectredesconegut["intensity",])
  #print(i)
  #print(spectredesconegut);print(rmatch)
  metadata <- resultIdent[resultIdent$UNKidSpectra==unkIDspectra & resultIdent$REFidSpectra==refIDspectra,]
  df1<-data.frame(x = rmatch[1,],y = rmatch[2,])
  p<-ggplot() +
    geom_linerange(data=df1,aes(x=x, ymax=y, ymin=0),colour="red")+
    coord_cartesian(ylim = c(max(df1$y), -max(df1$y)))
  df2<-data.frame(x = spectredesconegut[1,], y = spectredesconegut[2,])
  p <- p + geom_linerange(data=df2,aes(x=x, ymax=-y, ymin=0))
  text_x <- max(c(df1$x,df2$x))
  text_y <- 95
  p <- p + 
    annotate(geom="text", x=min(c(df1$x,df2$x)), y=-text_y, label=paste0("\nM: ", metadata$MATCHmassNum,"/",metadata$UNKmassNum,"\ncoSim: ", metadata$cossim), hjust = 0, vjust=1)+
    annotate(geom="text", x=min(c(df1$x,df2$x)), y=text_y, label=paste0("REFidmetab.: ", metadata$REFidMetabolite,"\nname: ", metadata$REFname,"\nMmi: ", round(metadata$REFMmi,3)),color="red",hjust = 0, vjust=0)+
    annotate(geom="text", x=text_x, y=text_y, label=paste0("REFidSpectra: ",refIDspectra,"\nadduct: ", metadata$REFprecAdduct,"\nprecMZ: ", round(metadata$REFprecMZ,3)),color="red",hjust = 1, vjust=0) + 
    annotate(geom="text", x=text_x, y=-text_y, label=paste0("UNKidSpectra: ",unkIDspectra,"\nadduct?: ", metadata$assmdUNKAdduct,"\nprecMZ: ", round(metadata$UNKprecMZ,3)),hjust = 1, vjust=1)
  p
  #print(ggplotly(p))
}  

binSpectra <- function(spectra, decimals2bin){
  keepRownames<-rownames(spectra[[1]])
  # round spectral masses and sum their intensities up when the rounded mass match
  spectra<-pblapply(spectra, function(x){
    x["mass-charge",]<-round(x["mass-charge",], decimals2bin)
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

createExcel <- function(data, filename){
  require(openxlsx)
  
  ## style for header
  headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#596e79", border="TopBottom", borderColour = "#4F81BD")
  
  ## style for body 
  bodyStyle <- createStyle(border="TopBottom", borderColour = "#4F81BD")
  m <- which(diff(data$UNKidSpectra)!=0)
  
  #Every precursor mass has alternate color
  StyleA1 <- createStyle(border="TopBottom", fgFill = "#ff9a3c")
  StyleA2 <- createStyle(border="TopBottom", fgFill = "#ffc93c")
  StyleB1 <- createStyle(border="TopBottom", fgFill = "#c8f4de")
  StyleB2 <- createStyle(border="TopBottom", fgFill = "#79a8a9")
  
  lengCol <- seq_len(ncol(data))
  cossimCol <- which(colnames(data)=="cossim")
  nameCol <- which(colnames(data)=="REFname")
  inchikeyCol <- which(colnames(data)=="REFinchikey")
  
  
  ## Create a new workbook
  wb <- createWorkbook(creator = "jmbadia", title="My name here")
  
  noInchikey <- is.na(data$REFinchikey)
  data$REFinchikey[!noInchikey] <- paste0('HYPERLINK("https://www.ncbi.nlm.nih.gov/pccompound?term=%22',data$REFinchikey[!noInchikey], '%22[InChIKey]", "',data$REFinchikey[!noInchikey],'")')
  
  
  for(idFile in unique(data$file)) {
    #select data
    tmp <- data[data$file==idFile,]
    
    #only first 30 characters of the filename
    idFile <- substr(idFile, start = 1, stop = min(30,nchar(idFile)))
    
    #for every precMass, which have the best cossim
    bst_cossim <- unlist(lapply(unique(tmp$UNKprecMZ), function(x) {
      maxCossim <- max(tmp[tmp$UNKprecMZ==x,"cossim"])
      which(tmp$cossim==maxCossim & tmp$UNKprecMZ==x)
    }))
    #for every precMass, which name repeats the most
    bst_name <- unlist(lapply(unique(tmp$UNKprecMZ), function(x) {
      mostRepName <- names(sort(table(tmp[tmp$UNKprecMZ==x,"REFname"]),decreasing=TRUE)[1])
      which(tmp$REFname==mostRepName & tmp$UNKprecMZ==x)
    }))
    
    ## Add a worksheet
    addWorksheet(wb, idFile, gridLines = FALSE) 
    
    ##write data to worksheet 1
    writeData(wb, sheet = idFile, tmp)
    #makeHyperlinkString(sheet = idFile, row = 1, col = 1, text = "NULL")
    colorA <- unique(tmp$UNKprecMZ)[c(TRUE, FALSE)]
    colorB <- unique(tmp$UNKprecMZ)[c(FALSE, TRUE)]
    
    color1 <- unique(tmp$UNKidSpectra)[c(TRUE, FALSE)]
    color2 <- unique(tmp$UNKidSpectra)[c(FALSE, TRUE)]
    
    addStyle(wb, sheet = idFile, headerStyle, rows = 1, cols = seq_len(ncol(tmp)), gridExpand = TRUE)
    addStyle(wb, sheet = idFile, StyleA1, rows = 1+which(tmp$UNKprecMZ %in% colorA & tmp$UNKidSpectra %in% color1), cols = lengCol, gridExpand = TRUE)
    addStyle(wb, sheet = idFile, StyleA2, rows = 1+which(tmp$UNKprecMZ %in% colorA & tmp$UNKidSpectra %in% color2), cols = lengCol, gridExpand = TRUE)
    addStyle(wb, sheet = idFile, StyleB1, rows = 1+which(tmp$UNKprecMZ %in% colorB & tmp$UNKidSpectra %in% color1), cols = lengCol, gridExpand = TRUE)
    addStyle(wb, sheet = idFile, StyleB2, rows = 1+which(tmp$UNKprecMZ %in% colorB & tmp$UNKidSpectra %in% color2), cols = lengCol, gridExpand = TRUE)
    
    addStyle(wb, sheet = idFile, createStyle(textDecoration = "bold", fontColour = "#d01257"), rows = 1+bst_cossim, cols = cossimCol, gridExpand = TRUE, stack = TRUE)  
    
    addStyle(wb, sheet = idFile, createStyle(textDecoration = "bold"), rows = 1+bst_name, cols = nameCol, gridExpand = TRUE, stack = TRUE)  
    
    #add link
    writeFormula(wb, sheet =idFile, startRow = 2, startCol = inchikeyCol, x = tmp$REFinchikey)
    
  }
  saveWorkbook(wb, paste0(filename,".xlsx"), overwrite = TRUE)
}