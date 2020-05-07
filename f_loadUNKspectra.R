#+++ ABOUT ------------------------------
#(based on structuredCode_Template.R)

#OBJECTIVE: Load Spectra scans (unknown MS2 spectra in mzml files). 

#SUMMARY: 

#COMMENTS:

#ABBREVIATIONS
#UNK => unknown. Refers to data belonging the spectra we want to identify
#REF => Reference. Refers to data belonging to MS2ID

#TO IMPROVE:

#+++ FUNCTION RELATED ------------------------------
# PARAMETERS:
#=> dirpath: directory containing the mzml sample files. By default, is the working directory
#=> scans: dataframe with filename and acquisition number of the scans we want to acquire. Default value is all the scans 
# Lets choose the unknown spectra to identify. One file per row, defining its filename & all the acquisitionNum we want from that file
#    #e.g.
  #scans <- data.frame(filename="NNK_200_TRAM2.mzML", acquisitionNum = I(list(c(283500:283900))))#acquisitionNum 283500 TO 283900
  #scans <- data.frame(filename="NNK_200_TRAM2.mzML", acquisitionNum = NA)
  #in case we want TO ADD other files spectra
  #scans <- rbind(scans, data.frame(filename="NNK_200_TRAM2.mzML",acquisitionNum = I(list(c(1879,1880)))))

#RESULT: list with 2 items
  #UNK$Metadata => Dataframe with Metadata with an idUNKSpectra
  #UNK$Spectra => Proper spectra with the same idUNKSpectra

#+++ CODE ------------------------------
loadUNKspectra <- function(dirpath=".", scans= NULL){
  if(is.null(scans)){
    mzml_files <- dir(path=dirpath, pattern="*\\.mzml$", ignore.case=TRUE, full.names = FALSE)
  } else{
    scans <- optimizeScansList(scans)
    mzml_files <- scans$filename
  }
    
  # Take and merge metadata & spectralMatrix from MS2 spectra
  UNKMetadata <- data.frame()
  UNKSpectra <- list()
  #Capture MS2 spectra info from all mzmlfiles
  for(id_file in seq_along(mzml_files)){
    arxiu_mzML <- openMSfile(file.path(dirpath, mzml_files[id_file]))
    temp <- header(arxiu_mzML)
    pos2Catch <- temp$msLevel==2# spectra MS2?
    
    if(!is.null(scans)){
      if(!is.na(scans$acquisitionNum[id_file])){
        acqNumb <- as.integer(unlist(scans$acquisitionNum[id_file]))
        pos2Catch <- pos2Catch & (temp$acquisitionNum %in% acqNumb)
      }
    }
    #load metadata & filename
    temp <- cbind(file= mzml_files[id_file], temp[pos2Catch,])#append filename to spectra metadata
    
    # load spectra matrix
    temp_UNKSpectra <- lapply(which(pos2Catch), function(x) {
      m <- t(peaks(arxiu_mzML,x))
      rownames(m) <- c("mass-charge","intensity")
      return(m)
    })
    
    #merge info
    UNKSpectra <- c(temp_UNKSpectra, UNKSpectra)
    UNKMetadata <- rbind(temp,UNKMetadata)
    close(arxiu_mzML)
  } 
  
  # Apply same ID to metadata spectra and metadata matrix
  UNKMetadata <- cbind(idUNKSpectra= seq_len(nrow(UNKMetadata)),UNKMetadata)
  UNKSpectra <- list(idUNKSpectra=seq_len(nrow(UNKMetadata)),spectra = UNKSpectra)
  result <- list(Metadata=UNKMetadata,Spectra=UNKSpectra)
  return(result)
}