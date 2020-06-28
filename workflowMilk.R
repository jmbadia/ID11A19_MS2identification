#+++ ABOUT ------------------------------
#(based on structuredCode_Template.R)

#OBJECTIVE: MS2 identification Workflow 

#SUMMARY: Using a list of MS2 scans to identify, load them from mzml sample files, processed them (sort, bin & filter) and apply cossim against a local spectra library in order to identify its metabolite

#COMMENTS:

#ABBREVIATIONS
#UNK => unknown. Refers to data belonging the spectra we want to identify
#REF => Reference. Refers to data belonging to MS2ID

#TO IMPROVE:
## Apply superspectra Roger's concept (see CliqueMSMS.R)
# Parse formula on HMDB & MoNA databases. Calculate all REF neutral mass using formula variable.
# Use experimental mass range to subset the cossim mass range application   
# Consider precursor CHARGE on the analysis of mz
# Add other REF library filters, like ionization or fragmentation type
# Remove cosinesim package dependency


#+++ PACKGS & SOURCES ------------------------------
require(dplyr)
source("f_ID11A19.R")
source("f_loadUNKspectra.R")
source("f_identify.R")
load("resources/adductsTable.RData")#load adducts table (enviPat) and change positive or negative polarity with usual terms (1 or 0)
adducts[adducts$Ion_mode=="positive","Ion_mode"] <- 1
adducts[adducts$Ion_mode=="negative","Ion_mode"] <- 0

#+++ THISFILE VARIABLES ------------------------------

# working variables. MUST COINCIDE WITH MS2ID filter and binning conditions
decimal_bin=2# decimal to bin with. A 2 value means a 0.01 binning i.e. 78.04,78.05,78.06.....
minNoiseAllowed <- 0.01# minimum intensity allowed in the noise filtering

# Path tho Big Data
base <- file.path("/home/jmbadia/esborra_mesuresperR")
metabolites_database_directory<-file.path(base,"metabolites_database")
bigdata_folder<-file.path(base,"ID11A19_identificarEspectresReals")


  #+++ CODE ------------------------------

#.------------------------------
# 1. LOAD SAMPLES ------------------------------

#_1.1 Load scans desiderata ----
# File with list of scans to identify (column n.1 (file) and column n.5  (acquisitionNum))
#load(file.path(bigdata_folder,"rawdataCarla/3_MET_LC-MS_MARCH20.Rdata"))

##_1.2 Load UNK spectra ----
# (acording scans desiderata)
samplefolder <- file.path("/home/jmbadia/Desktop","mzMLMilk")
UNK <- loadUNKspectra(dirpath = samplefolder)

#.------------------------------
# 2. PREPARE SAMPLES ------------------------------

##_2.1 Filter spectral noise ----
#lets remove nonconventional matrices 
validSpectra <- vapply(UNK$Spectra$spectra, function(x) {nrow(x)==2 & ncol(x)>0}, FUN.VALUE = T)
noValidIdSpectra <- UNK$Spectra$idUNKSpectra[!validSpectra]
UNK$Spectra$idUNKSpectra <- UNK$Spectra$idUNKSpectra[which(validSpectra)] 
UNK$Spectra$spectra <- UNK$Spectra$spectra[which(validSpectra)]
UNK$Metadata <-UNK$Metadata[!UNK$Metadata$idUNKSpectra %in% noValidIdSpectra,] 

# Remove fragments with intensity < 1% base peak (minNoiseAllowed=0.01)
UNK$Spectra$spectra <- lapply(UNK$Spectra$spectra, function(x) {
  x[,x["intensity",] > minNoiseAllowed * max(x["intensity",]), drop=FALSE]
})

##_2.2 Bin UNK spectra ----
UNK$Spectra$spectra <- binSpectra(spectra=UNK$Spectra$spectra, decimals2bin)
  
#.------------------------------
# 3. LOAD SPECTRAL DB ------------------------------
# LOAD precooked MS2ID DATABASE (i.e. already spectra mass-sorted (essential), filtered (minNoiseAllowed=0.01) and binned (binning 0.01))
DB <- readRDS(file=file.path(metabolites_database_directory,"MS2ID_binFiltered_20200605_102226.rds"))
  
#in order to avoid confussion with the UNK variables we modify the REF variables' names, adding REF to df_spectra and df_metabolite variables 
colnames(DB$df_spectra)[!colnames(DB$df_spectra) %in% c("ID_spectra","primalIdSpectra")] <- paste("REF",colnames(DB$df_spectra)[!colnames(DB$df_spectra) %in% c("ID_spectra","primalIdSpectra")],sep="")
colnames(DB$df_metabolite)[!colnames(DB$df_metabolite) %in% c("ID_metabolite","primalIdMetab")] <- paste("REF",colnames(DB$df_metabolite)[!colnames(DB$df_metabolite) %in% c("ID_metabolite","primalIdMetab")],sep="")
  
#.------------------------------
# 4. IDENTIFICATION ------------------------------
# see identifySpectra() documentation
result <- identifySpectra(UNKdata = UNK, stick2Adducts=T, cosSimTh=0.8, massErrorAllowed=10, samePrecMass=F, candidatesNumber=20, cl=2L, minCommonMassNum=2, topMassesNum = 5)

#FERHO F i F i després a l'excel q es vegi els T T, per veure si en podem treure algo de la resta de considerats
#saveRDS(result, file="resultMilk_stck2Add.rds")
#.------------------------------
# 5. POST-ID TOOLS ------------------------------
  
#SUMMARY: For every UNK spectra and REFmetabolite, ONLY THE BEST cossim
summary_result <- result %>% group_by(UNKidSpectra, REFidMetabolite) %>% top_n(1, cossim) %>% distinct(UNKidSpectra, REFidMetabolite, .keep_all = T) %>% arrange(UNKidSpectra)
  
plotSpectra(unk=59607, ref=471297, result)
infoSpectra()


#https://www.ncbi.nlm.nih.gov/pccompound?term=%22BSYNRYMUTXBXSQ-UHFFFAOYSA-N%22[InChIKey]
#makeHyperlinkString(sheet, row = 1, col = 1, text = NULL, file = NULL)
a <- summary_result
require(openxlsx)

## create and add a style to the column headers
headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#596e79", border="TopBottom", borderColour = "#4F81BD")

## style for body 
bodyStyle <- createStyle(border="TopBottom", borderColour = "#4F81BD")
m <- which(diff(a$UNKidSpectra)!=0)

#Every precursor mass has alternate color
StyleA1 <- createStyle(border="TopBottom", fgFill = "#ff9a3c")
StyleA2 <- createStyle(border="TopBottom", fgFill = "#ffc93c")
StyleB1 <- createStyle(border="TopBottom", fgFill = "#c8f4de")
StyleB2 <- createStyle(border="TopBottom", fgFill = "#79a8a9")

lengCol <- seq_len(ncol(a))
cossimCol <- which(colnames(a)=="cossim")
nameCol <- which(colnames(a)=="REFname")
inchikeyCol <- which(colnames(a)=="REFinchikey")

#we cannot use - in a excel hyperlink, so we replace them by its homonim %2D (https://www.w3schools.com/tags/ref_urlencode.asp)
#alink <- gsub("-","%2D",a$REFinchikey)
    

## Create a new workbook
wb <- createWorkbook(creator = "jmbadia", title="My name here")

noInchikey <- is.na(a$REFinchikey)
a$REFinchikey[!noInchikey] <- paste0('HYPERLINK("https://www.ncbi.nlm.nih.gov/pccompound?term=%22',a$REFinchikey[!noInchikey], '%22[InChIKey]", "',a$REFinchikey[!noInchikey],'")')


for(idFile in unique(a$file)) {
  #only first 30 charcaters of the filename
  idFile <- substr(idFile, start = 1, stop = min(30,nchar(idFile)))
  
  #select data
  tmp <- a[a$file==idFile,]
  
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
saveWorkbook(wb, "Clairo_ALLDB.xlsx", overwrite = TRUE)



