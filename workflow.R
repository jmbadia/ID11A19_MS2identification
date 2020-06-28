#+++ WORKING VARIABLES ------------------------------
workVariables <- list(
#Project variables
  projectName = "Milk_stck2Addct_CONS",
  identDate = format(Sys.time(),"%Y%m%d_%H%M"),
  mail2SendResult = "josepmaria.badia@gmail.com",
  projectDirectory = file.path("invisible2Git/MS2ID_Identifications/20200608_Milk"),

#Consensus spectra variables
  useCons = FALSE,
  dec2binCons= 2, #decimals to bin precMZ and decide whose are under the same bin so must resume their spectra in a consensus spectra,
  mzdiff = 0.02, #max mz difference among fragments (same spectra or not) in order to consider them as the same fragment
  minProp = 2/3, #min presence of a fragment in order to be on the consensus spectra

#Fragment bin & filter variables
  dec2binFrag = 2, # decimal to bin with. A 2 value means a 0.01 binning i.e. 78.04,78.05,78.06.....
  minNoiseAllowed = 0.01, # minimum intensity allowed in the noise filtering,

#Identification variables
  stck2Add = TRUE,
  samePrecMass = FALSE,
  cosSimTh = 0.8,
  massErrorTh = 10,
  candNum=20,
  cl=2L,
  minCommMassNum=2, 
  topMassNum = 5,
  REFfilt_db = NULL
)

#remove useless variables
if(!workVariables$useCons) workVariables$dec2binCons <- workVariables$mzdiff  <- 
    workVariables$minProp <- NULL




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
require(googledrive)
require(dplyr)
require(glue)
require(MSnbase)
require(blastula)
require(keyring)
require(enviPat)
require(zip)

source("f_ID11A19.R")
source("f_loadUNKspectra.R")
source("f_identify.R")
adducts <- readRDS("resources/adductsTable.rds")#load adducts table (enviPat) and change positive or negative polarity with usual terms (1 or 0)
adducts[adducts$Ion_mode=="positive","Ion_mode"] <- 1
adducts[adducts$Ion_mode=="negative","Ion_mode"] <- 0


#+++ CODE ------------------------------

#.------------------------------
# 1. LOAD SAMPLES ------------------------------

#_1.1 Load scans desiderata ----
# File with list of scans to identify (column n.1 (file) and column n.5  (acquisitionNum))
#load(file.path(bigdata_folder,"rawdataCarla/3_MET_LC-MS_MARCH20.Rdata"))

##_1.2 Load UNK spectra ----
# (acording scans desiderata)
samplefolder <- file.path("invisible2Git/MS2ID_Identifications/20200608_Milk")
UNK <- loadUNKspectra(dirpath = file.path(workVariables$projectDirectory,"rawData"))

#.------------------------------


#.------------------------------
# 2. PREPARE SAMPLES ------------------------------

##_2.0 Clean Spectra ----
#lets remove nonconventional matrices 
validSpectra <- vapply(UNK$Spectra$spectra, function(x) {nrow(x)==2 & ncol(x)>0}, FUN.VALUE = T)
noValidIdSpectra <- UNK$Spectra$idUNKSpectra[!validSpectra]
UNK$Spectra$idUNKSpectra <- UNK$Spectra$idUNKSpectra[which(validSpectra)] 
UNK$Spectra$spectra <- UNK$Spectra$spectra[which(validSpectra)]
UNK$Metadata <-UNK$Metadata[!UNK$Metadata$idUNKSpectra %in% noValidIdSpectra,] 

# 1BIS. Consensus Spectra ------------------------------
# Summarize different scans with the same precursor mass
#TO DO: A) WHICH RT? B)IF ONLY 2 SCANS C)IF SCANS SAME MZ PREC in different RT?
if(workVariables$useCons){
  CONS<-list()
  CONS$Metadata <- data.frame()
  CONS$Spectra <- list(idUNKSpectra=vector(),spectra=list())
  
  #file by file
  for (idfile in unique(UNK$Metadata$file)){
    #idfile<-"X20191015_013.mzML"
    TEMP_Metadata<-UNK$Metadata[UNK$Metadata$file==idfile,]
    TEMP_spectra<- UNK$Spectra$spectra[which(UNK$Spectra$idUNKSpectra %in% TEMP_Metadata$idUNKSpectra)]
    TEMP_idUNKSpectra<- UNK$Spectra$idUNKSpectra[which(UNK$Spectra$idUNKSpectra %in% TEMP_Metadata$idUNKSpectra)]
    
    #rounded precMZ
    rounded<-round(TEMP_Metadata$precursorMZ, workVariables$dec2binCons)
    #positions grouped considering unique(rounded)
    precMZ_grpd<-lapply(unique(rounded), function(precMZ_flag) which(rounded==precMZ_flag))
    #realprecMZ_grpd<-lapply(precMZ_grpd, function(x) TEMP_Metadata$precursorMZ[x])  
    #calculate consensus (and new rt and precMZ)
    spectra_grpd<-lapply(precMZ_grpd, function(x) unlist(lapply(x, function(y) which(TEMP_idUNKSpectra==TEMP_Metadata$idUNKSpectra[y]))))
    
    newSpectra <- lapply(spectra_grpd, function(x){
      #create a list of spectra
      Splist <- lapply(x, function(idsp){
        metadata<-TEMP_Metadata[TEMP_Metadata$idUNKSpectra==TEMP_idUNKSpectra[idsp],c("retentionTime","precursorMZ")]
        new("Spectrum2", rt = metadata[1,1], precursorMz = metadata[1,2], 
            mz = TEMP_spectra[[idsp]]["mass-charge",],
            intensity = TEMP_spectra[[idsp]]["intensity",])
      })
      consSp <- consensusSpectrum(Splist, mzd = workVariables$mzdiff, minProp = workVariables$minProp)
      #only when results a spectra
      if(peaksCount(consSp)==0){
        spectrMatrix<-NULL
      } else {
        spectrMatrix<-matrix(data=c(mz(consSp),intensity(consSp)),nrow=2, byrow = T,dimnames = list( c("mass-charge","intensity"), NULL))
      }
      list(matrix=spectrMatrix, rt=rtime(consSp), precMZ=precursorMz(consSp))
    })
    
    #We keep as REPRESENTATIVE metadata the metadata of the FIRST spectra
    keep_acqNum<-vapply(precMZ_grpd, function(x) paste(TEMP_Metadata$acquisitionNum[x], collapse=","), FUN.VALUE = "ko")
    TEMP_Metadata<-TEMP_Metadata[vapply(precMZ_grpd, `[[`, 1, FUN.VALUE = 1), colnames(TEMP_Metadata)%in% c("idUNKSpectra", "file", "msLevel", "polarity", "collisionEnergy", "ionisationEnergy","centroided")]
    TEMP_Metadata$acquisitionNum<- keep_acqNum
    TEMP_Metadata$retentionTime <- vapply(newSpectra, `[[`, 2, FUN.VALUE = 1.2)
    TEMP_Metadata$precursorMZ <- vapply(newSpectra, `[[`, 3, FUN.VALUE = 1.2)
    TEMP_idUNKSpectra<-TEMP_idUNKSpectra[vapply(precMZ_grpd, `[[`, 1, FUN.VALUE = 1)]
    TEMP_spectra<-lapply(newSpectra, `[[`, 1)
    CONS$Metadata<-rbind(CONS$Metadata,TEMP_Metadata)
    CONS$Spectra$idUNKSpectra<-c(CONS$Spectra$idUNKSpectra,TEMP_idUNKSpectra)
    CONS$Spectra$spectra<-c(CONS$Spectra$spectra,TEMP_spectra)
  }
  
  #Remove results with NO achieved consensus
  Cons<-!vapply(CONS$Spectra$spectra, is.null, FUN.VALUE = T)
  CONS$Spectra$spectra <- CONS$Spectra$spectra[which(Cons)]
  idCons<-CONS$Spectra$idUNKSpectra[which(Cons)]
  CONS$Spectra$idUNKSpectra<-CONS$Spectra$idUNKSpectra[which(Cons)]
  CONS$Metadata<-CONS$Metadata[CONS$Metadata$idUNKSpectra%in%idCons,]
  
  UNK<-CONS
}

##_2.1 Filter spectral noise ----
# Remove fragments with intensity < 1% base peak (minNoiseAllowed=0.01)
UNK$Spectra$spectra <- lapply(UNK$Spectra$spectra, function(x) {
  x[,x["intensity",] > workVariables$minNoiseAllowed * max(x["intensity",]), drop=FALSE]
})

##_2.2 Bin UNK spectra ----
UNK$Spectra$spectra <- binSpectra(spectra=UNK$Spectra$spectra, workVariables$dec2binFrag)
  
#.------------------------------
# 3. LOAD SPECTRAL DB ------------------------------
# LOAD precooked MS2ID DATABASE (i.e. already spectra mass-sorted (essential), filtered (minNoiseAllowed=0.01) and binned (binning 0.01))
DB <- readRDS(file.path("invisible2Git/rawData","MS2ID_binFiltered_20200605_102226.rds"))
  
#in order to avoid confussion with the UNK variables we modify the REF variables' names, adding REF to df_spectra and df_metabolite variables 
colnames(DB$df_spectra)[!colnames(DB$df_spectra) %in% c("ID_spectra","primalIdSpectra")] <- paste("REF",colnames(DB$df_spectra)[!colnames(DB$df_spectra) %in% c("ID_spectra","primalIdSpectra")],sep="")
colnames(DB$df_metabolite)[!colnames(DB$df_metabolite) %in% c("ID_metabolite","primalIdMetab")] <- paste("REF",colnames(DB$df_metabolite)[!colnames(DB$df_metabolite) %in% c("ID_metabolite","primalIdMetab")],sep="")

#.------------------------------
# 4. IDENTIFICATION ------------------------------
# see identifySpectra() documentation
result <- identifySpectra(UNKdata = UNK,
                          db = REFfilt_db,
                stick2Adducts = workVariables$stck2Add,
                cosSimTh = workVariables$cosSimTh,
                massErrorAllowed = workVariables$massErrorTh,
                samePrecMass = workVariables$samePrecMass,
                candidatesNumber = workVariables$candNum,
                cl = workVariables$cl,
                minCommonMassNum = workVariables$minCommMassNum,
                topMassesNum = workVariables$topMassNum)


#saveRDS(result, file=file.path(workVariables$projectDirectory, paste0(workVariables$projectName,"-",workVariables$identDate,".rds")))
result <- readRDS(file=file.path(workVariables$projectDirectory,"resultMilk_stck2Add.rds"))

#.------------------------------
# 5. POST-ID TOOLS ------------------------------

#plot one
#plotSpectra(unk=8071, ref=106216, result)

#obtain info
#infoSpectra(8071, 106216)

#dir to save files to zip
tmpdir <- tempdir()
filenameBase<- paste0(workVariables$projectName,"_",workVariables$identDate)
sub_tmpdir <- file.path(tmpdir,filenameBase)
dir.create(sub_tmpdir, showWarnings = FALSE)

#create EXCEL file with the summary
#summary: For every UNK spectra and REFmetabolite, ONLY THE BEST cossim
summary_result <- result %>% group_by(UNKidSpectra, REFidMetabolite) %>% top_n(1, cossim) %>% distinct(UNKidSpectra, REFidMetabolite, .keep_all = T) %>% arrange(UNKprecMZ, UNKidSpectra)
createExcel(data= summary_result, filename=file.path(sub_tmpdir,filenameBase))

#prepare LOT FILE (file with all the info final user needs to analize the identification result)
REFidspectra_involved <-  which(DB$list_fragments$ID_spectra %in% result$REFidSpectra)
REFfragments <- DB$list_fragments
REFfragments$ID_spectra <- REFfragments$ID_spectra[REFidspectra_involved] 
REFfragments$spectra <- REFfragments$spectra[REFidspectra_involved] 
UNKfragments <- UNK$Spectra
lot <- list(metadata=result, UNKfragments=UNKfragments, REFfragments=REFfragments, IDvariables=workVariables) 
saveRDS(lot, file.path(sub_tmpdir,paste0("lot_",filenameBase,".rds")))

zipFile<-file.path(workVariables$projectDirectory,paste0(filenameBase,".zip"))
zipr(zipfile=zipFile, sub_tmpdir, include_directories = F)

#save zip to gDrive
drive_auth(email = "josepmaria.badia@gmail.com")
chicken <- drive_upload(
  zipFile,
  basename(zipFile)
)
#share the file
chicken <- chicken %>% drive_share(role = "reader", type = "anyone")
link2share <- paste0("https://drive.google.com/uc?export=download&id=",chicken$id)
#https://drive.google.com/uc?export=download&id=1GNffEuoy2fcOJ-pi2921eXYuEbCohM4_

#to save gmail credentials first time
#sudo apt-get install libsecret-1-dev
#create_smtp_creds_key(id = "gmail_creds_forR",  provider = "gmail",  user = "josepmaria.badia@gmail.com")
vars2Email <- paste(names(workVariables), workVariables, sep=" = ", collapse=", ")
email <- compose_email(
  body = md(glue(
    "
### {workVariables$projectName} project identification
Please, download <a href='{link2share}' target='_blank'>here</a> a zip with the result of identification. It contains:  
* **A lot file** with all the resulting ID data. Use the <a href='https://jmbadia.shinyapps.io/MS2IDbrowser/' target='_blank'>MS2ID browser</a> in order to visualize and compare sample MS2 spectra (UNKnown) with their identification spectra (REFerence library spectra)
.  
*  **An Excel file** with the results distributed by samples. Use it to annotate the comments infered usng the MS2ID browser.  

The following are the variables used on the identification workflow. Please consult <a href='https://jmbadia.github.io/ID11A19_MS2identification/IDvariables.html' target='_blank'>here</a> their meaning if necessary.  
>  {vars2Email}
"   
  ))
)
email %>%
  smtp_send(
    from = "josepmaria.badia@gmail.com",
    to = workVariables$mail2SendResult,
    subject = paste(workVariables$projectName, "ID results"),
    credentials = creds_key("gmail_creds_forR")
  )
