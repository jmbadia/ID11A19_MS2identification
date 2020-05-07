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
source("f_ID11A19.R")
source("f_loadUNKspectra.R")
source("f_identify.R")
load("resources/adductsTable.RData")#load adducts table (enviPat) and change positive or negative polarity with usual terms (1 or 0)
adducts[adducts$Ion_mode=="positive","Ion_mode"] <- 1
adducts[adducts$Ion_mode=="negative","Ion_mode"] <- 0

#+++ THISFILE VARIABLES ------------------------------
# working variables
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
load(file.path(bigdata_folder,"rawdataCarla/3_MET_LC-MS_MARCH20.Rdata"))

# structure scans list as a dataframe (see loadUNKspectra() documentation)
scans <- data.frame()  
for(idRow in seq_len(nrow(MET))){
  scans <- rbind(scans,
                 data.frame(filename=MET$file[idRow],
                            acquisitionNum=I((MET$acquisitionNum[idRow]))))
}

##_1.2 Load UNK spectra ----
# (acording scans desiderata)
samplefolder <- file.path(bigdata_folder,"rawdataCarla")
UNK <- loadUNKspectra(dirpath = samplefolder, scans = scans)


#.------------------------------
# 2. PREPARE SAMPLES ------------------------------

##_2.1 Filter spectral noise ----
# Remove fragments with intensity < 1% base peak (minNoiseAllowed=0.01)
UNK$Spectra$spectra <- lapply(UNK$Spectra$spectra, function(x) {
  x[,x["intensity",] > minNoiseAllowed * max(x["intensity",]), drop=FALSE]
})

##_2.2 Bin UNK spectra ----
UNK$Spectra$spectra <- binSpectra(spectra=UNK$Spectra$spectra, decimals2bin)
  
#.------------------------------
# 3. LOAD SPECTRAL DB ------------------------------
# LOAD precooked MS2ID DATABASE (i.e. already spectra mass-sorted (essential), filtered (minNoiseAllowed=0.01) and binned (binning 0.01))
  load(file=file.path(metabolites_database_directory,"sp_MS2ID_filtradaiBinejada_20200424_142048.RData"))
  
# Procedures to incorporate in the precooking algorithm asap
  #Remove spectra with no metabolite associated (HMDB's fault)
  idspectraNoidmetabolite <- df_EspectreMetabolit[is.na(df_EspectreMetabolit$idmetabolit),"idespectre"]
  df_EspectreMetabolit <- df_EspectreMetabolit[!is.na(df_EspectreMetabolit$idmetabolit),]
  df_metaespectres <- df_metaespectres[!df_metaespectres$idespectre%in%idspectraNoidmetabolite,]
  spectraIndex <- which(!list_fragments$idespectre%in%idspectraNoidmetabolite)
  list_fragments$idespectre <- list_fragments$idespectre[spectraIndex]
  list_fragments$espectre <- list_fragments$espectre[spectraIndex]
  
  #in order to avoid confussion with the UNK variables we modify the REF variables' names, adding REF to df_metaespectres and df_metametabolits variables
  colnames(df_metaespectres)[!colnames(df_metaespectres) %in% c("idespectre","BDespectreoriginal")] <- paste("REF",colnames(df_metaespectres)[!colnames(df_metaespectres) %in% c("idespectre","BDespectreoriginal")],sep="")
  colnames(df_metametabolits)[!colnames(df_metametabolits) %in% c("idmetabolit","BDmetabolitoriginal")] <- paste("REF",colnames(df_metametabolits)[!colnames(df_metametabolits) %in% c("idmetabolit","BDmetabolitoriginal")],sep="")
  
  #Obtain metabolite monoisotopicMW using formula
  # because original monoisotopicMW contains a 5% NA.
  # Avoiding molecular ion formulas, that gets an error. In that case we dont modify the default metabolite monoisotopicMW (NA or not). 
  # PROBLEM: 95% metabolites have no FORMULA (HMDB & MoNA). Not useful until we reparse both db
    require(Rdisop)
    YAformula <- !is.na(df_metametabolits$REFformula)
    calc_monoiMW <- lapply(df_metametabolits$REFformula[YAformula], function(x) tryCatch({getMolecule(x)$exactmass}, error=identity))
    noerror <- !vapply(calc_monoiMW, is, logical(1), "error")
    df_metametabolits$REFmonoisotopic_molecular_weight[YAformula][noerror] <- unlist(calc_monoiMW[noerror])
  
    #order df_metaMetabolites by neutral mass. Need it to identify() (more precisely, in order to find faster the REF spectra with the neutral mass we want)
  df_metametabolits <- df_metametabolits[order(df_metametabolits$REFmonoisotopic_molecular_weight, decreasing=FALSE ),]
 
   #order by precursor mass the df_metaespectres. Need it to identify() (more precisely, in order to find faster the REF spectra with the precurso mass we want)
  df_metaespectres <- df_metaespectres[order(df_metaespectres$REFprecursor_mass, decreasing=FALSE ),]
  
#.------------------------------
# 4. IDENTIFICATION ------------------------------
# see identifySpectra() documentation
result <- identifySpectra(UNKdata = UNK, stick2Adducts=T, cosSimTh=0.8, massErrorAllowed=10, samePrecMass=F, candidatesNumber=20, cl=NULL, minCommonMassNum=2, topMassesNum = 5)

#.------------------------------
# 5. POST-ID TOOLS ------------------------------
  
#SUMMARY: For every UNK spectra and REFmetabolite, ONLY THE BEST cossim
summary_result <- result[,c("UNKidspectra", "REFidspectra","REFidmetabolite", "cossim","assmdUNKAdduct","REFprecAdduct","UNKprecMZ","REFmonoisotMW","REFinchikey")] %>% group_by(UNKidspectra, REFidmetabolite) %>% top_n(1,cossim) %>% distinct(UNKidspectra, REFidmetabolite, .keep_all = T)
  
plotSpectra(unk=305, ref=126709, result)
infoSpectra()
#write.csv(result, file="selectedScansAllDB.csv")
