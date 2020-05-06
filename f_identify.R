#+++ ABOUT ------------------------------
#(based on structuredCode_Template.R)

#OBJECTIVE: Given an unknown MS2 spectrum (or UNK spectrum), to obtain its most similar spectra (and its metadata) from a spectra library (REF spectrum). We use the variable "score" as a score variable (range 0-1); it is based on cosine similarity comparison and number of mz coincidences REF - UNK spectra. 
# (optional previous subsetting) It is possible to subset the REF spectra library using the arguments "db" (original data base), "stick2Adducts" (only REF metabolites with neutral mass compatible with possible adducts of UNK spectra), "samePrecMass" (indicating if REF and UNK precursor mz value should be equal) and topMassesNum (Nth most intensity fragments, both REF and UNK, which must have at list a common fragment). Subsetting arguments have and AND relation (REF spectra must accomplish ALL the conditions, e.g db="metlin" AND stick2Adducts=TRUE AND samePrecMass=TRUE) 
# The REF spectra (subsetted or not) is then filtered by fragment. This step is about finding REF spectra with, at list, a minimum number (minCommonMassNum argument) of REF - UNK common fragments, becoming REF spectra applicant to cos sim algorithm; the score algorithm will be applied only to such applicant REF spectra. 

# ARGUMENTS:
# =>UNKdata: list with unknown spectra and its metadata

# =>db: REF library subsetting argument. Vector with the db name we want to use as REF spectra library. (e.g. db=(db=c("metlin","nist")), default => no db restrictions)

#=>samePrecMass: REF library subsetting argument. Boolean indicating if REF precursor mass and UNK precursor mass must be equals. NA's values are always non excluyent (default samePrecMass=FALSE) 

#=>stick2Adducts: REF library subsetting argument. Boolean indicating if REF neutral mass must be compatible with UNK adducts (considering UNK precursor). We considere adducts the enviPat adducts list with UNK spectra polarization (positive or negative). NA's values are always non excluyent.

#=>topMassesNum: REF library subsetting argument. Only REF spectra with spectra whose top N (topMassesNum value) fragments share at list one mz with top N (topMassesNum value) UNK fragments (default, no subsetting)

#=>minCommonMassNum: Filter by fragment argument. Resulting REF spectra must have at least a number N (minCommonMassNum value) of common fragments with UNK spectra (default minCommonMassNum=2)

# =>cosSimTh: cossine similarity threshold. REF spectra applicants with lower values than cosSimTh are eliminated (e.g. cosSimTh=0.9, default 0.8)

# =>candidatesNumber: Lastly, number of best cossim REF spectra we keep (per UNK spectra) on the final result (default candidatesNumber=100)

#=>massErrorAllowed: mass  ppm error contemplated on stick2Adducts and samePrecMass arguments (default massErrorAllowed=5 ppm)

#=>cl: Processors to paralelize with the function. ONLY LINUX VERSION, otherwise, serial use (default value cl=NULL) should be used

#RETURN
#Dataframe with the proposal identifications and all its metadata (UNK and REf related)

#TIME BENCHMARKING
#Common Variables: intel core i7 4 cores, 16GRam, no SSD disk, cosSimTh=0.8, massErrorAllowed=5, candidatesNumber=20, cl=NULL, minCommonMassNum=2, REF library with 579390 spectra
# =>0.13 sec/scan. Variables ionMassInc=c(1.007276466812), samePrecMass=TRUE, topMassesNum = 5
#=>0.094 sec/scan: Variables ionMassInc=c(1.007276466812), samePrecMass=TRUE, topMassesNum = NULL    
#=>0.83 sec/scan: Variables ionMassInc=NULL, samePrecMass=TRUE, topMassesNum = 5    
#=>1.80 sec/scan: Variables ionMassInc=c(1.007276466812), samePrecMass=TRUE, topMassesNum = NULL    
#=>2.86 sec/scan: Variables ionMassInc=NULL, samePrecMass=TRUE, topMassesNum = NULL
#=>5.26 sec/scan: Variables ionMassInc=NULL, samePrecMass=NULL, topMassesNum = NULL

#TO IMPROVE: 

#COMMENTS

identifySpectra <- function(UNKdata , cl=NULL, db=NULL, topMassesNum=NULL, ...){
  
  #STEP1: SUBSET REF DataBase acording to db argument
  listfragOfInterest <- list()
  if(!is.null(db)){
    id <- df_metaespectres[tolower(df_metaespectres$BDespectreoriginal)%in%tolower(db),"idespectre"]
    filtered_positions <- which(list_fragments$idespectre %in% id)
    listfragOfInterest$idespectre <- list_fragments$idespectre[filtered_positions]
    listfragOfInterest$espectre <- list_fragments$espectre[filtered_positions]
  }else {
    listfragOfInterest <- list_fragments    
  } 
  
  #obtain values we'll need on the next function df_metametabolits & df_metaespectresTrimm without monoisotopic_molecular_weight or precursor_mass respectivaly. We'll need it on the identifyONEspectrum() function
  precalc <- list()
  # df_metaespectres with precursor Mass
  precalc$df_metaspectra_noNAprecMass <- df_metaespectres[!is.na(df_metaespectres$REFprecursor_mass),c("idespectre","REFprecursor_mass")]
  # idspectra with no monoisotopic_molecular_weight
  precalc$REFidspectraNAprecMass <- df_metaespectres[is.na(df_metaespectres$REFprecursor_mass),"idespectre"]
  #df_metametabolite with monoisotopic_molecular_weight
  precalc$df_metametabolite_noNAneutralMass <- df_metametabolits[!is.na(df_metametabolits$REFmonoisotopic_molecular_weight), c("idmetabolit","REFmonoisotopic_molecular_weight")]
  REFidmetaboliteNAmonoIW <- df_metametabolits[is.na(df_metametabolits$REFmonoisotopic_molecular_weight),"idmetabolit"]
  #idspectra with NA monoisotopic_molecular_weight
  precalc$REFidspectraNAmonoIW <- df_EspectreMetabolit[df_EspectreMetabolit$idmetabolit%in%REFidmetaboliteNAmonoIW,"idespectre"]
  #in case we use only top N REF masses argument, preselect top N most intense fragments of every REF spectra
  if(!is.null(topMassesNum)){
    precalc$REFmaxIntMZ<-lapply(listfragOfInterest$espectre , function(x) x["mass-charge",order(x[2,],decreasing = TRUE)[seq_len(min(topMassesNum,ncol(x)))]])
  }
  
  #identify every spectra using lapply
  result <- pblapply(seq_along(UNKdata$Spectra$idUNKSpectra), function(x)
    { 
    UNKidspect <- UNKdata$Spectra$idUNKSpectra[x]
    identifyONEspectrum(UNKidspectra=UNKidspect,UNKspectra=UNKdata$Spectra$spectra[[x]],precalc=precalc,UNKmetadata= UNKdata$Metadata[UNKdata$Metadata$idUNKSpectra==UNKidspect,] ,listfragOfInterest=listfragOfInterest, topMassesNum=topMassesNum, ...)
    }, cl=cl)
  #remove spectra with no results
  result <- result[!is.na(result)]
  
  # list to dataframe
  finaldf <- do.call(rbind,result) 
  rownames(finaldf) <- NULL
  
    #ADD REF spectral metadata to the results
    finaldf <- cbind(finaldf,
                      df_metaespectres[match(finaldf$id_espectre,df_metaespectres$idespectre), names(df_metaespectres)!="idespectre"])
    
    # ADD REF metabolite metadata
    # First, check if idspectra is related to more than one metabolite. In that case replicate n times the row with that idspectra in order to assign a idmetabolite to each row 
    metabolitsperespectre <- lapply(finaldf$id_espectre, function(x) df_EspectreMetabolit[df_EspectreMetabolit$idespectre==x,"idmetabolit"])
    #copy n times the row where n is the number of idmetabolites per idspectra
    finaldf <- finaldf[rep(seq_along(metabolitsperespectre),vapply(metabolitsperespectre, length, FUN.VALUE = 1)),]
    #add idmetabolites
    finaldf$idmetabolit <- unlist(metabolitsperespectre)
    
    #finally, add REF spectral metadata to the results
    finaldf <- cbind(finaldf,
                      df_metametabolits[match(finaldf$idmetabolit,df_metametabolits$idmetabolit), names(df_metametabolits)!="idmetabolit"])
    
    #obtain score
    finaldf$score <- round((finaldf$cossim+log(1+7*finaldf$num_masses_coincidents/finaldf$UNKmassNum))/(1+log(8)),2)
    #round cosine similarity
    finaldf$cossim <- round(finaldf$cossim,2)
    
    # add UNK metadata
    finaldf <- cbind(finaldf,
                    UNK$Metadata[match(finaldf$idUNKSpectra,UNK$Metadata$idUNKSpectra), c("file","acquisitionNum","msLevel","polarity","retentionTime","collisionEnergy","precursorMZ","precursorCharge","precursorIntensity")])
  
  #rename some column names
  colnames(finaldf)[match(c("num_masses_coincidents","id_espectre","REFenergia","REFpolaritat","REFnaturalesa","BDespectreoriginal","REFprecursor_mass","REFionsource","idmetabolit","REFnommetabolit","REFmonoisotopic_molecular_weight","REFcasNum","idUNKSpectra","msLevel","polarity","retentionTime","collisionEnergy","precursorMZ","precursorCharge","precursorIntensity","REFinstrument","REFinchikey","REFformula","acquisitionNum","REFadduct"),colnames(finaldf))] <- c("MATCHmassNum", "REFidspectra","REF_CE","REFpolarity","REFnature","REF_DBspectra","REFprecMZ","REFionSource","REFidmetabolite","REFmetaboliteName","REFmonoisotMW","REFcasNum","UNKidspectra","UNKmsLevel","UNKpolarity","UNKrt","UNK_CE","UNKprecMZ","UNKprecCharge","UNKprecInt","REFinstrument","REFinchikey","REFformula","UNKacqNum","REFprecAdduct")
  #order & subset columns
  finaldf<-finaldf[order(finaldf$UNKacqNum, finaldf$UNKprecMZ, -finaldf$score, decreasing = FALSE), c("UNKacqNum","file","UNKprecMZ","REFmonoisotMW","assmdUNKAdduct","REFprecAdduct","REFprecMZ","score","cossim","REFmetaboliteName","REFformula","REFinchikey","REFcasNum","MATCHmassNum","UNKmassNum","UNK_CE","REF_CE","UNKpolarity","REFpolarity","UNKrt","UNKprecCharge","UNKprecInt","REFnature","REFinstrument","REFionSource","UNKmsLevel","UNKidspectra","REFidspectra","REFidmetabolite","REF_DBspectra")]
}  
  
identifyONEspectrum <- function(UNKidspectra, UNKspectra, UNKmetadata, candidatesNumber=20, listfragOfInterest, cosSimTh=NULL, massErrorAllowed=5, samePrecMass=FALSE, precalc, minCommonMassNum=2, topMassesNum=NULL, stick2Adducts=T, range=FALSE){
  #Variables defined in the parent scope are not visible, but globally-defined variables are visible. If the parent scope is the same as the global scope – those variables will be visible!. Arguments are immutable – if you change the value of an argument, what you are actually doing is creating a new variable and changing it
  #listfragOfInterest variable is passed as argument BUT NOT MODIFIED, IN ORDER TO AVOID REPLICATION (ITS TOO BIG)
  #STEP0: Obtain possible assumed adducts of UNK precursor mass 
  ##Given a UNK precursor mass, obtain its adducts masses. We'll use it to attach this info on results with every resulting REF metabolite. Also it is a necessary info when applyin stick2Adduct=T
  #if we have the ionization polarity we use it to filter adducts candidates 
  if(UNKmetadata$polarity%in%c(1,0)){
    polarityMatch <- as.integer(adducts$Ion_mode)==UNKmetadata$polarity
  }else {
    polarityMatch <- rep(T,length(adducts$Ion_mode))
  }
  #neutral mass candidates
  NeutralMassCand <- (UNKmetadata$precursorMZ-adducts$Mass)*abs(adducts$Charge)/adducts$Mult
  subsettedAdducts <- data.frame(AdductName=adducts[polarityMatch,c("Name")],neutralmass=NeutralMassCand[polarityMatch])
  
  #range according massErrorAllowed
  subsettedAdducts$minNeutralMass <- subsettedAdducts$neutralmass*(1-massErrorAllowed/10^6)
  subsettedAdducts$maxNeutralMass <- subsettedAdducts$neutralmass*(1+massErrorAllowed/10^6)
  #Begin and End positions on REF neutral MAss 
  subsettedAdducts$rangeBegin <- 1+findInterval(subsettedAdducts$minNeutralMass, precalc$df_metametabolite_noNAneutralMass$REFmonoisotopic_molecular_weight, left.open = T)
  subsettedAdducts$rangeEnd <- findInterval(subsettedAdducts$maxNeutralMass, precalc$df_metametabolite_noNAneutralMass$REFmonoisotopic_molecular_weight)
  #only accept adducts with mass  
  subsettedAdducts <- subsettedAdducts[subsettedAdducts$rangeBegin<=subsettedAdducts$rangeEnd,]
  
  #STEP1: FILTER REF DB according filtering arguments such as ionMassInc or samePrecMass
  lst_idspectraFiltered <- list()#to keep which REF idspectra to use on identifiaction. Every items correspons to a condition. Eventually we will only keep idspectra present on all the conditions (using Reduce())  
  
  #if(UNKidspectra==36) browser()
    
  #REF DB subset according samePrecMass condition
  if(samePrecMass){
    minREFprecMass <- UNKmetadata$precursorMZ*(1-2*massErrorAllowed/10^6)
    maxREFprecMass <- UNKmetadata$precursorMZ*(1+2*massErrorAllowed/10^6)
    lst_idspectraFiltered$samePrecMass <- c(precalc$REFidspectraNAprecMass, precalc$df_metaspectra_noNAprecMass[findInterval(minREFprecMass, precalc$df_metaspectra_noNAprecMass$REFprecursor_mass):findInterval(maxREFprecMass, precalc$df_metaspectra_noNAprecMass$REFprecursor_mass),"idespectre"])
  }
  #REF DB subset according stick2Adducts condition
  if(stick2Adducts){
    #only if there are REF metabolites wuth neutral masses that match UNK precursor adducts
    lst_idspectraFiltered$stick2Adducts <- vector()
    if(nrow(subsettedAdducts)>0){
      idmetabCand <- lapply(seq_len(nrow(subsettedAdducts)), function(x) precalc$df_metametabolite_noNAneutralMass[subsettedAdducts$rangeBegin[x]:subsettedAdducts$rangeEnd[x],"idmetabolit"])
      lst_idspectraFiltered$stick2Adducts <- df_EspectreMetabolit[df_EspectreMetabolit$idmetabolit%in%unlist(idmetabCand),"idespectre"]
    }
    # Also add idspectra of metabolites with no neutral mass (monoIW) information 
    lst_idspectraFiltered$stick2Adducts <- c(precalc$REFidspectraNAmonoIW,lst_idspectraFiltered$stick2Adducts)
  }
                     
  #REF spectra subset according selected idspectra on previous filters
  if(length(lst_idspectraFiltered)>0){
    idspectraFiltered <- Reduce(intersect, lst_idspectraFiltered)
  }else{
    idspectraFiltered <- listfragOfInterest$idespectre
  }
  #IN case no REF spectra remains, EXIT
  if(length(idspectraFiltered)==0) return(NA)
  
  #STEP2: Loof after REFspectra with UNKspectra fragments  
  if(!is.null(topMassesNum)){
    #select only REF spectra whose top topMassesNum most intense masses contains at list one of the top topMassesNum most intense UNK masses
    #order by mz UNK spectra
    unkMasses<-UNKspectra["mass-charge",order(UNKspectra[2,], decreasing = TRUE)[seq_len(min(topMassesNum,ncol(UNKspectra)))]]
    
    newTrimm<-vapply(precalc$REFmaxIntMZ[listfragOfInterest$idespectre%in%idspectraFiltered] , function(x) any(x %in% unkMasses), FUN.VALUE = TRUE)
    #update idspectraFiltered with the results
    idspectraFiltered <- listfragOfInterest$idespectre[listfragOfInterest$idespectre%in%idspectraFiltered][newTrimm]
  }else{
    unkMasses<-UNKspectra["mass-charge",]
  }
  
  #filter by fragment
  newTrimm<-vapply(listfragOfInterest$espectre[listfragOfInterest$idespectre%in%idspectraFiltered] , function(x) sum(x["mass-charge",] %in% unkMasses)>=minCommonMassNum, FUN.VALUE = TRUE)
  idspectraFiltered <- listfragOfInterest$idespectre[listfragOfInterest$idespectre%in%idspectraFiltered][newTrimm]
    
  #IN case no fragmented mass coincidence, EXIT
  if(length(idspectraFiltered)==0) return(NA)
  
  #STEP3: IDENTIFY
  #obtain UNK range 
  UNKrange <- round(c(UNKmetadata$scanWindowLowerLimit,UNKmetadata$scanWindowUpperLimit))
  #print(UNKidspectra)
  #if(UNKidspectra==206) browser()
  #if(UNKidspectra==322) browser()
  if(range){
    cosineSim<-vapply(idspectraFiltered, function(id_spct) {
      REFrange <- df_metaespectres[df_metaespectres$idespectre==id_spct, c("REFmzRangeIni","REFmzRangeEnd")]  
      
      REFspectra <- listfragOfInterest$espectre[listfragOfInterest$idespectre==id_spct][[1]]
      #subset spectra against oposite spectra's range
      REFspectra <- REFspectra[,(1+findInterval(UNKrange[1], REFspectra["mass-charge",], left.open = TRUE)):findInterval(UNKrange[2],REFspectra["mass-charge",]),drop=F]
      UNKspectraSubset <- UNKspectra[,(1+findInterval(REFrange[1],UNKspectra["mass-charge",], left.open = TRUE)):findInterval(REFrange[2],UNKspectra["mass-charge",]),drop=F]
      
      #find common masses. Twice to speed up further code
      massesrepetides_desc<-UNKspectraSubset["mass-charge",] %in% REFspectra["mass-charge",]
      massesrepetides_x<-REFspectra["mass-charge",] %in% UNKspectraSubset["mass-charge",]
      #fem les masses coincidents
      matriu<-rbind(UNKspectraSubset[2,massesrepetides_desc],REFspectra[2,massesrepetides_x])
      #afegim les no coincidents DE LA INTERSECCIO DELS RANGS, de cada matriu 
        if(!all(massesrepetides_desc)) {
          matriu1<-rbind(UNKspectraSubset[2,!massesrepetides_desc,drop=F],0)
          matriu<-cbind(matriu,matriu1)  
        }
        if(!all(massesrepetides_x)) {
          matriu2<-rbind(0,REFspectra[2,!massesrepetides_x])
          matriu<-cbind(matriu,matriu2)  
        }
      c(cosine(matriu[1,],matriu[2,]),sum(massesrepetides_desc))
    }, FUN.VALUE=c(3.2,3))
  }else{
    cosineSim<-vapply(listfragOfInterest$espectre[listfragOfInterest$idespectre%in%idspectraFiltered], function(x) {
      #Ho fem dues vegades per fer un rbind rapi despres
      massesrepetides_desc<-UNKspectra["mass-charge",] %in% x["mass-charge",]
      massesrepetides_x<-x["mass-charge",] %in% UNKspectra["mass-charge",]
      #fem les masses coincidents
      matriu<-rbind(UNKspectra[1,massesrepetides_desc],UNKspectra[2,massesrepetides_desc],x[2,massesrepetides_x])
      #afegim les no coincidents DE LA INTERSECCIO DELS RANGS, de cada matriu 
        if(!all(massesrepetides_desc)) {
          matriu1<-rbind(UNKspectra[,!massesrepetides_desc,drop=F],0)
          matriu<-cbind(matriu,matriu1)  
        }
        if(!all(massesrepetides_x)) {
          matriu2<-rbind(x[1,!massesrepetides_x,drop=F],0,x[2,!massesrepetides_x])
          matriu<-cbind(matriu,matriu2)  
        }
      c(cosine(matriu[2,],matriu[3,]),sum(massesrepetides_desc))
    }, FUN.VALUE=c(3.2,3))
  }
  
    #Remove results with cosinus similarity < cosSimTh
  if(!is.null(cosSimTh)) acceptedCos <- cosineSim[1,]>=cosSimTh
  # NO cossim acceptable, EXIT
  if(!any(acceptedCos)) return(NA)
  
  reslt <- data.frame(cossim=cosineSim[1,acceptedCos], num_masses_coincidents=cosineSim[2,acceptedCos],id_espectre=listfragOfInterest$idespectre[listfragOfInterest$idespectre%in%idspectraFiltered][acceptedCos])
  
  #order results by cossim
  reslt<-reslt[order(reslt$cossim, decreasing = TRUE),]
  #subset the n results with best cossim
  reslt<-reslt[seq_len(min(candidatesNumber, nrow(reslt))),]
 
  #Results' neutralMasses 
  ResltNeutralMasses <- vapply(df_EspectreMetabolit[match(reslt$id_espectre,df_EspectreMetabolit$idespectre),"idmetabolit"], function(x) df_metametabolits[df_metametabolits$idmetabolit==x,"REFmonoisotopic_molecular_weight"], FUN.VALUE=3.2)
  
  #print(UNKidspectra)
  #adducts related to those REF neutral masses
  reslt$assmdUNKAdduct <- vapply(ResltNeutralMasses, function(NM){
    #assumed Adduct with that neutral mass related
    addctCnd <- NA_character_
    if(!is.na(NM)) {
      addMatch <- subsettedAdducts$minNeutralMass<=NM & subsettedAdducts$maxNeutralMass>=NM
      if(any(addMatch)) addctCnd <- paste(subsettedAdducts[addMatch,"AdductName"])
    }
    return(addctCnd)
  }, FUN.VALUE = "nameadduct") 
  
  #add common variables
  reslt$idUNKSpectra <- UNKidspectra 
  reslt$UNKmassNum <- ncol(UNKspectra)
  
  return(reslt)
}
