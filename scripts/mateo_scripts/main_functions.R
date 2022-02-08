chrominfo <- readRDS("data/mateo_hg19_chrom_sizes.rds")

locationOfThisScript = function() # Function LocationOfThisScript returns the location of this .R script (may be needed to source other files in same dir)
{
    this.file = NULL
    # This file may be 'sourced'
    for (i in -(1:sys.nframe())) {
        if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
    }
    
    if (!is.null(this.file)) return(dirname(this.file))
    
    # But it may also be called from the command line
    cmd.args = commandArgs(trailingOnly = FALSE)
    cmd.args.trailing = commandArgs(trailingOnly = TRUE)
    cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
    res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)
    
    # If multiple --file arguments are given, R uses the last one
    res = tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))
    
    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
this_path<-locationOfThisScript()
source(paste(this_path,"helper_functions.R",sep="/"))

quantifySignatures<-function(sample_by_component,component_by_signature=NULL)
{
    if(is.null(component_by_signature))
    {
        component_by_signature<-readRDS(paste(this_path,"data/feat_sig_mat.rds",sep="/"))
    }
    signature_by_sample<-YAPSA::LCD(t(sample_by_component), YAPSA:::normalize_df_per_dim(component_by_signature,2))
    signature_by_sample<-normaliseMatrix(signature_by_sample)
    signature_by_sample
}

generateSignatures<-function(sample_by_component,nsig,seed=77777,nmfalg="brunet", cores=1)
{
    NMF::nmf(t(sample_by_component),nsig,seed=seed,nrun=1000,method=nmfalg,.opt = paste0("p", cores) )
}

chooseNumberSignatures<-function(sample_by_component, outfile="numSigs.pdf", min_sig=3, max_sig=12, iter=100, cores=1)
{

    nmfalg<-"brunet"
    seed<-77777
    
    estim.r <- NMF::nmfEstimateRank(t(sample_by_component), min_sig:max_sig,seed = seed,nrun=iter,
                               verbose=FALSE, method=nmfalg, .opt=list(shared.memory=FALSE, paste0("p", cores) ) )

    V.random <- NMF::randomize(t(sample_by_component))
    estim.r.random <- NMF::nmfEstimateRank(V.random, min_sig:max_sig, seed =seed,nrun=iter,
                                      verbose=FALSE, method=nmfalg, .opt=list(shared.memory=FALSE, paste0("p", cores) ) )
    
    p<-NMF::plot(estim.r,estim.r.random, 
            what = c("cophenetic", "dispersion","sparseness", "silhouette"),
            xname="Observed",yname="Randomised",main="")
    pdf(file=outfile, width=10, height=10 )
    p
    dev.off()

    return(p)

}

extractCopynumberFeatures<-function(CN_data, cores = 1)
{
    # # #get chromosome lengths
    # chrlen<-read.table(paste(this_path,"data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]
    # chrlen[which(chrlen$V1 == "chrX"),]$V1 <- "chr23"
    # chrlen[which(chrlen$V1 == "chrY"),]$V1 <- "chr24"
    # 
    # # #get centromere locations
    # gaps<-read.table(paste(this_path,"data/gap_hg19.txt",sep="/"),sep="\t",header=F,stringsAsFactors = F)
    # centromeres<-gaps[gaps[,8]=="centromere",]
  
    
    
    
    if(cores > 1) {
        require(foreach)
        doMC::registerDoMC(cores)

        temp_list = foreach::foreach(i=1:11) %dopar% {
            if(i == 1){
                list(segsize = getSegsize(CN_data) )
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data,chrominfo) )
            } else if (i == 3) {
                list(osCN = lation(CN_data,chrominfo) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(CN_data,chrominfo) )
            } else if (i == 5) {
                list(changepoint = getChangepointCN(CN_data) )
            } else if (i == 6){
                list(copynumber = getCN(CN_data) )
            } else if (i == 7){
               list(distToCentromere = getDistanceToCentromere(CN_data, chrominfo))
            } else if( i == 8){
              list(distToTelomere = getDistanceToTelomere(CN_data, chrominfo))
            } else if( i == 9){
              list(lohFractionPerArm = getLOHFractionPerArm(CN_data, chrominfo))
            } else if( i == 10){
              list(lohBAF = getLOH_BAF(CN_data))
            } else if (i == 11) {
              #list(sizeOfDiploidSeg = getSizeOfDiploidSeg(CN_data, chrominfo))
            }
        
        }
        unlist( temp_list, recursive = FALSE )
    } else {  
        
        segsize<-getSegsize(CN_data)
        bp10MB<-getBPnum(CN_data,chrominfo)
        osCN<-getOscilation(CN_data)
        bpchrarm<-getCentromereDistCounts(CN_data,chrominfo)
        changepoint<-getChangepointCN(CN_data)
        copynumber<-getCN(CN_data)
        distToCentromere <- getDistanceToCentromere(CN_data, chrominfo)
        distToTelomere <- getDistanceToTelomere(CN_data, chrominfo)
        lohFractionPerArm <- getLOHFractionPerArm(CN_data, chrominfo)
        lohBAF <- getLOH_BAF(CN_data)
        #sizeOfDiploidSeg <- getSizeOfDiploidSeg(CN_data, chrominfo)
        
        list(segsize=segsize,
             bp10MB=bp10MB,
             osCN=osCN,
             bpchrarm=bpchrarm,
             changepoint=changepoint,
             copynumber=copynumber, 
             distToCentromere=distToCentromere,
             distToTelomere = distToTelomere,
             lohFractionPerArm = lohFractionPerArm,
             lohBAF = lohBAF
             #sizeOfDiploidSeg = sizeOfDiploidSeg
             )
    }

}

fitMixtureModels<-function(CN_features, seed=77777, min_comp=2, max_comp=8, min_prior=0.01, model_selection="BIC",
                            nrep=5, niter=1000, cores = 1, featsToFit = seq(1, 11), em_tolerance = 0.001, classify_method = "hard")
{

    if(cores > 1) {
        require(foreach)

        doMC::registerDoMC(cores)

        temp_list = foreach(i=1:11) %dopar% {

            if(i == 1 & i %in% featsToFit ){
                cat("Fitting components of segsize...\n")
                dat<-as.numeric(CN_features[["segsize"]][,2])
                list( segsize = fitComponent(dat,seed=seed,model_selection=model_selection,
                                             min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp,em_tolerance = em_tolerance, classify_method = "hard") )
            
            } else if (i == 2 & i %in% featsToFit ) {
                cat("Fitting components of bp10mb...\n")
                dat<-as.numeric(CN_features[["bp10MB"]][,2])
                list( bp10MB = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                            min_prior=0.01,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp,em_tolerance = em_tolerance, classify_method = "hard") )

            } else if (i == 3 & i %in% featsToFit ) {
                cat("Fitting components of osCN...\n")
                dat<-as.numeric(CN_features[["osCN"]][,2])
                list( osCN = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                          min_prior=0.01,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp,em_tolerance = em_tolerance, classify_method = "auto") )
            
            } else if (i == 4 & i %in% featsToFit ) {
                cat("Fitting components of bpchrarm...\n")
                dat<-as.numeric(CN_features[["bpchrarm"]][,2])
                list( bpchrarm = fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                              min_prior=0.1,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp,em_tolerance = em_tolerance, classify_method = "hard") )
            
            } else if (i == 5 & i %in% featsToFit ) {
            
                cat("Fitting components of changepoint...\n")
                dat<-as.numeric(CN_features[["changepoint"]][,2])
                changepoint_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                             min_prior=0.001,niter=1500,nrep=nrep,min_comp=min_comp,max_comp=max_comp, em_tolerance = 0.01, classify_method = "hard")
              
            } else if (i == 6 & i %in% featsToFit) {
                cat("Fitting components of copynumber...  \n")
                dat<-as.numeric(CN_features[["copynumber"]][,2])
                list( copynumber = fitComponent(dat,seed=seed,model_selection=model_selection,
                                                nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.1,niter=niter,em_tolerance = 0.1, classify_method = "hard") )

            } else if (i == 7 & i %in% featsToFit) {
              cat("Fitting components of Dist TO centromere  \n")
              dat<-as.numeric(CN_features[["distToCentromere"]][,2])
              list( distToCentromere = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") )
              
            } else if(i == 8 & i %in% featsToFit ){
              cat("Fitting components of dist to Telomere  \n")
              dat<-as.numeric(CN_features[["distToTelomere"]][,2])
              list( distToTelomere = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") )
              
            } else if(i == 9 & i %in% featsToFit ){
             cat("Fitting components of Loh Fraction Per Arm  \n")
             dat<-as.numeric(CN_features[["lohFractionPerArm"]][,2])
             list( lohFractionPerArm = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") )
            } else if(i == 10 & i %in% featsToFit ){
              cat("Fitting components of LOH BAFs \n")
              dat<-as.numeric(CN_features[["lohBAF"]][,2])
              list( lohBAF = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") )
            } else if(i == 11 & i %in% featsToFit ){
              # cat("Fitting components of size Of Diploid Seg \n")
              # dat<-as.numeric(CN_features[["sizeOfDiploidSeg"]][,2])
              # list( sizeOfDiploidSeg = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") )
            }

        }
        unlist( temp_list, recursive = FALSE ) 
    } else {
        cat("Fitting components of segsize...\n")
        dat<-as.numeric(CN_features[["segsize"]][,2])
        segsize_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                 min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp,em_tolerance = em_tolerance, classify_method = "hard")

        cat("Fitting components of bp10mb\n")
        dat<-as.numeric(CN_features[["bp10MB"]][,2])
        bp10MB_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                                min_prior=0.01,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp,em_tolerance = em_tolerance, classify_method = "hard")

        cat("Fitting components of osCN\n")
        dat<-as.numeric(CN_features[["osCN"]][,2])
        osCN_mm<-fitComponent(dat,dist="pois",seed=seed,model_selection=model_selection,
                              min_prior=0.1,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=max_comp,em_tolerance = 0.0001, classify_method = "hard")

        cat("Fitting components of bpchrarm\n")
        dat<-as.numeric(CN_features[["bpchrarm"]][,2])
        bpchrarm_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                  min_prior=min_prior,niter=niter,nrep=nrep,min_comp=min_comp,max_comp=3,em_tolerance = em_tolerance, classify_method = "hard")

        cat("Fitting components of changepoint\n")
        dat<-as.numeric(CN_features[["changepoint"]][,2])
        changepoint_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                     min_prior=0.01,niter=1500,nrep=nrep,min_comp=min_comp,max_comp=max_comp, em_tolerance = 0.07, classify_method = "hard")

        cat("Fitting components of copynumber  \n")
        dat<-as.numeric(CN_features[["copynumber"]][,2])
        copynumber_mm<-fitComponent(dat,seed=seed,model_selection=model_selection,
                                nrep=nrep,min_comp=min_comp,max_comp=3,min_prior=min_prior,niter=niter,em_tolerance = 0.1, classify_method = "hard")

        cat("Fitting components of Dist To centromere  \n")
        dat<-as.numeric(CN_features[["distToCentromere"]][,2])
        distToCentromere = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=3,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard")
        
        cat("Fitting components of Dist To Telomere  \n")
        dat<-as.numeric(CN_features[["distToTelomere"]][,2])
        distToTelomere = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") 

        cat("Fitting components of Loh Fraction Per Arm  \n")
        dat<-as.numeric(CN_features[["lohFractionPerArm"]][,2])
        lohFractionPerArm = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=0.001,niter=niter, em_tolerance = 0.01, classify_method = "hard")
        
        cat("Fitting components of Loh Fraction Per Arm  \n")
        dat<-as.numeric(CN_features[["lohBAF"]][,2])
        lohBAF = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=2,max_comp=3,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") 
        
        # cat("Fitting components of size Of Diploid Seg \n")
        # dat<-as.numeric(CN_features[["sizeOfDiploidSeg"]][,2])
        # sizeOfDiploidSeg = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "auto")

        
        list(
             segsize=segsize_mm,
             bp10MB=bp10MB_mm,
             osCN=osCN_mm,
             bpchrarm=bpchrarm_mm,
             changepoint=changepoint_mm,
             copynumber=copynumber_mm,
             distToCentromere=distToCentromere,
             distToTelomere=distToTelomere,
             lohFractionPerArm = lohFractionPerArm,
             lohBAF = lohBAF
             #sizeOfDiploidSeg = sizeOfDiploidSeg
             )
    }
}

generateSampleByComponentMatrix<-function(CN_features, all_components=NULL, cores = 1, rowIter = 1000, subcores = 2)
{
    if(is.null(all_components))
    {
        all_components<-readRDS(paste(this_path,"data/component_parameters.rds",sep="/"))
    }

    if(cores > 1){
        require(foreach)

        feats = c( "segsize", "bp10MB", "osCN", "changepoint", "copynumber", "bpchrarm", "distToCentromere", "distToTelomere","lohBAF")
        doMC::registerDoMC(cores)

        full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
            calculateSumOfPosteriors(CN_features[[feat]],all_components[[feat]], 
                feat, rowIter = rowIter, cores = subcores)
        }
    } else {
        full_mat<-cbind(
        calculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
        calculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
        calculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
        calculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
        calculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"),
        calculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"),
        calculateSumOfPosteriors(CN_features[["distToCentromere"]],all_components[["distToCentromere"]],"distToCentromere"),
        calculateSumOfPosteriors(CN_features[["distToTelomere"]],all_components[["distToTelomere"]],"distToTelomere"),
        calculateSumOfPosteriors(CN_features[["lohFractionPerArm"]],all_components[["lohFractionPerArm"]],"lohFractionPerArm"),
        calculateSumOfPosteriors(CN_features[["lohBAF"]],all_components[["lohBAF"]],"lohBAF"),
        #calculateSumOfPosteriors(CN_features[["sizeOfDiploidSeg"]],all_components[["sizeOfDiploidSeg"]],"sizeOfDiploidSeg")
        )
    }

    rownames(full_mat)<-unique(CN_features[["segsize"]][,1])
    full_mat[is.na(full_mat)]<-0
    full_mat
}
