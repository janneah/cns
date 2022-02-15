chrominfo <- readRDS("data/mateo_hg19_chrom_sizes.rds")

# Function LocationOfThisScript returns the location of this
#.R script (may be needed to source other files in same dir)
location_of_this_script <- function() {
    this_file <- NULL
    # This file may be 'sourced'
    for (i in - (1:sys.nframe())) {
        if (identical(sys.function(i), base::source))
        this_file <- (normalizePath(sys.frame(i)$ofile))
    }
    if (!is.null(this_file)) return(dirname(this_file))
    # But it may also be called from the command line
    cmd_args <- commandArgs(trailingOnly = FALSE)
    cmd_args_trailing <- commandArgs(trailingOnly = TRUE)
    cmd_args <- cmd_args[seq.int(
        from = 1,
        length_out = length(cmd_args) - length(cmd_args_trailing))]
    res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmd_args)
    # If multiple --file arguments are given, R uses the last one
    res <- tail(res[res != ""], 1)
    if (0 < length(res)) return(dirname(res))
    # Both are not the case. Maybe we are in an R GUI?
    return(NULL)
}
this_path <- locationOfThisScript()
source(paste(this_path, "helper_functions.R", sep = "/"))

quantify_signatures <- function(
    sample_by_component,
    component_by_signature = NULL) {
    if (is.null(component_by_signature)) {
        component_by_signature <- readRDS(
            paste(this_path, "data/feat_sig_mat.rds", sep = "/"))
    }
    signature_by_sample <- YAPSA::LCD(
        t(sample_by_component),
        YAPSA:::normalize_df_per_dim(component_by_signature, 2))
    signature_by_sample <- normaliseMatrix(
        signature_by_sample)
    signature_by_sample
}

generate_signatures <- function(
    sample_by_component,
    nsig,
    seed = 77777,
    nmfalg = "brunet",
    cores = 1) {
    NMF::nmf(
        t(sample_by_component),
        nsig,
        seed = seed,
        nrun = 1000,
        method = nmfalg,
        .opt = paste0("p", cores))
}

choose_number_signatures <- function(
    sample_by_component,
    outfile = "numSigs.pdf",
    min_sig = 3,
    max_sig = 12,
    iter = 100,
    cores = 1) {

    nmfalg <- "brunet"
    seed <- 77777
    estim_r <- NMF::nmfEstimateRank(
        t(sample_by_component),
        min_sig:max_sig,
        seed = seed,
        nrun = iter,
        verbose = FALSE,
        method = nmfalg,
        .opt = list(shared.memory = FALSE, paste0("p", cores)))

    v_random <- NMF::randomize(t(sample_by_component))
    estim_r_random <- NMF::nmfEstimateRank(
        v_random,
        min_sig:max_sig,
        seed = seed,
        nrun = iter,
        verbose = FALSE,
        method = nmfalg,
        .opt = list(shared.memory = FALSE, paste0("p", cores)))
    p <- NMF::plot(estim_r, estim_r_random,
            what = c("cophenetic", "dispersion", "sparseness", "silhouette"),
            xname = "Observed", yname = "Randomised", main = "")
    pdf(file = outfile, width = 10, height = 10)
    p
    dev.off()

    return(p)

}

extract_copynumber_features <- function(cn_data, cores = 1) {
    # # #get chromosome lengths
    # chrlen<-read.table(
        # paste(this_path,"data/hg19.chrom.sizes.txt", sep="/"),
        # sep = "\t", stringsAsFactors = F)[1:24,]
    # chrlen[which(chrlen$V1 == "chrX"),]$V1 <- "chr23" # nolint
    # chrlen[which(chrlen$V1 == "chrY"),]$V1 <- "chr24" # nolint
    # # #get centromere locations
    # gaps<-read.table(
        # paste(this_path, "data/gap_hg19.txt",
        # sep = "/"),
        # sep="\t",
        # header=F,
        # stringsAsFactors = F)
    # centromeres <- gaps[gaps[,8] == "centromere", ] # nolint
    if (cores > 1) {
        require(foreach)
        doMC::registerDoMC(cores)

        temp_list <- foreach::foreach(i = 1:11) %dopar% {
            if (i == 1) {
                list(segsize = getSegsize(
                    cn_data))
            } else if (i == 2) {
                list(bp10mb = getBPnum(
                    cn_data, chrominfo))
            } else if (i == 3) {
                list(oscn = lation(
                    cn_data, chrominfo))
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(
                    cn_data, chrominfo))
            } else if (i == 5) {
                list(changepoint = getChangepointCN(
                    cn_data))
            } else if (i == 6) {
                list(copynumber = getCN(
                    cn_data))
            } else if (i == 7) {
               list(dist_to_centromere = getDistanceToCentromere(
                   cn_data, chrominfo))
            } else if (i == 8) {
              list(dist_to_telomere = getDistanceToTelomere(
                  cn_data, chrominfo))
            } else if (i == 9) {
              list(loh_fraction_per_arm = getLOHFractionPerArm(
                  cn_data, chrominfo))
            } else if (i == 10) {
              list(lohBAF = getLOH_BAF(
                  cn_data))
            } else if (i == 11) {
              #list(sizeOfDiploidSeg = getSizeOfDiploidSeg(cn_data, chrominfo)) # nolint
            }
        }
        unlist(temp_list, recursive = FALSE)
    } else {
        segsize <- getSegsize(
            cn_data)
        bp10mb <- getBPnum(
            cn_data, chrominfo)
        oscn <- getOscilation(
            cn_data)
        bpchrarm <- getCentromereDistCounts(
            cn_data, chrominfo)
        changepoint <- getChangepointCN(
            cn_data)
        copynumber <- getCN(
            cn_data)
        dist_to_centromere <- getDistanceToCentromere(
            cn_data, chrominfo)
        dist_to_telomere <- getDistanceToTelomere(
            cn_data, chrominfo)
        loh_fraction_per_arm <- getLOHFractionPerArm(
            cn_data, chrominfo)
        loh_baf <- getLOH_BAF(
            cn_data)
        #size_of_diploid_seg <- getSizeOfDiploidSeg(cn_data, chrominfo) # nolint
        list(segsize = segsize,
             bp10mb = bp10mb,
             oscn = oscn,
             bpchrarm = bpchrarm,
             changepoint = changepoint,
             copynumber = copynumber,
             dist_to_centromere = dist_to_centromere,
             dist_to_telomere = dist_to_telomere,
             loh_fraction_per_arm = loh_fraction_per_arm,
             loh_baf = loh_baf
             #size_of_diploid_seg = size_of_diploid_seg # nolint
             )
    }

}

fit_mixture_models <- function(
    cn_features, seed=77777, min_comp=2, max_comp=8,
    min_prior=0.01, model_selection="BIC", nrep=5,
    niter=1000, cores = 1, feats_to_fit = seq(1, 11),
    em_tolerance = 0.001, classify_method = "hard") {
    if (cores > 1) {
        require(foreach)

        doMC::registerDoMC(cores)

        temp_list <- foreach(i = 1:11) %dopar% {

            if (i == 1 & i %in% featsToFit) {
                cat("Fitting components of segsize...\n")
                dat <- as.numeric(cn_features[["segsize"]][, 2])
                list(segsize = fitComponent(
                    dat,
                    seed = seed,
                    model_selection = model_selection,
                    min_prior = min_prior,
                    niter = niter,
                    nrep = nrep,
                    min_comp = min_comp,
                    max_comp = max_comp,
                    em_tolerance = em_tolerance,
                    classify_method = "hard"))
            } else if (i == 2 & i %in% featsToFit) {
                cat("Fitting components of bp10mb...\n")
                dat <- as.numeric(cn_features[["bp10MB"]][, 2])
                list(bp10MB = fitComponent(
                    dat,
                    dist = "pois",
                    seed = seed,
                    model_selection = model_selection,
                    min_prior = 0.01,
                    niter = niter,
                    nrep = nrep,
                    min_comp = min_comp,
                    max_comp = max_comp,
                    em_tolerance = em_tolerance,
                    classify_method = "hard"))

            } else if (i == 3 & i %in% featsToFit) {
                cat("Fitting components of osCN...\n")
                dat <- as.numeric(cn_features[["osCN"]][, 2])
                list(osCN = fitComponent(
                    dat,
                    dist = "pois",
                    seed = seed,
                    model_selection = model_selection,
                    min_prior = 0.01,
                    niter = niter,
                    nrep = nrep,
                    min_comp = min_comp,
                    max_comp = max_comp,
                    em_tolerance = em_tolerance,
                    classify_method = "auto"))
            } else if (i == 4 & i %in% featsToFit) {
                cat("Fitting components of bpchrarm...\n")
                dat <- as.numeric(cn_features[["bpchrarm"]][, 2])
                list(bpchrarm = fitComponent(
                    dat,
                    dist = "pois",
                    seed = seed,
                    model_selection = model_selection,
                    min_prior = 0.1,
                    niter = niter,
                    nrep = nrep,
                    min_comp = min_comp,
                    max_comp = max_comp,
                    em_tolerance = em_tolerance,
                    classify_method = "hard"))
            } else if (i == 5 & i %in% featsToFit) {
                cat("Fitting components of changepoint...\n")
                dat <- as.numeric(cn_features[["changepoint"]][, 2])
                changepoint_mm <- fitComponent(
                    dat, seed = seed, model_selection = model_selection,
                    min_prior = 0.001, niter = 1500, nrep = nrep,
                    min_comp = min_comp, max_comp = max_comp,
                    em_tolerance = 0.01, classify_method = "hard")
            } else if (i == 6 & i %in% featsToFit) {
                cat("Fitting components of copynumber...  \n")
                dat <- as.numeric(cn_features[["copynumber"]][, 2])
                list(copynumber = fitComponent(
                    dat, seed = seed, model_selection = model_selection,
                    nrep = nrep, min_comp = min_comp, max_comp = max_comp,
                    min_prior = 0.1, niter = niter, em_tolerance = 0.1,
                    classify_method = "hard"))

            } else if (i == 7 & i %in% featsToFit) {
              cat("Fitting components of Dist TO centromere  \n")
              dat <- as.numeric(cn_features[["distToCentromere"]][, 2])
              list(distToCentromere = fitComponent(
                  dat, seed = seed, model_selection = model_selection,
                  nrep = nrep, min_comp = min_comp, max_comp = max_comp,
                  min_prior = min_prior, niter = niter,
                  em_tolerance = em_tolerance, classify_method = "hard"))
            } else if (i == 8 & i %in% featsToFit) {
              cat("Fitting components of dist to Telomere  \n")
              dat <- as.numeric(cn_features[["distToTelomere"]][, 2])
              list(distToTelomere = fitComponent(
                  dat, seed = seed, model_selection = model_selection,
                  nrep = nrep, min_comp = min_comp, max_comp = max_comp,
                  min_prior = min_prior, niter = niter,
                  em_tolerance = em_tolerance, classify_method = "hard"))
            } else if (i == 9 & i %in% featsToFit) {
             cat("Fitting components of Loh Fraction Per Arm  \n")
             dat <- as.numeric(cn_features[["lohFractionPerArm"]][, 2])
             list(lohFractionPerArm = fitComponent(
                 dat, seed = seed, model_selection = model_selection,
                 nrep = nrep, min_comp = min_comp, max_comp = max_comp,
                 min_prior = min_prior, niter = niter,
                 em_tolerance = em_tolerance, classify_method = "hard"))
            } else if (i == 10 & i %in% featsToFit) {
              cat("Fitting components of LOH BAFs \n")
              dat <- as.numeric(cn_features[["lohBAF"]][, 2])
              list(lohBAF = fitComponent(
                  dat, seed = seed, model_selection = model_selection,
                  nrep = nrep, min_comp = min_comp, max_comp = max_comp,
                  min_prior = min_prior, niter = niter,
                  em_tolerance = em_tolerance, classify_method = "hard"))
            } else if (i == 11 & i %in% featsToFit) {
              # cat("Fitting components of size Of Diploid Seg \n") # nolint
              # dat<-as.numeric(cn_features[["sizeOfDiploidSeg"]][,2]) # nolint
              # list( sizeOfDiploidSeg = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "hard") ) # nolint
            }

        }
        unlist(temp_list, recursive = FALSE)
    } else {
        cat("Fitting components of segsize...\n")
        dat <- as.numeric(cn_features[["segsize"]][, 2])
        segsize_mm <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            min_prior = min_prior, niter = niter, nrep = nrep,
            min_comp = min_comp, max_comp = max_comp,
            em_tolerance = em_tolerance, classify_method = "hard")

        cat("Fitting components of bp10mb\n")
        dat <- as.numeric(cn_features[["bp10MB"]][, 2])
        bp10mb_mm <- fitComponent(
            dat, dist = "pois", seed = seed, model_selection = model_selection,
            min_prior = 0.01, niter = niter, nrep = nrep, min_comp = min_comp,
            max_comp = max_comp,
            em_tolerance = em_tolerance, classify_method = "hard")

        cat("Fitting components of osCN\n")
        dat <- as.numeric(cn_features[["osCN"]][, 2])
        os_cn_mm <- fitComponent(
            dat, dist = "pois", seed = seed, model_selection = model_selection,
            min_prior = 0.1, niter = niter, nrep = nrep,
            min_comp = min_comp, max_comp = max_comp,
            em_tolerance = 0.0001, classify_method = "hard")

        cat("Fitting components of bpchrarm\n")
        dat <- as.numeric(cn_features[["bpchrarm"]][, 2])
        bpchrarm_mm <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            min_prior = min_prior, niter = niter,
            nrep = nrep, min_comp = min_comp, max_comp = 3,
            em_tolerance = em_tolerance, classify_method = "hard")

        cat("Fitting components of changepoint\n")
        dat <- as.numeric(cn_features[["changepoint"]][, 2])
        changepoint_mm <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            min_prior = 0.01, niter = 1500, nrep = nrep,
            min_comp = min_comp, max_comp = max_comp,
            em_tolerance = 0.07, classify_method = "hard")

        cat("Fitting components of copynumber  \n")
        dat <- as.numeric(cn_features[["copynumber"]][, 2])
        copynumber_mm <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            nrep = nrep, min_comp = min_comp, max_comp = 3,
            min_prior = min_prior, niter = niter,
            em_tolerance = 0.1, classify_method = "hard")

        cat("Fitting components of Dist To centromere  \n")
        dat <- as.numeric(cn_features[["distToCentromere"]][, 2])
        dist_to_centromere <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            nrep = nrep, min_comp = min_comp, max_comp = 3,
            min_prior = min_prior, niter = niter,
            em_tolerance = em_tolerance, classify_method = "hard")
        cat("Fitting components of Dist To Telomere  \n")
        dat <- as.numeric(cn_features[["distToTelomere"]][, 2])
        dist_to_telomere <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            nrep = nrep, min_comp = min_comp, max_comp = max_comp,
            min_prior = min_prior, niter = niter,
            em_tolerance = em_tolerance, classify_method = "hard")

        cat("Fitting components of Loh Fraction Per Arm  \n")
        dat <- as.numeric(cn_features[["lohFractionPerArm"]][, 2])
        lohfraction_per_arm <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            nrep = nrep, min_comp = min_comp, max_comp = max_comp,
            min_prior = 0.001, niter = niter,
            em_tolerance = 0.01, classify_method = "hard")
        cat("Fitting components of Loh Fraction Per Arm  \n")
        dat <- as.numeric(cn_features[["lohBAF"]][, 2])
        loh_baf <- fitComponent(
            dat, seed = seed, model_selection = model_selection,
            nrep = nrep, min_comp = 2, max_comp = 3, min_prior = min_prior,
            niter = niter, em_tolerance = em_tolerance,
            classify_method = "hard")
        # cat("Fitting components of size Of Diploid Seg \n") # nolint
        # dat<-as.numeric(cn_features[["sizeOfDiploidSeg"]][,2]) # nolint
        # sizeOfDiploidSeg = fitComponent(dat,seed=seed,model_selection=model_selection,nrep=nrep,min_comp=min_comp,max_comp=max_comp,min_prior=min_prior,niter=niter, em_tolerance = em_tolerance, classify_method = "auto") # nolint
        list(
             segsize = segsize_mm,
             bp10MB = bp10MB_mm,
             osCN = osCN_mm,
             bpchrarm = bpchrarm_mm,
             changepoint = changepoint_mm,
             copynumber = copynumber_mm,
             distToCentromere = distToCentromere,
             distToTelomere = distToTelomere,
             lohFractionPerArm = lohFractionPerArm,
             lohBAF = lohBAF
             #sizeOfDiploidSeg = sizeOfDiploidSeg # nolint
             )
    }
}

generate_s_by_c_matrix <- function(
    cn_features, all_components=NULL,
    cores = 1, row_iter = 1000, subcores = 2) {
    if (is.null(all_components)) {
        all_components <- readRDS(
            paste(this_path,
            "data/component_parameters.rds",
            sep = "/"))
    }

    if (cores > 1) {
        require(foreach)

        feats <- c(
            "segsize", "bp10MB", "osCN", "changepoint",
            "copynumber", "bpchrarm", "distToCentromere",
            "distToTelomere", "lohBAF")
        doMC::registerDoMC(cores)

        full_mat <- foreach(feat = feats, .combine = cbind) %dopar% {
            calculateSumOfPosteriors(
                cn_features[[feat]], all_components[[feat]],
                feat, rowIter = rowIter, cores = subcores)
        }
    } else {
        full_mat <- cbind(
        calculateSumOfPosteriors(
            cn_features[["segsize"]], all_components[["segsize"]], "segsize"),
        calculateSumOfPosteriors(
            cn_features[["bp10MB"]], all_components[["bp10MB"]], "bp10MB"),
        calculateSumOfPosteriors(
            cn_features[["osCN"]], all_components[["osCN"]], "osCN"),
        calculateSumOfPosteriors(
            cn_features[["changepoint"]], all_components[["changepoint"]],
            "changepoint"),
        calculateSumOfPosteriors(
            cn_features[["copynumber"]], all_components[["copynumber"]],
            "copynumber"),
        calculateSumOfPosteriors(
            cn_features[["bpchrarm"]], all_components[["bpchrarm"]],
            "bpchrarm"),
        calculateSumOfPosteriors(
            cn_features[["distToCentromere"]],
            all_components[["distToCentromere"]],
            "distToCentromere"),
        calculateSumOfPosteriors(
            cn_features[["distToTelomere"]],
            all_components[["distToTelomere"]],
            "distToTelomere"),
        calculateSumOfPosteriors(
            cn_features[["lohFractionPerArm"]],
            all_components[["lohFractionPerArm"]],
            "lohFractionPerArm"),
        calculateSumOfPosteriors(
            cn_features[["lohBAF"]],
            all_components[["lohBAF"]],
            "lohBAF"),
        #calculateSumOfPosteriors(CN_features[["sizeOfDiploidSeg"]],all_components[["sizeOfDiploidSeg"]],"sizeOfDiploidSeg") # nolint
        )
    }

    rownames(full_mat) <- unique(cn_features[["segsize"]][, 1])
    full_mat[is.na(full_mat)] <- 0
    full_mat
}
