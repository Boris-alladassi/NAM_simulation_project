nam_simulation <- function(architecture = "A", nfamilies = 10, nfounders = 100, nChrom = 10, totalNqtn = 200,
                           traitMean = 300, traitVar = 64, traitH2 = 1, traitName = "Trait1"){
  #Architecture: A = Additive, AE = Additive + Epistasis
  #traitH2 = Heritability of the trait
  #suggestion: maker nfounders = 10 x nfamilies
  
  ### Loading the local source of simplePHENOTYPES
  path_codes <- "./simplePHENOTYPES-master/R"
  codes <- list.files(path = path_codes, pattern = "*.R", full.names = TRUE)
  for (i in 1:length(codes)) {
    source(paste0(codes[i]))
  }
  require(AlphaSimR)
  source("./get.me.my.SNPs.in.hapmap.format.R")
  
  # AlphaSimR creates the founders below.
  founders <- runMacs(nInd = nfounders, nChr = nChrom, segSites = 2*totalNqtn, species = "GENERIC", inbred = T)
  SP <<- SimParam$new(founders)
  
  if(architecture == "A"){
    ###### Additive ---------------- ####
    ###################### Founders
    SP$addTraitA(nQtl = floor(totalNqtn/nChrom), mean = traitMean, var = traitVar, name = traitName)
    
    ### Source population of parents
    pop <- newPop(founders, simParam = SP)
    
    ### Pulling out the genotypic/genomic data. This Does not work unless you have added a trait.
    founder_qtl_map <- getQtlMap(trait = 1, simParam = SP)
    founder_qtls <- pullQtlGeno(pop, trait = 1, simParam = SP)
    
    #The function by Alex transforms the snp data into hapmap format needed for simplePHENOTYPES
    founder_qtls_hapmap <- get.me.my.SNPs.in.hapmap.format(these.SNPs = founder_qtls, 
                                                           this.physical.map = founder_qtl_map)
    
    ### Effs from AlphaSimR: simulate founders phenotype in simplePHENOTYPES using Alphasim effects
    ## Additive effects and list of QTNs
    add_effs <- (SP$traits[[1]])@addEff
    QTN_list <- list()
    QTN_list$add <- list(trait1 = founder_qtls_hapmap$snp)
    
    founder_gv_alpha <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = traitMean,,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      QTN_list = QTN_list, #list(add = list(trait1 =marker_names)),
      add_effect = list(trait1 = add_effs),
      model = "A",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 1000)
    
    founder_pheno_alpha <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      QTN_list = QTN_list,
      add_effect = list(trait1 = add_effs),
      model = "A",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 1000)
    
    founder_alpha <- pop
    
    ### DIRECTLY in simplePHENOTYPES: Now, let us simulate phenotype DIRECTLY in simplePHENOTYPES for the founders
    founder_gv <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = traitMean,
      rep = 1,
      add_QTN_num = totalNqtn,
      big_add_QTN_effect = sqrt(traitVar)/2,
      add_effect = -0.9,
      sim_method = "geometric",
      mode = "A",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 2000)
    
    founder_pheno <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      add_QTN_num = totalNqtn,
      big_add_QTN_effect = sqrt(traitVar)/2,
      add_effect = -0.9,
      sim_method = "geometric",
      mode = "A",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 2000)
    
    founder <- pop
    
    #### NAM RILs -------------------####
    #### Clustering using PAM
    pop_geno <- pullQtlGeno(pop)
    clust <- cluster::pam(pop_geno, k = nfamilies+1, nstart = 100)
    parents <- pop[clust$id.med]
    
    ## Creating a mating design of parent 1 (recurrent) crossed to parents 2-6.
    crossplan <- cbind(rep(1, nfamilies), c(2:(nfamilies+1)))
    
    ## Make the crosses and generate the hybrids
    pop_NAM <- makeCross(parents, crossplan, nProgeny = 250, simParam = SP)
    
    ### Self and obtain NAM RILs using single seed descent
    nam_rils <- pop_NAM
    for (g in 1:8) {
      nam_rils <- self(nam_rils, nProgeny = 1)
    }
    
    ### Again, let us use simplePHENOTYPES to simulate the pheno of the NAM RILs
    nam_qtl_map <- getQtlMap(trait = 1, simParam = SP)
    nam_qtls <- pullQtlGeno(pop=nam_rils, trait = 1, simParam = SP)
    nam_qtls_hapmap <- get.me.my.SNPs.in.hapmap.format(these.SNPs = nam_qtls, 
                                                       this.physical.map = nam_qtl_map)
    
    ### AlphaSimR Effects: simulate RILs phenotype in simplePHENOTYPES
    nam_rils_gv_alpha <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = traitMean,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      QTN_list = QTN_list,
      add_effect = list(trait1 = add_effs),
      model = "A",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 3000)
    
    nam_rils_pheno_alpha <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      QTN_list = QTN_list, #list(add = list(trait1 =marker_names)),
      add_effect = list(trait1 = add_effs),
      model = "A",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 3000)
    
    nam_rils_alpha <- nam_rils
    
    
    ### Directly in simplePHENOTYPES: simulate NAM RILs phenotype in simplePHENOTYPES
    nam_rils_gv <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = 300,
      rep = 1,
      add_QTN_num = 100,
      big_add_QTN_effect = sqrt(traitVar)/2,
      add_effect = -0.9,
      sim_method = "geometric",
      mode = "A",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 4000)
    
    nam_rils_pheno <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      add_QTN_num = totalNqtn,
      big_add_QTN_effect = sqrt(traitVar)/2,
      add_effect = -0.9,
      sim_method = "geometric",
      mode = "A",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 4000)
    
    colnames(founder_gv_alpha) <- c("Founder_ID", paste0(traitName, "_GV"))
    colnames(founder_pheno_alpha) <- c("Founder_ID", paste0(traitName, "_Pheno"))
    colnames(founder_gv) <- c("Founder_ID", paste0(traitName, "_GV"))
    colnames(founder_pheno) <- c("Founder_ID", paste0(traitName, "_Pheno"))
    
    colnames(nam_rils_gv_alpha) <- c("RIL_ID", paste0(traitName, "_GV"))
    colnames(nam_rils_pheno_alpha) <- c("RIL_ID", paste0(traitName, "_Pheno"))
    colnames(nam_rils_gv) <- c("RIL_ID", paste0(traitName, "_GV"))
    colnames(nam_rils_pheno) <- c("RIL_ID", paste0(traitName, "_Pheno"))
    
  }else if(architecture == "AE"){
    ###################### Additive + Epistasis ------------------------ ####
    SP$addTraitAE(nQtl = totalNqtn/nChrom, mean = traitMean, var = traitVar, name = traitName, relAA = 0.5) 
    
    ### Source population of parents
    pop <- newPop(founders, simParam = SP)
    
    ### Pulling out the genotypic/genomic data. This Does not work unless you have added a trait.
    founder_qtl_map <- getQtlMap(trait = 1, simParam = SP)
    founder_qtls <- pullQtlGeno(pop, trait = 1, simParam = SP)
    
    founder_qtls_hapmap <- get.me.my.SNPs.in.hapmap.format(these.SNPs = founder_qtls, 
                                                           this.physical.map = founder_qtl_map)
    
    
    ### Effs from AlphaSimR: simulate founders phenotype in simplePHENOTYPES using Alphasim effects
    ## Additive effects and list of QTNs
    add_effs <- (SP$traits[[1]])@addEff
    epi_effs <- (SP$traits[[1]])@epiEff[,3]
    
    epi_qtns2 <- c()
    for (i in 1:length(epi_effs)) {
      index <- c((SP$traits[[1]])@epiEff[i,1], (SP$traits[[1]])@epiEff[i,2])
      epi_qtns2 <- c(epi_qtns2, colnames(founder_qtls)[index])
    }
    QTN_list <- list()
    QTN_list$add[[1]] <- founder_qtls_hapmap$snp
    QTN_list$epi[[1]] <- epi_qtns2
    
    founder_gv_alpha <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = traitMean,,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      epi_QTN_num = length(epi_qtns2)/2,
      QTN_list = QTN_list,
      add_effect = list(trait1 = add_effs),
      epi_effect = list(trait1 = epi_effs),
      model = "AE",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 4000)
    
    founder_pheno_alpha <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      epi_QTN_num = length(epi_qtns2)/2,
      QTN_list = QTN_list,
      add_effect = list(trait1 = add_effs),
      epi_effect = list(trait1 = epi_effs),
      model = "AE",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 4000)
    
    founder_alpha <- pop
    
    ### DIRECTLY in simplePHENOTYPES: simulate founders phenotype DIRECTLY in simplePHENOTYPES
    founder_gv <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = traitMean,
      rep = 1,
      big_add_QTN_effect = 5,
      add_QTN_num = floor(totalNqtn*2/3),
      epi_QTN_num = floor(totalNqtn/3),
      add_effect = -0.9,
      epi_effect = -0.5,
      epi_interaction = 2,
      sim_method = "geometric",
      mode = "AE",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 4000)
    
    founder_pheno <- create_phenotypes(
      geno_obj = founder_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      big_add_QTN_effect = 5,
      add_QTN_num = floor(totalNqtn*2/3),
      epi_QTN_num = floor(totalNqtn/3),
      add_effect = -0.9,
      epi_effect = -0.5,
      epi_interaction = 2,
      sim_method = "geometric",
      mode = "E",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 4000)
    
    founder <- pop
    
    #### NAM RILs -------------------####
    #### Clustering using PAM
    pop_geno <- pullQtlGeno(pop)
    clust <- cluster::pam(pop_geno, k = nfamilies+1, nstart = 100)
    parents <- pop[clust$id.med]
    
    ## Creating a mating design of parent 1 (recurrent) crossed to parents 2-6.
    crossplan <- cbind(rep(1, nfamilies), c(2:(nfamilies+1)))
    
    ## Make the crosses and generate the hybrids
    pop_NAM <- makeCross(parents, crossplan, nProgeny = 250, simParam = SP)
    
    ### Self and obtain NAM RILs using single seed descent
    nam_rils <- pop_NAM
    for (g in 1:8) {
      nam_rils <- self(nam_rils, nProgeny = 1)
    }
    
    ### Again, let us use simplePHENOTYPES to simulate the pheno of the NAM RILs
    nam_qtl_map <- getQtlMap(trait = 1, simParam = SP)
    nam_qtls <- pullQtlGeno(pop=nam_rils, trait = 1, simParam = SP)
    nam_qtls_hapmap <- get.me.my.SNPs.in.hapmap.format(these.SNPs = nam_qtls, 
                                                       this.physical.map = nam_qtl_map)
    
    ### AlphaSimR Effects: simulate RILs phenotype in simplePHENOTYPES
    nam_rils_gv_alpha <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = traitMean,,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      epi_QTN_num = length(epi_qtns2)/2,
      QTN_list = QTN_list,
      add_effect = list(trait1 = add_effs),
      epi_effect = list(trait1 = epi_effs),
      model = "AE",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 4000)
    
    nam_rils_pheno_alpha <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      sim_method = "custom",
      add_QTN_num = totalNqtn,
      epi_QTN_num = length(epi_qtns2)/2,
      QTN_list = QTN_list,
      add_effect = list(trait1 = add_effs),
      epi_effect = list(trait1 = epi_effs),
      model = "AE",
      home_dir = getwd(),
      to_r = T,
      vary_QTN = FALSE,
      quiet = T,
      seed = 4000)
    
    nam_rils_alpha <- nam_rils
    
    ### Directly in simplePHENOTYPES: simulate NAM RILs phenotype in simplePHENOTYPES
    nam_rils_gv <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = 1,
      mean = traitMean,
      rep = 1,
      big_add_QTN_effect = 5,
      add_QTN_num = floor(totalNqtn*2/3),
      epi_QTN_num = floor(totalNqtn/3),
      add_effect = -0.9,
      epi_effect = -0.5,
      epi_interaction = 2,
      sim_method = "geometric",
      mode = "AE",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 4000)
    
    nam_rils_pheno <- create_phenotypes(
      geno_obj = nam_qtls_hapmap,
      ntraits = 1,
      h2 = traitH2,
      mean = traitMean,
      rep = 1,
      big_add_QTN_effect = 5,
      add_QTN_num = floor(totalNqtn*2/3),
      epi_QTN_num = floor(totalNqtn/3),
      add_effect = -0.9,
      epi_effect = -0.5,
      epi_interaction = 2,
      sim_method = "geometric",
      mode = "AE",
      home_dir = getwd(),
      to_r = TRUE,
      quiet = T,
      seed = 4000)
    
    colnames(founder_gv_alpha) <- c("Founder_ID", paste0(traitName, "_GV"))
    colnames(founder_pheno_alpha) <- c("Founder_ID", paste0(traitName, "_Pheno"))
    colnames(founder_gv) <- c("Founder_ID", paste0(traitName, "_GV"))
    colnames(founder_pheno) <- c("Founder_ID", paste0(traitName, "_Pheno"))
    
    colnames(nam_rils_gv_alpha) <- c("RIL_ID", paste0(traitName, "_GV"))
    colnames(nam_rils_pheno_alpha) <- c("RIL_ID", paste0(traitName, "_Pheno"))
    colnames(nam_rils_gv) <- c("RIL_ID", paste0(traitName, "_GV"))
    colnames(nam_rils_pheno) <- c("RIL_ID", paste0(traitName, "_Pheno"))
  }
  
  ########## Create the return list ---------------------- ####
  out_list <- list(founder_list <- list(founder_gv_alphasim = founder_gv_alpha, founder_pheno_alphasim = founder_pheno_alpha,
                                        founder_gv_simplePheno = founder_gv, founder_pheno_simplePheno = founder_pheno),
                   nam_rils_list <- list(nam_rils_gv_alphasim = nam_rils_gv_alpha, nam_rils_pheno_alphasim = nam_rils_pheno_alpha,
                                         nam_rils_gv_simplePheno = nam_rils_gv, nam_rils_pheno_simplePheno = nam_rils_pheno))
  return(out_list)
} ################################################# END OF THE FUNCTION ##########################################################
