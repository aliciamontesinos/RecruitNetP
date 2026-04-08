###ALL INTERNAL FUNCTIONS

remove_no_cover_UNI <- function(int_data=NULL, cover_data=NULL) {

  # Find species present in RN but that lack data on cover.
  cover_list <- sort(unique(cover_data$Canopy))
  RN_list <- sort(unique(c(int_data$Canopy, int_data$Recruit)))
  lack_cover <- setdiff(RN_list, cover_list)

  # Remove species lacking cover from RN
  if (length(lack_cover) == 0) {
    df <- int_data
  } else {
    df <- int_data[-which(int_data$Recruit %in% lack_cover), ]
  }

  if (length(lack_cover) == 0) {
    df
  } else {
    if(length(which(df$Canopy %in% lack_cover))==0) {
      df
    } else {
      df <- df[-which(df$Canopy %in% lack_cover), ]
    }
  }
  return(df)
}

# -------------------------------------------------------

# aggr_RN_UNI

aggr_RN_UNI <- function(int_data) {

  if(length(int_data[int_data$Canopy=="Open","Frequency"])==0){

    int_data<-rbind(int_data, int_data[nrow(int_data), ])

    int_data[nrow(int_data),"Canopy"]<-"Open"
    int_data[nrow(int_data),"Standardized_Canopy"]<-"Open"
    int_data[nrow(int_data),"Family_Canopy"]<-"Open"
    int_data[nrow(int_data),"Frequency"]<-0
    int_data[nrow(int_data),"LifeHabit_Canopy"]<-NA

  }else{

    int_data<-int_data}

  # Sum the number of recruits per interaction across plots
  RN <- stats::aggregate(Frequency ~ Canopy*Recruit, data = int_data, FUN = sum)
  colnames(RN) <- c("Canopy", "Recruit", "Fcr")
  RN$Icr <- stats::aggregate(Frequency ~ Canopy*Recruit, data=int_data, FUN = NROW)[[3]]
  RN$Pcr <- ifelse(RN$Icr==0,0,1)
  RN$Canopy <- gsub("[[:space:]]", "_", RN$Canopy)
  RN$Recruit <- gsub("[[:space:]]", "_", RN$Recruit)

  RN[RN$Fcr==0,"Icr"]<-0
  RN[RN$Fcr==0,"Pcr"]<-0

  # Incorporate the unobserved interactions
  species_list <- unique(c(RN$Canopy, RN$Recruit))
  df <- expand.grid(Canopy = species_list, Recruit = species_list)
  RN <- merge(df, RN, all = TRUE)
  RN[is.na(RN)] <- 0
  RN$Canopy <- as.character(RN$Canopy)
  RN$Recruit <- as.character(RN$Recruit)
  # RN <- RN[which(RN$Recruit!="Open"),] # Remove Open from the Recruit species
  return(RN)

}

# -------------------------------------------------------

# aggr_cover_UNI

aggr_cover_UNI <- function(cover_data) {
  # Calculations
  cover_data$Ac <- (cover_data$Cover / 100) * cover_data$Sampled_distance_or_area
  cover <- stats::aggregate(Ac ~ Canopy, data = cover_data, sum)
  return(cover)
}

# -------------------------------------------------------

# comm_to_RN_UNI

comm_to_RN_UNI <- function(int_data, cover_data) {

  #add a row of canopy Open if it does not exist

  if(length(int_data[int_data$Canopy=="Open","Frequency"])==0){

    int_data<-rbind(int_data, int_data[nrow(int_data), ])

    int_data[nrow(int_data),"Canopy"]<-"Open"
    int_data[nrow(int_data),"Standardized_Canopy"]<-"Open"
    int_data[nrow(int_data),"Family_Canopy"]<-"Open"
    int_data[nrow(int_data),"Frequency"]<-0
    int_data[nrow(int_data),"LifeHabit_Canopy"]<-NA


  } else{

    int_data<-int_data}

  # Aggregate data sets across plots.
  int_df <- aggr_RN_UNI(int_data)
  cover_df <- aggr_cover_UNI(cover_data)

  # Find species present in RN but that lack data on cover.
  cover_list <- sort(unique(cover_df$Canopy))
  RN_list <- sort(unique(c(int_df$Canopy, int_df$Recruit)))

  lack_cover <- setdiff(RN_list, cover_list)

  # Remove species lacking cover from RN
  if (length(lack_cover) == 0) {
    RNc <- int_df
  } else {
    RNc <- int_df[-which(int_df$Recruit %in% lack_cover), ]
  }

  if (length(lack_cover) == 0) {
    RNc
  } else {
    if(length(which(RNc$Canopy %in% lack_cover))==0) {
      RNc
    } else {
      RNc <- RNc[-which(RNc$Canopy %in% lack_cover), ]
    }
  }

  # Add variables with the cover of the canopy (Ac) and recruit (Ar) species
  RNc$Ac <- RNc$Canopy
  RNc$Ar <- RNc$Recruit
  for (i in 1:dim(RNc)[1]) {
    RNc$Ac[i] <- as.numeric(replace(
      RNc$Ac[i],
      match(RN_list, RNc$Canopy[i]),
      cover_df$Ac[match(RNc$Canopy[i], cover_df$Canopy)]
    ))
  }

  for (i in 1:dim(RNc)[1]) {
    RNc$Ar[i] <- as.numeric(replace(
      RNc$Ar[i],
      match(RN_list, RNc$Recruit[i]),
      cover_df$Ac[match(RNc$Recruit[i], cover_df$Canopy)]
    ))
  }

  RNc <- utils::type.convert(RNc, as.is = TRUE)
  RNc <- RNc[order(RNc$Canopy, RNc$Recruit),]

  return(RNc)

}

# -------------------------------------------------------

# RN_to_matrix_UNI()

RN_to_matrix_UNI <- function(int_data=NULL, weight = "Fcr"){

  # Check column names
  if ("Canopy" %in% names(int_data) == FALSE) warning("ERROR: your data lacks a column named: Canopy")
  if ("Recruit" %in% names(int_data) == FALSE) warning("ERROR: your data lacks a column named: Recruit")

  data <- int_data
  # Formatting
  list_Canopy <- sort(unique(data$Canopy))
  list_Recruit <- sort(unique(data$Recruit))
  Num_canopy <- length(list_Canopy)
  Num_recruit <- length(list_Recruit)
  RNmat <- data[[weight]]
  dim(RNmat) <- c(Num_recruit, Num_canopy)
  colnames(RNmat) <- list_Canopy
  rownames(RNmat) <- list_Recruit
  return(RNmat)

}

# -------------------------------------------------------

# pre_associndex_UNISITE_UNI()

pre_associndex_UNISITE_UNI <- function(int_data = NULL) {
  # Formatting
  mydata <- int_data[,c("Canopy", "Recruit", "Fcr","Ac")]

  # Frequency in open
  Fcr_matrix_df <- as.data.frame(RN_to_matrix_UNI(mydata))
  Open_matrix <- rep(Fcr_matrix_df$Open, length(Fcr_matrix_df))
  dim(Open_matrix) <- c(length(Fcr_matrix_df), length(Fcr_matrix_df))
  mydata$Fro <- unlist(as.vector(Open_matrix))

  # Cover of Open
  cov_matrix_df <- as.data.frame(RN_to_matrix_UNI(mydata, "Ac"))
  Open_cov_matrix <- rep(cov_matrix_df$Open, length(cov_matrix_df))
  dim(Open_cov_matrix) <- c(length(cov_matrix_df), length(cov_matrix_df))
  mydata$Ao <- unlist(as.vector(Open_cov_matrix))




  # Remove absent interactions
  mydata <- mydata[mydata$Fcr+mydata$Fro > 0,]
  # Remove "Open" from the list of canopy species
  if(nrow(mydata[mydata$Canopy=="Open",])>0){
    mydata <- mydata[-which(mydata$Canopy=="Open"),]
  }else{
    mydata <- mydata
  }
  return(mydata)

}

# -------------------------------------------------------

# associndex_UNISITE_UNI()

associndex_UNISITE_UNI <- function(int_data = NULL, threshold_density=NULL) {

  if (!"Open" %in% int_data$Canopy) stop("tests cannot be conducted because your data does not contain a node named Open or it is spelled differently.")

  thr <- threshold_density

  # Assemble the data
  db_inter <- pre_associndex_UNISITE_UNI(int_data)

  # Incorporate density of recruitment (recruits/m2) under each canopy species and in open.
  db_inter$Dcr <- db_inter$Fcr/db_inter$Ac
  db_inter$Dro <- db_inter$Fro/db_inter$Ao

  # Retain the interactions with estimated density below the threshold.
  db_inter <- db_inter[which(db_inter$Dcr<thr & db_inter$Dro<thr), ]

  #Obtain the maximum recruitment density for each recruit under the canopy species or in open.
  db_inter$Max_Recr_Density <- pmax(db_inter$Dcr,db_inter$Dro)

  db_inter <- utils::type.convert(db_inter, as.is = TRUE)

  max_rd <- stats::aggregate(Max_Recr_Density ~ Recruit, data = db_inter, FUN = "max")

  # Add a variable max_Recr to each pair indicating the maximum recruitment density of the recruit species in the study site
  Recr_list <- sort(unique(c(db_inter$Recruit)))
  Dens_list <- sort(unique(max_rd$Recruit))
  lack_dens <- setdiff(Recr_list, Dens_list)

  db_inter$max_Recr <- db_inter$Recruit
  for (i in 1:(dim(db_inter)[1])) {
    db_inter$max_Recr[i] <- replace(
      db_inter$max_Recr[i],
      match(Recr_list, db_inter$max_Recr[i]),
      max_rd$Max_Recr_Density[match(db_inter$max_Recr[i], max_rd$Recruit)]
    )
  }

  db_inter <- utils::type.convert(db_inter, as.is = TRUE)

  # Calculate indices Ns, NintC, NintA and RII
  db_inter$Ns <- (db_inter$Dcr - db_inter$Dro)/db_inter$max_Recr
  db_inter$NintC <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dcr + db_inter$Dro)+abs(db_inter$Dcr-db_inter$Dro))
  db_inter$NintA <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dro) + abs(db_inter$Dcr-db_inter$Dro))
  db_inter$RII <- (db_inter$Dcr - db_inter$Dro)/(db_inter$Dcr + db_inter$Dro)

  removed <- names(db_inter) %in% c("Frequency", "Max_Recr_Density")
  db_inter <- db_inter[!removed]
  return(db_inter)

}


# -------------------------------------------------------

# int_significance_UNI()

int_significance_UNI <- function(int_data){

  if (!"Open" %in% int_data$Canopy) stop("tests cannot be conducted because your data does not contain a node named Open or it is spelled differently.")

  df <- pre_associndex_UNISITE_UNI(int_data)
  n_tests <- dim(df)[1]
  df$exp_p <- df$Ac/(df$Ac+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)

  # Testability through Binomial test

  df$Ftot <- df$Fcr+df$Fro

  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p

  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- stats::binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability

  # Binomial (or Chi square) Test Significance

  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Significance[i] <- stats::binom.test(df$Fcr[i], df$Fcr[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- stats::chisq.test(c(df$Fcr[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance

  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }
  df$Test_type <- Test_type
  #  if(length(unique(df$Test_type))>1) warning("Different tests were used for different canopy-recruit pairs. Check column Test_type")

  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fcr[i]/df$Ac[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Enhancing",
                         Effect_int[i] <- "Depressing")
           )
    )
  }

  df$Effect_int <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]

  return(df)
}

# -------------------------------------------------------

# bipartite_RNs_UNI()

bipartite_RNs_UNI <- function(int_data, effect_var = int_data$Effect_int, int_type = "Enhancing", frequency_var = int_data$Fcr) {

  # Subset the data set to retain only the indicated type of interactions
  type_df <- int_data
  type_df$typeCR <- frequency_var
  type_df$effect_var <- effect_var
  for(i in 1:dim(type_df)[1]) {
    type_df$typeCR[i] <- ifelse(type_df$effect_var[i] == int_type, type_df$typeCR[i], 0)
  }

  type_df <- type_df[which(type_df$typeCR >0),]
  type_df <- data.frame(cbind(type_df$Canopy, type_df$Recruit, type_df$typeCR))
  colnames(type_df) <- c("Canopy", "Recruit", "Fcr")
  type_df$Fcr <- as.numeric(type_df$Fcr)
  type_df$net_type <- rep(int_type, dim(type_df)[1])
  type_df <- type_df[, c(1,2,4,3)]

  # Make the interactions matrix
  bipartite_RN <- bipartite::frame2webs(type_df, varnames = c("Canopy", "Recruit", "net_type", "Fcr"))

  return(bipartite_RN[[1]])
}

# -------------------------------------------------------


# canopy_service_test_UNI()

canopy_service_test_UNI <- function(int_data){

  df <- pre_associndex_UNISITE_UNI(int_data)
  sp_Fc <- stats::aggregate(Fcr ~ Canopy, data = df, FUN = sum)
  sp_Ac <- stats::aggregate(Ac ~ Canopy, data = df, FUN = max)
  sp_Fro <- stats::aggregate(Fro ~ Canopy, data = df, FUN = sum)
  sp_Ao <- stats::aggregate(Ao ~ Canopy, data = df, FUN = max)
  n_tests <- dim(sp_Fc)[1]
  df <- data.frame(c(sp_Fc, sp_Ac, sp_Fro, sp_Ao))
  myvars <- names(df) %in% c("Canopy.1", "Canopy.2", "Canopy.3")
  df <- df[!myvars]
  colnames(df) <- c("Canopy", "Fc", "Ac", "Fro", "Ao")
  df$exp_p <- df$Ac/(df$Ac+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)

  # Testability through Binomial test

  df$Ftot <- df$Fc+df$Fro

  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p

  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- stats::binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability

  # Binomial (or Chi square) Test Significance

  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fc[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fc[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Significance[i] <- stats::binom.test(df$Fc[i], df$Fc[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- stats::chisq.test(c(df$Fc[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance

  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fc[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fc[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }
  df$Test_type <- Test_type

  if(length(unique(df$Test_type))>1) message("Different tests were used for different canopy-recruit pairs. Check column Test_type")

  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fc[i]/df$Ac[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Facilitative",
                         Effect_int[i] <- "Depressive")
           )
    )
  }

  df$Canopy_effect <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]
  return(df)
}

# -------------------------------------------------------

# recruitment_niche_test_UNI()

recruitment_niche_test_UNI <- function(int_data){

  df <- pre_associndex_UNISITE_UNI(int_data)
  sp_Fr <- stats::aggregate(Fcr ~ Recruit, data = df, FUN = sum)
  sp_Av <- stats::aggregate(Ac ~ Recruit, data = df, FUN = sum)
  sp_Fro <- stats::aggregate(Fro ~ Recruit, data = df, FUN = max)
  sp_Ao <- stats::aggregate(Ao ~ Recruit, data = df, FUN = max)
  n_tests <- dim(sp_Fr)[1]
  df <- data.frame(c(sp_Fr, sp_Av, sp_Fro, sp_Ao))
  myvars <- names(df) %in% c("Recruit.1", "Recruit.2", "Recruit.3")
  df <- df[!myvars]
  colnames(df) <- c("Recruit", "Fr", "Av", "Fro", "Ao")
  df$exp_p <- df$Av/(df$Av+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)

  # Testability through Binomial test

  df$Ftot <- df$Fr+df$Fro

  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p

  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- stats::binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability

  # Binomial (or Chi square) Test Significance

  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fr[i]+df$Fro[i])*(df$Av[i]/(df$Av[i]+df$Ao[i]))<=5 | (df$Fr[i]+df$Fro[i])*(df$Ao[i]/(df$Av[i]+df$Ao[i]))<=5),
           Significance[i] <- stats::binom.test(df$Fr[i], df$Fr[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- stats::chisq.test(c(df$Fr[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance

  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fr[i]+df$Fro[i])*(df$Av[i]/(df$Av[i]+df$Ao[i]))<=5 | (df$Fr[i]+df$Fro[i])*(df$Ao[i]/(df$Av[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }
  df$Test_type <- Test_type

  if(length(unique(df$Test_type))>1) warning("Different tests were used for different canopy-recruit pairs. Check column Test_type")

  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fr[i]/df$Av[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Facilitated",
                         Effect_int[i] <- "Depressed")
           )
    )
  }

  df$Veg_effect <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]
  return(df)
}

# -------------------------------------------------------

# node_degrees_UNI()

node_degrees_UNI <- function(int_data) {

  RNc <- int_data

  matrix_Fcr <- RN_to_matrix_UNI(RNc, weight = "Fcr")
  matrix_Pcr <- RN_to_matrix_UNI(RNc, weight = "Pcr")
  matrix_Ac <- RN_to_matrix_UNI(RNc, weight = "Ac")

  # Degree properties
  canopy_service_width <- colSums(matrix_Pcr)
  canopy_contribution <- colSums(matrix_Fcr)
  recruit_niche_width <- rowSums(matrix_Pcr)
  recruit_bank_abundance <- rowSums(matrix_Fcr)
  node_abund <- matrix_Ac[1,]

  # Specialization
  effective_canopy_service <- exp(-1*rowSums((t(matrix_Fcr)/canopy_contribution)*log(t(matrix_Fcr)/canopy_contribution), na.rm=TRUE))
  effective_recruit_niche <- exp(-1*rowSums((matrix_Fcr/recruit_bank_abundance)*log(matrix_Fcr/recruit_bank_abundance), na.rm=TRUE))

  node_deg <- data.frame(cbind(names(canopy_service_width), node_abund, canopy_service_width, canopy_contribution, effective_canopy_service, recruit_niche_width, recruit_bank_abundance, effective_recruit_niche))
  colnames(node_deg) <- c("Node", "Ac", "canopy_service_width", "canopy_contribution", "effective_canopy_service", "recruitment_niche_width", "recruit_bank_abundance", "effective_recruitment_niche")
  node_deg <- utils::type.convert(node_deg, as.is = TRUE)
  rownames(node_deg) <- NULL
  return(node_deg)
}

# -------------------------------------------------------

# partial_RNs_UNI()

partial_RNs_UNI <- function(int_data, k) {

  if (!"Plot" %in% names(int_data)) stop("your interactions data lacks a column named Plots. This function requires data assembled in plots.")

  # Prepare the data
  int_data$Pcr <- ifelse(int_data$Frequency==0,0,1)
  netRaw <- data.frame(cbind(int_data$Plot, int_data$Canopy, int_data$Recruit, int_data$Pcr))
  colnames(netRaw) <- c("Plot", "Canopy", "Recruit", "Pcr")
  netRaw <- transform(netRaw, Pcr = as.numeric(Pcr))
  nPlots <- length(unique(netRaw$Plot))

  if (nPlots == 1)
    stop("Data must be structured in multiple plots.")

  # Make the network of each plot
  plot_RNs <- c()
  for(i in 1:nPlots) {
    #plot_RNs[[i]] <- igraph::graph_from_data_frame(netRaw[netRaw$Plot == i, 2:4], directed = TRUE)
    plot_RNs[[i]] <- igraph::graph_from_data_frame(netRaw[netRaw$Plot == unique(netRaw$Plot)[i], 2:4], directed = TRUE)
  }

  # Combine the networks of n plots k times.
  union_RNs <- list()
  for(i in 1:nPlots) {

    #add simplify = FALSE
    union_RNs[[i]] <- replicate(k, do.call(igraph::union, sample(plot_RNs, i)),simplify = FALSE)
  }


  # Add the incidence of each interaction across plots as an edge weight property to each partial network.
  for(i in 1:nPlots){
    for(j in 1:k){
      Icr <- rowSums(do.call(cbind.data.frame, igraph::edge.attributes(union_RNs[[i]][[j]])), na.rm=TRUE)
      union_RNs[[i]][[j]]<-igraph::set_edge_attr(union_RNs[[i]][[j]], "Icr", value = Icr)
    }}

  return(union_RNs)
}

# -------------------------------------------------------

# RN_dims_UNI()

RN_dims_UNI <- function(int_data, P_int){
  #TEST
  if(is.null(P_int)) stop("P_int column not found in dataframe.")

  # FUNCTION
  df <- int_data
  n_nodes <- length(unique(c(df$Canopy, df$Recruit)))
  n_links <- sum(P_int)
  connectance <- n_links/(n_nodes^2 - n_nodes)

  out <- data.frame(c(n_nodes, n_links, connectance))
  colnames(out) <- c("Value")
  rownames(out) <- c("Num. Nodes", "Num. Links", "Connectance")
  return(out)
}

# -------------------------------------------------------

# node_topol_UNI()

node_topol_UNI <- function(int_data) {

  RN_igraph <- igraph::graph_from_adjacency_matrix(t(RN_to_matrix_UNI(int_data, weight = "Pcr")), mode = "directed")
  eigen_cent <- igraph::eigen_centrality(RN_igraph, directed=TRUE, scale=FALSE, options = list(which="LR"))$vector
  out_neigh <- igraph::neighborhood_size(RN_igraph, order=gorder(RN_igraph), mode="out", mindist=1)
  in_neigh <- igraph::neighborhood_size(RN_igraph, order=gorder(RN_igraph), mode="in", mindist=1)
  df <- data.frame(eigen_cent, out_neigh,in_neigh)
  df[, 1] <- round(df[, 1], digits = 4)
  colnames(df) <- c("Eigenvector centrality", "Extended canopy service", "Extended recruitment niche")
  return(df)
}

# -------------------------------------------------------

# funtopol_UNI()

funtopol_UNI <- function(int_data){

  if (!"Open" %in% int_data$Canopy) stop("your data does not contain a node named Open or it is spelled differently.")

  int_data <- int_data[which(int_data$Fcr!=0), c("Canopy", "Recruit")]
  g <- igraph::graph_from_data_frame(int_data, directed = TRUE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = FALSE)

  if (length(which(int_data$Canopy=="Open"))==0) {
    warning("Open is included as a node in the network, even though no recruits are associated with Open in this community")

    g <- igraph::add_vertices(g, 1, name = "Open")

  }


  NEdges <- igraph::gsize(g)
  NNodes <- igraph::gorder(g)
  SCCs <- igraph::components(g, mode = "strong")

  if(max(SCCs$csize)>1){

    numSCCs <- SCCs$no
    numNTSCCs <- sum(SCCs$csize > 1)
    coreSize <- max(SCCs$csize)
    SCC_memb <- SCCs$membership
    SCC_memb <- as.data.frame(SCC_memb)
    SCC_subgraphs <- igraph::decompose(g, mode = "strong") # Makes a subgraph of each SCC
    IDcore <- match(coreSize, SCCs$csize) # locates the position of the core in the list of SCCs
    MembersCore <- igraph::V(SCC_subgraphs[[IDcore]])$name # List of the species in the core
    IDOpen <- SCC_memb$SCC_memb[match("Open", row.names(SCC_memb))] # Locate the position of the "open" node in the list of SCCs
    outReachFromOpen <- names(igraph::subcomponent(g, "Open", "out")[-1]) # List of nodes reachable from the "open"
    outReachFromCore <- vector("list", coreSize) # List of nodes reachable from core nodes
    for (i in 1:coreSize) {
      outReachFromCore[[i]] <- igraph::subcomponent(g, MembersCore[i], mode = "out")
    }
    a <- unlist(outReachFromCore)
    a <- unique(names(a))
    MembersSatellites <- setdiff(a, MembersCore)
    MembersTransients <- setdiff(igraph::V(g)$name,c(MembersCore,MembersSatellites))
    MembersTransients <- MembersTransients[!MembersTransients == "Open"]
    MembersStrictTransients <- setdiff(MembersTransients, outReachFromOpen)
    MembersDdTransients <- setdiff(MembersTransients, MembersStrictTransients)
    numSat <- length(MembersSatellites)
    numTransAll <- length(MembersTransients)
    numDdTrans <- length(MembersDdTransients)
    numStrictTrans <- length(MembersStrictTransients)
    propCore <- coreSize/(NNodes - 1)
    propSat <- numSat/(NNodes - 1)
    propTrans <- numTransAll/(NNodes - 1)
    propStrTrans <- numStrictTrans/(NNodes - 1)
    propDdTrans <- numDdTrans/(NNodes - 1)
    persistence <- propCore + propSat

    # Function output

    df <- data.frame(
      c(NNodes,
        NEdges,
        numNTSCCs,
        coreSize,
        propCore,
        numSat,
        propSat,
        numDdTrans,
        propDdTrans,
        numStrictTrans,
        propStrTrans,
        persistence)
    )
    colnames(df) <- c("Value")
    rownames(df) <- c(
      "Num. nodes",
      "Num. edges",
      "Num. non-trivial SCCs",
      "Num. core species",
      "Prop. core species",
      "Num. satellite species",
      "Prop. satellite species",
      "Num. disturbance-dependent transients",
      "Prop. disturbance-dependent transients",
      "Num. strict transients",
      "Prop. strict transients",
      "Qualitative Persistence")

    classif <- list(
      MembersSatellites,
      MembersCore,
      MembersStrictTransients,
      MembersDdTransients)
    classif <- stats::setNames(classif,
                               c("Satellites",
                                 "Core",
                                 "Strict_transients",
                                 "Disturbance_dependent_transients")
    )

    df0 <- classif
    df_Sat <- data.frame(df0$Satellites, rep("Satellite", length(df0$Satellites)))
    colnames(df_Sat) <- c("id", "group")
    df_Core <- data.frame(df0$Core, rep("Core", length(df0$Core)))
    colnames(df_Core) <- c("id", "group")
    df_Str <- data.frame(df0$Strict_transients, rep("Strict_transients", length(df0$Strict_transients)))
    colnames(df_Str) <- c("id", "group")
    df_Ddtr <- data.frame(df0$Disturbance_dependent_transients, rep("Disturbance_dependent_transients", length(df0$Disturbance_dependent_transients)))
    colnames(df_Ddtr) <- c("id", "group")
    df0 <- rbind(df_Sat, df_Core, df_Str, df_Ddtr)
    df0 <- df0[order(df0$id),]

    outputs <- list("Descriptors" = df, "Functional_classification" = df0)


  }else{

    warning("This network does not have a Core, and thus the species membership to different roles can not be defined")

    coreSize<-0

    df <- data.frame(
      c(NNodes,
        NEdges,
        coreSize))


    colnames(df) <- c("Value")
    rownames(df) <- c(
      "Num. nodes",
      "Num. edges",
      "Num. core species")

    df0<-list(
      Satellites = character(0),
      Core = character(0),
      Strict_transients = character(0),
      Disturbance_dependent_transients = character(0)
    )

    outputs <-list("Descriptors" = df, "Functional_classification" = df0)

  }

  return(outputs)

}

# -------------------------------------------------------

# visu_funtopol_UNI()

visu_funtopol_UNI <- function(int_data){



  if (!"Open" %in% int_data$Canopy) stop("your data does not contain a node named Open or it is spelled differently.")

  nodes_list <- funtopol_UNI(int_data)$Functional_classification

  if(max(SCCs$csize)>1){

    open_df <- c("Open", "Open")
    nodes_list <- rbind(nodes_list, open_df)
    nodes_list$label <- nodes_list$id
    int_data <- int_data[which(int_data$Fcr!=0), c("Canopy", "Recruit")]
    g <- igraph::graph_from_data_frame(int_data, directed = TRUE)
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = FALSE)
    g<- igraph::as_data_frame(g, what = "both")
    edges_list <- g$edges

    # nodes data.frame for legend
    lnodes <- data.frame(label = c("Open", "Core", "Satellite", "Strict transient", "Disturbance-dependent transient"),
                         shape = c( "dot"), color = c("#F0E442", "#009E73", "#0072B2", "#D55E00", "#CC79A7"),
                         title = "Functional types", id = 1:5)

    # Network visualization and export to html
    network <- visNetwork::visNetwork(nodes_list, edges_list) |>
      visNetwork::visIgraphLayout(layout = "layout_with_fr") |>
      visNetwork::visEdges(arrows ="to") |>
      visNetwork::visGroups(groupname = "Open", color = "#F0E442") |>
      visNetwork::visGroups(groupname = "Core", color = "#009E73") |>
      visNetwork::visGroups(groupname = "Satellite", color = "#0072B2") |>
      visNetwork::visGroups(groupname = "Strict_transients", color = "#D55E00") |>
      visNetwork::visGroups(groupname = "Disturbance_dependent_transients", color = "#CC79A7") |>
      visNetwork::visOptions(nodesIdSelection = TRUE) |>
      visNetwork::visLegend(addNodes = lnodes, useGroups = FALSE)

    visNetwork::visSave(network, file = "network.html") # Save the html version of the network

    return(network)

  }else{stop("This network does not have a Core, and thus the functional topology can not be visualized")
  }

}

# -------------------------------------------------------

# RN_heatmap_UNI

RN_heatmap_UNI <- function(int_data, weight_var = c("Fcr", "Dcr", "Icr", "Pcr"), scale_top = 1) {

  # manually set node order
  canopy_order <- unique(int_data$Canopy)
  canopy_order <- canopy_order[!canopy_order %in% c('Open')]
  canopy_order <- c("Open", canopy_order)
  int_data$Canopy2 <- factor(int_data$Canopy, levels = canopy_order)
  recruit_order <- sort(unique(int_data$Canopy), decreasing = TRUE)
  recruit_order <- recruit_order[!recruit_order %in% c('Open')]
  recruit_order <- c(recruit_order, "Open")
  int_data$Recruit2 <- factor(int_data$Recruit, levels = recruit_order)

  # Add recruitment density as another weighting variable
  int_data$Dcr <- int_data$Fcr/int_data$Ac

  # Make weight variable
  int_data$weight <- int_data[weight_var]

  # Lowest (non-zero) and highest values of the weighting variable
  highest_W <- max(int_data$weight)
  lowest_W <- min(int_data$weight[int_data$weight>0])

  # Plot the heatmap
  ggplot2::ggplot(int_data, ggplot2::aes(Canopy2, Recruit2, fill= Dcr)) +
    ggplot2::geom_tile(colour="gray", size=0.25, aes(height = 1)) +
    ggplot2::scale_fill_gradientn(colours = c("#F5F5F5", "#E69F00","#0072B2"), values = c(0,lowest_W, scale_top*highest_W)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=0))
}

# -------------------------------------------------------

##### NUEVAS FUNCIONES

# -------------------------------------------------------

# remove_no_cover ONLY CANOPY SPECIES

remove_no_cover_BI <- function(int_data=NULL, cover_data=NULL) {

  # Find species present in RN but that lack data on cover.
  cover_list <- sort(unique(cover_data$Canopy))
  RN_list <- sort(unique(int_data$Canopy))
  lack_cover <- setdiff(RN_list, cover_list)

  # Remove species lacking cover from RN
  if (length(lack_cover) == 0) {
    df <- int_data
  } else {
    df <- int_data[-which(int_data$Canopy %in% lack_cover), ]
  }

  return(df)
}

# -------------------------------------------------------


# aggr_RN WITHOUT EXPANDING THE UNOBSERVED INTERACTIONS WITH ZEROS

aggr_RN_BI <- function(int_data) {

  # Sum the number of recruits per interaction across plots
  RN <- stats::aggregate(Frequency ~ Canopy*Recruit, data = int_data, FUN = sum)
  colnames(RN) <- c("Canopy", "Recruit", "Fcr")
  RN$Icr <- stats::aggregate(Frequency ~ Canopy*Recruit, data=int_data, FUN = NROW)[[3]]
  RN$Pcr <- ifelse(RN$Icr==0,0,1)
  RN$Canopy <- gsub("[[:space:]]", "_", RN$Canopy)
  RN$Recruit <- gsub("[[:space:]]", "_", RN$Recruit)

  # RN <- RN[which(RN$Recruit!="Open"),] # Remove Open from the Recruit species
  return(RN)

}

# -------------------------------------------------------


#add a row of canopy Open if it does not exist


comm_to_RN_BI <- function(int_data, cover_data) {


  if(length(int_data[int_data$Canopy=="Open","Frequency"])==0){

    int_data<-rbind(int_data, int_data[nrow(int_data), ])

    int_data[nrow(int_data),"Canopy"]<-"Open"
    int_data[nrow(int_data),"Standardized_Canopy"]<-"Open"
    int_data[nrow(int_data),"Family_Canopy"]<-"Open"
    int_data[nrow(int_data),"Frequency"]<-0
    int_data[nrow(int_data),"LifeHabit_Canopy"]<-NA


  } else{

    int_data<-int_data}


  # Aggregate data sets across plots.
  int_df <- aggr_RN_BI(int_data)
  cover_df <- aggr_cover_UNI(cover_data)


  # Find species present in RN but that lack data on cover.
  cover_list <- sort(unique(cover_df$Canopy))
  RN_list <- sort(unique(int_df$Canopy))
  lack_cover <- setdiff(RN_list, cover_list)

  # Remove species lacking cover from RN
  if (length(lack_cover) == 0) {
    RNc <- int_df
  } else {
    RNc <- int_df[-which(int_df$Canopy %in% lack_cover), ]
  }

  # Add variables with the cover of the canopy (Ac) and recruit (Ar) species
  RNc$Ac <- RNc$Canopy

  for (i in 1:dim(RNc)[1]) {
    RNc$Ac[i] <- as.numeric(replace(
      RNc$Ac[i],
      match(RN_list, RNc$Canopy[i]),
      cover_df$Ac[match(RNc$Canopy[i], cover_df$Canopy)]
    ))
  }


  RNc <- utils::type.convert(RNc, as.is = TRUE)
  RNc <- RNc[order(RNc$Canopy, RNc$Recruit),]

  return(RNc)

}


# -------------------------------------------------------

# comm_to_RN_UNI_COMP: expands non-observed interactions and removes only canopy species
# without cover data (and maintaining recruit species without cover data)

comm_to_RN_UNI_COMP <- function(int_data, cover_data) {

  if(length(int_data[int_data$Canopy=="Open","Frequency"])==0){

    int_data<-rbind(int_data, int_data[nrow(int_data), ])

    int_data[nrow(int_data),"Canopy"]<-"Open"
    int_data[nrow(int_data),"Standardized_Canopy"]<-"Open"
    int_data[nrow(int_data),"Family_Canopy"]<-"Open"
    int_data[nrow(int_data),"Frequency"]<-0
    int_data[nrow(int_data),"LifeHabit_Canopy"]<-NA


  } else{

    int_data<-int_data}


  # Aggregate data sets across plots.
  int_df <- aggr_RN_UNI(int_data)
  cover_df <- aggr_cover_UNI(cover_data)

  # Find species present in RN but that lack data on cover.
  cover_list <- sort(unique(cover_df$Canopy))
  RN_list <- sort(unique(int_df$Canopy))
  lack_cover <- setdiff(RN_list, cover_list)

  # Remove species lacking cover from RN
  if (length(lack_cover) == 0) {
    RNc <- int_df
  } else {
    RNc <- int_df[-which(int_df$Canopy %in% lack_cover), ]
  }

  # Add variables with the cover of the canopy (Ac) species
  RNc$Ac <- RNc$Canopy

  for (i in 1:dim(RNc)[1]) {
    RNc$Ac[i] <- as.numeric(replace(
      RNc$Ac[i],
      match(RN_list, RNc$Canopy[i]),
      cover_df$Ac[match(RNc$Canopy[i], cover_df$Canopy)]
    ))
  }



  RNc <- utils::type.convert(RNc, as.is = TRUE)
  RNc <- RNc[order(RNc$Canopy, RNc$Recruit),]

  return(RNc)

}



# -------------------------------------------------------


# RN_to_matrix_BI() make a non-squared matrix non expanding the unobserved interactions

RN_to_matrix_BI <- function(int_data=NULL, weight = "Fcr"){

  # Check column names
  if ("Canopy" %in% names(int_data) == FALSE) warning("your data lacks a column named: Canopy")
  if ("Recruit" %in% names(int_data) == FALSE) warning("your data lacks a column named: Recruit")

  data <- int_data
  # Formatting

  RNmat<-as.matrix(tibble::column_to_rownames(reshape2::dcast(data, Recruit ~ Canopy, value.var = weight, fill = 0), "Recruit"))


  return(RNmat)

}


# -------------------------------------------------------

# pre_associndex_UNISITE_BI() replace RN_to_matrix by RN_to_matrix_non_expand and build the matrix non assuming that is square

pre_associndex_UNISITE_BI <- function(int_data = NULL) {
  # Formatting
  mydata <- int_data[,c("Canopy", "Recruit", "Fcr","Ac")]

  # Frequency in open
  Fcr_matrix_df <- as.data.frame(RN_to_matrix_BI(mydata))


  Open<-mydata[mydata$Canopy=="Open",]
  Open$Ao<-Open$Ac
  Open$Fro<-Open$Fcr
  mydata<-merge(mydata[mydata$Canopy!="Open",],Open[,c("Recruit","Fro","Ao")], by="Recruit", all.x=T)

  mydata[which(is.na(mydata$Fro)), "Fro"]<-0
  mydata[which(is.na(mydata$Ao)), "Ao"]<-min(mydata$Ao, na.rm=T)
  mydata <- mydata[mydata$Fcr+mydata$Fro > 0,]
  return(mydata)

}


# -------------------------------------------------------


# pre_associndex_UNISITE_BI_COMP() uses RN_to_matrix_UNI and build the matrix non assuming that is square

pre_associndex_UNISITE_BI_COMP <- function(int_data = NULL) {
  # Formatting
  mydata <- int_data[,c("Canopy", "Recruit", "Fcr","Ac")]

  # Frequency in open
  Fcr_matrix_df <- as.data.frame(RN_to_matrix_UNI(mydata))


  Open<-mydata[mydata$Canopy=="Open",]
  Open$Ao<-Open$Ac
  Open$Fro<-Open$Fcr
  mydata<-merge(mydata[mydata$Canopy!="Open",],Open[,c("Recruit","Fro","Ao")], by="Recruit", all.x=T)

  mydata[which(is.na(mydata$Fro)), "Fro"]<-0
  mydata[which(is.na(mydata$Ao)), "Ao"]<-min(mydata$Ao, na.rm=T)
  mydata <- mydata[mydata$Fcr+mydata$Fro > 0,]
  return(mydata)

}

# -------------------------------------------------------

# associndex_p() use pre_associndex_non_expand instead of pre_associndex and
#remove line to remove a column named "Frequency" that does not exist

associndex_UNISITE_BI <- function(int_data = NULL,
                                  threshold_density = 100) {

  if (!"Open" %in% int_data$Canopy) stop("tests cannot be conducted because your data does not contain a node named Open or it is spelled differently.")

  thr <- threshold_density

  # Assemble the data
  db_inter <- pre_associndex_UNISITE_BI(int_data)

  # Incorporate density of recruitment (recruits/m2) under each canopy species and in open.
  db_inter$Dcr <- db_inter$Fcr/db_inter$Ac
  db_inter$Dro <- db_inter$Fro/db_inter$Ao

  # Retain the interactions with estimated density below the threshold.
  db_inter <- db_inter[which(db_inter$Dcr<thr & db_inter$Dro<thr), ]

  #Obtain the maximum recruitment density for each recruit under the canopy species or in open.
  db_inter$Max_Recr_Density <- pmax(db_inter$Dcr,db_inter$Dro)

  db_inter <- utils::type.convert(db_inter, as.is = TRUE)

  max_rd <- stats::aggregate(Max_Recr_Density ~ Recruit, data = db_inter, FUN = "max")

  # Add a variable max_Recr to each pair indicating the maximum recruitment density of the recruit species in the study site
  Recr_list <- sort(unique(c(db_inter$Recruit)))
  Dens_list <- sort(unique(max_rd$Recruit))
  lack_dens <- setdiff(Recr_list, Dens_list)

  db_inter$max_Recr <- db_inter$Recruit
  for (i in 1:(dim(db_inter)[1])) {
    db_inter$max_Recr[i] <- replace(
      db_inter$max_Recr[i],
      match(Recr_list, db_inter$max_Recr[i]),
      max_rd$Max_Recr_Density[match(db_inter$max_Recr[i], max_rd$Recruit)]
    )
  }

  db_inter <- utils::type.convert(db_inter, as.is = TRUE)

  # Calculate indices Ns, NintC, NintA and RII
  db_inter$Ns <- (db_inter$Dcr - db_inter$Dro)/db_inter$max_Recr
  db_inter$NintC <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dcr + db_inter$Dro)+abs(db_inter$Dcr-db_inter$Dro))
  db_inter$NintA <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dro) + abs(db_inter$Dcr-db_inter$Dro))
  db_inter$RII <- (db_inter$Dcr - db_inter$Dro)/(db_inter$Dcr + db_inter$Dro)

  removed <- names(db_inter) %in% c("max_Recr", "Max_Recr_Density")
  db_inter <- db_inter[!removed]
  return(db_inter)

}


# -------------------------------------------------------

# associndex_p() use pre_associndex_non_expand instead of pre_associndex and
#remove line to remove a column named "Frequency" that does not exist

associndex_UNISITE_BI_COMP <- function(int_data = NULL,
                                       threshold_density = 100) {

  if (!"Open" %in% int_data$Canopy) stop("tests cannot be conducted because your data does not contain a node named Open or it is spelled differently.")

  thr <- threshold_density

  # Assemble the data
  db_inter <- pre_associndex_UNISITE_BI_COMP(int_data)

  # Incorporate density of recruitment (recruits/m2) under each canopy species and in open.
  db_inter$Dcr <- db_inter$Fcr/db_inter$Ac
  db_inter$Dro <- db_inter$Fro/db_inter$Ao

  # Retain the interactions with estimated density below the threshold.
  db_inter <- db_inter[which(db_inter$Dcr<thr & db_inter$Dro<thr), ]

  #Obtain the maximum recruitment density for each recruit under the canopy species or in open.
  db_inter$Max_Recr_Density <- pmax(db_inter$Dcr,db_inter$Dro)

  db_inter <- utils::type.convert(db_inter, as.is = TRUE)

  max_rd <- stats::aggregate(Max_Recr_Density ~ Recruit, data = db_inter, FUN = "max")

  # Add a variable max_Recr to each pair indicating the maximum recruitment density of the recruit species in the study site
  Recr_list <- sort(unique(c(db_inter$Recruit)))
  Dens_list <- sort(unique(max_rd$Recruit))
  lack_dens <- setdiff(Recr_list, Dens_list)

  db_inter$max_Recr <- db_inter$Recruit
  for (i in 1:(dim(db_inter)[1])) {
    db_inter$max_Recr[i] <- replace(
      db_inter$max_Recr[i],
      match(Recr_list, db_inter$max_Recr[i]),
      max_rd$Max_Recr_Density[match(db_inter$max_Recr[i], max_rd$Recruit)]
    )
  }

  db_inter <- utils::type.convert(db_inter, as.is = TRUE)

  # Calculate indices Ns, NintC, NintA and RII
  db_inter$Ns <- (db_inter$Dcr - db_inter$Dro)/db_inter$max_Recr
  db_inter$NintC <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dcr + db_inter$Dro)+abs(db_inter$Dcr-db_inter$Dro))
  db_inter$NintA <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dro) + abs(db_inter$Dcr-db_inter$Dro))
  db_inter$RII <- (db_inter$Dcr - db_inter$Dro)/(db_inter$Dcr + db_inter$Dro)

  removed <- names(db_inter) %in% c("max_Recr", "Max_Recr_Density")
  db_inter <- db_inter[!removed]
  return(db_inter)

}


# -------------------------------------------------------

# int_significance() replace pre_associndex by pre_associndex_UNISITE_BI


int_significance_BI <- function(int_data){

  if (!"Open" %in% int_data$Canopy)
    stop("tests cannot be conducted because your data does not contain a node
         named Open or it is spelled differently.")

  df <- pre_associndex_UNISITE_BI(int_data)
  n_tests <- dim(df)[1]
  df$exp_p <- df$Ac/(df$Ac+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)

  # Testability through Binomial test

  df$Ftot <- df$Fcr+df$Fro

  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p

  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- stats::binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability

  # Binomial (or Chi square) Test Significance

  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Significance[i] <- stats::binom.test(df$Fcr[i], df$Fcr[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- stats::chisq.test(c(df$Fcr[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance

  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }
  df$Test_type <- Test_type
  #  if(length(unique(df$Test_type))>1) warning("Different tests were used for different canopy-recruit pairs. Check column Test_type")

  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fcr[i]/df$Ac[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Enhancing",
                         Effect_int[i] <- "Depressing")
           )
    )
  }

  df$Effect_int <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]

  return(df)
}


# -------------------------------------------------------

# int_significance() replace pre_associndex by pre_associndex_UNISITE_BI_COMP


int_significance_BI_COMP <- function(int_data){

  if (!"Open" %in% int_data$Canopy)
    stop("tests cannot be conducted because your data does not contain a node
         named Open or it is spelled differently.")

  df <- pre_associndex_UNISITE_BI_COMP(int_data)
  n_tests <- dim(df)[1]
  df$exp_p <- df$Ac/(df$Ac+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)

  # Testability through Binomial test

  df$Ftot <- df$Fcr+df$Fro

  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p

  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- stats::binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability

  # Binomial (or Chi square) Test Significance

  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Significance[i] <- stats::binom.test(df$Fcr[i], df$Fcr[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- stats::chisq.test(c(df$Fcr[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance

  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }

  df$Test_type <- Test_type
  #  if(length(unique(df$Test_type))>1) warning("Different tests were used for different canopy-recruit pairs. Check column Test_type")

  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fcr[i]/df$Ac[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Enhancing",
                         Effect_int[i] <- "Depressing")
           )
    )
  }

  df$Effect_int <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]

  return(df)
}


# -------------------------------------------------------

node_degrees_BI <- function(int_data,cover_data){


  matrix_Fcr<-RN_to_matrix(int_data,cover_data, int_type="fac", weight="Fcr")
  matrix_Pcr <-ifelse(matrix_Fcr>0,1,0)
  p<-associndex(int_data,cover_data,expand="no",rm_sp_no_cover="onlycanopy", threshold_density = NULL)
  node_abund <-unique(p[,c("Canopy","Ac")])
  node_abund <-node_abund[node_abund$Canopy%in%colnames(matrix_Fcr),]
  rownames(node_abund)<-node_abund$Canopy

  N_open <-unique(p[,c("Recruit","Fro")])
  N_can<-stats::aggregate(Fcr ~ Recruit, data = p, FUN = sum, na.rm = TRUE)

  node_abund_f <-merge(N_can, N_open, by="Recruit")
  node_abund_f$N_ind <-with(node_abund_f, Fcr+Fro)

  node_abund_f <-node_abund_f[node_abund_f$Recruit%in%rownames(matrix_Fcr),]
  rownames(node_abund_f)<-node_abund_f$Recruit


  # Degree properties
  canopy_service_width <- colSums(matrix_Pcr)
  recruit_niche_width <- rowSums(matrix_Pcr)

  node_abund <-node_abund[names(canopy_service_width),]
  node_abund_f <-node_abund_f[names(recruit_niche_width),]

  node_deg <-list()

  node_deg_can <- data.frame(cbind(names(canopy_service_width), node_abund$Ac, canopy_service_width))
  colnames(node_deg_can) <- c("Nurse_sp", "Ac", "N_enhanced_recruit_sp")
  node_deg_can <- utils::type.convert(node_deg_can, as.is = TRUE)
  rownames(node_deg_can) <- NULL

  node_deg_rec <- data.frame(cbind(names(recruit_niche_width),node_abund_f$N_ind, recruit_niche_width))
  colnames(node_deg_rec) <- c("Facilitated_sp","N_ind", "N_enhancing_canopy_sp")
  node_deg_rec <- utils::type.convert(node_deg_rec, as.is = TRUE)
  rownames(node_deg_rec) <- NULL


  node_deg[[1]]<-node_deg_can
  node_deg[[2]]<-node_deg_rec

  names(node_deg) <- c("Canopy", "Recruit")


  return(node_deg)

}

#--------------------------------------------------------
node_degrees_BI_COMP <- function(int_data,cover_data){


  matrix_Fcr<-RN_to_matrix(int_data,cover_data, int_type="comp", weight="RII")
  matrix_Pcr <-ifelse(matrix_Fcr<0,1,0)
  p<-associndex(int_data,cover_data,expand="yes",rm_sp_no_cover="onlycanopy", threshold_density = NULL)
  node_abund <-unique(p[,c("Canopy","Ac")])
  node_abund <-node_abund[node_abund$Canopy%in%colnames(matrix_Fcr),]
  rownames(node_abund)<-node_abund$Canopy

  N_open <-unique(p[,c("Recruit","Fro")])
  N_can<-stats::aggregate(Fcr ~ Recruit, data = p, FUN = sum, na.rm = TRUE)

  node_abund_f <-merge(N_can, N_open, by="Recruit")
  node_abund_f$N_ind <-with(node_abund_f, Fcr+Fro)

  node_abund_f <-node_abund_f[node_abund_f$Recruit%in%rownames(matrix_Fcr),]
  rownames(node_abund_f)<-node_abund_f$Recruit

  # Degree properties
  canopy_service_width <- colSums(matrix_Pcr)
  recruit_niche_width <- rowSums(matrix_Pcr)

  node_abund <-node_abund[names(canopy_service_width),]
  node_abund_f <-node_abund_f[names(recruit_niche_width),]

  node_deg <-list()

  node_deg_can <- data.frame(cbind(names(canopy_service_width), node_abund$Ac, canopy_service_width))
  colnames(node_deg_can) <- c("Canopy_sp", "Ac", "N_depressed_recruit_sp")
  node_deg_can <- utils::type.convert(node_deg_can, as.is = TRUE)
  rownames(node_deg_can) <- NULL

  node_deg_rec <- data.frame(cbind(names(recruit_niche_width),node_abund_f$N_ind, recruit_niche_width))
  colnames(node_deg_rec) <- c("Recruit_sp","N_ind", "N_depressing_canopy_sp")
  node_deg_rec <- utils::type.convert(node_deg_rec, as.is = TRUE)
  rownames(node_deg_rec) <- NULL


  node_deg[[1]]<-node_deg_can
  node_deg[[2]]<-node_deg_rec

  names(node_deg) <- c("Canopy", "Recruit")

  return(node_deg)
}
#--------------------------------------------------------
