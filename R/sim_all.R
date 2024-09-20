#' A wrap script
#' This package aims to simulate large gene and protein expression time profiles

#####################
### Load dependencies
library(tidyverse)
library(ggplot2)


#####################
### set up data frame


#' The wrapper function which generates the study design with the given parametrization.
#'
#' @param n Number of patients in each group. (Need >50 individuals to have a reasonable estimated covariance structure.)
#' @param m Number of repeated visits (time points).
#' @param q Number of groups (eg., control and multiple treatment groups).
#' @param SingleGroup Whether the individuals are in single group or multiple groups.
#'
#' @return Model matrix with desired design.
#' @export
#'
#' @examples
wrap_design <- function(n,m,q, SingleGroup = FALSE, eps = FALSE) {
  ### design2
  if (SingleGroup) {
    message("Individuals are different in different groups. n is the sample size in each group.")
    if (n < 50) {
      message("Need >50 individuals to have a reasonable estimated covariance structure.")
    }

    if (eps){
      message("n*eps samples are from identical subjects in different groups.")

      p <- round(runif(n,m,m)) # balanced
      # visit <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
      visit <- paste("T", rep(1:m, n, each=TRUE), sep="")
      ## each patient measured in only one group
      dat <- data.frame()
      nRep = ceiling(n*eps)
      for (n_trt in 1:q) {
        dat_trt <- data.frame(ID = paste("P", (rep(1:n, times=p) + (n_trt-1)*n), sep=""),
                              ID_num = rep(1:n, times=p) + (n_trt-1)*n,
                              VISIT=visit,
                              GROUP = rep(paste0("Group", n_trt), m*n))
        dat <- rbind(dat, dat_trt)
      }
      for (s in 1:nRep) {
        dat$ID_num[dat$ID_num %% n == s] <- s
      }
      dat$ID <-  paste("P", dat$ID_num, sep="")
      rownames(dat) <- paste("Sample", rownames(dat), sep="_")
    } else {
      p <- round(runif(n,m,m)) # balanced
      # visit <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
      visit <- paste("T", rep(1:m, n, each=TRUE), sep="")
      ## each patient measured in only one group
      dat <- data.frame()
      for (n_trt in 1:q) {
        dat_trt <- data.frame(ID = paste("P", (rep(1:n, times=p) + (n_trt-1)*n), sep=""),
                              VISIT=visit,
                              GROUP = rep(paste0("Group", n_trt), m*n))
        dat <- rbind(dat, dat_trt)
      }
      rownames(dat) <- paste("Sample", rownames(dat), sep="_")
    }
  } else {
    message("Samples from the same individual are analysed in different groups. n is the number of individuals.")
    if (n < 50) {
      message("Need >50 individuals to have a reasonable estimated covariance structure.")
    }
    ### simulate number of visits for each individual
    p <- round(runif(n,m,m)) # balanced
    # visit <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
    visit <- paste("T", rep(1:m, n, each=TRUE), sep="")
    ## each patient measured in both control group and treatment group
    ###### More than one treatment!!! (good to write in groups)
    dat <- data.frame()
    for (n_trt in 1:q) {
      dat_trt <- data.frame(ID = paste("P", rep(1:n, times=p), sep=""),
                            VISIT=visit,
                            GROUP = rep(paste0("Group", n_trt), m*n))
      dat <- rbind(dat, dat_trt)
    }
    rownames(dat) <- paste("Sample", rownames(dat), sep="_")
  }

  return(dat)
}

# dat <- wrap_design(n = 10, m = 4, q=3)
# dat <- wrap_design(n = 10, m = 4, q=3, SingleGroup=TRUE)
#
# dat[1:10, ] %>%  knitr::kable()
# dat[90:100, ] %>%  knitr::kable()

#####################
### Simulation
#' TKT: GPT version
#' This function simulates data based on a design matrix, with options to specify various parameters for the simulation. It allows for the simulation of correlated random effects for intercepts and slopes, and can handle single or multiple groups.
#'
#' @param data A data frame (`data.frame`) containing the design matrix for the simulation. It must include columns for 'GROUP', 'VISIT', and 'ID'.
#' @param SingleGroup A logical value (`logical`) indicating whether samples from the same individual are analyzed in a single group (TRUE) or multiple groups (FALSE).
#' @param beta0 The intercept for the simulation (`numeric`). If NULL, it is set to 0 by default.
#' @param beta1 The slope coefficients for the group effect (`numeric` vector). If NULL, a default sequence is created.
#' @param beta2 The slope coefficients for the time effect (`numeric` vector). If NULL, a default logarithmic sequence is created.
#' @param lambda The interaction matrix between GROUP and VISIT (`matrix`). If NULL, it is set to a matrix of zeros by default.
#' @param sigma The standard deviation of the error terms (`numeric`). If NULL, it is set to 1.5 by default.
#' @param mu The mean vector for the random effects (`numeric` vector). If NULL, it is set to c(0, 0) by default.
#' @param tau The vector containing the diagonal and lower triangle elements of the variance matrix for the random effects (`numeric` vector). If NULL, it is set to c(1.2, 1.6, 0.3) by default.
#' @param act.cor The autocorrelation coefficient for the AR(1) error structure (`numeric`). If NULL, it is set to 0.6 by default.
#' @return A data frame (`data.frame`) with the simulated data, including the simulated outcome variable 'Exprs'.
#' @export
#' @note The function assumes that the input data frame has the correct structure and necessary columns. It is the responsibility of the user to ensure that the data frame is properly formatted before calling this function. The default values for parameters are chosen based on typical simulation scenarios and may need to be adjusted based on specific study designs.
#' @examples
#' # This is an example of how to use the wrap_simulation function.
#' # You need to provide a valid 'data' data frame with the necessary structure.
#' # wrap_simulation(data, SingleGroup = TRUE)
wrap_simulation <- function(data, SingleGroup,
                            beta0 = NULL,
                            beta1 = NULL,
                            beta2 = NULL,
                            lambda = NULL,
                            sigma = NULL,
                            mu = NULL,
                            tau = NULL,
                            act.cor = NULL) {
  # Function body remains the same
}
#' TKT: OLDER version
#' Simulate the log-scale expression of one gene based on the given design and model parameters.
#' Wrap Simulation Function
#'
#' @param data Design matrix in long format.
#' @param beta0 Mean of the global effect.
#' @param beta1 Mean of group effect.
#' @param beta2 Mean of time effect.
#' @param lambda Coefficient of interaction terms.
#' @param sigma Standard deviation of the simulated AR(1) sequence.
#' @param mu Mean vector of random intercept and slope.
#' @param tau Elements of variance-covariance matrix of random intercept and slope.
#' @param act.cor Correlation parameter of the AR(1) structure.
#'
#' @return Data matrix in long data frame.
#' @export
#'
#' @examples
wrap_simulation <- function(data, SingleGroup,
                            beta0 = NULL,
                            beta1 = NULL,
                            beta2 = NULL,
                            lambda = NULL,
                            sigma = NULL,
                            mu = NULL,
                            tau = NULL,
                            act.cor = NULL) {
  if (!is.data.frame(data)) {
    message("Please define a design matrix for simulation.")
  } else{
    ### read the design matrix
    q = length(unique(data$GROUP))
    m = length(unique(data$VISIT))

    if (SingleGroup == FALSE) {
      # if ( length(data$ID) > m*length(unique(data$ID)) ){
      message("Samples from the same individual are analysed in different groups.")
      n = length(unique(data$ID))
      nRep = 0
      SGFLAG = 0
    } else {
      message("Individuals are different in different groups.")
      SGFLAG = 1
      GroupInd <- data$GROUP[1]
      n = length(unique(data[data$GROUP == GroupInd, ]$ID))
      nRep = (n*q - length(unique(data$ID)))/(max(q-1, 1))
      eps = nRep/n
      nUnq = n - nRep
      cat("\n", nRep, "subjects appear repeatedly in different groups. \n")
    }

    ### simulate (correlated) random effects for intercepts and slopes
    if (is.null(beta0)) {
      beta0 = 0
      message("Using the default beta0 = 0.")
    } else{
      if (length(beta0) != 1) {
        message("Please define beta0 as a constant of dimention 1.")
      }
    }

    if (is.null(beta1)) {
      b0 = 4 ## set up an effect size
      beta1 = seq(0, b0*(q-1), by=b0)
      message("Using the default beta1 with group effect size = 4.")
    } else {
      if (length(beta1) != q) {
        message("Please define beta1 as a vector of dimention q.")
      }
    }


    if (is.null(beta2)) {
      beta2 = log(c(1:m))  ## dim = m
      message("Using the default beta2 reflecting a logrithm dependence of time point.")
    } else {
      if (length(beta2) != m) {
        message("Please define beta2 as a vector of dimention m.")
      }
    }

    if (is.null(lambda)) {
      lambda = matrix(0, byrow = TRUE, ncol = m, nrow=q)
      message("Using the default lambda reflecting no interaction between GROUP and VISIT.")
    } else {
      if (!((is.matrix(lambda)) & (nrow(lambda) == q) & (ncol(lambda) == m))) {
        message("Please define lambda as a matrix of dimention (q, m).")
      }
    }


    ### random effect
    if (is.null(mu)) {
      mu = c(0,0)
      message("Using the default mu = (0, 0).")
    } else {
      if (length(mu) != 2) {
        message("Please define mu as a vector of dimention 2.")
      }
    }

    if (is.null(tau)) {
      tau =   c(1.2, 1.6, 0.3)
      message("Using the default tau0 = (1.2, 1.6) and tau01 = 0.3 to generate the variance matrix of the random intercept and slope.")
    } else {
      if (length(tau) != 3) {
        message("Please define tau as a vector of dimention 3 containing the diagonal and lower triangle element of the variance matrix.")
      }
    }

    ## n_tot*random effect (intercept and slope), alpha is dim=2
    tau0 <- tau[1:2]
    tau01 = tau[3]
    S   = matrix(c(1, tau01, tau01, 1), nrow=2)
    S   = diag(tau0) %*% S %*% diag(tau0)

    ## total number of individuals in different scenarios
    n_tot <- length(unique(data$ID))
    U = mvrnorm((n_tot + nRep*(q-1)), mu=mu, Sigma=S)
    ## Work out the formula  of the sigma^2 within and between subjects!!!

    ### correlation structure
    ### TODO: give the values and plot!!!!
    if (is.null(sigma)) {
      sigma = 1.5
      message("Using the default sigma = 1.5.")
    } else {
      if (length(sigma) != 1) {
        message("Please define sigma as a constant of dimention 1.")
      }
    }

    if (is.null(act.cor)) {
      act.cor = 0.6
      message("Using the default act.cor = 0.6 for the AR correlation structure.")
    } else {
      if (length(act.cor) != 1) {
        message("Please define act.cor as a constant of dimention 1.")
      }
    }

    ### simulate AR(1) errors and then the actual outcomes
    ### true autocorrelation
    ar.val <- act.cor
    # sqrt(1-ar.val^2)
    # ### ARIMASIM <- arima.sim(model=list(ar=ar.val), n=m) * sqrt(1-ar.val^2)* sigma

    if (SGFLAG==1) {
      p <- round(runif((n_tot + nRep*(q-1)),m,m))
      VecEij <- sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2)* sigma )
      if (nRep >0){
        for (nSub in 1:(n_tot + nRep*(q-1)) ) {
          if ((nSub %% n) %in% (1:nRep)){
            VecEij[, nSub] <- VecEij[, (nSub %% n)]
          }
        }
      }
      data$eij <- c(unlist(VecEij))
      ### model with random intercept and slope
      # data$Exprs <- (beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p)) * GROUP *log(data$VISIT) + data$eij
      y_cell <- array(0, c(m,n,q))
      for (i in 1:q) {
        for (j in 1:m) {
          for (k in 1:n){
            k_tot = n*(i-1) + k ### real individual index
            y_cell[j, k, i] <- beta0 + U[k_tot, 1]  + (beta1[i] + U[k_tot, 2] )*( beta2[j] ) +  lambda[i, j]
          }
        }
      }
    } else {
      p <- round(runif(n,m,m))
      Group_seq <- unique(data$GROUP)
      for (G_info in Group_seq) {
        data$eij[data$GROUP == G_info] <- c(unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2)* sigma )))
      }
      ### model with random intercept and slope
      # data$Exprs <- (beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p)) * GROUP *log(data$VISIT) + data$eij
      y_cell <- array(0, c(m,n,q))
      for (i in 1:q) {
        for (j in 1:m) {
          for (k in 1:n){
            y_cell[j, k, i] <- beta0 + U[k, 1]  + (beta1[i] + U[k, 2] )*( beta2[j] ) +  lambda[i, j]
          }
        }
      }


    }

    data$Exprs <- data$eij + c(y_cell)

    return(data)
  }
}



########################################
## Simulate more than one gene
########################################
### Simulate DGE List
#' Simulate Differential Gene Expression Data
#'
#' Simulates gene expression data for a given number of genes and samples, and returns a list
#' containing the normalized counts matrix, sample phenotype data, and feature (gene) data.
#'
#' @param n Integer, the number of biological replicates per group.
#' @param m Integer, the number of time points.
#' @param q Integer, the number of groups.
#' @param N_gene Integer, the number of genes to simulate.
#' @param global_eff Numeric matrix or NULL, global effects for each gene defined in each column.
#' @param trt_eff Numeric vector, treatment effects for each gene.
#' @param trt_slope Numeric vector or NULL, treatment slope effects for each gene.
#' @param inter Numeric vector or NULL, interaction effects for each gene.
#' @param sigma Numeric vector or NULL, residual standard deviations for each gene.
#' @param mu Numeric vector or NULL, mean expression levels for each gene.
#' @param tau Numeric vector or NULL, overdispersion parameters for each gene (applicable if `distribution` is "NB").
#' @param act.cor Numeric vector or NULL, actual correlations for each gene.
#' @param distribution Character string or NULL, specifies the distribution to use for simulation. If "NB", a negative binomial distribution is used.
#' @param SingleGroup Logical, if TRUE, simulates data for a single group.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{norm_counts}: Transposed normalized counts matrix with genes in rows and samples in columns.
#'   \item \code{pheno_data}: Data frame containing phenotype data for the samples.
#'   \item \code{feature_data}: Data frame containing feature (gene) data.
#' }
#'
#' @examples
#' \dontrun{
#' ngene = 10000
#' simulated_DGEList <- wrap_sim_DGEList(n = 5,
#'                                       m = 3,
#'                                       q = 2,
#'                                       N_gene = ngene,
#'                                       trt_eff = matrix(rep(c(1,2), ngene), ncol=ngene),
#'                                       SingleGroup = FALSE)
#' }
#'
#' @note The function assumes that the necessary libraries (such as `dplyr`, `edgeR`, and any other required packages) are already loaded in the R session.
#' The actual simulation functions `wrap_simulation` and `wrap_simulation_NB` are not defined within this documentation and should be provided for the function to work.
#'
#' @author [Zhiwen Jiang and Tony Kam-Thong]
#' @seealso \code{\link[edgeR]{DGEList}}
#'
#' @export

wrap_sim_DGEList <- function(n,
                             m,
                             q,
                             N_gene,
                             global_eff = NULL,
                             trt_eff,
                             trt_slope = NULL,
                             inter = NULL,
                             sigma = NULL,
                             mu = NULL,
                             tau = NULL,
                             act.cor = NULL,
                             distribution = NULL,
                             SingleGroup = FALSE){

  data1 <- wrap_design(n, m, q)
  loop_sim_gene <- function(N,
                            data1,
                            global_genewise,
                            trt_genewise,
                            slope_genewise,
                            inter_genewise,
                            sigma_genewise,
                            mu_genewise,
                            tau_genewise,
                            act.cor_genewise,
                            distribution,
                            SingleGroup=FALSE){
    if (distribution == NULL) {
    norm_counts<- wrap_simulation(
      data = data1,
      beta0 = global_genewise,
      beta1 = trt_genewise,
      beta2 = slope_genewise,
      lambda = inter_genewise,
      sigma = sigma_genewise,
      mu = mu_genewise,
      tau = tau_genewise,
      act.cor = act.cor_genewise,
      SingleGroup) %>%
      # dplyr::mutate((!!(GeneNames[N])) = Exprs) %>%
      dplyr::mutate(GENE = Exprs) %>%
      rownames_to_column(., var = "SampleID") %>%
      dplyr::select(-eij, -Exprs)
    colnames(norm_counts) <- gsub("GENE", paste0("GENE_", N),
                                  colnames(norm_counts))
    return(norm_counts)
    } else if (distribution == "NB"){
      norm_counts<- wrap_simulation_NB(
        data = data1,
        beta0 = global_genewise,
        beta1 = trt_genewise,
        beta2 = slope_genewise,
        lambda = inter_genewise,
        sigma = sigma_genewise,
        mu = mu_genewise,
        tau = tau_genewise,
        act.cor = act.cor_genewise,
        SingleGroup) %>%
        # dplyr::mutate((!!(GeneNames[N])) = Exprs) %>%
        dplyr::mutate(GENE = Exprs) %>%
        rownames_to_column(., var = "SampleID") %>%
        dplyr::select(-eij, -Exprs)
      colnames(norm_counts) <- gsub("GENE", paste0("GENE_", N),
                                    colnames(norm_counts))
      return(norm_counts)
    }
  }
  # browser()
  norm_counts_global <- mclapply(1:N_gene,
                                 function(N) loop_sim_gene(
                                   N,
                                   data1,
                                   global_eff[, N],
                                   trt_eff[, N],
                                   trt_slope[, N],
                                   inter[[N]],
                                   sigma[, N],
                                   mu[[N]],
                                   tau[[N]],
                                   act.cor[, N],
                                   SingleGroup=FALSE),
                                 mc.cores = 2)

  norm_counts_Matrix <- norm_counts_global[[1]]
  for (N in 2:N_gene) {
    norm_counts_Matrix <- full_join(norm_counts_global[[N]], norm_counts_Matrix,
                                    by = c("SampleID", "ID", "VISIT", "GROUP"))
  }

  counts_Matrix <- norm_counts_Matrix %>%
    column_to_rownames(., "SampleID") %>%
    dplyr::select(-ID, -VISIT, -GROUP)

  sample_pheno <- norm_counts_Matrix %>%
    column_to_rownames(., "SampleID") %>%
    dplyr::select(ID, VISIT, GROUP)

  feature_gene <- norm_counts_Matrix %>%
    dplyr::select(-SampleID, -ID, -VISIT, -GROUP) %>%
    colnames() %>% as.data.frame() %>%
    dplyr::rename(GeneName = ".") %>%
    mutate(FeatureID = GeneName)

  DGEList_norm <- list(
    norm_counts = t(counts_Matrix),
    pheno_data = sample_pheno,
    feature_data = feature_gene
  )


  return(DGEList_norm)
}

#################################
###### Simulate with corr defined
#################################
### AR(1) corr
#' Title AR(1) correlation matrix
#'
#' @param rho Correlation coefficient
#' @param d Dimension of repeated measures
#'
#' @return AR(1) correlation matrix
#' @export
#'
#' @examples
fVarAR1 <- function(rho, d){
  VarAR1 <- matrix(NA, nrow = d, ncol = d)
  return(rho^abs(row(VarAR1)-col(VarAR1)))
}
# fVarAR1(0.8, 4)

### CS corr
#' Title CS correaltion matrix
#'
#' @param rho Correlation coefficient
#' @param d Dimension of repeated measures
#'
#' @return CS correlation matrix
#' @export
#'
#' @examples
fVarCS <- function(rho, d) {
  VarCS = matrix(rho, nrow = d, ncol = d)
  diag(VarCS) = rep(1, d)
  return(VarCS)
}
# fVarCS(0.5, 5)

#' Generate the mariginal variance matrix
#'
#' @param CorrSub Subject level correlation matrix, dim nq*nq
#' @param SigmaGroup Group level variance, length q
#' @param CorrRep Correlation matrix between repeated measures, dim m*m
#' @param SigmaRep Variance of repeated measures, length m
#'
#' @return Marginal covariance matrix cross samples, dim nqm*nqm
#' @export
#'
#' @examples
CovMarginal_sim <- function(CorrSub, SigmaGroup, CorrRep, SigmaRep, sigmaRE=NULL){
  if(is.null(sigmaRE)){
    ### OPTION 1: Same time-dependent correlation structure applied to between-subject correlation
    m = length(SigmaRep)
    q = length(SigmaGroup)
    n = nrow(CorrSub)/q
    sigmaGroup <- rep(SigmaGroup, each=n)
    veps_sigmaGroup <- matrix(sqrt(sigmaGroup),nrow=n*q,ncol=1) %*% matrix(sqrt(sigmaGroup),nrow=1,ncol=n*q)*CorrSub
    CovMatrix <- kronecker(veps_sigmaGroup,
                           matrix(sqrt(SigmaRep),nrow=m,ncol=1) %*% matrix(sqrt(SigmaRep),nrow=1,ncol=m)*CorrRep)
  }  else {
    ### OPTION 2: No time-dependent correlation structure between subjects. Only CS correlation plugged in.
    m = length(SigmaRep)
    q = length(SigmaGroup)
    n = nrow(CorrSub)/q
    sigmaGroup <- rep(SigmaGroup, each=n)
    veps_sigmaGroup <- matrix(sqrt(sigmaGroup),nrow=nrow(CorrSub),ncol=1) %*%
      matrix(sqrt(sigmaGroup),nrow=1,ncol=nrow(CorrSub)) *
      diag(diag(CorrSub), nrow=nrow(CorrSub))
    CovMatrix <- kronecker(veps_sigmaGroup,
                           matrix(sqrt(SigmaRep),nrow=m,ncol=1) %*%
                             matrix(sqrt(SigmaRep),nrow=1,ncol=m)*CorrRep) +
      sigmaRE^2*kronecker(ceiling(abs(CorrSub)), fVarCS(1, m))
  }
  return(CovMatrix)
}

#' Title Plot the variance-covariance matrix in heatmap
#'
#' @param CMatrix A squared matrix with numeric elements
#'
#' @return
#' @export
#'
#' @examples
CovHeatPlot <- function(CMatrix){
  Cov_heatmap <- ggplot(
    data = reshape2::melt(CMatrix),
    aes(reorder(Var2,  -dplyr::desc(factor(Var2))),
        reorder(Var1,  dplyr::desc(factor(Var1))),
        fill = value)) +
    geom_tile() +
    xlab("") +
    ylab("") +
    coord_fixed(ratio = 1) +
    theme(legend.position = "right")  +
    scale_fill_gradient2(low="green", high="red", mid="white",
                         breaks=seq(-1,1,0.1))
  print(Cov_heatmap)
}

###############################################
#' Title Simulate DGE List
#'
#' @param data Design matrix in long format.
#' @param N Number of simulated genes in total.
#' @param eps_Signal Proportion of signal genes.
#' @param SingleGroup
#' @param coef.beta0 Mean of the global effect.
#' @param coef.beta1 Mean of group effect.
#' @param coef.beta2 Mean of time effect.
#' @param coef.int Coefficient of interaction terms.
#' @param sigma Standard deviation of the simulated AR(1) sequence.
#' @param mu Mean vector of random intercept and slope.
#' @param tau Elements of variance-covariance matrix of random intercept and slope.
#' @param CovMarginal Marginal covariance.
#' @param act.cor Correlation coefficient of the time effect.
#'
#' @return DGE list
#' @export
#'
#' @examples
wrap_sim_corr <- function(data, N, eps_Signal,
                          SingleGroup,
                          coef.beta0 = NULL,
                          coef.beta1 = NULL,
                          coef.beta2 = NULL,
                          coef.int = NULL,
                          sigma = NULL,
                          mu = NULL,
                          tau = NULL,
                          CovMarginal = NULL, ### top level
                          act.cor = NULL){
  ####
  if (!is.data.frame(data)) {
    message("Please define a design matrix for simulation.")
  } else{
    ### read the design matrix
    q = length(unique(data$GROUP))
    m = length(unique(data$VISIT))

    if (SingleGroup == FALSE) {
      # if ( length(data$ID) > m*length(unique(data$ID)) ){
      message("Samples from the same individual are analysed in different groups.")
      n = length(unique(data$ID))
      nRep = 0
      SGFLAG = 0
    } else {
      message("Individuals are different in different groups.")
      SGFLAG = 1
      GroupInd <- data$GROUP[1]
      n = length(unique(data[data$GROUP == GroupInd, ]$ID))
      nRep = (n*q - length(unique(data$ID)))/(max(q-1, 1))
      eps = nRep/n
      nUnq = n - nRep
      cat("\n", nRep, "subjects appear repeatedly in different groups. \n")
    }

    ## default values of params
    G_Signal <- ceiling(N*eps_Signal)
    G_Null <- N - G_Signal

    if (is.null(coef.beta0)) {
      coef.beta0_null = 0
      coef.beta0_Signal = 0.5
      coef.beta0 = c(rep(coef.beta0_Signal, G_Signal), rep(coef.beta0_null, G_Null))
      message("Using the default coef.beta0.")
    } else{
      if (length(coef.beta0) != N) {
        message("Please define coef.beta0 as a constant vector of length N.")
      }
    }

    if (is.null(coef.beta1)) {
      b0 = 4 ## set up an effect size
      coef.beta1_Signal <- lapply(1:G_Signal,
                                  function(x) {seq(0, b0*(q-1), length=q)})
      coef.beta1_null <- lapply(1:G_Null,
                                function(x) {seq(0, b0*(q-1)/2, length=q)})
      coef.beta1 <- c(coef.beta1_Signal, coef.beta1_null)
      message("Using the default beta1 for group effect size.")
    } else {
      if (length(coef.beta1) != N ) {
        message("Please define coef.beta1 as a list of length N with vector of dimention q.")
      }
    }


    if (is.null(coef.beta2)) {
      coef.beta2_Signal <- lapply(1:G_Signal,
                                  function(x) {log(c(1:m))})
      coef.beta2_null <- lapply(1:G_Null,
                                function(x) {log(c(1:m))/2})
      coef.beta2 <- c(coef.beta2_Signal, coef.beta2_null)
      message("Using the default coef.beta2 reflecting a logrithm dependence of time point.")
    } else {
      if (length(coef.beta2) != N) {
        message("Please define coef.beta2 as a list of length N with vector of dimention m.")
      }
    }

    if (is.null(coef.int)) {
      coef.int_Signal <- lapply(1:G_Signal,
                                function(x) {matrix(0, byrow = TRUE, ncol = m, nrow=q)})
      coef.int_null <- lapply(1:G_Null,
                              function(x) {matrix(0, byrow = TRUE, ncol = m, nrow=q)})
      coef.int <- c(coef.int_Signal, coef.int_null)
      message("Using the default coef.int reflecting no interaction between GROUP and VISIT.")
    } else {
      if (!((is.matrix(coef.int[[1]])) & (nrow(coef.int[[1]]) == q) & (ncol(coef.int[[1]]) == m))) {
        message("Please define coef.int as a list of lenth N with matrix of dimention (q, m).")
      }
    }


    ### random effect
    if (is.null(mu)) {
      mu = c(0,0)
      message("Using the default mu = (0, 0).")
    } else {
      if (length(mu) != 2) {
        message("Please define mu as a vector of dimention 2.")
      }
    }

    if (is.null(tau)) {
      tau =   c(1.2, 1.6, 0.3)
      message("Using the default tau0 = (1.2, 1.6) and tau01 = 0.3 to generate the variance matrix of the random intercept and slope.")
    } else {
      if (length(tau) != 3) {
        message("Please define tau as a vector of dimention 3 containing the diagonal and lower triangle element of the variance matrix.")
      }
    }

    if (is.null(sigma)) {
      # default marginal variance
      sigma <- seq(1, m, length = m)
      message("Using the default marginal variance.")
    } else {
      if (length(sigma) != m) {
        message("Please define sigma as a constant of dimention m.")
      }
    }

    if (is.null(act.cor)) {
      act.cor_Signal <- rep(0.6, G_Signal)
      act.cor_null <- rep(0, G_Null)
      act.cor <- c(act.cor_Signal, act.cor_null)
      message("Using the default act.cor = 0.6 for the AR correlation structure.")
    } else {
      if (length(act.cor) != N) {
        message("Please define act.cor as a vector of length N with constant of dimention 1.")
      }
    }

    ######
    ### Design 1: single-treatment on each individual
    # Random intercept
    n_tot <- length(unique(data$ID))
    coef.alpha <- MASS::mvrnorm(1, rep(0, n_tot), tau[1]*diag(n_tot))
    # browser()


    G.tot <- N
    SampleNum <- m*n*q
    YG <- matrix(0, G.tot, SampleNum)
    for (g in 1:G.tot) {
      beta0 <- coef.beta0[g]
      beta1 <- coef.beta1[[g]]
      beta2 <- coef.beta2[[g]]
      int <- coef.int[[g]]
      sigmasqrt<-sqrt(sigma)

      if (is.null(CovMarginal)) {
        CovMarginal = kronecker(diag(n*q),
                                matrix(sqrt(sigma),nrow=m,ncol=1) %*% matrix(sqrt(sigma),nrow=1,ncol=m)*fVarAR1(act.cor[g], m) )
        # message("Using the default AR1 correlation structure.") ##
      } else {
        if (nrow(CovMarginal) != m*n*q ) {
          message("Please define the marginal covariance matrix.")
        }
      }
      epsg <- MASS::mvrnorm(1, rep(0, SampleNum), CovMarginal)

      Yg <- NULL
      for (i in (1:SampleNum)) {
        Xgi1 <- ceiling(i/(m*n))
        Xgi2 <- (i-1) %% m+1
        Zgi <- ceiling(i/m)
        Yg[i] <-  beta0 + beta1[Xgi1] + beta2[Xgi2] + coef.alpha[Zgi] + int[Xgi1,Xgi2] + epsg[i]
        # Yg_count[i] <- ceiling(2^(Yg[i] + log2(1e+07)-log2(1e+06)))
      }
      YG[g, ] <- Yg
    }
    # browser()

    rownames(YG) <- rownames(YG, do.NULL = FALSE, prefix = "GENE_")
    colnames(YG) <- colnames(YG, do.NULL = FALSE, prefix = "Sample_")
    counts_Matrix <- as.data.frame(YG)
    sample_pheno <- data ## %>% mutate(RandomInt = ) ####
    feature_gene <- YG %>%
      rownames() %>% as.data.frame() %>%
      dplyr::rename(GeneName = ".") %>%
      mutate(FeatureID = GeneName) ##  %>%

    DGEList_norm <- list(
      norm_counts = counts_Matrix,
      pheno_data = sample_pheno,
      feature_data = feature_gene
    )

    return(DGEList_norm)
  }

}






################## Plot profiles
#' Title Plot gene profiles
#'
#' @param data Simulated data frame of gene profiles
#'
#' @return
#' @export
#'
#' @examples
plot_prof <- function(data) {
  n_Group <- length(unique(data$GROUP))
  n_Rep <- data %>% group_by(VISIT, ID, GeneName) %>%
    dplyr::summarise(N = length(ID)) %>% ungroup()
  ## length(ID) to count repeated patients
  n_Rep <- n_Rep %>% dplyr::select(ID, N) %>% unique() %>%
    mutate(N = paste0("Freq_", N))
  data <- data %>% left_join(n_Rep, by=c("ID" = "ID"))

  prof = ggplot(data, aes(x=VISIT, y=Exprs), color=ID) +
    geom_line(aes(color=ID, group=ID)) +
    geom_smooth(aes(group=GROUP), method = "lm") +
    facet_grid(cols = vars(GROUP), rows = vars(GeneName), scales="free_y") +
    theme(legend.position = "None") +
    ggtitle("") +
    ylab("Expression") +
    xlab("Visit")
  print(prof)
}

# prof_diff <- function(data) {
# data.wide <- data %>% dplyr::select(-eij) %>%
#   tidyr::pivot_wider(names_from = GROUP, values_from = Exprs, names_prefix ="GROUP") %>%
#   mutate(diff = GROUPTreatment - GROUPControl)
# return(data.wide)
# }


# plot_diff <- function(datawide) {
# ggplot(datawide, aes(x=VISIT, y=diff, group = ID), color=ID) +
#     geom_line(aes(color=ID))  +
#     theme(legend.position = "none") +
# # scale_color_gradientn(colours = rainbow(n)) +
#     ggtitle("Difference between treatments") +
#     ylab("Differential expressions") +
#     xlab("Visit")
# }


## eij
#' Title The correlation and covariance
#'
#' @param data Simulated data frame
#' @param condition A specific group level
#' @param plot.line Plot the line graph (TRUE or FALSE)
#' @param plot.heat Plot the heatmap (TRUE or FALSE)
#'
#' @return
#' @export
#'
#' @examples
CovCor_err <- function(data, condition,
                       plot.line = FALSE,
                       plot.heat = FALSE) {
  data_table_simu.eij <- reshape2::dcast(data, ID + GROUP  ~ VISIT, value.var = "eij")
  ### Focus on one treatment group
  data_table_simu.1 <- data_table_simu.eij %>% filter(GROUP == condition)
  VisitTime <- sort(unique(data$VISIT))
  cor_struc <- cor(data_table_simu.1[, VisitTime], use="pairwise.complete.obs",method="pearson")
  cov_struc <- cov(data_table_simu.1[, VisitTime], use="pairwise.complete.obs",method="pearson")
  ### Count time lag
  LagTime <- seq(1, length(VisitTime), by=1)
  rownames(cor_struc) <- as.character(LagTime)
  colnames(cor_struc) <- as.character(LagTime)
  rownames(cov_struc) <- as.character(LagTime)
  colnames(cov_struc) <- as.character(LagTime)

  Covplot <- CovStructRM::plotCovar(CovStructRM::vecSym(cov_struc),
                                    x.var='lag', col="black", xlab="Lag",
                                    ylab="Empirical Covariance", main = "Covariance")
  Corplot <- CovStructRM::plotCovar(CovStructRM::vecSym(cor_struc),
                                    x.var = 'lag', col="black", xlab="Lag",
                                    # ylim=c(0.5,1.05),
                                    ylab="Empirical Correlation", main = "Correlation")

  CorHeatplot <- ggplot(
    data = reshape2::melt(cor_struc),
    aes(reorder(Var2,  -dplyr::desc(factor(Var2))),
        reorder(Var1,  dplyr::desc(factor(Var1))),
        fill = value)) +
    geom_tile() +
    xlab("Visit") +
    ylab("Visit") +
    coord_fixed(ratio = 1) +
    theme(legend.position = "right")  +
    scale_fill_gradient2(low="green", high="red", mid="white",
                         breaks=seq(-1,1,0.1)) +
    ggtitle("Empirical correlation")

  if (plot.line) {
    print(Covplot, split=c(1,1,2,1), more=TRUE)
    print(Corplot, split=c(2,1,2,1), more=FALSE)
  }


  if (plot.heat) {
    print(CorHeatplot)
  }
  return(list(ECor = cor_struc,
              ECov = cov_struc))
}














garsim_theta <- function (n, phi, X = matrix(0, nrow = n), beta = as.matrix(0),
                          sd = 1, family = "gaussian", transform.Xbeta = "identity",
                          link = "identity", minimum = 0, zero.correction = "zq1",
                          c = 1, theta = 0)
{
  y <- matrix(0, nrow = n)
  m <- matrix(0, nrow = n)
  lambda <- matrix(0, nrow = n)
  y.transformed <- matrix(0, nrow = n)
  phi <- as.matrix(phi)
  p <- nrow(phi)
  if (transform.Xbeta == "identity") {
    Xbeta <- as.matrix(as.matrix(X) %*% as.vector(beta))
  }
  if (transform.Xbeta == "exponential") {
    Xbeta <- exp(as.matrix(as.matrix(X) %*% as.vector(beta)))
  }
  if (link == "identity") {
    if (family == "gaussian") {
      for (t in 1:n) {
        if (t < (p + 1)) {
          m[t] <- Xbeta[t]
          y[t] <- rnorm(1, m[t], sd = sd)
        }
        else {
          yphi <- t(y[(t - p):(t - 1)]) %*% rev(phi)
          Xbetaphi <- t(-Xbeta[(t - p):(t - 1)]) %*%
            rev(phi)
          m[t] <- yphi + Xbetaphi + Xbeta[t]
          y[t] <- rnorm(1, m[t], sd = sd)
        }
      }
    }
    if (family == "poisson") {
      for (t in 1:n) {
        if (t < (p + 1)) {
          m[t] <- Xbeta[t]
          y[t] <- rpois(1, m[t])
        }
        else {
          yphi <- t(y[(t - p):(t - 1)]) %*% rev(phi)
          Xbetaphi <- t(-Xbeta[(t - p):(t - 1)]) %*%
            rev(phi)
          m[t] <- max(yphi + Xbetaphi + Xbeta[t], minimum)
          y[t] <- rpois(1, m[t])
        }
      }
    }
    if (family == "negative.binomial") {
      for (t in 1:n) {
        if (t < (p + 1)) {
          m[t] <- Xbeta[t]
          y[t] <- rnegbin(1, m[t], theta[t])
        }
        else {
          yphi <- t(y[(t - p):(t - 1)]) %*% rev(phi)
          Xbetaphi <- t(-Xbeta[(t - p):(t - 1)]) %*%
            rev(phi)
          m[t] <- max(yphi + Xbetaphi + Xbeta[t], minimum)
          y[t] <- rnegbin(1, m[t], theta[t])
        }
      }
    }
  }
  if (link == "log") {
    if (family == "gaussian") {
      stop(cat("Gaussian model with log link is not implemented, "))
    }
    if (family == "poisson") {
      if (zero.correction == "zq1") {
        for (t in 1:n) {
          if (t < (p + 1)) {
            m[t] <- Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rpois(1, lambda[t])
            y.transformed[t] <- max(c, y[t])
          }
          else {
            yphi <- t(log(y.transformed[(t - p):(t -
                                                   1)])) %*% rev(phi)
            Xbetaphi <- t(-Xbeta[(t - p):(t - 1)]) %*%
              rev(phi)
            m[t] <- yphi + Xbetaphi + Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rpois(1, lambda[t])
            y.transformed[t] <- max(c, y[t])
          }
        }
      }
      if (zero.correction == "zq2") {
        for (t in 1:n) {
          if (t < (p + 1)) {
            m[t] <- Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rpois(1, lambda[t])
            y.transformed[t] <- y[t] + c
          }
          else {
            yphi <- t(log(y.transformed[(t - p):(t -
                                                   1)])) %*% rev(phi)
            Xbetaphi <- t(-log(exp(Xbeta[(t - p):(t -
                                                    1)]) + c)) %*% rev(phi)
            m[t] <- yphi + Xbetaphi + Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rpois(1, lambda[t])
            y.transformed[t] <- y[t] + c
          }
        }
      }
    }
    if (family == "negative.binomial") {
      if (zero.correction == "zq1") {
        for (t in 1:n) {
          if (t < (p + 1)) {
            m[t] <- Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rnegbin(1, lambda[t], theta[t])
            y.transformed[t] <- max(c, y[t])
          }
          else {
            yphi <- t(log(y.transformed[(t - p):(t -
                                                   1)])) %*% rev(phi)
            Xbetaphi <- t(-Xbeta[(t - p):(t - 1)]) %*%
              rev(phi)
            m[t] <- yphi + Xbetaphi + Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rnegbin(1, lambda[t], theta[t])
            y.transformed[t] <- max(c, y[t])
          }
        }
      }
      if (zero.correction == "zq2") {
        for (t in 1:n) {
          if (t < (p + 1)) {
            m[t] <- Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rnegbin(1, lambda[t], theta[t])
            y.transformed[t] <- y[t] + c
          }
          else {
            yphi <- t(log(y.transformed[(t - p):(t -
                                                   1)])) %*% rev(phi)
            Xbetaphi <- t(-log(exp(Xbeta[(t - p):(t -
                                                    1)]) + c)) %*% rev(phi)
            m[t] <- yphi + Xbetaphi + Xbeta[t]
            lambda[t] <- exp(m[t])
            y[t] <- rnegbin(1, lambda[t], theta[t])
            y.transformed[t] <- y[t] + c
          }
        }
      }
    }
  }
  return(y)
}


wrap_simulation_NB <- function(data, SingleGroup,
                               beta0 = NULL,
                               beta1 = NULL,
                               beta2 = NULL,
                               lambda = NULL,
                               sigma = NULL,
                               theta = NULL,
                               mu = NULL,
                               tau = NULL,
                               act.cor = NULL) {
  if (!is.data.frame(data)) {
    message("Please define a design matrix for simulation.")
  } else{
    ### read the design matrix
    q = length(unique(data$GROUP))
    m = length(unique(data$VISIT))

    if (SingleGroup == FALSE) {
      # if ( length(data$ID) > m*length(unique(data$ID)) ){
      message("Samples from the same individual are analysed in different groups.")
      n = length(unique(data$ID))
      nRep = 0
      SGFLAG = 0
    } else {
      message("Individuals are different in different groups.")
      SGFLAG = 1
      GroupInd <- data$GROUP[1]
      n = length(unique(data[data$GROUP == GroupInd, ]$ID))
      nRep = (n*q - length(unique(data$ID)))/(max(q-1, 1))
      eps = nRep/n
      nUnq = n - nRep
      cat("\n", nRep, "subjects appear repeatedly in different groups. \n")
    }

    ### simulate (correlated) random effects for intercepts and slopes
    if (is.null(beta0)) {
      beta0 = 0
      message("Using the default beta0 = 0.")
    } else{
      if (length(beta0) != 1) {
        message("Please define beta0 as a constant of dimention 1.")
      }
    }

    if (is.null(beta1)) {
      b0 = 100 ## set up an effect size
      beta1 = seq(0, b0*(q-1), by=b0)
      message("Using the default beta1 with group effect size = 100.")
    } else {
      if (length(beta1) != q) {
        message("Please define beta1 as a vector of dimention q.")
      }
    }


    if (is.null(beta2)) {
      # beta2 = log(c(1:m))  ## dim = m
      beta2 = 100*(c(1:m))
      message("Using the default beta2 reflecting a linear dependence of time point.")
    } else {
      if (length(beta2) != m) {
        message("Please define beta2 as a vector of dimention m.")
      }
    }

    if (is.null(lambda)) {
      lambda = matrix(0, byrow = TRUE, ncol = m, nrow=q)
      message("Using the default lambda reflecting no interaction between GROUP and VISIT.")
    } else {
      if (!((is.matrix(lambda)) & (nrow(lambda) == q) & (ncol(lambda) == m))) {
        message("Please define lambda as a matrix of dimention (q, m).")
      }
    }


    ### random effect
    if (is.null(mu)) {
      mu = c(0,0)
      message("Using the default mu = (0, 0).")
    } else {
      if (length(mu) != 2) {
        message("Please define mu as a vector of dimention 2.")
      }
    }

    if (is.null(tau)) {
      tau =   c(1.2, 1.6, 0.3)
      message("Using the default tau0 = (1.2, 1.6) and tau01 = 0.3 to generate the variance matrix of the random intercept and slope.")
    } else {
      if (length(tau) != 3) {
        message("Please define tau as a vector of dimention 3 containing the diagonal and lower triangle element of the variance matrix.")
      }
    }

    ## n_tot*random effect (intercept and slope), alpha is dim=2
    tau0 <- tau[1:2]
    tau01 = tau[3]
    S   = matrix(c(1, tau01, tau01, 1), nrow=2)
    S   = diag(tau0) %*% S %*% diag(tau0)

    ## total number of individuals in different scenarios
    n_tot <- length(unique(data$ID))
    U = mvrnorm((n_tot + nRep*(q-1)), mu=mu, Sigma=S)
    ## Work out the formula  of the sigma^2 within and between subjects!!!

    ### correlation structure
    ### TODO: give the values and plot!!!!
    if (is.null(sigma)) {
      sigma = 1.5
      message("Using the default sigma = 1.5.")
    } else {
      if (length(sigma) != 1) {
        message("Please define sigma as a constant of dimention 1.")
      }
    }

    if (is.null(act.cor)) {
      act.cor = 0.6
      message("Using the default act.cor = 0.6 for the AR correlation structure.")
    } else {
      if (length(act.cor) != 1) {
        message("Please define act.cor as a constant of dimention 1.")
      }
    }

    ### simulate AR(1) errors and then the actual outcomes
    ### true autocorrelation
    ar.val <- act.cor
    # sqrt(1-ar.val^2)
    # ### ARIMASIM <- arima.sim(model=list(ar=ar.val), n=m) * sqrt(1-ar.val^2)* sigma

    if (SGFLAG==1) {
      p <- round(runif((n_tot + nRep*(q-1)),m,m))
      # X <- matrix(1, nrow = m, ncol = 1)
      X = 1
      sigma <- sqrt(beta2)+100
      theta <- beta2^2/(sigma^2-beta2)
      #VecEij <- sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2)* sigma )  ### ~ N(0, sigma^2)
      VecEij <- sapply(p, function(x) garsim_theta(n=x,
                                                   phi=ar.val,
                                                   X = X,
                                                   beta = beta2,
                                                   theta = theta,
                                                   family= "negative.binomial",
                                                   zero.correction = "zq1")) ##
      ## X is needed and associated with beta
      ## X is of dimension n by m, where n is number of samples, m is number of repeated measures.
      ## phi is the AR coef
      ## beta is the mean vector
      ## ~ NB() ## remove sigma and calculate theta using (beta, sigma^2)
      ## theta = beta^2/(sigma^2-beta)
      ## sqrt(1-ar.val^2) is needed or not?
      if (nRep>0){
        for (nSub in 1:(n_tot + nRep*(q-1)) ) {
          if ((nSub %% n) %in% (1:nRep)){
            VecEij[, nSub] <- VecEij[, (nSub %% n)]
          }
        }
      }
      # data$eij <- c(unlist(VecEij))
      data$Exprs <- c(unlist(VecEij))
      y_cell <- array(0, c(m,n,q))
      for (i in 1:q) {
        for (j in 1:m) {
          for (k in 1:n){
            k_tot = n*(i-1) + k ### real individual index
            # y_cell[j, k, i] <- beta0 + U[k_tot, 1]  + (beta1[i] + U[k_tot, 2] )*( beta2[j] ) +  lambda[i, j]
            y_cell[j, k, i] <- beta0 + beta1[i]
          }
        }
      }
      data$Exprs <- data$Exprs + c(y_cell)

      ### model with random intercept and slope
      ### data$Exprs <- (beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p)) * GROUP *log(data$VISIT) + data$eij
      # y_cell <- array(0, c(m,n,q))
      # for (i in 1:q) {
      #  for (j in 1:m) {
      #   for (k in 1:n){
      #     k_tot = n*(i-1) + k ### real individual index
      #     # y_cell[j, k, i] <- beta0 + U[k_tot, 1]  + (beta1[i] + U[k_tot, 2] )*( beta2[j] ) +  lambda[i, j]
      #     y_cell[j, k, i] <- beta0 + U[k_tot, 1]  + (beta1[i] + U[k_tot, 2] )*( beta2[j] ) +  lambda[i, j]
      #   }
      #  }
      # }
    } else {
      p <- round(runif(n,m,m))
      # X <- matrix(1, nrow = m, ncol = 1)
      X = 1
      sigma <- sqrt(beta2)+100
      theta <- beta2^2/(sigma^2-beta2)
      # theta <- theta[m]
      Group_seq <- unique(data$GROUP)
      for (nG in 1:(length(Group_seq))) {
        # browser()
        G_info <- Group_seq[nG]
        #data$eij[data$GROUP == G_info] <- c(unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2)* sigma )))
        data$Exprs[data$GROUP == G_info] <- beta0 + beta1[nG] + c(
          unlist(sapply(p, function(x) garsim_theta(n=x,
                                                    phi=ar.val,
                                                    X = X,
                                                    beta = beta2,
                                                    theta = theta,
                                                    family= "negative.binomial",
                                                    zero.correction = "zq1"))))
      }

    }

    # data$Exprs <- data$eij + c(y_cell)

    return(data)
  }
}
