

#### ---- function 1   ---- 
#### ---- function 1   ---- 

# input: data imported with `geteegdata()`
# output: list, with one element for each trial
#           \_ each element is a list of voltage data matrices (multivariate time serie, for each sensor) for each subject

crea_lista_vltmtx_TRIAL <- function(trial, dataset) {
  dati <- na.omit(dataset[dataset$trial == trial, ])
  
  split_sogg <- split(dati, dati$subject)
  
  res <- lapply(split_sogg, function(a) {
    a$channel <- factor(a$channel)
    
    vltgmtx <- matrix(
      a$voltage,
      nrow = nlevels(a$channel),
      byrow = TRUE
    )
    
    rn <- levels(a$channel)
    rownames(vltgmtx) <- rn
    colnames(vltgmtx) <- 0:255
    
    vltgmtx[!rn %in% c("X", "Y", "nd"), , drop = FALSE]
  })
  
  res
}


  
#### ---- function 2   ---- 
#### ---- function 2   ---- 
  
Wgg.stattest.2smpls <- function(l1, l2) {
  A1 <- simplify2array(l1)
  A2 <- simplify2array(l2)
  
  diff_mean <- apply(A1 - A2, c(1, 2), mean)
  Wgg <- sum(abs(diff_mean)[upper.tri(diff_mean)])
  
  structure(Wgg, names = "Wgg")
}



#### ---- function 3   ---- 
#### ---- function 3   ---- 
  
bootsXeeg <- function(set.grafi1, set.grafi2, numero.repliche) {
  B <- numero.repliche
  congiunto <- c(set.grafi1, set.grafi2)
  n1 <- length(set.grafi1)
  n2 <- length(set.grafi2)
  
  bb <- vapply(seq_len(B), function(i) {
    Wgg.stattest.2smpls(
      sample(congiunto, n1, replace = TRUE),
      sample(congiunto, n2, replace = TRUE)
    )
  }, numeric(1))
  
  Wgg.oss <- Wgg.stattest.2smpls(set.grafi1, set.grafi2)
  pval <- mean(bb > Wgg.oss)
  
  structure(pval, names = "nonparam.pvalue")
}



#### ---- function 4   ---- 
#### ---- function 4   ---- 

crea.Lfin.Lsogg <- function(Lsogg, quante.finestre) {
  require(stringi)
  
  Lfin <- vector("list", quante.finestre)
  names(Lfin) <- stri_c("fin_", seq_len(quante.finestre))
  
  for (j in seq_along(Lsogg)) {
    M <- Lsogg[[j]]
    cuts <- floor(seq(0, ncol(M), length.out = quante.finestre + 1))
    
    for (k in seq_len(quante.finestre)) {
      Lfin[[k]][[j]] <- M[, (cuts[k] + 1):cuts[k + 1], drop = FALSE]
      names(Lfin[[k]])[j] <- names(Lsogg)[j]
    }
  }
  
  Lfin
}

#### ---- function 5   ---- 
#### ---- function 5   ---- 

mediacamp.matrice.V4 <- function(Lsogg) {
  Reduce(`+`, Lsogg) / length(Lsogg)
}

#### ---- function 6   ---- 
#### ---- function 6   ---- 

listaDIliste.mt.cor2adj.V4.opzno05.VERS2.BIS <- function(
    L.obj.corr,
    tipo.di.soglia,
    valore.soglia.fissa = 0.5,
    solo.signif,
    soglia.pval
) {
  R <- L.obj.corr$r
  P <- L.obj.corr$p
  
  ut <- upper.tri(R)
  q <- quantile(R[ut], c(.25, .75))
  
  soglie <- if (tipo.di.soglia == "fissa") {
    c(-valore.soglia.fissa, valore.soglia.fissa)
  } else {
    c(min(-valore.soglia.fissa, q[1]),
      max(valore.soglia.fissa, q[2]))
  }
  
  adj <- (R < soglie[1] | R > soglie[2])
  
  if (solo.signif) {
    adj <- adj & (P < soglia.pval)
  }
  
  storage.mode(adj) <- "numeric"
  dimnames(adj) <- dimnames(R)
  
  adj
}



#### ---- MAIN function ----
#### ---- MAIN function ----
#  was processo.2st!

fn_MAIN <- function(Lsogg.A,
                    Lsogg.C,
                    nr.win = 10,
                    meth = "spearman",
                    correction = "holm",
                    cutoff_type = "entrambi",
                    value_fixed_cutoff = .5,
                    only_signif_TF = TRUE,
                    signif_pvalue = .05,
                    nr_replicas = 1000) {

    
  # # 1) Cut time series in non-overlapping windows:
  
  Lfin.Lsogg.C <- crea.Lfin.Lsogg(Lsogg = Lsogg.C, quante.finestre = nr.win)
  Lfin.Lsogg.A <- crea.Lfin.Lsogg(Lsogg = Lsogg.A, quante.finestre = nr.win)
  
  
  # # 2) Preprocess - for now empty
  
  
  # # 3) Mean matrices of the series
  
  # 3.a - calculate mean
  
  Lfin.Mcamp.C <- lapply(X = Lfin.Lsogg.C, FUN = mediacamp.matrice.V4)
  Lfin.Mcamp.A <- lapply(X = Lfin.Lsogg.A, FUN = mediacamp.matrice.V4)
  
  
  # 3.b - Transpose
  
  T.Lfin.Mcamp.C <- lapply(X = Lfin.Mcamp.C, FUN = t)
  T.Lfin.Mcamp.A <- lapply(X = Lfin.Mcamp.A, FUN = t)
  
  
  # # 4) From mean matrices to correlations (w/ multiple tests correction)
  
  
  Lfin.corr.C <- lapply(X = T.Lfin.Mcamp.C, FUN =  
                          corr.test, method = meth, adjust = correction, ci = FALSE)
  Lfin.corr.A <- lapply(X = T.Lfin.Mcamp.A, FUN = 
                          corr.test, method = meth, adjust = correction, ci = FALSE)
  
  
  # # 5) from correlatio matrices to adjacency matrices
  
  
  Lfin.adj.C <- lapply(X = Lfin.corr.C, FUN = listaDIliste.mt.cor2adj.V4.opzno05.VERS2.BIS,
                       tipo.di.soglia = cutoff_type, valore.soglia.fissa = value_fixed_cutoff,
                       solo.signif = only_signif_TF, soglia.pval = signif_pvalue)
  Lfin.adj.A <- lapply(X = Lfin.corr.A, FUN = listaDIliste.mt.cor2adj.V4.opzno05.VERS2.BIS, 
                       tipo.di.soglia = cutoff_type, valore.soglia.fissa = value_fixed_cutoff,
                       solo.signif = only_signif_TF, soglia.pval = signif_pvalue)
  
  
  # # 6) Calculate stat-test W(g,g')
  
  Wgg.oss <- Wgg.stattest.2smpls(l1 = Lfin.adj.C, l2 = Lfin.adj.A)
  
  # # # --- 
  
  OUTPUT <- vector(mode = 'list', length = 2)
  names(OUTPUT) <- c("Wgg_obs", "estimated_pvalue")
  OUTPUT[[1]] <- Wgg.oss
  
  
  # # 7) P-value bootstrap estimation
  
  stima.p.value <- bootsXeeg(set.grafi1 = Lfin.adj.C, set.grafi2 = Lfin.adj.A, numero.repliche = nr_replicas)
  OUTPUT[[2]] <- stima.p.value
  print( '=============' )
  print( OUTPUT )
  

  return(OUTPUT)
}