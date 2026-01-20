
#### ---- Set wd, load libraries ----
#### ---- Set wd, load libraries ----

# setwd(dir = "/media/beoebdue/DATA/R_WD_survival")
setwd(dir = "/media/beoebdue/DATA/R_WD_stage_eeg_2026")  # # <~~ Move everything to here!
require(eegkit)
require(psych) 
require(stringi)

source("stage_eeg_2026_functions.R")

#  https://archive.ics.uci.edu/dataset/121/eeg+database must be downloaded



#### ---- Import csv small dataset ---- 
#### ---- Import csv small dataset ---- 

# create csv format <~~ Do this!
geteegdata(indir = "/media/beoebdue/DATA/R_WD_stage_eeg_2026/eeg+database/SMNI_CMI_TRAIN/",
                           outdir = "/media/beoebdue/DATA/R_WD_stage_eeg_2026/",
                           cond = "S2n",
                           filename = "small_eeg",
                           filetype = ".csv")

small_eeg <- read.table("small_eeg.csv", header = T, sep = ",")



#### ---- Exploratory analysis data small ---- 
#### ---- Exploratory analysis data small ---- 

sbj_tr <- (table(small_eeg$subject, small_eeg$trial) > 0) * 1
sbj_tr_A <- sbj_tr[1:10, ]
sbj_tr_C <- sbj_tr[11:20, ]

# we keep only trial 5 and 19, because we have enough subjects that have performed it both for group A (alcoholic) and group C (controls)
pr_min_trials <- cbind(5, 19)

to_keep <- (sbj_tr[, "5"] * sbj_tr[, "19"]) %>% `>`(0) %>% 
  which %>% names

sm_eeg_5_19 <- do.call(rbind, lapply(to_keep, function(s) small_eeg[small_eeg$subject == s, ]))
sm_eeg_5_19_A <- sm_eeg_5_19[sm_eeg_5_19$group == "a", ]
sm_eeg_5_19_C <- sm_eeg_5_19[sm_eeg_5_19$group == "c", ]



#### ---- Produce list (one elem per trial) of lists (one elem per subject) of voltage data matrices ----
#### ---- Produce list (one elem per trial) of lists (one elem per subject) of voltage data matrices ----

l_trials_l_sogg_A <- apply(X = pr_min_trials, MARGIN = 2, FUN = crea_lista_vltmtx_TRIAL, dataset = sm_eeg_5_19_A)
l_trials_l_sogg_C <- apply(X = pr_min_trials, MARGIN = 2, FUN = crea_lista_vltmtx_TRIAL, dataset = sm_eeg_5_19_C)
names(l_trials_l_sogg_A) <- as.character(pr_min_trials)
names(l_trials_l_sogg_C) <- as.character(pr_min_trials)

# ... but this example focuses on a single trial (to simplify): 
l_sogg_A <- l_trials_l_sogg_A$`5`
l_sogg_C <- l_trials_l_sogg_C$`5`



#### ---- Parameters list ---- 
#### ---- Parameters list ---- 

# PARAMETERS of the MAIN function: Set as desired:
  # Parameter used in each step: 
  # # # #  1)
  #   · nr.win                     # integer - nr of non-overlapping windows in which to cut the data
  # # # #  2) <empty>
  # # # #  3) None
  # # # #  4) 
  #   · meth                       # {"pearson", "spearman", "kendall"} - type of correlations
  #   · correction                 # { "holm", "hochberg", "hommel"} - multiple test correction method for correlations
  # # # #  5)
  #   · cutoff_type                # {"fissa", "entrambi"} - from COR to ADJ, do you want only a fixed cutoff or chosen between a fixed cutoff and a quartile cutoff?
  #   · value_fixed_cutoff         # [0,1] - value of the fixed cutoff
  #   · only_signif_TF             # T/F - filter non-significant correlations?
  #   · signif_pvalue              # [0,1] - cutoff significance of the p-value
  # # # #  6) None
  # # # #  7) 
  #   · nr_replicas                # integer - nr bootstrap replication 



#### --- Execute MAIN function ---- 
#### --- Execute MAIN function ---- 

set.seed(321)
  
fn_MAIN(
    Lsogg.A = l_sogg_A,
    Lsogg.C = l_sogg_C,
    nr.win = 10,
    meth = "spearman",
    correction = "holm",
    cutoff_type = "entrambi",
    value_fixed_cutoff = .5,
    only_signif_TF = TRUE,
    signif_pvalue = .05,
    nr_replicas = 1000)
  


