# read data from the synthetic NCDB file
data <- readLines("NCDBPUF_Rectum_SYNTHETIC.dat")

# find data labels and points
lookup <- read.csv("NCDB PUF QuickStart Data Structure 2023.csv")
#names(lookup)

# only uses necessary columns and make a lookup table and remove empty rows
lookup <- lookup[!is.na(lookup$End), c("PUF.Item.Name","Start","End")]

# makes start and end numeric type
lookup$Start <- as.numeric(lookup$Start)
lookup$End <- as.numeric(lookup$End)

ncdb <- as.data.frame(
  lapply(seq_len(nrow(lookup)), function(i) {
    substr(
      data,
      lookup$Start[i],
      lookup$End[i]
    )
  }),
  stringsAsFactors = FALSE
)

names(ncdb) <- lookup$PUF.Item.Name


# Age (above 18) AGE

ncdb_age <- ncdb[as.integer(ncdb$"AGE") >= 18 & as.integer(ncdb$"AGE") < 150, ]

# Patient survived for more than 30 days after surgery
ncdb_mortality <- ncdb_age[ncdb_age$"PUF_30_DAY_MORT_CD" == "0", ]

# Patient had surgery
ncdb_surg <- ncdb_mortality[ncdb_mortality$"REASON_FOR_NO_SURGERY" == "0", ]

# Patient has APR, LAR, or TME (not TAE)
codes_old <- c(30, 50)
codes_new <- c("A300", "A500")
ncdb_surgtype <- ncdb_surg[ncdb_surg$"RX_SUMM_SURG_PRIM_SITE" %in% codes_old |
                             ncdb_surg$"RX_SUMM_SURG_PRIM_SITE_2023" %in% codes_new, ]

# Patient did not have neoadjuvant therapy of any kind
ncdb_adjuvant <- ncdb_surgtype[ncdb_surgtype$"RX_SUMM_SYSTEMIC_SUR_SEQ" %in% c("0", "3"), ]

# Patient did not have neoadjuvant radiation
ncdb_adjrad <- ncdb_adjuvant[ncdb_adjuvant$"RX_SUMM_SURGRAD_SEQ" %in% c("0", "3"), ]


# Patient did not have any metastisis at diagnosis
ncdb_mets <- ncdb_adjrad[(ncdb_adjrad$"METS_AT_DX_BONE" == 0 &
                            ncdb_adjrad$"METS_AT_DX_BRAIN" == 0 &
                            ncdb_adjrad$"METS_AT_DX_DISTANT_LN" == 0 &
                            ncdb_adjrad$"METS_AT_DX_LIVER" == 0 &
                            ncdb_adjrad$"METS_AT_DX_LUNG" == 0 &
                            ncdb_adjrad$"METS_AT_DX_OTHER" == 0) |
                           (ncdb_adjrad$"CS_METS_DX_BONE" == 0 &
                              ncdb_adjrad$"CS_METS_DX_BRAIN" == 0 &
                              ncdb_adjrad$"CS_METS_DX_LIVER" == 0 &
                              ncdb_adjrad$"CS_METS_DX_LUNG" == 0) |
                           (ncdb_adjrad$"CS_METS_AT_DX" == "00"), ]

# Patient was clinically staged T0/T1/T2)
cT <- c(unique(ncdb$"TNM_CLIN_T"), unique(ncdb$"AJCC_TNM_CLIN_T"))
cT_012 <- cT[grepl("^\\s*cT?[012]", cT)]
ncdb_cT <- ncdb_mets[ncdb_mets$"TNM_CLIN_T" %in% cT_012 |
                       ncdb_mets$"AJCC_TNM_CLIN_T" %in% cT_012, ]

# Patient was clinically staged N0
cN <- c(unique(ncdb$"TNM_CLIN_N"), unique(ncdb$"AJCC_TNM_CLIN_N"))
cN0 <- cN[grepl("^\\s*cN?0", cN)]
ncdb_cN <- ncdb_cT[ncdb_cT$"TNM_CLIN_N" %in% cN0 |
                     ncdb_cT$"AJCC_TNM_CLIN_N" %in% cN0, ]

# Filter out unknown pathology for T and N stages
ncdb_tstage <- ncdb_cN[ncdb_cN$"TNM_PATH_T" != "   pX" &
                         ncdb_cN$"AJCC_TNM_PATH_T" != "pTX            ", ]
ncdb_nstage <- ncdb_tstage[ncdb_tstage$"TNM_PATH_N" != "   pX" &
                             ncdb_tstage$"AJCC_TNM_PATH_N" != "pNX            ", ]

# Patient N stage was upstaged pathologically

# Create pathology column to simplify upstaging categorization
ncdb_nstage$"path_N_numeric" <- NA_integer_
pN <- c(unique(ncdb$"TNM_PATH_N"), unique(ncdb$"AJCC_TNM_PATH_N"))

pN0 <- pN[grepl("^\\s*pN?0", pN)]
ncdb_nstage$"path_N_numeric"[ncdb_nstage$"TNM_PATH_N" %in% pN0 |
                               ncdb_nstage$"AJCC_TNM_PATH_N" %in% pN0] <- 0

pN1 <- pN[grepl("^\\s*pN?1", pN)]
ncdb_nstage$"path_N_numeric"[ncdb_nstage$"TNM_PATH_N" %in% pN1 |
                               ncdb_nstage$"AJCC_TNM_PATH_N" %in% pN1] <- 1

pN2 <- pN[grepl("^\\s*pN?2", pN)]
ncdb_nstage$"path_N_numeric"[ncdb_nstage$"TNM_PATH_N" %in% pN2 |
                               ncdb_nstage$"AJCC_TNM_PATH_N" %in% pN2] <- 2



ncdb_nstage$"path_T_numeric" <- NA_integer_
pT <- c(unique(ncdb$"TNM_PATH_T"), unique(ncdb$"AJCC_TNM_PATH_T"))

pT0 <- pT[grepl("^\\s*pT?0", pT)]
ncdb_nstage$"path_T_numeric"[ncdb_nstage$"TNM_PATH_T" %in% pT0 |
                               ncdb_nstage$"AJCC_TNM_PATH_T" %in% pT0] <- 0

pT1 <- pT[grepl("^\\s*pT?1", pT)]
ncdb_nstage$"path_T_numeric"[ncdb_nstage$"TNM_PATH_T" %in% pT1 |
                               ncdb_nstage$"AJCC_TNM_PATH_T" %in% pT1] <- 1

pT2 <- pT[grepl("^\\s*pT?2", pT)]
ncdb_nstage$"path_T_numeric"[ncdb_nstage$"TNM_PATH_T" %in% pT2 |
                               ncdb_nstage$"AJCC_TNM_PATH_T" %in% pT2] <- 2

pT3 <- pT[grepl("^\\s*pT?3", pT)]
ncdb_nstage$"path_T_numeric"[ncdb_nstage$"TNM_PATH_T" %in% pT3 |
                               ncdb_nstage$"AJCC_TNM_PATH_T" %in% pT3] <- 3

pT4 <- pT[grepl("^\\s*pT?4", pT)]
ncdb_nstage$"path_T_numeric"[ncdb_nstage$"TNM_PATH_T" %in% pT4 |
                               ncdb_nstage$"AJCC_TNM_PATH_T" %in% pT4] <- 4




# Create clinical column to simplify upstaging categorization
ncdb_nstage$"clin_N_numeric" <- NA_integer_

ncdb_nstage$"clin_N_numeric"[ncdb_nstage$"TNM_CLIN_N" %in% cN0 |
                               ncdb_nstage$"AJCC_TNM_CLIN_N" %in% cN0] <- 0


ncdb_nstage$"clin_T_numeric" <- NA_integer_

cT0 <- cT[grepl("^\\s*cT?0", cT)]
ncdb_nstage$"clin_T_numeric"[ncdb_nstage$"TNM_CLIN_T" %in% cT0 |
                               ncdb_nstage$"AJCC_TNM_CLIN_T" %in% cT0] <- 0

cT1 <- cT[grepl("^\\s*cT?1", cT)]
ncdb_nstage$"clin_T_numeric"[ncdb_nstage$"TNM_CLIN_T" %in% cT1 |
                               ncdb_nstage$"AJCC_TNM_CLIN_T" %in% cT1] <- 1

cT2 <- cT[grepl("^\\s*cT?2", cT)]
ncdb_nstage$"clin_T_numeric"[ncdb_nstage$"TNM_CLIN_T" %in% cT2 |
                               ncdb_nstage$"AJCC_TNM_CLIN_T" %in% cT2] <- 2



# Identify which patient's were upstaged from cN0 to pN>0 and from cT1/2 to pT>2
ncdb_nstage$"upstaged" <- ncdb_nstage$"path_N_numeric" > ncdb_nstage$"clin_N_numeric" |
  (ncdb_nstage$"path_T_numeric" > ncdb_nstage$"clin_T_numeric" &
     ncdb_nstage$"path_T_numeric" > 2)

ncdb_upstaged <- ncdb_nstage[ncdb_nstage$"upstaged" %in% TRUE, ]


#testing groups
ncdb_treatment <- ncdb_upstaged
chemo  <- trimws(ncdb_treatment$RX_SUMM_CHEMO)          # 01-03 given; 00, 82-87 not given; 88/99 unknown
immuno <- trimws(ncdb_treatment$RX_SUMM_IMMUNOTHERAPY)  # 01 given; 00, 82-87 not given
rad_after <- trimws(ncdb_treatment$RX_SUMM_SURGRAD_SEQ) == "3"
no_rad    <- trimws(ncdb_treatment$RX_SUMM_SURGRAD_SEQ) == "0"
sys_after <- trimws(ncdb_treatment$RX_SUMM_SYSTEMIC_SUR_SEQ) == "3"
got_chemo <- chemo %in% c("01","02","03");  no_chemo  <- chemo %in% c("00","82","85","86","87")
got_imm   <- immuno == "01";                no_imm    <- immuno %in% c("00","82","85","86","87")

ncdb_treatment$adj_group <- NA_character_
ncdb_treatment$adj_group[no_rad & no_chemo & no_imm]                <- "Surgery alone"
ncdb_treatment$adj_group[got_chemo & sys_after & no_rad & no_imm]   <- "Chemotherapy"
ncdb_treatment$adj_group[rad_after & no_chemo & no_imm]             <- "Radiotherapy"
ncdb_treatment$adj_group[got_chemo & sys_after & rad_after & no_imm] <- "Chemoradiotherapy"
ncdb_treatment$adj_group[got_imm & sys_after]                       <- "Immunotherapy (any)"
table(ncdb_treatment$adj_group, useNA = "ifany")   # NA = unknown or inconsistent; exclude
