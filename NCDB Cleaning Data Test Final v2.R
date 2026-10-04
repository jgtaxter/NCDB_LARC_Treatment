# set working directory
setwd("/Users/juliana.taxter/Downloads/NCDB Data")

# read data from the NCDB file
data <- readLines("NCDBPUF_Rectum.0.2023.0.dat")

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


test_data <- data[1:1000]

ncdb_test <- as.data.frame(
  lapply(seq_len(nrow(lookup)), function(i) {
    substr(
      test_data,
      lookup$Start[i],
      lookup$End[i]
    )
  }),
  stringsAsFactors = FALSE
)

names(ncdb_test) <- lookup$PUF.Item.Name

# Age (above 18) AGE

ncdb_test_age <- ncdb_test[as.integer(ncdb_test$"AGE") >= 18 & as.integer(ncdb_test$"AGE") < 150, ]

# Patient was diagnosed with adenocarcinoma
histology <- c(8480, 8481, 8255, 8262, 8263, 8574, 8140, 8144, 8210, 8211, 8244) # review
ncdb_test_hist <- ncdb_test_age[ncdb_test_age$"HISTOLOGY" %in% histology, ]

# Patient survived for more than 30 days after surgery
ncdb_test_mortality <- ncdb_test_hist[ncdb_test_hist$"PUF_30_DAY_MORT_CD" == "0", ]

# Patient had surgery
ncdb_test_surg <- ncdb_test_mortality[ncdb_test_mortality$"REASON_FOR_NO_SURGERY" == "0", ]

# Patient has APR, LAR, or TME (not TAE)
codes_old <- c(30, 50)
codes_new <- c("A300", "A500")
ncdb_test_surgtype <- ncdb_test_surg[ncdb_test_surg$"RX_SUMM_SURG_PRIM_SITE" %in% codes_old |
                                       ncdb_test_surg$"RX_SUMM_SURG_PRIM_SITE_2023" %in% codes_new, ]

# Patient did not have neoadjuvant therapy of any kind
ncdb_test_adjuvant <- ncdb_test_surgtype[ncdb_test_surgtype$"RX_SUMM_SYSTEMIC_SUR_SEQ" %in% c("0", "3"), ]

# Patient did not have neoadjuvant radiation
ncdb_test_adjrad <- ncdb_test_adjuvant[ncdb_test_adjuvant$"RX_SUMM_SURGRAD_SEQ" %in% c("0", "3"), ]


# Patient did not have any metastisis at diagnosis
cM <- c(unique(ncdb$"TNM_CLIN_M"), unique(ncdb$"AJCC_TNM_CLIN_M"))
cM0 <- cM[!grepl("^\\s*cM?0", cM)]
pM <- c(unique(ncdb$"TNM_PATH_M"), unique(ncdb$"AJCC_TNM_PATH_M"))
pM0 <- pM[!grepl("^\\s*pM?0", pM)]
ncdb_test_mets <- ncdb_test_adjrad[(ncdb_test_adjrad$"TNM_CLIN_M" %in% cM0 |
                                 ncdb_test_adjrad$"AJCC_TNM_CLIN_M" %in% cM0) &
                                 (ncdb_test_adjrad$"TNM_PATH_M" %in% pM0 |
                                 ncdb_test_adjrad$"AJCC_TNM_PATH_M" %in% pM0), ]

# Patient was clinically staged T0/T1/T2)
cT <- c(unique(ncdb$"TNM_CLIN_T"), unique(ncdb$"AJCC_TNM_CLIN_T"))
cT_012 <- cT[grepl("^\\s*cT?[012]", cT)]
ncdb_test_cT <- ncdb_test_mets[ncdb_test_mets$"TNM_CLIN_T" %in% cT_012 |
                                 ncdb_test_mets$"AJCC_TNM_CLIN_T" %in% cT_012, ]

# Patient was clinically staged N0
cN <- c(unique(ncdb$"TNM_CLIN_N"), unique(ncdb$"AJCC_TNM_CLIN_N"))
cN0 <- cN[grepl("^\\s*cN?0", cN)]
ncdb_test_cN <- ncdb_test_cT[ncdb_test_cT$"TNM_CLIN_N" %in% cN0 |
                               ncdb_test_cT$"AJCC_TNM_CLIN_N" %in% cN0, ]

# Filter out unknown pathology for T and N stages
ncdb_test_tstage <- ncdb_test_cN[ncdb_test_cN$"TNM_PATH_T" != "   pX" &
                                   ncdb_test_cN$"AJCC_TNM_PATH_T" != "pTX            ", ]
ncdb_test_nstage <- ncdb_test_tstage[ncdb_test_tstage$"TNM_PATH_N" != "   pX" &
                                       ncdb_test_tstage$"AJCC_TNM_PATH_N" != "pNX            ", ]

# Patient N stage was upstaged pathologically

# Create pathology column to simplify upstaging categorization
ncdb_test_nstage$"path_N_numeric" <- NA_integer_
pN <- c(unique(ncdb$"TNM_PATH_N"), unique(ncdb$"AJCC_TNM_PATH_N"))

pN0 <- pN[grepl("^\\s*pN?0", pN)]
ncdb_test_nstage$"path_N_numeric"[ncdb_test_nstage$"TNM_PATH_N" %in% pN0 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_N" %in% pN0] <- 0

pN1 <- pN[grepl("^\\s*pN?1", pN)]
ncdb_test_nstage$"path_N_numeric"[ncdb_test_nstage$"TNM_PATH_N" %in% pN1 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_N" %in% pN1] <- 1

pN2 <- pN[grepl("^\\s*pN?2", pN)]
ncdb_test_nstage$"path_N_numeric"[ncdb_test_nstage$"TNM_PATH_N" %in% pN2 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_N" %in% pN2] <- 2



ncdb_test_nstage$"path_T_numeric" <- NA_integer_
pT <- c(unique(ncdb$"TNM_PATH_T"), unique(ncdb$"AJCC_TNM_PATH_T"))

pT0 <- pT[grepl("^\\s*pT?0", pT)]
ncdb_test_nstage$"path_T_numeric"[ncdb_test_nstage$"TNM_PATH_T" %in% pT0 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_T" %in% pT0] <- 0

pT1 <- pT[grepl("^\\s*pT?1", pT)]
ncdb_test_nstage$"path_T_numeric"[ncdb_test_nstage$"TNM_PATH_T" %in% pT1 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_T" %in% pT1] <- 1

pT2 <- pT[grepl("^\\s*pT?2", pT)]
ncdb_test_nstage$"path_T_numeric"[ncdb_test_nstage$"TNM_PATH_T" %in% pT2 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_T" %in% pT2] <- 2

pT3 <- pT[grepl("^\\s*pT?3", pT)]
ncdb_test_nstage$"path_T_numeric"[ncdb_test_nstage$"TNM_PATH_T" %in% pT3 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_T" %in% pT3] <- 3

pT4 <- pT[grepl("^\\s*pT?4", pT)]
ncdb_test_nstage$"path_T_numeric"[ncdb_test_nstage$"TNM_PATH_T" %in% pT4 |
                                    ncdb_test_nstage$"AJCC_TNM_PATH_T" %in% pT4] <- 4




# Create clinical column to simplify upstaging categorization
ncdb_test_nstage$"clin_N_numeric" <- NA_integer_

ncdb_test_nstage$"clin_N_numeric"[ncdb_test_nstage$"TNM_CLIN_N" %in% cN0 |
                                    ncdb_test_nstage$"AJCC_TNM_CLIN_N" %in% cN0] <- 0


ncdb_test_nstage$"clin_T_numeric" <- NA_integer_

cT0 <- cT[grepl("^\\s*cT?0", cT)]
ncdb_test_nstage$"clin_T_numeric"[ncdb_test_nstage$"TNM_CLIN_T" %in% cT0 |
                                    ncdb_test_nstage$"AJCC_TNM_CLIN_T" %in% cT0] <- 0

cT1 <- cT[grepl("^\\s*cT?1", cT)]
ncdb_test_nstage$"clin_T_numeric"[ncdb_test_nstage$"TNM_CLIN_T" %in% cT1 |
                                    ncdb_test_nstage$"AJCC_TNM_CLIN_T" %in% cT1] <- 1

cT2 <- cT[grepl("^\\s*cT?2", cT)]
ncdb_test_nstage$"clin_T_numeric"[ncdb_test_nstage$"TNM_CLIN_T" %in% cT2 |
                                    ncdb_test_nstage$"AJCC_TNM_CLIN_T" %in% cT2] <- 2



# Identify which patient's were upstaged from cN0 to cN>0
ncdb_test_nstage$"upstaged" <- ncdb_test_nstage$"path_N_numeric" > ncdb_test_nstage$"clin_N_numeric" |
  (ncdb_test_nstage$"path_T_numeric" > ncdb_test_nstage$"clin_T_numeric" &
     ncdb_test_nstage$"path_T_numeric" > 2)

ncdb_test_N_upstaged <- ncdb_test_nstage[ncdb_test_nstage$"upstaged" == TRUE, ]

