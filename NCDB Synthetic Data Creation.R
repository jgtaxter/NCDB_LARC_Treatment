# ==================================================================
# make_synthetic_ncdb.R
# Builds a FAKE NCDB-style fixed-width .dat file with the same column layout
# as the real PUF, so the cohort and analysis code can run without the real
# data. Every value is randomly generated. No real patient data is read or
# used, only the column positions from the QuickStart layout file.
#
# All values are random, so any results from this file are meaningless.
# It is only for checking that code runs.
# ==================================================================

# ---- Settings -----------------------------------------------------
lookup_file <- "NCDB PUF QuickStart Data Structure 2023.csv"   # layout file
out_dir     <- "synthetic"
out_file    <- file.path(out_dir, "NCDBPUF_Rectum_SYNTHETIC.dat")
n           <- 5000
set.seed(2026)

# ---- Column layout ------------------------------------------------
lookup <- read.csv(lookup_file)
lookup <- lookup[!is.na(lookup$End), c("PUF.Item.Name", "Start", "End")]
lookup$Start <- as.numeric(lookup$Start)
lookup$End   <- as.numeric(lookup$End)

pick <- function(x, prob = NULL) sample(x, n, replace = TRUE, prob = prob)
blank_if <- function(cond, x) ifelse(cond, "", x)

# ---- Fake patients ------------------------------------------------
year     <- pick(2004:2023)
new_tnm  <- year >= 2018                 # AJCC 8 fields; older years use TNM_* fields
site_era <- year >= 2016                 # METS_AT_DX_* site fields
cs_site  <- year >= 2010 & year <= 2015  # CS_METS_DX_* site fields
cs_era   <- year <= 2015                 # CS_METS_AT_DX

# Clinical and pathologic stage (mix of eligible and ineligible patients)
ct <- pick(c("1", "2", "3", "4"), c(.25, .35, .30, .10))
cn <- pick(c("0", "1", "2"),      c(.75, .20, .05))
cm <- pick(c("0", "1"),           c(.93, .07))
pt <- ifelse(ct %in% c("1", "2"),
             pick(c("1", "2", "3", "4A"), c(.20, .40, .35, .05)),
             pick(c("3", "4A", "4B"),     c(.70, .20, .10)))
pn <- pick(c("0", "1A", "1B", "1C", "2A", "2B"), c(.55, .15, .12, .03, .10, .05))

# Node counts consistent with pN (N1c = tumor deposits, 0 positive nodes)
pos  <- ifelse(pn %in% c("0", "1C"), 0L,
               ifelse(pn == "1A", 1L,
                      ifelse(pn == "1B", pick(2:3), pick(4:10))))
exam <- pmax(pos, pick(6:30))

# Treatment
neo       <- pick(c(TRUE, FALSE), c(.40, .60))   # some neoadjuvant (should be excluded)
got_chemo <- pick(c(TRUE, FALSE), c(.55, .45))
got_rad   <- pick(c(TRUE, FALSE), c(.35, .65))
got_imm   <- pick(c(TRUE, FALSE), c(.04, .96))

# Survival
months <- pmin(rexp(n, 1 / 70), 180)
dead   <- months < 180 & runif(n) < 0.4

fields <- list(
  PRIMARY_SITE          = rep("C209", n),
  YEAR_OF_DIAGNOSIS     = as.character(year),
  AGE                   = sprintf("%03d", pick(25:90)),
  SEX                   = pick(c("1", "2")),
  CDCC_TOTAL_BEST       = pick(c("0", "1", "2", "3"), c(.70, .20, .07, .03)),
  HISTOLOGY             = pick(c("8140", "8480", "8210", "8490", "8246"),
                               c(.80, .08, .05, .04, .03)),
  PUF_30_DAY_MORT_CD    = pick(c("0", "1", "9"), c(.95, .03, .02)),
  REASON_FOR_NO_SURGERY = pick(c("0", "1"), c(.95, .05)),
  RX_SUMM_SURG_PRIM_SITE      = blank_if(year == 2023,
                                         pick(c("27", "30", "40", "50", "60", "70", "80"),
                                              c(.15, .45, .08, .20, .04, .04, .04))),
  RX_SUMM_SURG_PRIM_SITE_2023 = blank_if(year != 2023,
                                         pick(c("A270", "A300", "A400", "A500"),
                                              c(.15, .55, .10, .20))),
  RX_SUMM_SURG_OTH_REGDIS     = ifelse(cm == "1", pick(c("0", "4")), "0"),
  
  # TNM: old fields right-justified with c/p prefix, AJCC 8 left-justified
  TNM_CLIN_T      = blank_if(new_tnm, paste0("c", ct)),
  TNM_CLIN_N      = blank_if(new_tnm, paste0("c", cn)),
  TNM_CLIN_M      = blank_if(new_tnm, paste0("c", cm)),
  TNM_PATH_T      = blank_if(new_tnm, paste0("p", pt)),
  TNM_PATH_N      = blank_if(new_tnm, paste0("p", pn)),
  TNM_PATH_M      = blank_if(new_tnm, "88"),
  AJCC_TNM_CLIN_T = blank_if(!new_tnm, paste0("cT", ct)),
  AJCC_TNM_CLIN_N = blank_if(!new_tnm, paste0("cN", cn)),
  AJCC_TNM_CLIN_M = blank_if(!new_tnm, paste0("cM", cm)),
  AJCC_TNM_PATH_T = blank_if(!new_tnm, paste0("pT", tolower(pt))),
  AJCC_TNM_PATH_N = blank_if(!new_tnm, paste0("pN", tolower(pn))),
  AJCC_TNM_PATH_M = blank_if(!new_tnm, "88"),
  ANALYTIC_STAGE_GROUP = ifelse(cm == "1", "4",
                                ifelse(pn != "0", "3",
                                       ifelse(pt %in% c("3", "4A", "4B"), "2", "1"))),
  
  # Metastases (M1 patients get liver mets)
  METS_AT_DX_BONE       = blank_if(!site_era, "0"),
  METS_AT_DX_BRAIN      = blank_if(!site_era, "0"),
  METS_AT_DX_DISTANT_LN = blank_if(!site_era, "0"),
  METS_AT_DX_LIVER      = blank_if(!site_era, ifelse(cm == "1", "1", "0")),
  METS_AT_DX_LUNG       = blank_if(!site_era, "0"),
  METS_AT_DX_OTHER      = blank_if(!site_era, "0"),
  CS_METS_DX_BONE       = blank_if(!cs_site, "0"),
  CS_METS_DX_BRAIN      = blank_if(!cs_site, "0"),
  CS_METS_DX_DISTANT_LN = blank_if(!cs_site, "0"),
  CS_METS_DX_LIVER      = blank_if(!cs_site, ifelse(cm == "1", "1", "0")),
  CS_METS_DX_LUNG       = blank_if(!cs_site, "0"),
  CS_METS_DX_OTHER      = blank_if(!cs_site, "0"),
  CS_METS_AT_DX         = blank_if(!cs_era, ifelse(cm == "1", "40", "00")),
  
  REGIONAL_NODES_POSITIVE = sprintf("%02d", pos),
  REGIONAL_NODES_EXAMINED = sprintf("%02d", exam),
  
  # Treatment sequence: 0 = none, 2 = before surgery, 3 = after surgery
  RX_SUMM_SYSTEMIC_SUR_SEQ = ifelse(got_chemo | got_imm, ifelse(neo, "2", "3"), "0"),
  RX_SUMM_SURGRAD_SEQ      = ifelse(got_rad, ifelse(neo, "2", "3"), "0"),
  RX_SUMM_CHEMO            = ifelse(got_chemo, pick(c("01", "02", "03")),
                                    pick(c("00", "82", "87", "99"), c(.85, .05, .05, .05))),
  RX_SUMM_IMMUNOTHERAPY    = ifelse(got_imm, "01", "00"),
  DX_CHEMO_STARTED_DAYS    = ifelse(got_chemo, as.character(pick(30:120)), ""),
  DX_RAD_STARTED_DAYS      = ifelse(got_rad,   as.character(pick(40:150)), ""),
  DX_IMMUNO_STARTED_DAYS   = ifelse(got_imm,   as.character(pick(30:120)), ""),
  
  DX_LASTCONTACT_DEATH_MONTHS = sprintf("%.2f", months),
  PUF_VITAL_STATUS            = ifelse(dead, "0", "1")   # 0 = dead, 1 = alive
)

# ---- Write fixed-width lines --------------------------------------
skipped <- setdiff(names(fields), lookup$PUF.Item.Name)
if (length(skipped) > 0)
  message("Not in this layout file, left out: ", paste(skipped, collapse = ", "))

lines <- rep(strrep(" ", max(lookup$End)), n)    # every other column stays blank
for (i in seq_len(nrow(lookup))) {
  name <- lookup$PUF.Item.Name[i]
  if (!name %in% names(fields)) next
  w <- lookup$End[i] - lookup$Start[i] + 1
  v <- as.character(fields[[name]])
  v[is.na(v)] <- ""
  v <- substr(v, 1, w)
  # Old TNM_* fields are right-justified ("   c1"); everything else left-justified
  v <- if (startsWith(name, "TNM_")) formatC(v, width = w) else formatC(v, width = -w)
  substr(lines, lookup$Start[i], lookup$End[i]) <- v
}

dir.create(out_dir, showWarnings = FALSE)
writeLines(lines, out_file)
message("Wrote ", n, " synthetic records to ", out_file)

