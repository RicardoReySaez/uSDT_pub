# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : greenwald2017_extraction.R                                 ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 02-11-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Greenwald, A. G., & De Houwer, J. (2017). Unconscious conditioning:
# Demonstration of existence and difference from conscious conditioning.
# Journal of Experimental Psychology: General, 146(12), 1705-1721.
# https://doi.org/10.1037/xge0000371. Data: https://osf.io/ydeza/, saved in
# ydeza-osfstorage-archive/ with one zip per experiment folder.
#
# This script extracts the trial-level data of the 18 data sets of Series 1
# (unconscious conditioning, 14 experiments) and saves them in one long data
# frame, Data/Greenwald 2017/greenwald2017_long.csv, with one row per trial of
# the two tasks.
#
# The 18 data sets, with the IDs of the Supplement (Table A2 of the article):
#
#   Study  Article   Archive folder   CS ms (acq./test)   N total   N analysed
#   S1     Exp. 1a   UCC3 (2007)      75 / 25                 180           97
#   S2     -         UCC2 (2006)      75 / 25                  78           47
#   S3     Exp. 1b   UCC4 (2007)      75 / 25                  93           85
#   S4     Exp. 1b   UCC4 (2007)      25 / 75                  93           85
#   S5     Exp. 1b   UCC4 (2007)      25 / 25                  93           85
#   S6     -         UCC17b (2010)    33 / 33                  51           43
#   S7     -         UCC1 (2006)      42 / 42                 113           73
#   S8     -         UCC27 (2015)     42 / 42 (*)              75           67
#   S9     -         UCC27 (2015)     42 / 42                  74           69
#   S10    -         UCC18a (2010)    58 / 58                 104           92
#   S11    -         UCC23b (2011)    58 / 58                  56           53
#   S12    -         UCC23c (2012)    58 / 58                  89           82
#   S13    Exp. 1a   UCC3 (2007)      75 / 75                 180           96
#   S14    Exp. 1b   UCC4 (2007)      75 / 75                  93           85
#   S15    -         UCC17b (2010)    75 / 75                  64           51
#   S16    -         UCC24 (2012)     75 / 75                 174          145
#   S17    -         UCC24 (Belg)     75 / 75                  32           32
#   S18    Exp. 1c   UCC26 (2014)     75 / 75                 386          350
#
#   (*) In S8 a 33-ms mask preceded the CS; in S9 it followed it.
#   N total is Table A2's "Total subject N"; the archive differs for S8, S9
#   and S16 (Section 9).
#
# Columns of the long data frame:
#   study       ID of the data set (S1 to S18)
#   experiment  Experiment in the article ("1a", "1b", "1c") or "Supplement"
#   archive     Experiment folder of the OSF archive
#   cs_acq      CS duration in acquisition (ms)
#   cs_test     CS duration in the conditioning and visibility tests (ms)
#   subject     Participant, unique across data sets ("<archive ID>_<number>")
#   task        "indirect" (conditioning test) or "direct" (visibility test)
#   block       Block number in the Inquisit script
#   trial       Trial number within the block
#   test_block  Order of the block within the conditioning test. From 2011
#               on, the authors computed d' in each of its two blocks and
#               averaged them.
#   cs          Response key associated with the CS in acquisition ("L" or
#               "R"). "R" is the signal in the authors' d'.
#   us          Response key required by the US ("L" or "R")
#   response    Key pressed ("L" or "R"; NA for any other key)
#   latency     Response latency (ms)
#   exclusion   "none", or the authors' reason for excluding the trial
#
#
# Each experiment is processed as the authors did in their SPSS syntax.
# Section 9 checks the extraction. From the trials with exclusion == "none",
# it recomputes the authors' d' and compares, for each data set, the number
# of participants, the mean and SD of both d's, and the intercept and slope
# of the regression of conditioning on visibility with the values of the
# authors' SPSS outputs (also in Table A2 and Figure 4 of the article).
#
# Paths are relative to the project root. Open uSDT_pub.Rproj, or set the
# working directory to the folder containing it, before running. If the OSF
# archive is not in ydeza-osfstorage-archive/, the script downloads it
# (Section 3). It unzips the raw data files it needs into
# Data/Greenwald 2017/archive/.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(haven)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Stop if this script is run outside the project root.
if (!file.exists("uSDT_pub.Rproj")) {
  stop("Run this script from the project root (the folder holding ",
       "uSDT_pub.Rproj).", call. = FALSE)
}

# Download and verify the OSF archive if experiment zips are missing.
# osf.io/ydeza; _Archive.zip contains one zip per experiment.
osf_url <- "https://osf.io/download/n9273/"
osf_md5 <- "25cdd4616a90bf7afdfb656bccc6083c"
osf_dir <- "ydeza-osfstorage-archive"

find_zips <- function() {
  list.files(osf_dir, pattern = "^UCC3 \\(2007\\)\\.zip$", recursive = TRUE,
             full.names = TRUE)
}

if (length(find_zips()) == 0) {
  dir.create(osf_dir, showWarnings = FALSE)
  osf_zip <- file.path(osf_dir, "_Archive.zip")
  if (!file.exists(osf_zip) || unname(tools::md5sum(osf_zip)) != osf_md5) {
    options(timeout = max(3600, getOption("timeout")))
    download.file(osf_url, osf_zip, mode = "wb")
  }
  if (unname(tools::md5sum(osf_zip)) != osf_md5) {
    stop("The file downloaded from OSF is not the archive analysed here ",
         "(MD5 ", osf_md5, ").", call. = FALSE)
  }
  utils::unzip(osf_zip, exdir = file.path(osf_dir, "_Archive"))
}

# Paths for the experiment zips, extracted files, and output.
# Manual downloads nest under ydeza-osfstorage-archive/_Archive/_Archive.
zip_dir <- dirname(find_zips()[1])
raw_dir <- "Data/Greenwald 2017/archive"
out_dir <- "Data/Greenwald 2017"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# Extract .dat files and `extra`, skipping files already extracted.
# On Windows, keep the project path under ~110 characters.
unzip_folder <- function(folder, extra = character(0)) {
  zip_file <- file.path(zip_dir, paste0(folder, ".zip"))
  inside   <- utils::unzip(zip_file, list = TRUE)$Name
  wanted   <- inside[grepl("\\.dat$", inside) | basename(inside) %in% extra]
  missing  <- wanted[!file.exists(file.path(raw_dir, wanted))]
  if (length(missing) > 0) {
    utils::unzip(zip_file, files = missing, exdir = raw_dir)
  }
  file.path(raw_dir, folder)
}

# Combine Inquisit files in SPSS order, recording station and file line.
# Remove nontrial rows (those without a numeric subject ID).
read_inquisit <- function(files, station) {
  lapply(seq_along(files), function(i) {
    read.delim(files[i], colClasses = "character", quote = "",
               comment.char = "", na.strings = "") |>
      mutate(station = station[i], line = row_number())
  }) |>
    bind_rows() |>
    filter(grepl("^[0-9]+$", subject)) |>
    transmute(date, time, subject = as.numeric(subject), station, line,
              block = as.numeric(blocknum), trial = as.numeric(trialnum),
              code = trialcode, key = as.numeric(response),
              latency = as.numeric(latency))
}

# Since 2014, SPSS sorted booth files and removed repeated rows.
# Source: 1.DATA.UcsConditioning27.9Apr2015.sps, lines 22-25.
read_inquisit_sorted <- function(file, station) {
  d <- read_inquisit(file, station) |> arrange(subject, block, trial)
  repeated <- d$subject == lag(d$subject) & d$trial == lag(d$trial) &
    d$latency == lag(d$latency) & d$code == lag(d$code, default = "")
  d[!(repeated %in% TRUE), ]
}

# SPSS assigned conditions by subject number (remainder 0 becomes m).
condition_of <- function(subject, m) {
  r <- subject %% m
  ifelse(r == 0, m, r)
}

# Scancodes 18 (E) and 23 (I) indicate left and right responses.
# Other keys count toward criteria but not d' (2010 SPSS, line 121).
key_side <- function(key) case_when(key == 18 ~ "L", key == 23 ~ "R")

# Compare adjacent rows, including participant/task boundaries.
# Apply after SELECT IF filtering to match SPSS's LAG.
add_alternation <- function(d) {
  d |> mutate(alternation = ifelse(key != lag(key), 1, NA))
}

# Behavioural criteria from all trials per participant/task (2007 SPSS).
# alternate: % key changes; mean_key: 20.5 indicates balanced use.
# pct_*: latency bins <250, 250-300, 600-800, and >800 ms.
behaviour_criteria <- function(d) {
  d |>
    group_by(subject, task) |>
    summarise(alternate     = 100 * sum(alternation, na.rm = TRUE) / n(),
              mean_lat      = mean(latency, na.rm = TRUE),
              mean_key      = mean(key, na.rm = TRUE),
              pct_very_fast = 100 * mean(latency < 250, na.rm = TRUE),
              pct_fast      = 100 * mean(latency >= 250 & latency <= 300, na.rm = TRUE),
              pct_slow      = 100 * mean(latency >= 600 & latency <= 800, na.rm = TRUE),
              pct_very_slow = 100 * mean(latency > 800, na.rm = TRUE),
              .groups = "drop")
}

# Pre-2011 task exclusions based on latency, alternation, and key use.
# Source: (2)Intercept.UcsConditioning.3.23Jan07.sps, lines 106-111.
flag_until_2010 <- function(crit) {
  crit |>
    mutate(flagged =
             (task == "indirect" & pct_very_fast + pct_very_slow > 10) |
             (task == "indirect" & (alternate <= 33 | alternate >= 67)) |
             (task == "indirect" & pct_fast + pct_slow > 20) |
             (task == "indirect" & (mean_key < 19.5 | mean_key > 21.5)) |
             (task == "direct" & (mean_key < 19 | mean_key > 22)))
}

# From 2011, the indirect task had 96 trials and relaxed criteria.
# Source: 2.Intercept.UcsConditioning23a.10Oct2011.sps, lines 72-84.
flag_from_2011 <- function(crit) {
  crit |>
    mutate(flagged =
             (task == "indirect" & pct_very_fast + pct_very_slow > 20) |
             (task == "indirect" & (alternate < 25 | alternate > 75)) |
             (task == "indirect" & pct_fast + pct_slow > 25) |
             (task == "indirect" & (mean_key < 19.25 | mean_key > 21.75)) |
             (task == "direct" & (mean_key < 19 | mean_key > 22)))
}

# From 2010, retain marginal tasks pending US-discrimination checks.
# Source: (2a)Intercept.UcsConditioning17b.17May2010.13Sep2014.sps,
# lines 94-102.
add_marginal <- function(crit) {
  crit |>
    mutate(marginal = flagged & mean_lat < 475 & pct_fast + pct_very_fast < 30)
}

# Exclude flagged nonmarginal tasks and indirect trials outside
# 250-800 ms (2007 SPSS syntax, lines 113-129).
exclude_trials <- function(d, crit) {
  if (!"marginal" %in% names(crit)) crit$marginal <- FALSE
  d |>
    left_join(select(crit, subject, task, flagged, marginal),
              by = c("subject", "task")) |>
    mutate(exclusion = case_when(
      !(flagged %in% FALSE) & !(marginal %in% TRUE) ~
        paste0("behaviour criteria (", task, " task)"),
      task == "indirect" & !(latency > 250 & latency < 800) ~
        "latency outside 250-800 ms",
      .default = "none")) |>
    select(-flagged, -marginal)
}

# Authors' edge correction: rates of 0 and 1 become .25/N and
# 1 - .25/N before applying the probit transformation.
dprime_correct <- function(hit, fa, n_signal, n_noise) {
  hit <- case_when(hit == 0 ~ .25 / n_signal, hit == 1 ~ 1 - .25 / n_signal,
                   .default = hit)
  fa  <- case_when(fa == 0 ~ .25 / n_noise, fa == 1 ~ 1 - .25 / n_noise,
                   .default = fa)
  qnorm(hit) - qnorm(fa)
}

# Compute d' for CS (main outcome) or US (exclusion checks).
# cells=TRUE: unweighted mean of CS x US cell rates (2011-2013).
# count_all=TRUE: use all CS trials for N, including other keys,
# even when computing US d' (2006-2007).
# Sources: 2.Intercept.UcsConditioning23a.10Oct2011.sps, lines 219-238;
# UCC3 syntax, lines 164-183.
dprime_by <- function(d, by, stimulus = "cs", cells = FALSE, count_all = FALSE) {
  if (!"response" %in% names(d)) d$response <- key_side(d$key)
  d |>
    mutate(signal = .data[[stimulus]] == "R",
           yes    = ifelse(is.na(response), NA, response == "R"),
           cell   = if (cells) paste(cs, us) else "all") |>
    group_by(across(all_of(c(by, "cell")))) |>
    summarise(n_signal = if (count_all) sum(cs == "R", na.rm = TRUE)
              else sum(signal & !is.na(yes), na.rm = TRUE),
              n_noise  = if (count_all) sum(cs == "L", na.rm = TRUE)
              else sum(!signal & !is.na(yes), na.rm = TRUE),
              hit      = mean(yes[signal], na.rm = TRUE),
              fa       = mean(yes[!signal], na.rm = TRUE),
              .groups  = "drop") |>
    group_by(across(all_of(by))) |>
    summarise(dprime = dprime_correct(mean(hit, na.rm = TRUE), mean(fa, na.rm = TRUE),
                                      sum(n_signal), sum(n_noise)),
              .groups = "drop") |>
    mutate(dprime = ifelse(is.nan(dprime), NA, dprime))
}

# Map study identifiers to article experiments.
experiment_of <- function(study) {
  case_when(study %in% c("S1", "S13") ~ "1a",
            study %in% c("S3", "S4", "S5", "S14") ~ "1b",
            study == "S18" ~ "1c",
            .default = "Supplement")
}

# Standard trial-level output columns (Section 1).
as_long <- function(d, archive, id) {
  d |>
    group_by(subject, task) |>
    mutate(test_block = ifelse(task == "indirect",
                               match(block, sort(unique(block))), NA)) |>
    ungroup() |>
    transmute(study, experiment = experiment_of(study), archive = archive,
              cs_acq, cs_test, subject = paste(id, subject, sep = "_"),
              task, block, trial, test_block, cs, us,
              response = key_side(key), latency, exclusion)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Experiments of 2006-2007 (UCC1, UCC2, UCC3, UCC4)
# ─────────────────────────────────────────────────────────────────────────────

# UCC1-UCC4 share the UCC3 analysis syntax
# ((2)Intercept.UcsConditioning.3.23Jan07.sps).
# Blocks 6-11 = indirect; 14-19 = direct; 12-13 = practice.
# Example "i25RLwB": indirect, 25-ms CS, CS=R, US=L, letters=w.
# DATA syntax: lines 111, 122-156.
decode_2006 <- function(d) {
  d |>
    filter(block >= 6, !block %in% c(12, 13, 20)) |>
    mutate(task    = ifelse(substr(code, 1, 1) == "i", "indirect", "direct"),
           cs_test = as.numeric(substr(code, 2, 3)),
           cs      = substr(code, 4, 4),
           us      = substr(code, 5, 5),
           letters = substr(code, 6, 6))
}

# Apply behavioural/latency exclusions (UCC3 syntax, lines 92-129).
# In UCC1-UCC3, a coding error linked CS and US in the direct task.
# Exclude participants responding to the US, using study-specific d'
# cutoffs; exclude those without a valid direct task. In UCC3, check
# 25- and 75-ms CSs separately (UCC3 syntax, lines 4-11, 224-237;
# UCC1/UCC2 syntax: (2b)Intercept.UcsConditioning.7Nov06.sps,
# lines 198-204).
exclusions_2006 <- function(d, us_cutoff = NULL) {
  crit <- d |> add_alternation() |> behaviour_criteria() |> flag_until_2010()
  d <- exclude_trials(d, crit)
  if (is.null(us_cutoff)) return(d)

  us_dprime <- d |>
    filter(task == "direct", exclusion == "none") |>
    dprime_by(by = c("subject", "study"), stimulus = "us", count_all = TRUE)
  d |>
    left_join(us_dprime, by = c("subject", "study")) |>
    mutate(exclusion = case_when(
      is.na(dprime) & task == "indirect" ~ "no valid visibility test",
      dprime > us_cutoff ~ "responded to the US in the visibility test",
      .default = exclusion)) |>
    select(-dprime)
}

# ── UCC1 (2006): S7 ───────────────────────────────────────────────────────────

# UCC1's DATA syntax reads UCC2 files; repeated subject IDs were
# renumbered without documentation. Use the authors' analysed .sav
# file in its original order ((2b)Intercept.UcsConditioning.3Jul06.sps,
# line 82).
path <- unzip_folder("UCC1 (2006)", extra = "UcsConditioning.3Jul06.sav")
ucc1 <- haven::read_sav(file.path(
  path, "Analyses", "UcsConditioning.3Jul06.sav")) |>
  transmute(subject = as.numeric(subject), station = as.character(station),
            block   = as.numeric(block), trial = as.numeric(trialnum),
            key     = as.numeric(response), latency = as.numeric(latency),
            task    = ifelse(task == 1, "indirect", "direct"),
            cs_test = as.numeric(primedur), cs = as.character(prime),
            us      = as.character(target),
            study   = "S7", cs_acq = 42) |>
  exclusions_2006(us_cutoff = 1.86) |>
  as_long(archive = "UCC1 (2006)", id = "UCC1")

# ── UCC2 (2006): S2 ───────────────────────────────────────────────────────────

# Follow the authors' booth order, omitting a file marked "don't use
# this" and excluding duplicate/incomplete IDs and later sessions
# (DATA.UcsConditioning.7Nov06.sps, lines 11, 54-103).
path   <- unzip_folder("UCC2 (2006)")
booths <- c("CW", "DA", "EW", "EA", "FW", "FA", "GW", "GA")
ucc2 <- read_inquisit(file.path(
  path, "Analyses",
  paste0(booths, ".UcsConditioning.7Nov06.dat")), booths) |>
  filter(!subject %in% c(501, 503, 509, 513, 514, 515, 517, 520, 521, 529, 534,
                         552, 562, 589),
         subject <= 593 | subject == 999) |>
  decode_2006() |>
  mutate(study = "S2", cs_acq = 75) |>
  exclusions_2006(us_cutoff = 1.86) |>
  as_long(archive = "UCC2 (2006)", id = "UCC2")

# ── UCC3 (2007): S1 and S13 (Experiment 1a) ───────────────────────────────────

# Acquisition: 75-ms CSs; tests: 25 ms (S1) or 75 ms (S13).
# Drop incomplete/test sessions and renumber duplicate IDs (DATA syntax,
# lines 79-95); also exclude subject 1 (analysis syntax, lines 88-90).
path   <- unzip_folder("UCC3 (2007)")
booths <- c("C", "D", "E", "F", "G")
ucc3 <- read_inquisit(file.path(
  path, "Analyses",
  paste0(booths, ".UcsConditioning.23Jan07.dat")), booths) |>
  filter(!subject %in% c(707, 721, 740, 741, 785, 816, 817, 820, 834, 863, 878,
                         890, 9998, 9999, 1)) |>
  mutate(subject = case_when(subject == 798 & station == "E" ~ 898,
                             subject == 897 & station == "G" ~ 997,
                             .default = subject)) |>
  decode_2006() |>
  mutate(study = ifelse(cs_test == 25, "S1", "S13"), cs_acq = 75) |>
  exclusions_2006(us_cutoff = 1.146) |>
  # Exclude booth C's 36 participants (faulty latency recording;
  # analysis syntax, lines 12-14, 224).
  mutate(exclusion = ifelse(
    station == "C", "apparatus failure (booth C)", exclusion)) |>
  as_long(archive = "UCC3 (2007)", id = "UCC3")

# ── UCC4 (2007): S3, S4, S5 and S14 (Experiment 1b) ───────────────────────────

# CS duration (25/75 ms) varied in acquisition and tests.
# Drop incomplete/test sessions; no US check is needed. Infer acquisition
# duration from condition and letters (DATA syntax, line 60; analysis
# syntax, lines 170-193).
path <- unzip_folder("UCC4 (2007)")
ucc4 <- read_inquisit(file.path(
  path, "Analyses",
  paste0(booths, ".UcsConditioning.4Apr07.dat")), booths) |>
  filter(!subject %in% c(105, 173, 999)) |>
  decode_2006() |>
  mutate(cond   = condition_of(subject, 8),
         cs_acq = case_when(
           cond %in% 1:2 ~ ifelse(letters %in% c("w", "x"), 25, 75),
           cond %in% 3:4 ~ ifelse(letters %in% c("y", "z"), 25, 75),
           cond %in% 5:6 ~ ifelse(letters %in% c("s", "t"), 25, 75),
           cond %in% 7:8 ~ ifelse(letters %in% c("u", "v"), 25, 75)),
         study  = case_when(cs_acq == 75 & cs_test == 25 ~ "S3",
                            cs_acq == 25 & cs_test == 75 ~ "S4",
                            cs_acq == 25 & cs_test == 25 ~ "S5",
                            cs_acq == 75 & cs_test == 75 ~ "S14")) |>
  exclusions_2006() |>
  as_long(archive = "UCC4 (2007)", id = "UCC4")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Experiments of 2010 (UCC17b, UCC18a)
# ─────────────────────────────────────────────────────────────────────────────

# 2010 codes: "i058LvtA" = indirect, 58-ms CS, US=L, letters=vt;
# "d58RLvtA" also specifies CS=R. In indirect trials, the CS side
# depends on condition and letter order (DATA.UcsConditioning17b.
# 17May2010.sps, lines 114-137).
decode_2010 <- function(d) {
  d |>
    mutate(task    = ifelse(substr(code, 1, 1) == "i", "indirect", "direct"),
           cs_test = as.numeric(ifelse(task == "indirect", substr(code, 2, 4),
                                       substr(code, 2, 3))),
           us      = substr(code, 5, 5),
           cs      = ifelse(task == "direct", substr(code, 4, 4), NA),
           letters = substr(code, 6, 7))
}

# Apply the earlier behavioural/latency rules, retaining marginal
# tasks initially. Base exclusion on the first remaining task (indirect
# first). For marginal cases, require indirect-task US d' >= .65.
# Source: (2a)Intercept.UcsConditioning17b.17May2010.13Sep2014.sps,
# lines 73-116, 267-301.
exclusions_2010 <- function(d) {
  crit <- d |> add_alternation() |> behaviour_criteria() |>
    flag_until_2010() |> add_marginal()
  d <- exclude_trials(d, crit)

  first_task <- d |>
    filter(exclusion == "none") |>
    distinct(subject, task) |>
    left_join(crit, by = c("subject", "task")) |>
    arrange(subject, task != "indirect") |>
    slice(1, .by = subject)
  to_check <- first_task |> filter(marginal %in% TRUE) |> pull(subject)
  good_us  <- d |>
    filter(subject %in% to_check, task == "indirect", exclusion == "none") |>
    dprime_by(by = "subject", stimulus = "us") |>
    filter(dprime >= .65) |>
    pull(subject)
  d |>
    mutate(exclusion = ifelse(subject %in% setdiff(to_check, good_us),
                              "marginal behaviour and poor US discrimination",
                              exclusion))
}

# ── UCC17b (2010): S6 and S15 ─────────────────────────────────────────────────

# CS duration: 33 ms (S6) or 75 ms (S15), in both phases.
# Renumber one duplicate ID, drop incomplete sessions and keep test
# blocks only (DATA syntax, lines 43-55, 74; analysis, lines 177-178).
path <- unzip_folder("UCC17b (2010)")
ucc17b <- read_inquisit(file.path(path, "Analyses",
                                  paste0(booths, ".UcsConditioning17b.17May2010.dat")), booths) |>
  mutate(subject = ifelse(subject == 134 & time == "15:28", 534, subject)) |>
  filter(!(subject %in% c(204, 227) & station == "E"),
         !subject %in% c(168, 176, 251)) |>
  filter(block >= 8, !block %in% c(12, 13, 18)) |>
  decode_2010() |>
  # Infer CS side from condition and letters (analysis, lines 153-161).
  mutate(cond   = condition_of(subject, 8),
         cs     = ifelse(task == "direct", cs, case_when(
           cond %in% c(1, 5) & letters == "vt" ~ "L", cond %in% c(1, 5) & letters == "tv" ~ "R",
           cond %in% c(2, 6) & letters == "tv" ~ "L", cond %in% c(2, 6) & letters == "vt" ~ "R",
           cond %in% c(3, 7) & letters == "zx" ~ "L", cond %in% c(3, 7) & letters == "xz" ~ "R",
           cond %in% c(4, 8) & letters == "xz" ~ "L", cond %in% c(4, 8) & letters == "zx" ~ "R")),
         cs_acq = ifelse(cond <= 4, 75, 33),
         study  = ifelse(cs_acq == 33, "S6", "S15")) |>
  exclusions_2010() |>
  as_long(archive = "UCC17b (2010)", id = "UCC17b")

# ── UCC18a (2010): S10 ────────────────────────────────────────────────────────

# 58-ms CSs. Exclude block 8 (CS exposure/contingency instructions)
# before computing criteria (analysis syntax, line 143).
# Keep indirect blocks 9-12 and direct blocks 15-18 (DATA, line 63).
path <- unzip_folder("UCC18a (2010)")
ucc18a <- read_inquisit(file.path(
  path, "Analyses",
  paste0(booths, ".UcsConditioning18a.23Jun2010.dat")), booths) |>
  filter(block >= 9, !block %in% c(13, 14, 19)) |>
  decode_2010() |>
  # Infer CS side in acquisition (analysis syntax, lines 222-225).
  mutate(cond   = condition_of(subject, 6),
         cs     = ifelse(task == "direct", cs, case_when(
           cond %in% 1:3 & letters == "vt" ~ "L", cond %in% 1:3 & letters == "tv" ~ "R",
           cond %in% 4:6 & letters == "zx" ~ "L", cond %in% 4:6 & letters == "xz" ~ "R")),
         cs_acq = 58, study = "S10") |>
  exclusions_2010() |>
  # Exclude subject 122 for responding to the US in visibility
  # (analysis syntax, lines 349-351).
  mutate(exclusion = ifelse(subject == 122, "responded to the US in the visibility test",
                            exclusion)) |>
  as_long(archive = "UCC18a (2010)", id = "UCC18a")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Experiments of 2011-2013 (UCC23b, UCC23c, UCC24, UCC24 Belgium)
# ─────────────────────────────────────────────────────────────────────────────

# Indirect tests use two 48-trial blocks. Codes encode CS and US details
# (e.g., "LA1unpB_358"); the US key is in position 7 for UCC24
# ("LA1unpRA_075"). Direct codes include both keys ("d58LRA2valB").
# CS mapping depends on the counterbalanced sub-alphabet order.
# Sources: 1rev.DATA.UcsConditioning24.20Feb2012.sps, lines 232-238;
# 1.DATA.UcsConditioning23a.10Oct2011.sps, lines 209-215, 255-261.
decode_2011 <- function(d, us_position) {
  d |>
    mutate(task       = case_when(block %in% c(7, 8) ~ "indirect",
                                  block %in% c(11, 12) ~ "direct",
                                  .default = "acquisition"),
           test_block = case_when(block == 7 ~ 1, block == 8 ~ 2),
           us_set     = ifelse(task == "direct", NA, substr(code, 4, 6)),
           us         = ifelse(task == "direct", substr(code, 5, 5),
                               substr(code, us_position, us_position)),
           cs_n       = substr(code, 3, 3),
           cs         = case_when(
             task == "direct" ~ substr(code, 4, 4),
             cs_n == "1" ~ ifelse(sub_alphabet == "12", "L", "R"),
             cs_n == "2" ~ ifelse(sub_alphabet == "12", "R", "L"))) |>
    # Match SPSS's US-category sort before computing alternation;
    # direct-test trials come first (SORT CASES BY acq_task, line 217).
    arrange(!is.na(us_set), !us_set %in% c("mal", "fem"))
}

# Apply the relaxed criteria across both indirect-test blocks; exclude
# flagged tasks unless marginal, and indirect latencies outside
# 250-800 ms.
# Marginal cases require mean US d' >= .65 across blocks. Alternation
# must already be computed in the authors' trial order.
# Source: 2.Intercept.UcsConditioning23a.10Oct2011.sps, lines 79-155,
# 296-341 (same rules in UCC23c, UCC24 and UCC24 Belgium).
exclusions_2011 <- function(d) {
  crit <- d |> behaviour_criteria() |> flag_from_2011() |> add_marginal()
  d <- exclude_trials(d, crit)

  to_check <- crit |> filter(marginal %in% TRUE) |> pull(subject)
  good_us  <- d |>
    filter(subject %in% to_check, task == "indirect", exclusion == "none") |>
    dprime_by(by = c("subject", "test_block"), stimulus = "us", cells = TRUE) |>
    group_by(subject) |>
    summarise(dprime = mean(dprime, na.rm = TRUE)) |>
    filter(dprime >= .65) |>
    pull(subject)
  d |>
    mutate(exclusion = ifelse(subject %in% setdiff(to_check, good_us),
                              "marginal behaviour and poor US discrimination",
                              exclusion))
}

# ── UCC23b (2011): S11 ────────────────────────────────────────────────────────

# Of 32 conditions varying CS duration, consistency, US category and
# sub-alphabet, S11 includes only 58-ms CSs with 100% consistency.
# Duplicate IDs are fixed; incomplete sessions and six wrong-acquisition
# participants are removed. Acquisition trials in 100% conditions lack
# task codes and are removed before computing alternation.
# Sources: 1.DATA.UcsConditioning23a.10Oct2011.sps, lines 72-86, 90-152,
# 178-191; 2.Intercept.UcsConditioning23a.10Oct2011.sps,
# lines 1-17, 37-39.
path <- unzip_folder("UCC23b (2011)")
ucc23b <- read_inquisit(file.path(path, "Analyses",
                                  paste0(booths, ".UcsConditioning23a.10Oct2011.dat")), booths) |>
  mutate(subject = case_when(
    subject == 1   & station == "D" ~ 385,
    subject == 1   & station == "G" ~ 353,
    subject == 143 & station == "D" ~ 463,
    subject == 283 & station == "E" ~ 603,
    subject == 286 & station == "E" ~ 606,
    (subject == 115 & station == "D") | (subject == 254 & station == "C") |
      (subject == 293 & station == "G") | (subject == 327 & station == "D") |
      (subject == 279 & time == "09:36") ~ 0,
    .default = subject)) |>
  filter(subject != 0, !subject %in% c(103, 104, 107, 108, 111, 112)) |>
  mutate(cond         = condition_of(subject, 32),
         consistency  = ifelse(trunc((cond - .5) / 2) %% 2 == 0, 100, 80),
         cs_acq       = ifelse(cond %% 2 == 0, 58, 358),
         sub_alphabet = ifelse(trunc((cond - .5) / 16) %% 2 == 0, "12", "21")) |>
  filter(block %in% c(1:3, 5:8, 11:12), block >= 7 | consistency == 80) |>
  decode_2011(us_position = 1) |>
  add_alternation() |>
  filter(task != "acquisition", cs_acq == 58, consistency == 100) |>
  mutate(cs_test = cs_acq, study = "S11") |>
  exclusions_2011() |>
  as_long(archive = "UCC23b (2011)", id = "UCC23b")

# ── UCC23c (2012): S12 ────────────────────────────────────────────────────────

# S12 includes only the 100%-consistent conditions (58-ms CSs).
# Incomplete sessions are removed, one ID renumbered, and UCC23b's
# block/coding rules reused.
# Sources: 1.DATA.UcsConditioning23c.11Jan2012.sps, lines 48-50, 55-97,
# 124-140; 2.Intercept.UcsConditioning23c.11Jan2012.sps, lines 4-7.
path <- unzip_folder("UCC23c (2012)")
ucc23c <- read_inquisit(file.path(path, "Analyses",
                                  paste0(booths, ".UcsConditioning23c.11Jan2012.dat")), booths) |>
  filter(!subject %in% c(1, 101, 284)) |>
  mutate(subject = ifelse(subject == 109 & station == "G", 429, subject)) |>
  mutate(cond         = condition_of(subject, 16),
         consistency  = ifelse(trunc(cond - .5) %% 2 == 0, 100, 75),
         sub_alphabet = ifelse(trunc((cond - .5) / 8) %% 2 == 0, "12", "21")) |>
  filter(block %in% c(1:3, 5:8, 11:12), block >= 7 | consistency == 75) |>
  decode_2011(us_position = 1) |>
  add_alternation() |>
  filter(task != "acquisition", consistency == 100) |>
  mutate(cs_acq = 58, cs_test = 58, study = "S12") |>
  exclusions_2011() |>
  as_long(archive = "UCC23c (2012)", id = "UCC23c")

# ── UCC24 (2012): S16 ─────────────────────────────────────────────────────────

# 48 conditions with 75-ms CSs: unchanged CS-US pairs, swapped US keys,
# or new USs. S16 includes only the unchanged-pair conditions.
# Incomplete sessions/duplicate IDs are handled, and a sub-alphabet
# coding error is corrected for conditions 25-48.
# Sources: 1rev.DATA.UcsConditioning24.20Feb2012.sps, lines 59-78,
# 156-177, 192-195, 212; 2.Intercept.UcsConditioning24.20Feb2012.sps,
# lines 4-6.
conditions_24 <- function(d) {
  d |>
    mutate(cond         = condition_of(subject, 48),
           first_half   = trunc((cond + 11.5) / 12) %% 2,
           first_half   = ifelse(cond >= 25, 1 - first_half, first_half),
           sub_alphabet = ifelse(first_half == 1, "12", "21"),
           treatment    = condition_of(trunc(cond + .5), 3))
}

path <- unzip_folder("UCC24 (2012)")
ucc24 <- read_inquisit(file.path(path, "Analyses",
                                 paste0(booths, ".UcsConditioning24.20Feb2012.dat")), booths) |>
  filter(!(subject %in% c(175, 269) & station == "D"),
         !subject %in% c(1, 118, 130, 134, 1036, 1237, 1276)) |>
  mutate(subject = case_when(subject == 158 & station == "E" ~ 638,
                             subject == 291 & station == "E" ~ 771,
                             subject == 339 & station == "D" ~ 819,
                             subject == 1209 & date == "112612" ~ 1689,
                             .default = subject)) |>
  conditions_24() |>
  filter(block %in% c(7, 8, 11, 12)) |>
  decode_2011(us_position = 7) |>
  add_alternation() |>
  filter(treatment == 1) |>
  mutate(cs_acq = 75, cs_test = 75, study = "S16") |>
  exclusions_2011() |>
  as_long(archive = "UCC24 (2012)", id = "UCC24")

# ── UCC24 Belgium (2013-2014): S17 ────────────────────────────────────────────

# Ghent replication of UCC24. Reused IDs are offset by 48 from the
# second participant set in each file. The January file receives +4800
# and is combined first, preserving the authors' alternation order.
# Sources: 1rev.DATA.UcsConditioning24.13Dec2013.sps, lines 10-20;
# 1rev.DATA.UcsConditioning24.Dutch_rev.4Mar2014.sps, lines 15-23;
# 2.Intercept.UcsConditioning24.Dutch_rev.combined.sps, lines 4-9.
path <- unzip_folder("UCC24 (Belg)")
prepare_belgium <- function(d) {
  d |>
    conditions_24() |>
    filter(block %in% c(7, 8, 11, 12)) |>
    decode_2011(us_position = 7)
}
ucc24_oct <- read_inquisit(file.path(path, "Analyses",
                                     "UcsConditioning24.Dutch_16Oct13_corrected.dat"), "B") |>
  mutate(subject = ifelse(line >= 10452, subject + 48, subject),
         subject = ifelse(subject == 78 & line < 21852, 30, subject),
         subject = ifelse(subject == 88 & line < 19477, 40, subject)) |>
  prepare_belgium()
ucc24_jan <- read_inquisit(file.path(path, "Analyses",
                                     "UcsConditioning24.Dutch_20Jan14_corrected(v2).dat"), "B") |>
  mutate(subject = ifelse(line >= 11345, subject + 48, subject)) |>
  prepare_belgium() |>
  mutate(subject = subject + 4800)
ucc24_belg <- bind_rows(ucc24_jan, ucc24_oct) |>
  add_alternation() |>
  filter(treatment == 1) |>
  mutate(cs_acq = 75, cs_test = 75, study = "S17") |>
  exclusions_2011() |>
  as_long(archive = "UCC24 (Belg)", id = "UCC24Belg")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Experiments of 2014-2015 (UCC27, UCC26)
# ─────────────────────────────────────────────────────────────────────────────

# Two counterbalanced acquisition-test sequences and one visibility test;
# only the first sequence is analysed to avoid carryover (1.DATA.UcsConditioning27.
# 9Apr2015.sps, lines 1-9). In "LB1femA_33+42", "B1" is the CS and "fem" the
# US category; female/unpleasant USs map to L, and CS keys are counterbalanced
# (same syntax, lines 214-230).
# Tests were analysed separately using 2011-2013 behavioural and latency rules.
# Marginal tasks were retained without a final US check
# (2.Intercept.UcsConditioning27.9Apr2015.sps, lines 28-159, 410-503, 656).
split_tests_2014 <- function(d, test_blocks, visibility_blocks, cs_position) {
  test <- d |>
    filter(block %in% test_blocks) |>
    mutate(task       = "indirect",
           test_block = match(block, test_blocks),
           us_set     = substr(code, 4, 6),
           us         = ifelse(us_set %in% c("fem", "unp"), "L", "R"),
           cs         = case_when(
             substr(code, 3, 3) == "1" ~ ifelse(sub_alphabet == "12", "L", "R"),
             substr(code, 3, 3) == "2" ~ ifelse(sub_alphabet == "12", "R", "L"))) |>
    # Sort by US category, gender first (line 251)
    arrange(!us_set %in% c("mal", "fem"))
  visibility <- d |>
    filter(block %in% visibility_blocks) |>
    mutate(task = "direct",
           cs   = substr(code, cs_position, cs_position),
           us   = substr(code, cs_position + 1, cs_position + 1))
  bind_rows(test, visibility)
}

# Compute alternation within each test, matching the separate files.
exclusions_2014 <- function(d) {
  crit <- d |>
    group_by(task) |>
    add_alternation() |>
    ungroup() |>
    behaviour_criteria() |>
    flag_from_2011() |>
    add_marginal()
  exclude_trials(d, crit)
}

# Shared UCC26/UCC27 counterbalancing: sub-alphabet order and first procedure.
conditions_2014 <- function(d, first_procedure) {
  d |>
    mutate(cond         = condition_of(subject, 16),
           sub_alphabet = ifelse(trunc((cond + 7.5) / 8) %% 2 == 0, "21", "12"),
           first        = ifelse(trunc((cond + 1.5) / 2) %% 2 == 0,
                                 first_procedure[1], first_procedure[2]))
}

# ── UCC27 (2015): S8 and S9 ───────────────────────────────────────────────────

# 42-ms CS: a 33-ms mask precedes it in S8 ("33+42") or follows it in S9
# ("42+33"); S8 has the only non-75-ms CS-US interval. Booth files are sorted
# and cleaned, excluding five sessions lost to a programming error. First
# test: blocks 6-7; visibility: 19-21 (1.DATA.UcsConditioning27.9Apr2015.sps,
# lines 16-94, 350-378).
path  <- unzip_folder("UCC27 (2015)")
ucc27 <- lapply(booths, function(b) {
  read_inquisit_sorted(file.path(path, "Analyses",
                                 paste0(b, ".UcsConditioning27.9Apr2015.dat")), b)
}) |>
  bind_rows() |>
  arrange(subject, block, trial) |>
  filter(subject >= 105, subject != 129) |>
  conditions_2014(first_procedure = c("33+42", "42+33")) |>
  split_tests_2014(test_blocks = 6:7, visibility_blocks = 19:21, cs_position = 7) |>
  mutate(cs_acq = 42, cs_test = 42, study = ifelse(first == "33+42", "S8", "S9")) |>
  exclusions_2014() |>
  as_long(archive = "UCC27 (2015)", id = "UCC27")

# ── UCC26 (2014): S18 (Experiment 1c) ─────────────────────────────────────────

# Three replications counterbalance 75-ms (masked) and 375-ms (visible)
# sequences. S18 (Experiment 1c) includes only participants with 75 ms first;
# 375 ms first belongs to Experiment 2b (Series 2). Archive folders 26a/b/c
# correspond to data files 26/26a/26b(2) (26 combined/3.Intercept.combine
# UCC26,26a,26b2.sps, lines 9-52).

# Replication 1 (26a): a logging error created one file per session. Combine
# by booth in filename order, without subject sorting; drop test/incomplete
# runs. First test: blocks 6-7; visibility: 19-21
# (26a/Analyses/1.DATA.UcsConditioning26.13Feb2014.sps, lines 24-38, 150, 292).
path   <- unzip_folder("UCC26 (2014)")
raw_26 <- file.path(path, "26a", "UCC26a Raw Data")
files_26 <- unlist(lapply(booths, function(b) {
  sort(list.files(file.path(raw_26, paste0(b, ".UcsConditioning26.13Feb2014")),
                  pattern = "^[C-G][.]UcsConditioning26[.]13Feb2014.*[.]dat$",
                  full.names = TRUE))
}))
ucc26_1 <- read_inquisit(files_26, substr(basename(files_26), 1, 1)) |>
  filter(subject > 100, subject != 205) |>
  conditions_2014(first_procedure = c(375, 75)) |>
  split_tests_2014(test_blocks = 6:7, visibility_blocks = 19:21, cs_position = 4) |>
  mutate(id = "UCC26a")

# Replication 2 (26b): sort/clean booth files; drop test runs, incomplete
# sessions and a false start (26b/Analyses/1.DATA.UcsConditioning26a.
# 9Apr2014.sps, lines 21-98). Same blocks as replication 1.
ucc26_2 <- lapply(booths, function(b) {
  read_inquisit_sorted(file.path(path, "26b", "Analyses",
                                 paste0(b, ".UcsConditioning26a.9Apr2014.dat")), b)
}) |>
  bind_rows() |>
  arrange(subject, block, trial) |>
  filter(subject >= 101, subject != 118, !(subject == 222 & station == "E"),
         !(subject == 291 & time %in% "11:39:06")) |>
  conditions_2014(first_procedure = c(375, 75)) |>
  split_tests_2014(test_blocks = 6:7, visibility_blocks = 19:21, cs_position = 4) |>
  mutate(id = "UCC26b")

# Replication 3 (26c): two runs; first test in blocks 16-17, visibility in
# 20-22 (26c/Analyses/1.DATA.UcsConditioning26b.4Mar.2015.sps, lines 127-140).
# Run 1: drop condition 1 (procedural error; 1.DATA.UcsConditioning26b.
# 8Nov2014.sps, lines 20-128). Run 2: remove incomplete/reused sessions and
# renumber a duplicate; booth G reuses Run 1's file (lines 22-104).
# Keep Run 1 before Run 2 in each test (2.Intercept.UcsConditioning26b+b2.
# 4Mar2015.sps, lines 8, 386).
run_1 <- lapply(booths, function(b) {
  read_inquisit_sorted(file.path(path, "26c", "Analyses", "Analysis.8Nov2014",
                                 paste0(b, ".UcsConditioning26b.6Oct2014.dat")), b)
}) |>
  bind_rows() |>
  arrange(subject, block, trial) |>
  filter(subject >= 105, subject != 9999) |>
  conditions_2014(first_procedure = c(375, 75)) |>
  filter(cond != 1)
run_2_files <- c(paste0(c("C", "D", "E", "F"), ".UcsConditioning26b2.15Nov2014.dat"),
                 "G.UcsConditioning26b.6Oct2014.dat")
run_2 <- lapply(seq_along(booths), function(i) {
  read_inquisit_sorted(file.path(path, "26c", "Analyses", run_2_files[i]), booths[i])
}) |>
  bind_rows() |>
  filter(!(station == "C" & subject %in% c(337, 352, 435, 615)),
         !(station == "F" & subject %in% c(276, 1436694)),
         !(station == "G" & subject %in% c(185, 198, 210, 219, 224, 236, 246,
                                           353, 449, 593, 468))) |>
  arrange(subject, block, trial) |>
  mutate(subject = ifelse(subject == 408 & station == "F", 16408, subject)) |>
  conditions_2014(first_procedure = c(375, 75))
ucc26_3 <- bind_rows(
  split_tests_2014(run_1, test_blocks = 16:17, visibility_blocks = 20:22, cs_position = 4),
  split_tests_2014(run_2, test_blocks = 16:17, visibility_blocks = 20:22, cs_position = 4)
) |>
  mutate(id = "UCC26c")

# Apply exclusions per replication; retain participants with 75 ms first.
ucc26 <- lapply(list(ucc26_1, ucc26_2, ucc26_3), function(r) {
  r |>
    exclusions_2014() |>
    filter(first == 75) |>
    mutate(cs_acq = 75, cs_test = 75, study = "S18") |>
    as_long(archive = "UCC26 (2014)", id = r$id[1])
}) |>
  bind_rows()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Long data frame
# ─────────────────────────────────────────────────────────────────────────────

greenwald2017 <- bind_rows(ucc1, ucc2, ucc3, ucc4, ucc17b, ucc18a, ucc23b,
                           ucc23c, ucc24, ucc24_belg, ucc27, ucc26) |>
  arrange(as.integer(sub("S", "", study)), subject, desc(task), block, trial)

write.csv(greenwald2017, file.path(out_dir, "greenwald2017_long.csv"),
          row.names = FALSE)

# Check participants, trials and exclusions.
greenwald2017 |>
  summarise(participants = n_distinct(subject), trials = n(), .by = study)
greenwald2017 |>
  count(exclusion, task) |>
  arrange(exclusion, task)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Verification against the authors' results
# ─────────────────────────────────────────────────────────────────────────────

# d' of each participant in each task, from the trials kept, with the
# authors' procedure in each period:
#   2006-2007  rates over all trials; .25/N counts all trials of each CS
#   2010       rates over all trials
#   2011-2013  unweighted mean of the four CS x US cells; conditioning d' is
#              the mean of its two blocks (MEAN, at least one block; e.g.,
#              2.Intercept.UcsConditioning23a.10Oct2011.sps, line 379)
#   2014-2015  rates over all trials of each block; conditioning d' is the
#              mean of its two blocks (MEAN.2, both blocks;
#              2.Intercept.UcsConditioning27.9Apr2015.sps, line 686)
period_of <- function(archive) {
  case_when(archive %in% c("UCC1 (2006)", "UCC2 (2006)", "UCC3 (2007)", "UCC4 (2007)") ~ "2006",
            archive %in% c("UCC17b (2010)", "UCC18a (2010)") ~ "2010",
            archive %in% c("UCC23b (2011)", "UCC23c (2012)", "UCC24 (2012)", "UCC24 (Belg)") ~ "2011",
            .default = "2014")
}

kept <- greenwald2017 |>
  filter(exclusion == "none") |>
  mutate(period = period_of(archive))

dprimes <- bind_rows(
  # Both tasks, d' over all their trials
  kept |>
    filter(period == "2006") |>
    dprime_by(by = c("study", "subject", "task"), count_all = TRUE),
  kept |>
    filter(period == "2010" | (period == "2014" & task == "direct")) |>
    dprime_by(by = c("study", "subject", "task")),
  kept |>
    filter(period == "2011", task == "direct") |>
    dprime_by(by = c("study", "subject", "task"), cells = TRUE),
  # Conditioning test from 2011 on: mean of the d' of its two blocks
  kept |>
    filter(period == "2011", task == "indirect") |>
    dprime_by(by = c("study", "subject", "task", "test_block"), cells = TRUE) |>
    summarise(dprime = mean(dprime, na.rm = TRUE), .by = c(study, subject, task)),
  kept |>
    filter(period == "2014", task == "indirect") |>
    dprime_by(by = c("study", "subject", "task", "test_block")) |>
    summarise(dprime = if (n() == 2) mean(dprime) else NA_real_,
              .by = c(study, subject, task))
) |>
  mutate(dprime = ifelse(is.nan(dprime), NA, dprime)) |>
  pivot_wider(names_from = task, values_from = dprime)

# The authors' d' of each participant, for the reanalysis of their regressions
write.csv(dprimes, file.path(out_dir, "greenwald2017_dprimes.csv"), row.names = FALSE)

# Our results per data set: participants with each d' (as in the authors'
# one-sample t tests), and the regression of conditioning on visibility on
# the participants with both
ours <- dprimes |>
  summarise(n          = sum(!is.na(indirect) & !is.na(direct)),
            intercept  = coef(lm(indirect ~ direct))[1],
            slope      = coef(lm(indirect ~ direct))[2],
            n_ind      = sum(!is.na(indirect)),
            mean_ind   = mean(indirect, na.rm = TRUE),
            sd_ind     = sd(indirect, na.rm = TRUE),
            n_dir      = sum(!is.na(direct)),
            mean_dir   = mean(direct, na.rm = TRUE),
            sd_dir     = sd(direct, na.rm = TRUE),
            .by = study)

# The authors' results, from their SPSS outputs ("(2)Intercept...spv",
# "2.Intercept...spv" and "26 combined/3.Intercept...spv" in each folder):
# REGRESSION of conditioning d' on visibility d' (N, intercept and slope;
# the N and intercepts are also in Table A2 and Figure 4 of the article)
# and T-TEST of each d' against zero (N, mean and SD).
reported <- tribble(
  ~study, ~n, ~intercept, ~slope,    ~n_ind, ~mean_ind, ~sd_ind,  ~n_dir, ~mean_dir, ~sd_dir,
  "S1",    97,  0.006108, -0.056517,  97,  0.010525,  0.204090, 112, -0.112987, 0.497278,
  "S2",    47,  0.013790, -0.027361,  47,  0.016124,  0.109976,  63, -0.059114, 0.522778,
  "S3",    85,  0.063068,  0.054696,  86,  0.062389,  0.231623,  92, -0.005435, 0.294156,
  "S4",    85, -0.007642,  0.034456,  86, -0.006537,  0.336975,  92, -0.055649, 0.312057,
  "S5",    85, -0.026943,  0.102326,  86, -0.021553,  0.243185,  92,  0.056482, 0.296985,
  "S6",    43,  0.030657, -0.056742,  43,  0.025515,  0.181006,  51,  0.079069, 0.187833,
  "S7",    73,  0.043694,  0.013422,  73,  0.043495,  0.155800,  81,  0.006405, 0.608641,
  "S8",    67,  0.045026,  0.017641,  67,  0.048294,  0.238923,  74,  0.181127, 0.572404,
  "S9",    69,  0.013116,  0.020696,  70,  0.015720,  0.227893,  74,  0.144097, 0.414364,
  "S10",   92,  0.131959, -0.033831,  92,  0.123939,  0.196960, 101,  0.221931, 0.479603,
  "S11",   53,  0.051262,  0.078956,  54,  0.076690,  0.233689,  55,  0.310956, 0.463936,
  "S12",   82,  0.100941,  0.130164,  83,  0.126865,  0.273282,  88,  0.179915, 0.344226,
  "S13",   96,  0.117318,  0.108333,  96,  0.114885,  0.245705, 112, -0.064507, 0.499386,
  "S14",   85,  0.124632,  0.024194,  86,  0.124090,  0.285137,  92, -0.023902, 0.315852,
  "S15",   51,  0.176055, -0.034341,  51,  0.162819,  0.177548,  62,  0.333620, 0.644718,
  "S16",  145,  0.176926, -0.002426, 146,  0.175899,  0.227062, 160,  0.654798, 1.059574,
  "S17",   32,  0.209244,  0.060480,  32,  0.271475,  0.266242,  32,  1.028950, 1.213443,
  "S18",  350,  0.103403, -0.011027, 354,  0.099257,  0.242277, 379,  0.614354, 0.966189
)

# Side by side, with the largest absolute difference per data set. With six
# decimals in the reported values, any difference below 5e-7 is rounding.
comparison <- ours |>
  pivot_longer(-study, names_to = "statistic", values_to = "ours") |>
  inner_join(reported |> pivot_longer(-study, names_to = "statistic",
                                      values_to = "authors"),
             by = c("study", "statistic")) |>
  mutate(difference = abs(ours - authors))

comparison |>
  summarise(max_difference = max(difference), .by = study)

if (all(comparison$difference < 5e-7)) {
  message("All 18 data sets reproduce the authors' N, means, SDs, ",
          "intercepts and slopes.")
} else {
  comparison |> filter(difference >= 5e-7)
}

# Participants per data set, against the "Total subject N" of Table A2. This
# count is informative only. Table A2 seems to swap the totals of S8 and S9:
# the authors' own file of UCC27 (UcsCond27.1st_full.9Apr2015.sav) has 74
# and 75 participants, as here, and the analysed Ns (67 and 69) match. For
# S16, Table A2 gives 174, but the archive holds 163 participants in that
# condition after the authors' data preparation (161 in their final file).
table_a2_total <- c(S1 = 180, S2 = 78, S3 = 93, S4 = 93, S5 = 93, S6 = 51,
                    S7 = 113, S8 = 75, S9 = 74, S10 = 104, S11 = 56, S12 = 89,
                    S13 = 180, S14 = 93, S15 = 64, S16 = 174, S17 = 32,
                    S18 = 386)
greenwald2017 |>
  summarise(participants = n_distinct(subject), .by = study) |>
  mutate(table_a2 = table_a2_total[study]) |>
  filter(participants != table_a2)

# ─────────────────────────────────────────────────────────────────────────────
