#!/usr/bin/env Rscript

# This script reimplements what IPS2FpGs.sh does:
# It filters tables from Pfam and Superfamily to identify R genes

library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)

# Args are:
# file1 (e.g. Pfam)
# file2 (e.g. Superfamily)
# outfile (Basename only)
# conf1
# conf2

## Alternative: read from conf1 and conf2 IPS2fpGs
# domains <- read_tsv(args[5], col_names = F) %>% set_colnames(c("Number", "Full", "Class", "Gene")) %>%
#  left_join(read_tsv(args[4], col_names = F) %>% set_colnames(c("Family", "Number")), by = "Family")

# Define domain - pfam or superfam relationships
## Currently, these use what is defined by HRP

domains <-
  data.frame(
    Class = c(
      "C",
      "T",
      "T",
      "T",
      "R",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "L",
      "N",
      "N"
    ),
    Family = c(
      "PF18052",
      "PF01582",
      "PF13676",
      "SSF52200",
      "PF05659",
      "SSF52047",
      "PF00560",
      "PF08263",
      "PF12799",
      "PF13516",
      "PF13855",
      "PF14580",
      "SSF52058",
      "SSF52540",
      "PF00931"
    )
  )

# Subset into individual classes

TIR <-  filter(domains, Class == "T")
COIL <- filter(domains, Class == "C")
RPW8 <- filter(domains, Class == "R")
LRR <-  filter(domains, Class == "L")
NB <-  filter(domains, Class == "N")

# Read in files

rbind(
  read_tsv(args[1], col_names = F, show_col_types = FALSE),
  read_tsv(args[2], col_names = F, show_col_types = FALSE)
) %>%
  set_colnames(
    c(
      "ID",
      "hash",
      "Length",
      "Analysis",
      "Family",
      "Description",
      "Start",
      "End",
      "Score",
      "Status",
      "Date",
      "InterPro_Accession",
      "InterPro_description"
    )
  ) %>%
  # Check TIR, COIL, RPW, NB and LRR presence
  group_by(ID) %>%
  mutate(
    has_TIR = case_when(Family %in% TIR$Family ~ TRUE,
                        TRUE ~ FALSE),
    has_TIR_evidence = case_when(Family %in% TIR$Family ~ Family,
                                 TRUE ~ NA_character_),
    has_COIL = case_when(Family %in% COIL$Family ~ TRUE,
                         TRUE ~ FALSE),
    has_COIL_evidence = case_when(Family %in% COIL$Family ~ Family,
                                  TRUE ~ NA_character_),
    has_RPW8 = case_when(Family %in% RPW8$Family ~ TRUE,
                         TRUE ~ FALSE),
    has_RPW8_evidence = case_when(Family %in% RPW8$Family ~ Family,
                                  TRUE ~ NA_character_),
    has_NB = case_when(Family %in% NB$Family ~ TRUE,
                       TRUE ~ FALSE),
    has_NB_evidence = case_when(Family %in% NB$Family ~ Family,
                                TRUE ~ NA_character_),
    has_LRR = case_when(Family %in% LRR$Family ~ TRUE,
                        TRUE ~ FALSE),
    has_LRR_evidence = case_when(Family %in% LRR$Family ~ Family,
                                 TRUE ~ NA_character_)
  ) %>%
  # Make uniques
  summarize(
    has_TIR = case_when(sum(has_TIR)  > 0 ~ TRUE, TRUE ~ FALSE),
    has_COIL = case_when(sum(has_COIL) > 0 ~ TRUE, TRUE ~ FALSE),
    has_RPW8 = case_when(sum(has_RPW8) > 0 ~ TRUE, TRUE ~ FALSE),
    has_NB = case_when(sum(has_NB) > 0 ~ TRUE, TRUE ~ FALSE),
    has_LRR = case_when(sum(has_LRR) > 0 ~ TRUE, TRUE ~ FALSE),
    TIR_evidence = str_flatten(has_TIR_evidence %>% {
      .[!is.na(.)]
    }, collapse = " "),
    COIL_evidence = str_flatten(has_COIL_evidence %>% {
      .[!is.na(.)]
    }, collapse = " "),
    RPW8_evidence = str_flatten(has_RPW8_evidence %>% {
      .[!is.na(.)]
    }, collapse = " "),
    NB_evidence = str_flatten(has_NB_evidence %>% {
      .[!is.na(.)]
    }, collapse = " "),
    LRR_evidence = str_flatten(has_LRR_evidence %>% {
      .[!is.na(.)]
    }, collapse = " ")
  ) %>%
  # Filter of NB-LRR
  filter(has_NB & has_LRR) %>%
  # Annotate
  mutate(Type = case_when(has_TIR ~ "TNL",
                          has_COIL ~ "CNL",
                          has_RPW8 ~ "RNL",
                          TRUE ~ "NL")) %>%
  ungroup() %T>%
  # Write larger table
  write_tsv(paste0(args[3], "_NLR_table.tsv")) %>%
  dplyr::select(ID, Type) %>%
  # Write relevant table
  write_tsv(paste0(args[3], "_NLR_genes.tsv"))
