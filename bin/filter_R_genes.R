# This script reimplements what IPS2FpGs.sh does
# It filters tables from Pfam and Superfamily to identify R genes

library(tidyverse)
library(magrittr)


args <- commandArgs(trailingOnly = TRUE)

# Args are:
# file1 (e.g. Pfam)
# file2 (e.g. Superfamily)
# outfile (Basename only)

domains <- read_tsv("assets/conf2.tsv", col_names = F) %>% set_colnames(c("Number", "Full", "Class", "Gene")) %>% 
  left_join(read_tsv("assets/conf1.tsv", col_names = F) %>% set_colnames(c("Family", "Number")))

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

TIR <-  filter(domains, Class == "T")
COIL <- filter(domains, Class == "C")
RPW8 <- filter(domains, Class == "R")
LRR <-  filter(domains, Class == "L")
NB <-  filter(domains, Class == "N")

rbind(
  read_tsv(args[1], col_names = F),
  read_tsv(args[2], col_names = F)
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
  group_by(ID) %>% 
  mutate(has_TIR = case_when(Family %in% TIR$Family ~ TRUE,
                             TRUE ~ FALSE),
         has_COIL = case_when(Family %in% COIL$Family ~ TRUE,
                              TRUE ~ FALSE),
         has_RPW8 = case_when(Family %in% RPW8$Family ~ TRUE,
                              TRUE ~ FALSE),
         has_NB = case_when(Family %in% NB$Family ~ TRUE,
                              TRUE ~ FALSE),
         has_LRR = case_when(Family %in% LRR$Family ~ TRUE,
                              TRUE ~ FALSE)) %>% 
  summarize(has_TIR = case_when(sum(has_TIR)  > 0 ~TRUE, TRUE ~ FALSE),
            has_COIL = case_when(sum(has_COIL) > 0 ~TRUE, TRUE ~ FALSE),
            has_RPW8 = case_when(sum(has_RPW8) > 0 ~TRUE, TRUE ~ FALSE),
            has_NB = case_when(sum(has_NB) > 0 ~TRUE, TRUE ~ FALSE),
            has_LRR = case_when(sum(has_LRR) > 0 ~TRUE, TRUE ~ FALSE),
            has_COIL = case_when(sum(has_COIL) > 0 ~TRUE, TRUE ~ FALSE)) %>% 
  #filter((has_TIR | has_COIL | has_RPW8) & has_NB & has_LRR)
  filter(has_NB & has_LRR) %>% 
  mutate(Type = case_when(
    has_TIR ~ "TNL",
    has_COIL ~ "CNL",
    has_RPW8 ~ "RNL",
    TRUE ~ "NL"
  )) %>% 
  ungroup() %T>%
  write_tsv(paste0(args[3],"_NLR_table.tsv")) %>% 
  dplyr::select(ID, Type) %>% 
  write_tsv(paste0(args[3],"NLR_genes.tsv"))
         