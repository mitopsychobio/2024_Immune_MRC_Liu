library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)

getwd()

#qPCR platemap----

qPCR_replicates <- readxl::read_xlsx("plate_maps/qPCR_platemap.xlsx", sheet = 3)%>%
  pivot_longer(cols = -...1,names_to = "Column", values_to = "replicate")%>%
  rename(Row = ...1)%>%
  mutate(across('Column', str_replace, 'Col', '')) %>%
  unite(Well, Row, Column, sep = "", remove = F) 
qPCR_platemap <- readxl::read_xlsx("plate_maps/qPCR_platemap.xlsx", sheet = 1)%>%
  pivot_longer(cols = -...1,names_to = "Column", values_to = "sample_no")%>%
  rename(Row = ...1)%>%
  mutate(across('Column', str_replace, 'Col', '')) %>%
  unite(Well, Row, Column, sep = "", remove = F) %>%
  full_join(., qPCR_replicates)

qPCR_96 <- readxl::read_xlsx("plate_maps/qPCR_platemap.xlsx", sheet = 2)%>%
  pivot_longer(cols = -...1,names_to = "Column", values_to = "sample_no")%>%
  rename(Row = ...1)%>%
  mutate(across('Column', str_replace, 'Col', '')) %>%
  unite(Well, Row, Column, sep = "", remove = F)


MiSBIE_lymp <- readxl::read_xlsx("plate_maps/MiSBIE_lymp.xlsx", sheet = 1)%>%
  rename(row_96 = Row, column_96 = Column, well_96 = Well)


#read files----

MiSBIE_lymp_nd1_b2m_A <- readxl::read_xls("raw_data/qPCR_data/lymphocytes/MiSBIE_MHI_lymphocytes_A_ND1_B2M.xls", sheet = 3)[44:812,] %>%
  row_to_names(row_number = 1) %>%
  select(`Well Position`,`Target Name`, CT) 
MiSBIE_lymp_nd1_b2m_B <- readxl::read_xls("raw_data/qPCR_data/lymphocytes/MiSBIE_MHI_lymphocytes_B_ND1_B2M.xls", sheet = 3)[44:812,] %>%
  row_to_names(row_number = 1) %>%
  select(`Well Position`,`Target Name`, CT) 

MiSBIE_lymp_cox1_rnasep_A <- readxl::read_xls("raw_data/qPCR_data/lymphocytes/MiSBIE_MHI_lymphocytes_A_COXI_RNaseP.xls", sheet = 3)[44:812,] %>%
  row_to_names(row_number = 1) %>%
  select(`Well Position`,`Target Name`, CT)
MiSBIE_lymp_cox1_rnasep_B <- readxl::read_xls("raw_data/qPCR_data/lymphocytes/MiSBIE_MHI_lymphocytes_B_COXI_RNaseP.xls", sheet = 3)[44:812,] %>%
  row_to_names(row_number = 1) %>%
  select(`Well Position`,`Target Name`, CT) 




cutoff_cv <- 0.005

# ### ND1/B2M ####

## nd1/b2m plate_A ----

#nd1_A
MiSBIE_lymp_nd1_A <- MiSBIE_lymp_nd1_b2m_A %>%
  filter(`Target Name` == "ND1") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_nd1_A, by = "Well") %>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(nd1_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(nd1_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(nd1_ct_raw_cv = nd1_ct_raw_sd/nd1_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(nd1_ct_raw_cv > cutoff_cv & abs(CT_r1-nd1_ct_raw_mean) > abs(CT_r2-nd1_ct_raw_mean) & abs(CT_r1-nd1_ct_raw_mean) > abs(CT_r3-nd1_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(nd1_ct_raw_cv > cutoff_cv & abs(CT_r2-nd1_ct_raw_mean) > abs(CT_r1-nd1_ct_raw_mean) & abs(CT_r2-nd1_ct_raw_mean) > abs(CT_r3-nd1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(nd1_ct_raw_cv > cutoff_cv & abs(CT_r3-nd1_ct_raw_mean) > abs(CT_r1-nd1_ct_raw_mean) & abs(CT_r3-nd1_ct_raw_mean) > abs(CT_r2-nd1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(nd1_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(nd1_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(nd1_ct_cleaned_cv = nd1_ct_cleaned_sd/nd1_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))

#b2m_A

MiSBIE_lymp_b2m_A <- MiSBIE_lymp_nd1_b2m_A %>%
  filter(`Target Name` == "B2M") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_b2m_A, by = "Well")%>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(b2m_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(b2m_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(b2m_ct_raw_cv = b2m_ct_raw_sd/b2m_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(b2m_ct_raw_cv > cutoff_cv & abs(CT_r1-b2m_ct_raw_mean) > abs(CT_r2-b2m_ct_raw_mean) & abs(CT_r1-b2m_ct_raw_mean) > abs(CT_r3-b2m_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(b2m_ct_raw_cv > cutoff_cv & abs(CT_r2-b2m_ct_raw_mean) > abs(CT_r1-b2m_ct_raw_mean) & abs(CT_r2-b2m_ct_raw_mean) > abs(CT_r3-b2m_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(b2m_ct_raw_cv > cutoff_cv & abs(CT_r3-b2m_ct_raw_mean) > abs(CT_r1-b2m_ct_raw_mean) & abs(CT_r3-b2m_ct_raw_mean) > abs(CT_r2-b2m_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(b2m_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(b2m_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(b2m_ct_cleaned_cv = b2m_ct_cleaned_sd/b2m_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))



## nd1/b2m plate_B ----

#nd1_B

MiSBIE_lymp_nd1_B <- MiSBIE_lymp_nd1_b2m_B %>%
  filter(`Target Name` == "ND1") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_nd1_B, by = "Well")%>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(nd1_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(nd1_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(nd1_ct_raw_cv = nd1_ct_raw_sd/nd1_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(nd1_ct_raw_cv > cutoff_cv & abs(CT_r1-nd1_ct_raw_mean) > abs(CT_r2-nd1_ct_raw_mean) & abs(CT_r1-nd1_ct_raw_mean) > abs(CT_r3-nd1_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(nd1_ct_raw_cv > cutoff_cv & abs(CT_r2-nd1_ct_raw_mean) > abs(CT_r1-nd1_ct_raw_mean) & abs(CT_r2-nd1_ct_raw_mean) > abs(CT_r3-nd1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(nd1_ct_raw_cv > cutoff_cv & abs(CT_r3-nd1_ct_raw_mean) > abs(CT_r1-nd1_ct_raw_mean) & abs(CT_r3-nd1_ct_raw_mean) > abs(CT_r2-nd1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(nd1_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(nd1_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(nd1_ct_cleaned_cv = nd1_ct_cleaned_sd/nd1_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))


#b2m_B

MiSBIE_lymp_b2m_B <- MiSBIE_lymp_nd1_b2m_B %>%
  filter(`Target Name` == "B2M") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_b2m_B, by = "Well")%>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(b2m_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(b2m_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(b2m_ct_raw_cv = b2m_ct_raw_sd/b2m_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(b2m_ct_raw_cv > cutoff_cv & abs(CT_r1-b2m_ct_raw_mean) > abs(CT_r2-b2m_ct_raw_mean) & abs(CT_r1-b2m_ct_raw_mean) > abs(CT_r3-b2m_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(b2m_ct_raw_cv > cutoff_cv & abs(CT_r2-b2m_ct_raw_mean) > abs(CT_r1-b2m_ct_raw_mean) & abs(CT_r2-b2m_ct_raw_mean) > abs(CT_r3-b2m_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(b2m_ct_raw_cv > cutoff_cv & abs(CT_r3-b2m_ct_raw_mean) > abs(CT_r1-b2m_ct_raw_mean) & abs(CT_r3-b2m_ct_raw_mean) > abs(CT_r2-b2m_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(b2m_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(b2m_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(b2m_ct_cleaned_cv = b2m_ct_cleaned_sd/b2m_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))


#average the duplicate plates
MiSBIE_lymp_nd1_A_ct <- MiSBIE_lymp_nd1_A %>%
  select(sampleName, sample_no, nd1_ct_cleaned_mean) %>%
  rename(nd1_A_ct = nd1_ct_cleaned_mean) %>%
  na.omit()
MiSBIE_lymp_nd1_B_ct <- MiSBIE_lymp_nd1_B %>%
  select(sampleName, sample_no, nd1_ct_cleaned_mean) %>%
  rename(nd1_B_ct = nd1_ct_cleaned_mean) %>%
  na.omit()
MiSBIE_lymp_b2m_A_ct <- MiSBIE_lymp_b2m_A %>%
  select(sampleName, sample_no, b2m_ct_cleaned_mean) %>%
  rename(b2m_A_ct = b2m_ct_cleaned_mean) %>%
  na.omit()
MiSBIE_lymp_b2m_B_ct <- MiSBIE_lymp_b2m_B %>%
  select(sampleName, sample_no, b2m_ct_cleaned_mean) %>%
  rename(b2m_B_ct = b2m_ct_cleaned_mean) %>%
  na.omit()

MiSBIE_lymp_nd1_b2m <- full_join(MiSBIE_lymp_nd1_A_ct, MiSBIE_lymp_nd1_B_ct, by = c("sampleName", "sample_no")) %>%
  full_join(MiSBIE_lymp_b2m_A_ct, by = c("sampleName", "sample_no"))%>%
  full_join(MiSBIE_lymp_b2m_B_ct, by = c("sampleName", "sample_no")) %>%
  arrange(by_group = sample_no) %>%
  mutate(nd1_av_ct = rowMeans(select(., starts_with("nd1_")), na.rm = TRUE)) %>%
  mutate(b2m_av_ct = rowMeans(select(., starts_with("b2m_")), na.rm = TRUE)) %>%
  mutate (delt_ct_1 = b2m_av_ct - nd1_av_ct) %>%
  mutate(mtDNAcn_1 = ifelse(delt_ct_1>2, 2*2^delt_ct_1, NA)) 


# ### COXI / RNaseP ####

## cox1/rnasep plate_A ----

#cox1_A
MiSBIE_lymp_cox1_A <- MiSBIE_lymp_cox1_rnasep_A %>%
  filter(`Target Name` == "COXI") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_cox1_A, by = "Well")%>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(cox1_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(cox1_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(cox1_ct_raw_cv = cox1_ct_raw_sd/cox1_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(cox1_ct_raw_cv > cutoff_cv & abs(CT_r1-cox1_ct_raw_mean) > abs(CT_r2-cox1_ct_raw_mean) & abs(CT_r1-cox1_ct_raw_mean) > abs(CT_r3-cox1_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(cox1_ct_raw_cv > cutoff_cv & abs(CT_r2-cox1_ct_raw_mean) > abs(CT_r1-cox1_ct_raw_mean) & abs(CT_r2-cox1_ct_raw_mean) > abs(CT_r3-cox1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(cox1_ct_raw_cv > cutoff_cv & abs(CT_r3-cox1_ct_raw_mean) > abs(CT_r1-cox1_ct_raw_mean) & abs(CT_r3-cox1_ct_raw_mean) > abs(CT_r2-cox1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(cox1_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(cox1_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(cox1_ct_cleaned_cv = cox1_ct_cleaned_sd/cox1_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))


#rnaseP_A

MiSBIE_lymp_rnasep_A <- MiSBIE_lymp_cox1_rnasep_A %>%
  filter(`Target Name` == "RNaseP") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_rnasep_A, by = "Well")%>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(rnasep_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(rnasep_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(rnasep_ct_raw_cv = rnasep_ct_raw_sd/rnasep_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(rnasep_ct_raw_cv > cutoff_cv & abs(CT_r1-rnasep_ct_raw_mean) > abs(CT_r2-rnasep_ct_raw_mean) & abs(CT_r1-rnasep_ct_raw_mean) > abs(CT_r3-rnasep_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(rnasep_ct_raw_cv > cutoff_cv & abs(CT_r2-rnasep_ct_raw_mean) > abs(CT_r1-rnasep_ct_raw_mean) & abs(CT_r2-rnasep_ct_raw_mean) > abs(CT_r3-rnasep_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(rnasep_ct_raw_cv > cutoff_cv & abs(CT_r3-rnasep_ct_raw_mean) > abs(CT_r1-rnasep_ct_raw_mean) & abs(CT_r3-rnasep_ct_raw_mean) > abs(CT_r2-rnasep_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(rnasep_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(rnasep_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(rnasep_ct_cleaned_cv = rnasep_ct_cleaned_sd/rnasep_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))



## cox1/rnasep plate_B ----

#cox1_B

MiSBIE_lymp_cox1_B <- MiSBIE_lymp_cox1_rnasep_B %>%
  filter(`Target Name` == "COXI") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_cox1_B, by = "Well")%>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(cox1_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(cox1_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(cox1_ct_raw_cv = cox1_ct_raw_sd/cox1_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(cox1_ct_raw_cv > cutoff_cv & abs(CT_r1-cox1_ct_raw_mean) > abs(CT_r2-cox1_ct_raw_mean) & abs(CT_r1-cox1_ct_raw_mean) > abs(CT_r3-cox1_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(cox1_ct_raw_cv > cutoff_cv & abs(CT_r2-cox1_ct_raw_mean) > abs(CT_r1-cox1_ct_raw_mean) & abs(CT_r2-cox1_ct_raw_mean) > abs(CT_r3-cox1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(cox1_ct_raw_cv > cutoff_cv & abs(CT_r3-cox1_ct_raw_mean) > abs(CT_r1-cox1_ct_raw_mean) & abs(CT_r3-cox1_ct_raw_mean) > abs(CT_r2-cox1_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(cox1_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(cox1_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(cox1_ct_cleaned_cv = cox1_ct_cleaned_sd/cox1_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))


#rnaseP_B

MiSBIE_lymp_rnasep_B <- MiSBIE_lymp_cox1_rnasep_B %>%
  filter(`Target Name` == "RNaseP") %>%
  rename(Well = `Well Position`) %>%
  full_join(qPCR_platemap, MiSBIE_lymp_rnasep_B, by = "Well")%>%
  select(sample_no, replicate, Well, CT) %>%
  pivot_wider(names_from = replicate, values_from = c(Well, CT)) %>%
  na.omit() %>% select(-c("Well_NA", "CT_NA")) %>%
  arrange(by_group = sample_no) %>%
  full_join(qPCR_96,., by = "sample_no") %>%
  rename(well_96 = Well) %>% rename (row_96 = Row) %>% rename (column_96 = Column) %>%
  full_join(MiSBIE_lymp,.) %>%
  mutate(across(starts_with("CT"), ~replace(., . == "Undetermined", NA))) %>%
  mutate(CT_r1 = as.numeric(CT_r1)) %>%
  mutate(CT_r2 = as.numeric(CT_r2)) %>%
  mutate(CT_r3 = as.numeric(CT_r3)) %>%
  mutate(rnasep_ct_raw_mean = rowMeans(select(., starts_with("CT_")), na.rm = TRUE)) %>%
  mutate(rnasep_ct_raw_sd = rowSds(as.matrix(select(., starts_with("CT"))), na.rm = TRUE)) %>%
  mutate(rnasep_ct_raw_cv = rnasep_ct_raw_sd/rnasep_ct_raw_mean) %>%
  mutate(CT_r1_cleaned = ifelse(rnasep_ct_raw_cv > cutoff_cv & abs(CT_r1-rnasep_ct_raw_mean) > abs(CT_r2-rnasep_ct_raw_mean) & abs(CT_r1-rnasep_ct_raw_mean) > abs(CT_r3-rnasep_ct_raw_mean) & !is.na(CT_r2) & !is.na(CT_r3),NA, CT_r1)) %>%
  mutate(CT_r2_cleaned = ifelse(rnasep_ct_raw_cv > cutoff_cv & abs(CT_r2-rnasep_ct_raw_mean) > abs(CT_r1-rnasep_ct_raw_mean) & abs(CT_r2-rnasep_ct_raw_mean) > abs(CT_r3-rnasep_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r3),NA, CT_r2)) %>%
  mutate(CT_r3_cleaned = ifelse(rnasep_ct_raw_cv > cutoff_cv & abs(CT_r3-rnasep_ct_raw_mean) > abs(CT_r1-rnasep_ct_raw_mean) & abs(CT_r3-rnasep_ct_raw_mean) > abs(CT_r2-rnasep_ct_raw_mean) & !is.na(CT_r1) & !is.na(CT_r2),NA, CT_r3)) %>%
  mutate(rnasep_ct_cleaned_mean = rowMeans(select(., starts_with("CT") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(rnasep_ct_cleaned_sd = rowSds(as.matrix(select(., starts_with("CT") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(rnasep_ct_cleaned_cv = rnasep_ct_cleaned_sd/rnasep_ct_cleaned_mean) %>%
  transform(sampleName = case_when(
    sample_no == "97" ~ "STD_1", sample_no == "98" ~ "STD_2",
    sample_no == "99" ~ "STD_3", sample_no == "100" ~ "STD_4",
    sample_no == "101" ~ "STD_5", sample_no == "102" ~ "STD_6",
    sample_no == "103" ~ "STD_7", sample_no == "104" ~ "STD_8",
    TRUE ~ sampleName
  ))


#average the duplicate plates----
MiSBIE_lymp_cox1_A_ct <- MiSBIE_lymp_cox1_A %>%
  select(sampleName, sample_no, cox1_ct_cleaned_mean) %>%
  rename(cox1_A_ct = cox1_ct_cleaned_mean) %>%
  na.omit()
MiSBIE_lymp_cox1_B_ct <- MiSBIE_lymp_cox1_B %>%
  select(sampleName, sample_no, cox1_ct_cleaned_mean) %>%
  rename(cox1_B_ct = cox1_ct_cleaned_mean) %>%
  na.omit()
MiSBIE_lymp_rnasep_A_ct <- MiSBIE_lymp_rnasep_A %>%
  select(sampleName, sample_no, rnasep_ct_cleaned_mean) %>%
  rename(rnasep_A_ct = rnasep_ct_cleaned_mean) %>%
  na.omit()
MiSBIE_lymp_rnasep_B_ct <- MiSBIE_lymp_rnasep_B %>%
  select(sampleName, sample_no, rnasep_ct_cleaned_mean) %>%
  rename(rnasep_B_ct = rnasep_ct_cleaned_mean) %>%
  na.omit()

MiSBIE_lymp_cox1_rnasep <- full_join(MiSBIE_lymp_cox1_A_ct, MiSBIE_lymp_cox1_B_ct, by = c("sampleName", "sample_no")) %>%
  full_join(MiSBIE_lymp_rnasep_A_ct, by = c("sampleName", "sample_no"))%>%
  full_join(MiSBIE_lymp_rnasep_B_ct, by = c("sampleName", "sample_no")) %>%
  arrange(by_group = sample_no) %>%
  mutate(cox1_av_ct = rowMeans(select(., starts_with("cox1_")), na.rm = TRUE)) %>%
  mutate(rnasep_av_ct = rowMeans(select(., starts_with("rnasep_")), na.rm = TRUE)) %>%
  mutate (delt_ct_2 = rnasep_av_ct - cox1_av_ct) %>%
  mutate(mtDNAcn_2 = ifelse(delt_ct_2>2, 2*2^delt_ct_2, NA)) 

#calculate mtDNAcn----


MiSBIE_lymp_mtDNAcn <- full_join(MiSBIE_lymp_nd1_b2m, MiSBIE_lymp_cox1_rnasep, by = c("sampleName", "sample_no")) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3")) %>%
  mutate(mtDNAcn_av = rowMeans(select(., starts_with("mtDNAcn")))) %>%
  mutate(mtDNAcn_sd = rowSds(as.matrix(select(., starts_with("mtDNAcn"))))) %>%
  mutate(mtDNAcn_cv = mtDNAcn_sd/mtDNAcn_av) %>%
  mutate (deviation = mtDNAcn_2 / mtDNAcn_1 -1) %>%
  transform(deviation = as.numeric(deviation))

average_deviation = mean(MiSBIE_lymp_mtDNAcn$deviation, na.rm = T)

mtDNAcn_cutoff_a <- average_deviation - 0.2
mtDNAcn_cutoff_b <- average_deviation + 0.2

MiSBIE_lymp_mtDNAcn <- MiSBIE_lymp_mtDNAcn %>%
  mutate (check = ifelse(deviation > mtDNAcn_cutoff_b | deviation < mtDNAcn_cutoff_a, "check", ""))

#address failures for mtDNAcn ---- 
#lin reg using plate 1 only. in MiSBIE_lymp samples used both plates together


MiSBIE_lymp_mtDNAcn_summ <- MiSBIE_lymp_mtDNAcn %>%
  select(sampleName, nd1_av_ct, b2m_av_ct, cox1_av_ct, rnasep_av_ct, mtDNAcn_1, mtDNAcn_2, mtDNAcn_av, mtDNAcn_cv, deviation, check)

# ENTER FAILED PRIMERS ####

ND1_failure = c()
COX1_failure = c("Mi017","Mi073","Mi005","Mi093")
B2M_failure = c()
RNaseP_failure = c()


mtDNAcn_lin_reg <- MiSBIE_lymp_mtDNAcn_summ %>%
  lm(mtDNAcn_2 ~ mtDNAcn_1, data = .)
mtDNAcn_lin_reg_slope <- coef(mtDNAcn_lin_reg)[2]
mtDNAcn_lin_reg_int <- coef(mtDNAcn_lin_reg)[1]  

MiSBIE_lymp_mtDNAcn_summ <- MiSBIE_lymp_mtDNAcn_summ %>%
  mutate(failed_primer = case_when(
    sampleName %in% ND1_failure ~ "ND1_failure",
    sampleName %in% COX1_failure ~ "COX1_failure",
    sampleName %in% B2M_failure ~ "B2M_failure",
    sampleName %in% RNaseP_failure ~ "RNaseP_failure")) %>%
  mutate(mtDNAcn_1_fail_checked = ifelse(
    sampleName %in% ND1_failure | sampleName %in% B2M_failure, 
    mtDNAcn_2/mtDNAcn_lin_reg_slope - mtDNAcn_lin_reg_int, mtDNAcn_1)) %>%
  mutate(mtDNAcn_2_fail_checked = ifelse(
    sampleName %in% COX1_failure |sampleName %in% RNaseP_failure, 
    mtDNAcn_1*mtDNAcn_lin_reg_slope + mtDNAcn_lin_reg_int, mtDNAcn_2)) %>%
  mutate(mtDNAcn_final = rowMeans(select(., starts_with("mtDNAcn") & ends_with("fail_checked"))))


#calculate nNORM ----

MiSBIE_lymp_nNORM_summ <- MiSBIE_lymp_mtDNAcn_summ %>%
  select(sampleName, b2m_av_ct, rnasep_av_ct) %>%
  mutate(nNORM_1 = 2^b2m_av_ct/100000) %>%
  mutate(nNORM_2 = 2^rnasep_av_ct/100000) %>%
  mutate(nNORM_av = rowMeans(select(., starts_with("nNORM"))))

nNORM_lin_reg <- MiSBIE_lymp_nNORM_summ %>%
  lm(nNORM_2 ~ nNORM_1, data = .)
nNORM_lin_reg_slope <- coef(nNORM_lin_reg)[2]
nNORM_lin_reg_int <- coef(nNORM_lin_reg)[1]  

MiSBIE_lymp_nNORM_summ <- MiSBIE_lymp_nNORM_summ %>%
  mutate(failed_primer = case_when(
    sampleName %in% B2M_failure ~ "B2M_failure",
    sampleName %in% RNaseP_failure ~ "RNaseP_failure")) %>%
  mutate(nNORM_1_fail_checked = ifelse(
    sampleName %in% B2M_failure, 
    nNORM_2/nNORM_lin_reg_slope - nNORM_lin_reg_int, nNORM_1)) %>%
  mutate(nNORM_2_fail_checked = ifelse(
    sampleName %in% RNaseP_failure, 
    nNORM_1*nNORM_lin_reg_slope + nNORM_lin_reg_int, nNORM_2)) %>%
  mutate(nNORM_final = rowMeans(select(., starts_with("nNORM") & ends_with("fail_checked"))))

MiSBIE_lymp_qPCR <- full_join(MiSBIE_lymp, MiSBIE_lymp_mtDNAcn_summ, by = c("sampleName")) %>%
  full_join(MiSBIE_lymp_nNORM_summ, by = c("sampleName")) %>%
  # select(well_96, sampleName, mtDNAcn_final, nNORM_final) %>%
  rename(Well = well_96, MiSBIE_lymp_mtDNAcn_final = mtDNAcn_final, MiSBIE_lymp_nNORM_final = nNORM_final) %>%
  transform(plate = case_when(plate == "lymphocyte_F.2" ~ plate, TRUE ~ "lymphocyte_F.2"))
write_xlsx(MiSBIE_lymp_qPCR,"interim_data/qPCR/MiSBIE_lymp_qPCR_all.xlsx")
