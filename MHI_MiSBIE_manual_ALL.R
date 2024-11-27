library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)

## read all files
getwd()

#PBMC_F.1
PBMC_F.1_cs <- readxl::read_xlsx("interim_data/manual/PBMC_F.1_cs.xlsx", sheet = 1) %>%
  rename(cs = PBMC_F.1_cs_spec_cleaned, cs_cell = PBMC_F.1_cs_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_F.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
PBMC_F.1_ci <- readxl::read_xlsx("interim_data/manual/PBMC_F.1_ci.xlsx", sheet = 1) %>%
  rename(ci = PBMC_F.1_ci_spec_cleaned, ci_cell = PBMC_F.1_ci_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_F.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
PBMC_F.1_sdh <- readxl::read_xlsx("interim_data/manual/PBMC_F.1_sdh.xlsx", sheet = 1) %>%
  rename(sdh = PBMC_F.1_sdh_spec_cleaned, sdh_cell = PBMC_F.1_sdh_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_F.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
PBMC_F.1_cox <- readxl::read_xlsx("interim_data/manual/PBMC_F.1_cox.xlsx", sheet = 1) %>%
  rename(cox = PBMC_F.1_cox_spec_cleaned, cox_cell = PBMC_F.1_cox_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_F.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))

PBMC_F.1_qpcr <- readxl::read_xlsx("interim_data/qPCR/MiSBIE_PBMC_F.1_qPCR.xlsx", sheet = 1) %>%
  rename(mtDNAcn = MiSBIE_PBMC_F.1_mtDNAcn_final, nNORM = MiSBIE_PBMC_F.1_nNORM_final)%>% 
  mutate (plate = "PBMC_F.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F)  %>%
  filter(!sampleName %in% c("EMPTY",NA)) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA))


#PBMC_S.1
PBMC_S.1_cs <- readxl::read_xlsx("interim_data/manual/PBMC_S.1_cs.xlsx", sheet = 1) %>%
  rename(cs = PBMC_S.1_cs_spec_cleaned, cs_cell = PBMC_S.1_cs_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_S.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
PBMC_S.1_ci <- readxl::read_xlsx("interim_data/manual/PBMC_S.1_ci.xlsx", sheet = 1) %>%
  rename(ci = PBMC_S.1_ci_spec_cleaned, ci_cell = PBMC_S.1_ci_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_S.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
PBMC_S.1_sdh <- readxl::read_xlsx("interim_data/manual/PBMC_S.1_sdh.xlsx", sheet = 1) %>%
  rename(sdh = PBMC_S.1_sdh_spec_cleaned, sdh_cell = PBMC_S.1_sdh_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_S.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
PBMC_S.1_cox <- readxl::read_xlsx("interim_data/manual/PBMC_S.1_cox.xlsx", sheet = 1) %>%
  rename(cox = PBMC_S.1_cox_spec_cleaned, cox_cell = PBMC_S.1_cox_spec_cleaned_cell) %>% 
  mutate (plate = "PBMC_S.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))

PBMC_S.1_qpcr <- readxl::read_xlsx("interim_data/qPCR/MiSBIE_PBMC_S.1_qPCR.xlsx", sheet = 1) %>%
  rename(mtDNAcn = MiSBIE_PBMC_S.1_mtDNAcn_final, nNORM = MiSBIE_PBMC_S.1_nNORM_final)%>% 
  mutate (plate = "PBMC_S.1") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F)  %>%
  filter(!sampleName %in% c("EMPTY",NA)) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA))

#lymp
lymp_cs <- readxl::read_xlsx("interim_data/manual/lymp_cs.xlsx", sheet = 1) %>%
  rename(cs = lymp_cs_spec_cleaned, cs_cell = lymp_cs_spec_cleaned_cell) %>% 
  mutate (plate = "lymp") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
lymp_ci <- readxl::read_xlsx("interim_data/manual/lymp_ci.xlsx", sheet = 1) %>%
  rename(ci = lymp_ci_spec_cleaned, ci_cell = lymp_ci_spec_cleaned_cell) %>% 
  mutate (plate = "lymp") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
lymp_sdh <- readxl::read_xlsx("interim_data/manual/lymp_sdh.xlsx", sheet = 1) %>%
  rename(sdh = lymp_sdh_spec_cleaned, sdh_cell = lymp_sdh_spec_cleaned_cell) %>% 
  mutate (plate = "lymp") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
lymp_cox <- readxl::read_xlsx("interim_data/manual/lymp_cox.xlsx", sheet = 1) %>%
  rename(cox = lymp_cox_spec_cleaned, cox_cell = lymp_cox_spec_cleaned_cell) %>% 
  mutate (plate = "lymp") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))

lymp_qpcr <- readxl::read_xlsx("interim_data/qPCR/MiSBIE_lymp_qPCR.xlsx", sheet = 1) %>%
  rename(mtDNAcn = MiSBIE_lymp_mtDNAcn_final, nNORM = MiSBIE_lymp_nNORM_final)%>% 
  mutate (plate = "lymp") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F)  %>%
  filter(!sampleName %in% c("EMPTY",NA)) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA))


#mono
mono_cs <- readxl::read_xlsx("interim_data/manual/mono_cs.xlsx", sheet = 1) %>%
  rename(cs = mono_cs_spec_cleaned, cs_cell = mono_cs_spec_cleaned_cell) %>% 
  mutate (plate = "mono") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
mono_ci <- readxl::read_xlsx("interim_data/manual/mono_ci.xlsx", sheet = 1) %>%
  rename(ci = mono_ci_spec_cleaned, ci_cell = mono_ci_spec_cleaned_cell) %>% 
  mutate (plate = "mono") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
mono_sdh <- readxl::read_xlsx("interim_data/manual/mono_sdh.xlsx", sheet = 1) %>%
  rename(sdh = mono_sdh_spec_cleaned, sdh_cell = mono_sdh_spec_cleaned_cell) %>% 
  mutate (plate = "mono") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
mono_cox <- readxl::read_xlsx("interim_data/manual/mono_cox.xlsx", sheet = 1) %>%
  rename(cox = mono_cox_spec_cleaned, cox_cell = mono_cox_spec_cleaned_cell) %>% 
  mutate (plate = "mono") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))

mono_qpcr <- readxl::read_xlsx("interim_data/qPCR/MiSBIE_mono_qPCR.xlsx", sheet = 1) %>%
  rename(mtDNAcn = MiSBIE_mono_mtDNAcn_final, nNORM = MiSBIE_mono_nNORM_final)%>% 
  mutate (plate = "mono") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F)  %>%
  filter(!sampleName %in% c("EMPTY",NA)) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA))


#neu
neu_cs <- readxl::read_xlsx("interim_data/manual/neu_cs.xlsx", sheet = 1) %>%
  rename(cs = neu_cs_spec_cleaned, cs_cell = neu_cs_spec_cleaned_cell) %>% 
  mutate (plate = "neu") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
neu_ci <- readxl::read_xlsx("interim_data/manual/neu_ci.xlsx", sheet = 1) %>%
  rename(ci = neu_ci_spec_cleaned, ci_cell = neu_ci_spec_cleaned_cell) %>% 
  mutate (plate = "neu") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
neu_sdh <- readxl::read_xlsx("interim_data/manual/neu_sdh_4.24.xlsx", sheet = 1) %>%
  rename(sdh = neu_sdh_spec_cleaned, sdh_cell = neu_sdh_spec_cleaned_cell) %>% 
  mutate (plate = "neu") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
neu_cox <- readxl::read_xlsx("interim_data/manual/neu_cox.xlsx", sheet = 1) %>%
  rename(cox = neu_cox_spec_cleaned, cox_cell = neu_cox_spec_cleaned_cell) %>% 
  mutate (plate = "neu") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))

neu_qpcr <- readxl::read_xlsx("interim_data/qPCR/MiSBIE_neu_qPCR.xlsx", sheet = 1) %>%
  rename(mtDNAcn = MiSBIE_neu_mtDNAcn_final, nNORM = MiSBIE_neu_nNORM_final)%>% 
  mutate (plate = "neu") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F)  %>%
  filter(!sampleName %in% c("EMPTY",NA)) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA))


#plt
plt_cs <- readxl::read_xlsx("interim_data/manual/plt_cs.xlsx", sheet = 1) %>%
  rename(cs = plt_cs_spec_cleaned, cs_cell = plt_cs_spec_cleaned_cell) %>% 
  mutate (plate = "plt") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
plt_ci <- readxl::read_xlsx("interim_data/manual/plt_ci.xlsx", sheet = 1) %>%
  rename(ci = plt_ci_spec_cleaned, ci_cell = plt_ci_spec_cleaned_cell) %>% 
  mutate (plate = "plt") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
plt_sdh <- readxl::read_xlsx("interim_data/manual/plt_sdh.xlsx", sheet = 1) %>%
  rename(sdh = plt_sdh_spec_cleaned, sdh_cell = plt_sdh_spec_cleaned_cell) %>% 
  mutate (plate = "plt") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))
plt_cox <- readxl::read_xlsx("interim_data/manual/plt_cox.xlsx", sheet = 1) %>%
  rename(cox = plt_cox_spec_cleaned, cox_cell = plt_cox_spec_cleaned_cell) %>% 
  mutate (plate = "plt") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) %>%
  filter(!sampleName %in% c("EMPTY",NA))

plt_qpcr <- readxl::read_xlsx("interim_data/qPCR/MiSBIE_plt_qPCR.xlsx", sheet = 1) %>%
  rename(mtDNAcn = MiSBIE_plt_mtDNAcn_final, nNORM = MiSBIE_plt_nNORM_final)%>% 
  mutate (plate = "plt") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F)  %>%
  filter(!sampleName %in% c("EMPTY",NA)) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA))


plt_qpcr_mtDNAct <- readxl::read_xlsx("interim_data/qPCR/MiSBIE_plt_qPCR_all.xlsx", sheet = 1) %>%
  select(Well, sampleName, nd1_av_ct, cox1_av_ct) %>%
  mutate (plate = "plt") %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F)  %>%
  filter(!sampleName %in% c("EMPTY",NA, "PBMC_S.1_Mi051", "neu_F.2_Mi068")) %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA, "STD_1","STD_2","STD_3","STD_4","STD_5","STD_6","STD_7","STD_8"))


## bind qpcr together


misbie_qpcr <- rbind(PBMC_F.1_qpcr, PBMC_S.1_qpcr, lymp_qpcr, mono_qpcr, neu_qpcr, plt_qpcr)


rm(PBMC_F.1_qpcr, PBMC_S.1_qpcr, lymp_qpcr, mono_qpcr, neu_qpcr, plt_qpcr)  

## bind all enzyme activities together

misbie_cs <- rbind(PBMC_F.1_cs, PBMC_S.1_cs, lymp_cs, mono_cs, neu_cs, plt_cs) %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) 
misbie_ci <- rbind(PBMC_F.1_ci, PBMC_S.1_ci, lymp_ci, mono_ci, neu_ci, plt_ci) %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) 
misbie_sdh <- rbind(PBMC_F.1_sdh, PBMC_S.1_sdh, lymp_sdh, mono_sdh, neu_sdh, plt_sdh) %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) 
misbie_cox <- rbind(PBMC_F.1_cox, PBMC_S.1_cox, lymp_cox, mono_cox, neu_cox, plt_cox) %>% 
  unite(plate_ID, plate, sampleName, sep = "_", remove = F) 

rm(PBMC_F.1_cs, PBMC_S.1_cs, lymp_cs, mono_cs, neu_cs, plt_cs,
   PBMC_F.1_ci, PBMC_S.1_ci, lymp_ci, mono_ci, neu_ci, plt_ci,
   PBMC_F.1_sdh, PBMC_S.1_sdh, lymp_sdh, mono_sdh, neu_sdh, plt_sdh,
   PBMC_F.1_cox, PBMC_S.1_cox, lymp_cox, mono_cox, neu_cox, plt_cox
)

#combine all activities activity in units mmol/min/10^6 cells

misbie_MHI_all <- full_join(misbie_cs, misbie_ci, by = c("plate_ID", "Well", "plate", "sampleName")) %>%
  select("plate_ID", "Well", "plate", "sampleName", "cs_cell", "ci_cell") %>%
  full_join(misbie_sdh) %>%
  select(-"sdh") %>%
  full_join(misbie_cox) %>%
  select(-"cox") %>%
    full_join (misbie_qpcr) %>%
  rename(cs = cs_cell, ci = ci_cell, sdh = sdh_cell, cox = cox_cell) %>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  ) %>%
  mutate(Column = as.numeric(Column)) %>%
  arrange(plate, Column) 

# data cleaning

#these are empty wells

cs_excl <- c()
ci_excl <- c()
sdh_excl <- c()
cox_excl <- c()
qPCR_excl <- c()

misbie_MHI_all <- misbie_MHI_all %>%
  mutate(cs = case_when(
    sampleName %in% cs_excl ~ NA,
    TRUE ~ cs
  )) %>%
  mutate(ci = case_when(
    sampleName %in% ci_excl ~ NA,
    TRUE ~ ci
  )) %>%
  mutate(sdh = case_when(
    sampleName %in% sdh_excl ~ NA,
    TRUE ~ sdh
  )) %>%
  mutate(cox = case_when(
    sampleName %in% cs_excl ~ NA,
    TRUE ~ cox
  )) %>%
  mutate(mtDNAcn = case_when(
    sampleName %in% qPCR_excl ~ NA,
    TRUE ~ mtDNAcn
  )) %>%
  mutate(nNORM = case_when(
    sampleName %in% qPCR_excl ~ NA,
    TRUE ~ nNORM
  )) 

# plate adjustment

ref_samples <-  misbie_MHI_all %>%
  filter(sampleName %in% c("REF_1", "REF_2", "REF_3")) %>%
  dplyr::select(plate, sampleName, cs, ci, sdh, cox)

cs_ref_av = mean(ref_samples$cs, na.rm = T)
ci_ref_av = mean(ref_samples$ci, na.rm = T)
sdh_ref_av = mean(ref_samples$sdh, na.rm = T)
cox_ref_av = mean(ref_samples$cox, na.rm = T)


ref_samples <- ref_samples %>%
  group_by(plate, sampleName) %>%
  mutate(average_cs = mean(cs, na.rm = T), 
         average_ci = mean(ci, na.rm = T), 
         average_sdh = mean(sdh, na.rm = T), 
         average_cox = mean(cox, na.rm = T)) %>%
  dplyr::select(sampleName, plate, average_cs, average_ci, average_sdh, average_cox) %>%
  unique() %>%
  mutate(cs_adj = average_cs/cs_ref_av, ci_adj = average_ci/ci_ref_av, sdh_adj = average_sdh/sdh_ref_av, cox_adj = average_cox/cox_ref_av)


ref_samples <- ref_samples %>%
  group_by(plate) %>%
  mutate(cs_adj_factor = mean(cs_adj, , na.rm = T), 
         ci_adj_factor = mean(ci_adj, na.rm = T), 
         sdh_adj_factor = mean(sdh_adj, na.rm = T), 
         cox_adj_factor = mean(cox_adj, na.rm = T)) %>%
  select(plate, cs_adj_factor, ci_adj_factor, sdh_adj_factor, cox_adj_factor) %>%
  unique()


misbie_MHI_all_adj <- misbie_MHI_all %>%
  # filter(sampleName != "empty") %>%
  # filter(!Well %in% c("A1", "D6", "H12")) %>%
  mutate(cs_plate_adj = case_when(str_detect(plate, "PBMC_F.1") ~ cs/ref_samples$cs_adj_factor[ref_samples$plate == "PBMC_F.1"],
                                  str_detect(plate, "PBMC_S.1") ~ cs/ref_samples$cs_adj_factor[ref_samples$plate == "PBMC_S.1"],
                                  str_detect(plate, "lymp") ~ cs/ref_samples$cs_adj_factor[ref_samples$plate == "lymp"],
                                  str_detect(plate, "mono") ~ cs/ref_samples$cs_adj_factor[ref_samples$plate == "mono"],
                                  str_detect(plate, "neu") ~ cs/ref_samples$cs_adj_factor[ref_samples$plate == "neu"],
                                  str_detect(plate, "plt") ~ cs/ref_samples$cs_adj_factor[ref_samples$plate == "plt"])) %>%
  mutate(ci_plate_adj = case_when(str_detect(plate, "PBMC_F.1") ~ ci/ref_samples$ci_adj_factor[ref_samples$plate == "PBMC_F.1"],
                                  str_detect(plate, "PBMC_S.1") ~ ci/ref_samples$ci_adj_factor[ref_samples$plate == "PBMC_S.1"],
                                  str_detect(plate, "lymp") ~ ci/ref_samples$ci_adj_factor[ref_samples$plate == "lymp"],
                                  str_detect(plate, "mono") ~ ci/ref_samples$ci_adj_factor[ref_samples$plate == "mono"],
                                  str_detect(plate, "neu") ~ ci/ref_samples$ci_adj_factor[ref_samples$plate == "neu"],
                                  str_detect(plate, "plt") ~ ci/ref_samples$ci_adj_factor[ref_samples$plate == "plt"])) %>%
  mutate(sdh_plate_adj = case_when(str_detect(plate, "PBMC_F.1") ~ sdh/ref_samples$sdh_adj_factor[ref_samples$plate == "PBMC_F.1"],
                                   str_detect(plate, "PBMC_S.1") ~ sdh/ref_samples$sdh_adj_factor[ref_samples$plate == "PBMC_S.1"],
                                   str_detect(plate, "lymp") ~ sdh/ref_samples$sdh_adj_factor[ref_samples$plate == "lymp"],
                                   str_detect(plate, "mono") ~ sdh/ref_samples$sdh_adj_factor[ref_samples$plate == "mono"],
                                   str_detect(plate, "neu") ~ sdh/ref_samples$sdh_adj_factor[ref_samples$plate == "neu"],
                                   str_detect(plate, "plt") ~ sdh/ref_samples$sdh_adj_factor[ref_samples$plate == "plt"])) %>%
  mutate(cox_plate_adj = case_when(str_detect(plate, "PBMC_F.1") ~ cox/ref_samples$cox_adj_factor[ref_samples$plate == "PBMC_F.1"],
                                   str_detect(plate, "PBMC_S.1") ~ cox/ref_samples$cox_adj_factor[ref_samples$plate == "PBMC_S.1"],
                                   str_detect(plate, "lymp") ~ cox/ref_samples$cox_adj_factor[ref_samples$plate == "lymp"],
                                   str_detect(plate, "mono") ~ cox/ref_samples$cox_adj_factor[ref_samples$plate == "mono"],
                                   str_detect(plate, "neu") ~ cox/ref_samples$cox_adj_factor[ref_samples$plate == "neu"],
                                   str_detect(plate, "plt") ~ cox/ref_samples$cox_adj_factor[ref_samples$plate == "plt"])) %>%
  transform(cs_plate_adj = as.numeric(cs_plate_adj)) %>%
  transform(cs_plate_adj = as.numeric(cs_plate_adj)) %>%
  transform(sdh_plate_adj = as.numeric(sdh_plate_adj)) %>%
  transform(cox_plate_adj = as.numeric(cox_plate_adj))



# n normalization  
misbie_MHI_all_adj_noplt <- misbie_MHI_all_adj %>%
  filter(!(plate=="plt"))
nDNA_av = mean(misbie_MHI_all_adj_noplt$nNORM, na.rm = T)


misbie_MHI_all_adj <-  misbie_MHI_all_adj %>%
  mutate(nNORM_adj = nNORM/nDNA_av)

misbie_MHI_all_adj <- misbie_MHI_all_adj %>%
  mutate(cs_nNORM = cs_plate_adj*nNORM_adj) %>%
  mutate(ci_nNORM = ci_plate_adj*nNORM_adj) %>%
  mutate(sdh_nNORM = sdh_plate_adj*nNORM_adj) %>%
  mutate(cox_nNORM = cox_plate_adj*nNORM_adj)
# 
# 
#cell count normalization

MiSBIE_mono_mhi_cellcount <- readxl:: read_xlsx ("raw_data/MiSBIE_mhi_cellcount.xlsx") %>%
  dplyr::select(!cell_type) %>%
  transform(cell_count = as.numeric(cell_count)) %>%
  mutate(cell_adj = 5/cell_count)

misbie_MHI_all_adj <- misbie_MHI_all_adj %>%
  full_join(MiSBIE_mono_mhi_cellcount, by = c("sampleName", "plate"))
# %>%
#   mutate(cs_cell_adj = cs_plate_adj*cell_adj) %>%
#   mutate(ci_cell_adj = ci_plate_adj*cell_adj) %>%
#   mutate(sdh_cell_adj = sdh_plate_adj*cell_adj) %>%
#   mutate(cox_cell_adj = cox_plate_adj*cell_adj)
# 

#protein normalization

MiSBIE_mono_mhi_protein <- readxl:: read_xlsx ("interim_data/protein_4.10.xlsx") %>%
  rename(average_conc = `average conc`, sampleName = `subject ID`) %>%
  mutate(prot_norm = average_conc/(mean(average_conc, na.rm = T))) %>%
  rename(prot_conc = average_conc)
  # dplyr::select(!(average_conc))
  
  plt_prot_average = MiSBIE_mono_mhi_protein %>%
    filter(plate == "plt") %>%
    filter(!sampleName %in% c("neu_F.2_Mi068", "PBMC_S.1_Mi051")) %>%
    summarise(plt_average = mean(prot_conc, na.rm = T)) %>%
    pull(plt_average)
  
MiSBIE_mono_mhi_protein <- MiSBIE_mono_mhi_protein %>%
  mutate(plt_prot_norm = case_when(plate == "plt" ~ prot_conc/plt_prot_average,
                                   TRUE ~ NA)) %>%
  transform(plt_prot_norm = case_when(sampleName %in% c("neu_F.2_Mi068", "PBMC_S.1_Mi051") ~ NA,
                                   TRUE ~ plt_prot_norm))
  
    

misbie_MHI_all_adj <- misbie_MHI_all_adj %>%
  full_join(MiSBIE_mono_mhi_protein, by = c("sampleName", "plate")) %>%
  mutate(cs_prot_adj = cs_plate_adj/prot_norm) %>%
  mutate(ci_prot_adj = ci_plate_adj/prot_norm) %>%
  mutate(sdh_prot_adj = sdh_plate_adj/prot_norm) %>%
  mutate(cox_prot_adj = cox_plate_adj/prot_norm) %>%
  # platelet protein norm
  mutate(cs_prot_adj_plt = cs_plate_adj/plt_prot_norm) %>%
  mutate(ci_prot_adj_plt = ci_plate_adj/plt_prot_norm) %>%
  mutate(sdh_prot_adj_plt = sdh_plate_adj/plt_prot_norm) %>%
  mutate(cox_prot_adj_plt = cox_plate_adj/plt_prot_norm)
  

#mtDNA per protein conc for platelets

plt_mtDNA_prot <- all_pivot %>%
  select(sampleName, prot_conc_plt) %>%
  full_join(plt_qpcr_mtDNAct) %>%
  mutate(mtDNA_av_ct = (nd1_av_ct+cox1_av_ct)/2) %>%
  mutate(mtDNA_av_ct_lin = 2^(mtDNA_av_ct)) %>%
  mutate(mtDNA_over_prot = mtDNA_av_ct_lin/prot_conc_plt) %>%
  select(plate_ID, sampleName, mtDNA_over_prot)

misbie_MHI_all_adj <- misbie_MHI_all_adj %>%
  full_join(plt_mtDNA_prot)
# write_xlsx(plt_mtDNA_prot, "final_data/plt_mtDNA_over_protein_conc_4.30.xlsx")


# rename samples that were run on separate plates
# wrote it as PBMC_S.1 in the platemap but it is actually PBMC_F.1 

misbie_MHI_all_adj <- misbie_MHI_all_adj %>%
  mutate(cell_type = case_when(
                plate_ID == "plt_PBMC_S.1_Mi051" ~ "PBMC_F.1",
                plate_ID == "plt_neu_F.2_Mi068" ~ "neu",
                TRUE ~ plate
              )) %>%
  transform(sampleName = case_when(
    plate_ID == "plt_PBMC_S.1_Mi051" ~ "Mi051",
    plate_ID == "plt_neu_F.2_Mi068" ~ "Mi068",
    TRUE ~ sampleName
              )) %>%
  transform(plate_ID = case_when(
    plate_ID == "plt_PBMC_S.1_Mi051" ~ "PBMC_F.1_Mi051",
    plate_ID == "plt_neu_F.2_Mi068" ~ "neu_Mi068",
    TRUE ~ plate_ID
              )) %>%
  filter(!plate_ID =="PBMC_F.1_1.28.2020") 

# outlier removal (OR)

misbie_MHI_all_adj_outlier_removed <- misbie_MHI_all_adj %>%
  group_by(cell_type) %>%
  #nnorm OR
  mutate(cs_nnorm_Lower = quantile(cs_nNORM, probs = .25, na.rm = T)- 3*(IQR(cs_nNORM, na.rm = T))) %>%
  mutate (ci_nnorm_Lower = quantile(ci_nNORM, probs = .25, na.rm = T) - 3*(IQR(ci_nNORM, na.rm = T))) %>%
  mutate(sdh_nnorm_Lower = quantile(sdh_nNORM, probs = .25, na.rm = T) - 3*(IQR(sdh_nNORM, na.rm = T))) %>%
  mutate (cox_nnorm_Lower = quantile(cox_nNORM, probs = .25, na.rm = T) - 3*(IQR(cox_nNORM, na.rm = T))) %>%
  mutate(cs_nnorm_Upper = quantile(cs_nNORM, probs = .75, na.rm = T) + 3*(IQR(cs_nNORM, na.rm = T))) %>%
  mutate (ci_nnorm_Upper = quantile(ci_nNORM, probs = .75, na.rm = T) + 3*(IQR(ci_nNORM, na.rm = T))) %>%
  mutate (sdh_nnorm_Upper = quantile(sdh_nNORM, probs = .75, na.rm = T) + 3*(IQR(sdh_nNORM, na.rm = T))) %>%
  mutate (cox_nnorm_Upper = quantile(cox_nNORM, probs = .75, na.rm = T) + 3*(IQR(cox_nNORM, na.rm = T))) %>%
  mutate(
    cs_nNORM_OR = case_when(
      (cs_nNORM < cs_nnorm_Lower | cs_nNORM > cs_nnorm_Upper) ~ NA,
      TRUE ~ cs_nNORM),
    ci_nNORM_OR = case_when(
      (ci_nNORM < ci_nnorm_Lower | ci_nNORM > ci_nnorm_Upper) ~ NA,
      TRUE ~ ci_nNORM),
    sdh_nNORM_OR = case_when(
      (sdh_nNORM < sdh_nnorm_Lower | sdh_nNORM > sdh_nnorm_Upper) ~ NA,
      TRUE ~ sdh_nNORM),
    cox_nNORM_OR = case_when(
      (cox_nNORM < cox_nnorm_Lower | cox_nNORM > cox_nnorm_Upper) ~ NA,
      TRUE ~ cox_nNORM)
  ) %>%
  #protein OR (normalized to average of all)
  mutate(cs_prot_adj_Lower = quantile(cs_prot_adj, probs = .25, na.rm = T)- 3*(IQR(cs_prot_adj, na.rm = T))) %>%
  mutate (ci_prot_adj_Lower = quantile(ci_prot_adj, probs = .25, na.rm = T) - 3*(IQR(ci_prot_adj, na.rm = T))) %>%
  mutate(sdh_prot_adj_Lower = quantile(sdh_prot_adj, probs = .25, na.rm = T) - 3*(IQR(sdh_prot_adj, na.rm = T))) %>%
  mutate (cox_prot_adj_Lower = quantile(cox_prot_adj, probs = .25, na.rm = T) - 3*(IQR(cox_prot_adj, na.rm = T))) %>%
  mutate(cs_prot_adj_Upper = quantile(cs_prot_adj, probs = .75, na.rm = T) + 3*(IQR(cs_prot_adj, na.rm = T))) %>%
  mutate (ci_prot_adj_Upper = quantile(ci_prot_adj, probs = .75, na.rm = T) + 3*(IQR(ci_prot_adj, na.rm = T))) %>%
  mutate (sdh_prot_adj_Upper = quantile(sdh_prot_adj, probs = .75, na.rm = T) + 3*(IQR(sdh_prot_adj, na.rm = T))) %>%
  mutate (cox_prot_adj_Upper = quantile(cox_prot_adj, probs = .75, na.rm = T) + 3*(IQR(cox_prot_adj, na.rm = T))) %>%
  mutate(
    cs_prot_adj_OR = case_when(
      (cs_prot_adj < cs_prot_adj_Lower | cs_prot_adj > cs_prot_adj_Upper) ~ NA,
      TRUE ~ cs_prot_adj),
    ci_prot_adj_OR = case_when(
      (ci_prot_adj < ci_prot_adj_Lower | ci_prot_adj > ci_prot_adj_Upper) ~ NA,
      TRUE ~ ci_prot_adj),
    sdh_prot_adj_OR = case_when(
      (sdh_prot_adj < sdh_prot_adj_Lower | sdh_prot_adj > sdh_prot_adj_Upper) ~ NA,
      TRUE ~ sdh_prot_adj),
    cox_prot_adj_OR = case_when(
      (cox_prot_adj < cox_prot_adj_Lower | cox_prot_adj > cox_prot_adj_Upper) ~ NA,
      TRUE ~ cox_prot_adj)
  ) %>%
  #protein OR (normalized to average of plateletes) --PLATELETS ONLY
  mutate(cs_prot_adj_plt_Lower = quantile(cs_prot_adj_plt, probs = .25, na.rm = T)- 3*(IQR(cs_prot_adj_plt, na.rm = T))) %>%
  mutate (ci_prot_adj_plt_Lower = quantile(ci_prot_adj_plt, probs = .25, na.rm = T) - 3*(IQR(ci_prot_adj_plt, na.rm = T))) %>%
  mutate(sdh_prot_adj_plt_Lower = quantile(sdh_prot_adj_plt, probs = .25, na.rm = T) - 3*(IQR(sdh_prot_adj_plt, na.rm = T))) %>%
  mutate (cox_prot_adj_plt_Lower = quantile(cox_prot_adj_plt, probs = .25, na.rm = T) - 3*(IQR(cox_prot_adj_plt, na.rm = T))) %>%
  mutate(cs_prot_adj_plt_Upper = quantile(cs_prot_adj_plt, probs = .75, na.rm = T) + 3*(IQR(cs_prot_adj_plt, na.rm = T))) %>%
  mutate (ci_prot_adj_plt_Upper = quantile(ci_prot_adj_plt, probs = .75, na.rm = T) + 3*(IQR(ci_prot_adj_plt, na.rm = T))) %>%
  mutate (sdh_prot_adj_plt_Upper = quantile(sdh_prot_adj_plt, probs = .75, na.rm = T) + 3*(IQR(sdh_prot_adj_plt, na.rm = T))) %>%
  mutate (cox_prot_adj_plt_Upper = quantile(cox_prot_adj_plt, probs = .75, na.rm = T) + 3*(IQR(cox_prot_adj_plt, na.rm = T))) %>%
  mutate(
    cs_prot_adj_plt_OR = case_when(
      (cs_prot_adj_plt < cs_prot_adj_plt_Lower | cs_prot_adj_plt > cs_prot_adj_plt_Upper) ~ NA,
      TRUE ~ cs_prot_adj_plt),
    ci_prot_adj_plt_OR = case_when(
      (ci_prot_adj_plt < ci_prot_adj_plt_Lower | ci_prot_adj_plt > ci_prot_adj_plt_Upper) ~ NA,
      TRUE ~ ci_prot_adj_plt),
    sdh_prot_adj_plt_OR = case_when(
      (sdh_prot_adj_plt < sdh_prot_adj_plt_Lower | sdh_prot_adj_plt > sdh_prot_adj_plt_Upper) ~ NA,
      TRUE ~ sdh_prot_adj_plt),
    cox_prot_adj_plt_OR = case_when(
      (cox_prot_adj_plt < cox_prot_adj_plt_Lower | cox_prot_adj_plt > cox_prot_adj_plt_Upper) ~ NA,
      TRUE ~ cox_prot_adj_plt)
  ) %>%
  #mtDNAcn OR
  mutate(mtDNAcn_Lower = quantile(mtDNAcn, probs = .25, na.rm = T)- 3*(IQR(mtDNAcn, na.rm = T))) %>%
  mutate(mtDNAcn_Upper = quantile(mtDNAcn, probs = .75, na.rm = T) + 3*(IQR(mtDNAcn, na.rm = T))) %>%
  mutate(
    mtDNAcn_OR = case_when(
      (mtDNAcn < mtDNAcn_Lower | mtDNAcn > mtDNAcn_Upper) ~ NA,
      TRUE ~ mtDNAcn)) %>%
  #mtDNA normalized to protein conc OR --PLATELETS ONLY
  mutate(mtDNA_over_prot_Lower = quantile(mtDNA_over_prot, probs = .25, na.rm = T)- 3*(IQR(mtDNA_over_prot, na.rm = T))) %>%
  mutate(mtDNA_over_prot_Upper = quantile(mtDNA_over_prot, probs = .75, na.rm = T) + 3*(IQR(mtDNA_over_prot, na.rm = T))) %>%
  mutate(
    mtDNA_over_prot_OR = case_when(
      (mtDNA_over_prot < mtDNA_over_prot_Lower | mtDNA_over_prot > mtDNA_over_prot_Upper) ~ NA,
      TRUE ~ mtDNA_over_prot))



# means for mhi
means_formhi <- misbie_MHI_all_adj_outlier_removed %>%
  filter(!cell_type == "plt")
cs_prot_adj_OR_mean <- mean(means_formhi$cs_prot_adj_OR, na.rm = T)
ci_prot_adj_OR_mean <- mean(means_formhi$ci_prot_adj_OR, na.rm = T)
sdh_prot_adj_OR_mean <- mean(means_formhi$sdh_prot_adj_OR, na.rm = T)
cox_prot_adj_OR_mean <- mean(means_formhi$cox_prot_adj_OR, na.rm = T)
cs_nNORM_OR_mean <- mean(means_formhi$cs_nNORM_OR, na.rm = T)
ci_nNORM_OR_mean <- mean(means_formhi$ci_nNORM_OR, na.rm = T)
sdh_nNORM_OR_mean <- mean(means_formhi$sdh_nNORM_OR, na.rm = T)
cox_nNORM_OR_mean <- mean(means_formhi$cox_nNORM_OR, na.rm = T)
mtDNAcn_OR_mean <- mean(means_formhi$mtDNAcn_OR, na.rm = T)

# mhi calculations  
misbie_MHI_all_adj_OR <- misbie_MHI_all_adj_outlier_removed %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA)) %>%
  # mutate(z_CS = cs / mean(cs, na.rm = T)) %>%
  # mutate(z_CI = ci / mean(ci, na.rm = T)) %>%
  # mutate(z_SDH = sdh / mean(sdh, na.rm = T)) %>%
  # mutate(z_COX = cox / mean(cox, na.rm = T)) %>%
  #z plate adjust
  # mutate(z_CS_plate_adj = cs_plate_adj / mean(cs_plate_adj, na.rm = T)) %>%
  # mutate(z_CI_plate_adj = ci_plate_adj / mean(ci_plate_adj, na.rm = T)) %>%
  # mutate(z_SDH_plate_adj = sdh_plate_adj / mean(sdh_plate_adj, na.rm = T)) %>%
  # mutate(z_COX_plate_adj = cox_plate_adj / mean(cox_plate_adj, na.rm = T)) %>%
  # #z cell adjust
  # mutate(z_CS_cell_adj = cs_cell_adj / mean(cs_cell_adj, na.rm = T)) %>%
  # mutate(z_CI_cell_adj = ci_cell_adj / mean(ci_cell_adj, na.rm = T)) %>%
  # mutate(z_SDH_cell_adj = sdh_cell_adj / mean(sdh_cell_adj, na.rm = T)) %>%
  # mutate(z_COX_cell_adj = cox_cell_adj / mean(cox_cell_adj, na.rm = T)) %>%
  #z protein norm
  mutate(z_CS_prot_adj = 
           case_when(cell_type =="plt" ~ (cs_prot_adj_plt_OR / mean(cs_prot_adj_plt_OR, na.rm = T)), 
                     TRUE ~ (cs_prot_adj_OR / cs_prot_adj_OR_mean))) %>%
  mutate(z_CI_prot_adj = 
           case_when(cell_type =="plt" ~ (ci_prot_adj_plt_OR / mean(ci_prot_adj_plt_OR, na.rm = T)), 
                     TRUE ~ ci_prot_adj_OR / ci_prot_adj_OR_mean)) %>%
  mutate(z_SDH_prot_adj = case_when(cell_type =="plt" ~ (sdh_prot_adj_plt_OR / mean(sdh_prot_adj_plt_OR, na.rm = T)), 
                                    TRUE ~ sdh_prot_adj_OR / sdh_prot_adj_OR_mean)) %>%
  mutate(z_COX_prot_adj = case_when(cell_type =="plt" ~ (cox_prot_adj_plt_OR / mean(cox_prot_adj_plt_OR, na.rm = T)), 
                                    TRUE ~ cox_prot_adj_OR / cox_prot_adj_OR_mean)) %>%
  #z nnorm
  mutate(z_nNORM_CS = cs_nNORM_OR / cs_nNORM_OR_mean)%>%
  mutate(z_nNORM_CI = ci_nNORM_OR / ci_nNORM_OR_mean) %>%
  mutate(z_nNORM_SDH = sdh_nNORM_OR / sdh_nNORM_OR_mean) %>%
  mutate(z_nNORM_COX = cox_nNORM_OR / cox_nNORM_OR_mean) %>%
  mutate(z_mtDNAcn = case_when(cell_type == "plt" ~ mtDNA_over_prot_OR / mean(mtDNA_over_prot_OR, na.rm = T),
                               TRUE ~ mtDNAcn_OR / mtDNAcn_OR_mean)) %>%
  #calculate mrc/mhi
  rowwise() %>%
  mutate(MHI_nNORM = 100*(z_nNORM_CI+z_nNORM_SDH+z_nNORM_COX)/(z_nNORM_CS+z_mtDNAcn+1))  %>%
  mutate(MRC_nNORM = 100*(mean(c(z_nNORM_CI^(1/2),z_nNORM_SDH^(1/2),z_nNORM_COX^(1/2)), na.rm = T))/(mean(c(z_nNORM_CS^(1/3),z_mtDNAcn^(1/3)), na.rm = T))) %>%
  mutate(MHI_prot_adj = 100*(z_CI_prot_adj+z_SDH_prot_adj+z_COX_prot_adj)/(z_CS_prot_adj+z_mtDNAcn+1)) %>%
  mutate(MRC_prot_adj = 100*(mean(c(z_CI_prot_adj^(1/2),z_SDH_prot_adj^(1/2),z_COX_prot_adj^(1/2)), na.rm = T))/(mean(c(z_CS_prot_adj^(1/3),z_mtDNAcn^(1/3)), na.rm = T))) 

# %>%
#  
#   
#   
#   
#    rowwise() %>%
#   mutate(TRC = mean(c_across(ends_with("sqrt")), na.rm = T))
#          
# 
# test <- misbie_MHI_all_adj_OR %>% 
#   summarise_all(funs(class(.)))
# 
#     
#   mutate(MHI_prot_adj = 100*(z_CI_prot_adj+z_SDH_prot_adj+z_COX_prot_adj)/(z_CS_prot_adj+z_mtDNAcn+1)) %>%
#   mutate(MRC_prot_adj = 100*(mean(c(z_CI_prot_adj^(1/2),z_SDH_prot_adj^(1/2),z_COX_prot_adj^(1/2)), na.rm = T))/(mean(c(z_CS_prot_adj^(1/3),z_mtDNAcn^(1/3)), na.rm = T))) 

# how many outliers based on 3x IQR method?
count_outliers <- misbie_MHI_all_adj_OR %>%
  select(cell_type, sampleName, starts_with(c("cs","ci","sdh","cox","mtDNA"))) %>%
  select(!ends_with(c("Lower","Upper"))) %>%
  select(cell_type, sampleName, ends_with(c("nNORM","prot_adj", "prot_adj_plt", "nNORM_OR","prot_adj_OR", "prot_adj_plt_OR")), starts_with("mtDNA")) %>%
  group_by(cell_type) %>%
  mutate (cs_outliers = case_when(cell_type == 'plt' ~ (sum(is.na(cs_prot_adj_plt_OR))-sum(is.na(cs_prot_adj_plt))),
                                  TRUE ~ (sum(is.na(cs_nNORM_OR))-sum(is.na(cs_nNORM)))),
          ci_outliers = case_when(cell_type == 'plt' ~ (sum(is.na(ci_prot_adj_plt_OR))-sum(is.na(ci_prot_adj_plt))),
                                  TRUE ~ (sum(is.na(ci_nNORM_OR))-sum(is.na(ci_nNORM)))),
          sdh_outliers = case_when(cell_type == 'plt' ~ (sum(is.na(sdh_prot_adj_plt_OR))-sum(is.na(sdh_prot_adj_plt))),
                                   TRUE ~ (sum(is.na(sdh_nNORM_OR))-sum(is.na(sdh_nNORM)))),
          cox_outliers = case_when(cell_type == 'plt' ~ (sum(is.na(cox_prot_adj_plt_OR))-sum(is.na(cox_prot_adj_plt))),
                                   TRUE ~ (sum(is.na(cox_nNORM_OR))-sum(is.na(cox_nNORM)))),
          mtDNAcn_outliers = case_when(cell_type == 'plt' ~ (sum(is.na(mtDNA_over_prot_OR))-sum(is.na(mtDNA_over_prot))),
                                       TRUE ~ (sum(is.na(mtDNAcn_OR))-sum(is.na(mtDNAcn))))
          ) %>%
  mutate( cs_total = sum(!is.na(cs_nNORM)),
          ci_total = sum(!is.na(ci_nNORM)),
          sdh_total = sum(!is.na(sdh_nNORM)),
          cox_total = sum(!is.na(cox_nNORM)),
          mtDNAcn_total = sum(!is.na(mtDNAcn))
          
          ) %>%
  mutate ( cs_outlier_percent = cs_outliers/cs_total,
           ci_outlier_percent = ci_outliers/ci_total,
           sdh_outlier_percent = sdh_outliers/sdh_total,
           cox_outlier_percent = cox_outliers/cox_total,
           mtDNAcn_outlier_percent = mtDNAcn_outliers/mtDNAcn_total) %>%
  select(cell_type, ends_with("outliers"), ends_with("total"), ends_with("percent")) %>%
  unique() %>%
  mutate(cs_remaining = cs_total - cs_outliers,
         ci_remaining = ci_total - ci_outliers,
         sdh_remaining = sdh_total - sdh_outliers,
         cox_remaining = cox_total - cox_outliers,
         mtDNAcn_remaining = mtDNAcn_total - mtDNAcn_outliers) %>%
  select(c(cell_type, starts_with("cs"), starts_with("ci"), starts_with("sdh"), starts_with("cox"), starts_with("mtDNA")))
write_xlsx(count_outliers, "final_data/count_outliers_4.30_2.xlsx")  

misbie_MHI <- misbie_MHI_all_adj_OR %>%
  select(!c(cs, ci, sdh, cox)) %>%
  select(!(ends_with(c("Upper","Lower")))) %>%
  select(!starts_with("z")) %>%
  select(!cell_adj)

# write_xlsx(misbie_MHI,"final_data/misbie_MHI_long_4.26.xlsx")

# pivot
all_pivot <- misbie_MHI %>%
  dplyr::select(cell_type, sampleName, cell_count, prot_conc, prot_norm, nNORM, nNORM_adj,
                mtDNAcn, mtDNAcn_OR,
                mtDNA_over_prot, mtDNA_over_prot_OR,
                cs_plate_adj, ci_plate_adj, sdh_plate_adj, cox_plate_adj,
                cs_nNORM, ci_nNORM, sdh_nNORM, cox_nNORM, 
                cs_prot_adj, ci_prot_adj, sdh_prot_adj, cox_prot_adj,
                cs_prot_adj_plt, ci_prot_adj_plt, sdh_prot_adj_plt, cox_prot_adj_plt,
                cs_nNORM_OR, ci_nNORM_OR, sdh_nNORM_OR, cox_nNORM_OR, MHI_nNORM, MRC_nNORM,
                cs_prot_adj_OR, ci_prot_adj_OR, sdh_prot_adj_OR, cox_prot_adj_OR,
                cs_prot_adj_plt_OR, ci_prot_adj_plt_OR, sdh_prot_adj_plt_OR, cox_prot_adj_plt_OR,
                MHI_prot_adj, MRC_prot_adj,
                # cs_nNORM, ci_nNORM, sdh_nNORM, cox_nNORM, MHI_nNORM, MRC_nNORM
                # ,cs, ci, sdh, cox, mtDNAcn, nNORM, MHI, MRC, 
                # cs_cell_adj, ci_cell_adj, sdh_cell_adj, cox_cell_adj, MHI_cell_adj, MRC_cell_adj, 
                # cs_prot_adj, ci_prot_adj, sdh_prot_adj, cox_prot_adj, MHI_prot_adj, MRC_prot_adj, 
               # MHI_plate_adj, MRC_plate_adj
               ) %>%
  pivot_wider(names_from = cell_type, values_from = c(cell_count, prot_conc, prot_norm, nNORM, nNORM_adj,
                                                      mtDNAcn, mtDNA_over_prot,
                                                      cs_plate_adj, ci_plate_adj, sdh_plate_adj, cox_plate_adj,
                                                      cs_nNORM, ci_nNORM, sdh_nNORM, cox_nNORM, 
                                                      cs_prot_adj, ci_prot_adj, sdh_prot_adj, cox_prot_adj, 
                                                      cs_prot_adj_plt, ci_prot_adj_plt, sdh_prot_adj_plt, cox_prot_adj_plt,
                                                      mtDNAcn_OR, mtDNA_over_prot_OR,
                                                      cs_nNORM_OR, ci_nNORM_OR, sdh_nNORM_OR, cox_nNORM_OR, MHI_nNORM, MRC_nNORM,
                                                      cs_prot_adj_OR, ci_prot_adj_OR, sdh_prot_adj_OR, cox_prot_adj_OR,
                                                      cs_prot_adj_plt_OR, ci_prot_adj_plt_OR, sdh_prot_adj_plt_OR, cox_prot_adj_plt_OR,
                                                      MHI_prot_adj, MRC_prot_adj,
  )) 

write_xlsx(misbie_MHI,"final_data/misbie_MHI_long_4.30_2.xlsx")
write_xlsx(all_pivot,"final_data/misbie_MHI_wide_4.30_2.xlsx")





# format for redcap


raw_data_pivot <- misbie_MHI %>%
  dplyr::select(cell_type, sampleName, cell_count, prot_conc, prot_norm, nNORM, nNORM_adj,
                mtDNAcn, mtDNA_over_prot,
                cs_plate_adj, ci_plate_adj, sdh_plate_adj, cox_plate_adj,
                cs_nNORM, ci_nNORM, sdh_nNORM, cox_nNORM, 
                cs_prot_adj, ci_prot_adj, sdh_prot_adj, cox_prot_adj,
                cs_prot_adj_plt, ci_prot_adj_plt, sdh_prot_adj_plt, cox_prot_adj_plt,
                  ) %>%
  rename_with(~paste0(., "_x"), ends_with("_plt")) %>%
  transform(cell_type = 
              case_when(
                cell_type =="PBMC_F.1" ~ "pbmc_f.1",
                cell_type == "PBMC_S.1" ~ "pbmc_s.1",
                TRUE ~ cell_type
              )) %>%
  rename(participant_id = sampleName, mtdnacn = mtDNAcn, mtdna_over_prot = mtDNA_over_prot, nnorm = nNORM, nnorm_adj = nNORM_adj) %>%
  rename(cs_nnorm = cs_nNORM, ci_nnorm = ci_nNORM, sdh_nnorm = sdh_nNORM, cox_nnorm = cox_nNORM) %>%
  pivot_wider(names_from = cell_type, values_from = c(cell_count, prot_conc, prot_norm, nnorm, nnorm_adj,
                                                      mtdnacn, mtdna_over_prot,
                                                      cs_plate_adj, ci_plate_adj, sdh_plate_adj, cox_plate_adj,
                                                      cs_nnorm, ci_nnorm, sdh_nnorm, cox_nnorm, 
                                                      cs_prot_adj, ci_prot_adj, sdh_prot_adj, cox_prot_adj, 
                                                      cs_prot_adj_plt_x, ci_prot_adj_plt_x, sdh_prot_adj_plt_x, cox_prot_adj_plt_x
                                                    )) %>%
  select(!c(cs_prot_adj_plt_x_pbmc_f.1,cs_prot_adj_plt_x_pbmc_s.1, cs_prot_adj_plt_x_lymp, cs_prot_adj_plt_x_mono, cs_prot_adj_plt_x_neu,
            ci_prot_adj_plt_x_pbmc_f.1, ci_prot_adj_plt_x_pbmc_s.1, ci_prot_adj_plt_x_lymp, ci_prot_adj_plt_x_mono, ci_prot_adj_plt_x_neu,
            sdh_prot_adj_plt_x_pbmc_f.1, sdh_prot_adj_plt_x_pbmc_s.1, sdh_prot_adj_plt_x_lymp, sdh_prot_adj_plt_x_mono, sdh_prot_adj_plt_x_neu,
            cox_prot_adj_plt_x_pbmc_f.1, cox_prot_adj_plt_x_pbmc_s.1, cox_prot_adj_plt_x_lymp, cox_prot_adj_plt_x_mono, cox_prot_adj_plt_x_neu,
            mtdna_over_prot_pbmc_f.1, mtdna_over_prot_pbmc_s.1, mtdna_over_prot_lymp, mtdna_over_prot_mono, mtdna_over_prot_neu))
    





final_data_pivot <- misbie_MHI %>%
  dplyr::select(cell_type, sampleName, 
                mtDNAcn_OR, mtDNA_over_prot_OR,
                cs_nNORM_OR, ci_nNORM_OR, sdh_nNORM_OR, cox_nNORM_OR, MHI_nNORM, MRC_nNORM,
                # cs_prot_adj_OR, ci_prot_adj_OR, sdh_prot_adj_OR, cox_prot_adj_OR, 
                cs_prot_adj_plt_OR, ci_prot_adj_plt_OR, sdh_prot_adj_plt_OR, cox_prot_adj_plt_OR,
                MHI_prot_adj, MRC_prot_adj
  ) %>%
  transform(cell_type = 
              case_when(
                cell_type =="PBMC_F.1" ~ "pbmc_f.1",
                cell_type == "PBMC_S.1" ~ "pbmc_s.1",
                TRUE ~ cell_type
              )) %>%
  mutate(final_cs = case_when(
    cell_type == 'plt' ~ cs_prot_adj_plt_OR,
    TRUE ~ cs_nNORM_OR,
  )) %>%
  mutate(final_ci = case_when(
    cell_type == 'plt' ~ ci_prot_adj_plt_OR,
    TRUE ~ ci_nNORM_OR,
  )) %>%
  mutate(final_sdh = case_when(
    cell_type == 'plt' ~ sdh_prot_adj_plt_OR,
    TRUE ~ sdh_nNORM_OR,
  )) %>%
  mutate(final_cox = case_when(
    cell_type == 'plt' ~ cox_prot_adj_plt_OR,
    TRUE ~ cox_nNORM_OR,
  )) %>%
  mutate(final_mtdnacn = case_when(
    cell_type == 'plt' ~ mtDNA_over_prot_OR,
    TRUE ~ mtDNAcn_OR,
  )) %>%
  mutate(final_mrc = case_when(
    cell_type == 'plt' ~ MRC_prot_adj,
    TRUE ~ MRC_nNORM,
  )) %>%
  mutate(final_mhi = case_when(
    cell_type == 'plt' ~ MHI_prot_adj,
    TRUE ~ MHI_nNORM,
  )) %>%
  select(cell_type, sampleName, starts_with("final")) %>%
  rename(participant_id = sampleName) %>%
  pivot_wider(names_from = cell_type, values_from = c(final_mtdnacn, final_cs, final_ci, final_sdh, final_cox, final_mrc, final_mhi)) %>%
  rename_with(~paste0(., "_c"), starts_with("final")) %>%
  rename_with(~str_remove(., "final_"), starts_with("final"))



write_xlsx (raw_data_pivot, "final_data/misbie_MHI_raw_4.30_2.xlsx")
write_xlsx (final_data_pivot, "final_data/misbie_MHI_final_4.30_withmhi_2.xlsx")













# {misbie_MHI} %>%
#   dplyr::group_by(sampleName, cell_type) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1L) 
# 
# length(misbie_MHI_all_adj$z_CS)

# write_xlsx(misbie_MHI_all_adj, "final_data/misbie_MHI_all_4.24.xlsx")

# 
# # biological duplicates
# 
# misbie_counts <- readxl::read_xlsx("plate_maps/misbie_counts.xlsx")
# 
# misbie_manifest_ID <- misbie_counts %>%
#   unite(id_all, subject_id, timepoints,box_no, box_loc, sep = "_", remove = F)%>%
#   mutate(id_all_short = str_replace_all(id_all, pattern = ("Cares mtDNA"), replacement = "")) %>%
#   mutate(sampleName = str_replace_all(id_all_short, pattern = ("01-817-"), replacement = "")) %>%
#   select(sampleName,subject_id,draw_date,timepoints,box_no, box_loc)
# 
# misbie_bio_dup <- misbie_MHI_all_adj %>%
#   filter(!Well %in% c("A1", "D6", "H12")) %>%
#   full_join(misbie_manifest_ID,.) %>%
#   group_by(subject_id, draw_date, timepoints) %>%
#   summarise(count = n(), across(c(box_no, box_loc, Well, sampleName, 
#                                   cs, ci, sdh, cox, mtDNAcn, nNORM, 
#                                   cs_plate_adj, ci_plate_adj, sdh_plate_adj, cox_plate_adj, 
#                                   cs_nNORM, ci_nNORM, sdh_nNORM, cox_nNORM, nNORM_adj, 
#                                   MHI, MHI_plate_adj, MHI_nNORM, MRC, MRC_plate_adj, MRC_nNORM), list)) %>%
#   unnest_wider(c(box_no, box_loc, Well,
#                  sampleName, 
#                  cs, ci, sdh, cox, mtDNAcn, nNORM, 
#                  cs_plate_adj, ci_plate_adj, sdh_plate_adj, cox_plate_adj, 
#                  cs_nNORM, ci_nNORM, sdh_nNORM, cox_nNORM, nNORM_adj, 
#                  MHI, MHI_plate_adj, MHI_nNORM, MRC, MRC_plate_adj, MRC_nNORM), 
#                names_sep = "_")
# 
# write_xlsx(misbie_bio_dup, "final_data/misbie_bio_dup.xlsx")


          



