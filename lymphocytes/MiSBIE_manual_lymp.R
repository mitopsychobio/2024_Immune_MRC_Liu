library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)

## read all files

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
  

#combine all activities activity in units mmol/min/10^6 cells

lymp_MHI_all <- full_join(lymp_cs, lymp_ci, by = c("plate_ID", "Well", "plate", "sampleName")) %>%
  # dplyr::select("plate_ID", "Well", "plate", "sampleName", "cs_cell", "ci_cell") %>%
  full_join(lymp_sdh) %>%
  # dplyr::select(-"sdh") %>%
  full_join(lymp_cox) %>%
  # dplyr::select(-"cox") %>%
  full_join (lymp_qpcr) %>%
  # rename(cs = cs_cell, ci = ci_cell, sdh = sdh_cell, cox = cox_cell) %>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  ) %>%
  mutate(Column = as.numeric(Column)) %>%
  arrange(plate, Column)

# data cleaning??

#these are empty wells

cs_excl <- c()
ci_excl <- c()
sdh_excl <- c()
cox_excl <- c()
qPCR_excl <- c()

lymp_MHI_all <- lymp_MHI_all %>%
  mutate(cs = case_when(
    plate_ID %in% cs_excl ~ NA,
    TRUE ~ cs
  )) %>%
  mutate(ci = case_when(
    plate_ID %in% ci_excl ~ NA,
    TRUE ~ ci
  )) %>%
  mutate(sdh = case_when(
    plate_ID %in% sdh_excl ~ NA,
    TRUE ~ sdh
  )) %>%
  mutate(cox = case_when(
    plate_ID %in% cs_excl ~ NA,
    TRUE ~ cox
  )) %>%
  mutate(mtDNAcn = case_when(
    plate_ID %in% qPCR_excl ~ NA,
    TRUE ~ mtDNAcn
  )) %>%
  mutate(nNORM = case_when(
    plate_ID %in% qPCR_excl ~ NA,
    TRUE ~ nNORM
  )) 

# plate adjustment

## from excel
sampleName <- c("REF_1","REF_2","REF_3")
cs_x <- c(36.01,	35.10,	30.69)
ci_x <- c(7.52,	5.12,	6.88)
sdh_x <- c(3.07,	1.87,	2.55)
cox_x <- c(27.18,	16.39,	20.54)
averages <- data.frame(sampleName, cs_x, ci_x, sdh_x, cox_x)


ref_samples <-  lymp_MHI_all %>%
  filter(sampleName %in% c("REF_1", "REF_2", "REF_3")) %>%
  dplyr::select(plate, sampleName, cs, ci, sdh, cox) %>%
  group_by(sampleName) %>%
  mutate(average_cs = mean(cs), average_ci = mean(ci), average_sdh = mean(sdh), average_cox = mean(cox)) %>%
  dplyr::select(sampleName, average_cs, average_ci, average_sdh, average_cox) %>%
  unique() %>%
  full_join(averages) %>%
  mutate(cs_adj = average_cs/cs_x, ci_adj = average_ci/ci_x, sdh_adj = average_sdh/sdh_x, cox_adj = average_cox/cox_x)

cs_adj_factor <- mean(ref_samples$cs_adj, na.rm = T)
ci_adj_factor <- mean(ref_samples$ci_adj, na.rm = T)
sdh_adj_factor <- mean(ref_samples$sdh_adj, na.rm = T)
cox_adj_factor <- mean(ref_samples$cox_adj, na.rm = T)

lymp_MHI_all <- lymp_MHI_all %>%
  dplyr::select("plate_ID", "Well", "plate", "sampleName", "cs_cell", "ci_cell", "sdh_cell","cox_cell", "mtDNAcn", "nNORM") %>%
  rename(cs = cs_cell, ci = ci_cell, sdh = sdh_cell, cox = cox_cell)

lymp_MHI_all <- lymp_MHI_all %>%
  mutate(cs_plate_adj = cs/cs_adj_factor, cs) %>%
  mutate(ci_plate_adj = ci/ci_adj_factor, ci) %>%
  mutate(sdh_plate_adj = sdh/sdh_adj_factor, sdh) %>%
  mutate(cox_plate_adj = cox/cox_adj_factor, cox)


# n normalization  
nDNA_av = mean(lymp_MHI_all$nNORM, na.rm = T)


lymp_MHI_all <-  lymp_MHI_all %>%
  mutate(nNORM_adj = nNORM/nDNA_av)

lymp_MHI_all <- lymp_MHI_all %>%
  mutate(cs_nNORM = cs_plate_adj*nNORM_adj) %>%
  mutate(ci_nNORM = ci_plate_adj*nNORM_adj) %>%
  mutate(sdh_nNORM = sdh_plate_adj*nNORM_adj) %>%
  mutate(cox_nNORM = cox_plate_adj*nNORM_adj)


#cell count normalization

MiSBIE_lymp_mhi_cellcount <- readxl:: read_xlsx ("raw_data/MiSBIE_mhi_cellcount.xlsx") %>%
  filter(plate == "lymp") %>%
  dplyr::select(!cell_type) %>%
  transform(cell_count = as.numeric(cell_count)) %>%
  mutate(cell_adj = 5/cell_count)

lymp_MHI_all <- lymp_MHI_all %>%
  full_join(MiSBIE_lymp_mhi_cellcount, by = c("sampleName", "plate")) %>%
  mutate(cs_cell_adj = cs_plate_adj*cell_adj) %>%
  mutate(ci_cell_adj = ci_plate_adj*cell_adj) %>%
  mutate(sdh_cell_adj = sdh_plate_adj*cell_adj) %>%
  mutate(cox_cell_adj = cox_plate_adj*cell_adj)

#protein normalization

MiSBIE_lymp_mhi_protein <- readxl:: read_xlsx ("interim_data/protein.xlsx") %>%
  filter(plate == "lymp") %>%
  rename(average_conc = `average conc`, sampleName = `subject ID`) %>%
  mutate(prot_norm = average_conc/mean(average_conc)) %>%
  dplyr::select(!c(plate, average_conc))

lymp_MHI_all <- lymp_MHI_all %>%
  full_join(MiSBIE_lymp_mhi_protein) %>%
  mutate(cs_prot_adj = cs_plate_adj/prot_norm) %>%
  mutate(ci_prot_adj = ci_plate_adj/prot_norm) %>%
  mutate(sdh_prot_adj = sdh_plate_adj/prot_norm) %>%
  mutate(cox_prot_adj = cox_plate_adj/prot_norm)

# mhi calculations  
lymp_MHI_all <- lymp_MHI_all %>%
  filter(!sampleName %in% c("REF_1", "REF_2", "REF_3",NA)) %>%
  mutate(z_CS = cs / mean(cs, na.rm = T)) %>%
  mutate(z_CI = ci / mean(ci, na.rm = T)) %>%
  mutate(z_SDH = sdh / mean(sdh, na.rm = T)) %>%
  mutate(z_COX = cox / mean(cox, na.rm = T)) %>%
  #z plate adjust
  mutate(z_CS_plate_adj = cs_plate_adj / mean(cs_plate_adj, na.rm = T)) %>%
  mutate(z_CI_plate_adj = ci_plate_adj / mean(ci_plate_adj, na.rm = T)) %>%
  mutate(z_SDH_plate_adj = sdh_plate_adj / mean(sdh_plate_adj, na.rm = T)) %>%
  mutate(z_COX_plate_adj = cox_plate_adj / mean(cox_plate_adj, na.rm = T)) %>%
  #z cell adjust
  mutate(z_CS_cell_adj = cs_cell_adj / mean(cs_cell_adj, na.rm = T)) %>%
  mutate(z_CI_cell_adj = ci_cell_adj / mean(ci_cell_adj, na.rm = T)) %>%
  mutate(z_SDH_cell_adj = sdh_cell_adj / mean(sdh_cell_adj, na.rm = T)) %>%
  mutate(z_COX_cell_adj = cox_cell_adj / mean(cox_cell_adj, na.rm = T)) %>%
  #z protein norm 
  mutate(z_CS_prot_adj = cs_prot_adj / mean(cs_prot_adj, na.rm = T)) %>%
  mutate(z_CI_prot_adj = ci_prot_adj / mean(ci_prot_adj, na.rm = T)) %>%
  mutate(z_SDH_prot_adj = sdh_prot_adj / mean(sdh_prot_adj, na.rm = T)) %>%
  mutate(z_COX_prot_adj = cox_prot_adj / mean(cox_prot_adj, na.rm = T)) %>%
  #z nnorm
  mutate(z_nNORM_CS = cs_nNORM / mean(cs_nNORM, na.rm = T))%>%
  mutate(z_nNORM_CI = ci_nNORM / mean(ci_nNORM, na.rm = T)) %>%
  mutate(z_nNORM_SDH = sdh_nNORM / mean(sdh_nNORM, na.rm = T)) %>%
  mutate(z_nNORM_COX = cox_nNORM / mean(cox_nNORM, na.rm = T)) %>%
  mutate(z_mtDNAcn = mtDNAcn / mean(mtDNAcn, na.rm = T)) %>%
  mutate(MHI = 100*(z_CI+z_SDH+z_COX)/(z_CS+z_mtDNAcn+1)) %>%
  mutate(MRC = (2/3)*100*(z_CI^(1/2)+z_SDH^(1/2)+z_COX^(1/2))/(z_CS^(1/3)+z_mtDNAcn^(1/3))) %>%
  mutate(MHI_plate_adj = 100*(z_CI_plate_adj+z_SDH_plate_adj+z_COX_plate_adj)/(z_CS_plate_adj+z_mtDNAcn+1)) %>%
  mutate(MRC_plate_adj = (2/3)*100*(z_CI_plate_adj^(1/2)+z_SDH_plate_adj^(1/2)+z_COX_plate_adj^(1/2))/(z_CS_plate_adj^(1/3)+z_mtDNAcn^(1/3))) %>%
  mutate(MHI_cell_adj = 100*(z_CI_cell_adj+z_SDH_cell_adj+z_COX_cell_adj)/(z_CS_cell_adj+z_mtDNAcn+1)) %>%
  mutate(MRC_cell_adj = (2/3)*100*(z_CI_cell_adj^(1/2)+z_SDH_cell_adj^(1/2)+z_COX_cell_adj^(1/2))/(z_CS_cell_adj^(1/3)+z_mtDNAcn^(1/3))) %>%
  mutate(MHI_nNORM = 100*(z_nNORM_CI+z_nNORM_SDH+z_nNORM_COX)/(z_nNORM_CS+z_mtDNAcn+1))  %>%
  mutate(MRC_nNORM = (2/3)*100*(z_nNORM_CI^(1/2)+z_nNORM_SDH^(1/2)+z_nNORM_COX^(1/2))/(z_nNORM_CS^(1/3)+z_mtDNAcn^(1/3))) %>%
  mutate(MRC_cell_adj = (2/3)*100*(z_CI_cell_adj^(1/2)+z_SDH_cell_adj^(1/2)+z_COX_cell_adj^(1/2))/(z_CS_cell_adj^(1/3)+z_mtDNAcn^(1/3))) %>%
  mutate(MHI_prot_adj = 100*(z_CI_prot_adj+z_SDH_prot_adj+z_COX_prot_adj)/(z_CS_prot_adj+z_mtDNAcn+1)) %>%
  mutate(MRC_prot_adj = (2/3)*100*(z_CI_prot_adj^(1/2)+z_SDH_prot_adj^(1/2)+z_COX_prot_adj^(1/2))/(z_CS_prot_adj^(1/3)+z_mtDNAcn^(1/3)))

write_xlsx(lymp_MHI_all, "final_data/manual/lymp_MHI_all_manual.xlsx")
