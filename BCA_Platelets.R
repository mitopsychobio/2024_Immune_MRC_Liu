library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)

# BCA platemaps


MiSBIE_BCA_plt <- readxl::read_xlsx("plate_maps/MHI_platemap_MiSBIE_BCA.xlsx", sheet = 4)%>%
  pivot_longer(cols = -...1,names_to = "Column", values_to = "sampleName")%>%
  rename(Row = ...1) %>%
  mutate(across('Row', str_replace, 'Row', '')) %>%
  mutate(across('Column', str_replace, 'C', '')) %>%
  unite(Well, Row, Column, sep = "", remove = F) %>%
  mutate(plate = "platelet_F.2")

# write_xlsx(MiSBIE_BCA_plt, "plate_maps/MiSBIE_BCA_plt.xlsx")

MiSBIE_BCA <- readxl::read_xlsx ("raw_data/MiSBIE_MHI_BCA_x.xlsx")


##plt data

plt_BCA_A <- MiSBIE_BCA[62:63, 3:98] %>%
  t()
plt_BCA_A <-  transform(as.data.frame(plt_BCA_A)) %>%
  rename(Well = V1, OD_A = V2) %>%
  full_join(MiSBIE_BCA_plt)


plt_BCA_B <- MiSBIE_BCA[67:68, 3:98] %>%
  t()
plt_BCA_B <-  transform(as.data.frame(plt_BCA_B)) %>%
  rename(Well = V1, OD_B = V2) %>%
  full_join(MiSBIE_BCA_plt)

plt_BCA <- full_join(plt_BCA_A, plt_BCA_B) %>%
  select(plate, Well, Row, Column, sampleName, OD_A, OD_B)

# blank adjustment 
plt_BCA_blank <- plt_BCA %>%
  mutate(OD_A = as.numeric(OD_A)) %>%
  mutate(OD_B = as.numeric(OD_B)) %>%
  filter(sampleName == "BLANK")

plt_OD_A_blank <- mean(plt_BCA_blank$OD_A)
plt_OD_B_blank <- mean(plt_BCA_blank$OD_B)
rm(plt_BCA_blank)

plt_BCA_bl_corr <- plt_BCA %>%
  mutate(OD_A = as.numeric(OD_A)) %>%
  mutate(OD_B = as.numeric(OD_B)) %>%
  mutate(OD_A_corr = OD_A - plt_OD_A_blank) %>%
  mutate(OD_B_corr = OD_B - plt_OD_B_blank)
rm(plt_OD_A_blank, plt_OD_B_blank)

write_xlsx(plt_BCA_bl_corr, "interim_data/BCA/plt_BCA.xlsx")
