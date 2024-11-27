library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)



#read file, get rid of temperature column ----

PBMC_F.1_cox <- read.delim ("raw_data/spectramax_data/PBMC_F.1/COX-PBMC_F.1.xls", header = TRUE, sep = "\t", fileEncoding = "UTF-16LE", skip = 2)[,-2]
PBMC_F.1 <- readxl::read_xlsx("plate_maps/MiSBIE_PBMC_F.1.xlsx", sheet = 1)



################### POSITIVE 1 ########################

#changing spectramax data into long format for 1 plate (p1 in this case)
PBMC_F.1_cox_p1 <- PBMC_F.1_cox[1:27, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_F.1_cox_p1 <- PBMC_F.1_cox_p1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_F.1_cox_p1$Column <- as.numeric(PBMC_F.1_cox_p1$Column)
PBMC_F.1_cox_p1$mOD <- as.numeric(PBMC_F.1_cox_p1$mOD)

PBMC_F.1_cox_p1$Time <- period_to_seconds(hms(PBMC_F.1_cox_p1$Time))

PBMC_F.1_cox_p1 <- arrange(PBMC_F.1_cox_p1, Row)
PBMC_F.1_cox_p1 <- arrange(PBMC_F.1_cox_p1, Column)


#plot all traces as interactive
p1 <- ggplot(data = PBMC_F.1_cox_p1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row)) +
  geom_point(aes(color = Row))

plotly::ggplotly(p1)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 598

t2a <- 0
t2b <- 299

PBMC_F.1_cox_p1_excl <- c()
PBMC_F.1_cox_p1_it_diff <- c("D6", "A1", "H1", "H12")

PBMC_F.1_cox_p1_it1 <- PBMC_F.1_cox_p1[PBMC_F.1_cox_p1$Time>=t1a & PBMC_F.1_cox_p1$Time<=t1b,] 
PBMC_F.1_cox_p1_it2 <- PBMC_F.1_cox_p1[PBMC_F.1_cox_p1$Time>=t2a & PBMC_F.1_cox_p1$Time<=t2b,] 

#time1--
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_p1_models <- PBMC_F.1_cox_p1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_p1_models <- PBMC_F.1_cox_p1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_p1_models2 <- PBMC_F.1_cox_p1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_p1_models2 <- PBMC_F.1_cox_p1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times- 
PBMC_F.1_cox_p1_it1_mod <- select(PBMC_F.1_cox_p1_models, Well, slope)
PBMC_F.1_cox_p1_it2_mod <- select(PBMC_F.1_cox_p1_models2, Well, slope)


PBMC_F.1_cox_p1_allmod <- full_join(PBMC_F.1_cox_p1_it1_mod, PBMC_F.1_cox_p1_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_F.1_cox_p1_excl ~ NA,
    Well%in%PBMC_F.1_cox_p1_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_F.1_cox_p1_allmod$final_it_act <- -1*1000*60*(PBMC_F.1_cox_p1_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_F.1_cox_p1_allmod$Well)
PBMC_F.1_cox_p1_allmod <- PBMC_F.1_cox_p1_allmod[well_order, ]

#extract activities per well
PBMC_F.1_cox_p1_act <- select(PBMC_F.1_cox_p1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)

rm(p1, t1a, t1b, t2a, t2b,
   PBMC_F.1_cox_p1_excl, PBMC_F.1_cox_p1_it_diff, 
   PBMC_F.1_cox_p1_it1, PBMC_F.1_cox_p1_it2, 
   PBMC_F.1_cox_p1_models, PBMC_F.1_cox_p1_models2,
   PBMC_F.1_cox_p1_it1_mod, PBMC_F.1_cox_p1_it2_mod,
   PBMC_F.1_cox_p1_allmod
)




# ################## POSITIVE 2 ########################

#changing spectramax data into long format for 1 plate (p2 in this case)
PBMC_F.1_cox_p2 <- PBMC_F.1_cox[31:57, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_F.1_cox_p2 <- PBMC_F.1_cox_p2[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_F.1_cox_p2$Column <- as.numeric(PBMC_F.1_cox_p2$Column)
PBMC_F.1_cox_p2$mOD <- as.numeric(PBMC_F.1_cox_p2$mOD)

PBMC_F.1_cox_p2$Time <- period_to_seconds(hms(PBMC_F.1_cox_p2$Time))

PBMC_F.1_cox_p2 <- arrange(PBMC_F.1_cox_p2, Row)
PBMC_F.1_cox_p2 <- arrange(PBMC_F.1_cox_p2, Column)


#plot all traces as interactive
p2 <- ggplot(data = PBMC_F.1_cox_p2, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p2)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 598

t2a <- 0
t2b <- 230

PBMC_F.1_cox_p2_excl <- c()
PBMC_F.1_cox_p2_it_diff <- c("D6", "A1", "H1", "H12")

PBMC_F.1_cox_p2_it1 <- PBMC_F.1_cox_p2[PBMC_F.1_cox_p2$Time>=t1a & PBMC_F.1_cox_p2$Time<=t1b,] 
PBMC_F.1_cox_p2_it2 <- PBMC_F.1_cox_p2[PBMC_F.1_cox_p2$Time>=t2a & PBMC_F.1_cox_p2$Time<=t2b,] 

#time1--
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_p2_models <- PBMC_F.1_cox_p2_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_p2_models <- PBMC_F.1_cox_p2_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_p2_models2 <- PBMC_F.1_cox_p2_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_p2_models2 <- PBMC_F.1_cox_p2_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times- 
PBMC_F.1_cox_p2_it1_mod <- select(PBMC_F.1_cox_p2_models, Well, slope)
PBMC_F.1_cox_p2_it2_mod <- select(PBMC_F.1_cox_p2_models2, Well, slope)

PBMC_F.1_cox_p2_allmod <- full_join(PBMC_F.1_cox_p2_it1_mod, PBMC_F.1_cox_p2_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_F.1_cox_p2_excl ~ NA,
    Well%in%PBMC_F.1_cox_p2_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_F.1_cox_p2_allmod$final_it_act <- -1*1000*60*(PBMC_F.1_cox_p2_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_F.1_cox_p2_allmod$Well)
PBMC_F.1_cox_p2_allmod <- PBMC_F.1_cox_p2_allmod[well_order, ]

#extract activities per well
PBMC_F.1_cox_p2_act <- select(PBMC_F.1_cox_p2_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(p2, t1a, t1b, t2a, t2b,
   PBMC_F.1_cox_p2_excl, PBMC_F.1_cox_p2_it_diff, 
   PBMC_F.1_cox_p2_it1, PBMC_F.1_cox_p2_it2, 
   PBMC_F.1_cox_p2_models, PBMC_F.1_cox_p2_models2,
   PBMC_F.1_cox_p2_it1_mod, PBMC_F.1_cox_p2_it2_mod,
   PBMC_F.1_cox_p2_allmod
)



# ################## POSITIVE 3 ########################

#changing spectramax data into long format for 1 plate (p3 in this case)
PBMC_F.1_cox_p3 <- PBMC_F.1_cox[61:87, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_F.1_cox_p3 <- PBMC_F.1_cox_p3[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_F.1_cox_p3$Column <- as.numeric(PBMC_F.1_cox_p3$Column)
PBMC_F.1_cox_p3$mOD <- as.numeric(PBMC_F.1_cox_p3$mOD)

PBMC_F.1_cox_p3$Time <- period_to_seconds(hms(PBMC_F.1_cox_p3$Time))

PBMC_F.1_cox_p3 <- arrange(PBMC_F.1_cox_p3, Row)
PBMC_F.1_cox_p3 <- arrange(PBMC_F.1_cox_p3, Column)


#plot all traces as interactive
p3 <- ggplot(data = PBMC_F.1_cox_p3, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p3)


##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 598

t2a <- 0
t2b <- 184

PBMC_F.1_cox_p3_excl <- c()
PBMC_F.1_cox_p3_it_diff <- c("D6", "A1", "H1", "H12")

PBMC_F.1_cox_p3_it1 <- PBMC_F.1_cox_p3[PBMC_F.1_cox_p3$Time>=t1a & PBMC_F.1_cox_p3$Time<=t1b,] 
PBMC_F.1_cox_p3_it2 <- PBMC_F.1_cox_p3[PBMC_F.1_cox_p3$Time>=t2a & PBMC_F.1_cox_p3$Time<=t2b,] 

#time1-- 
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_p3_models <- PBMC_F.1_cox_p3_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_p3_models <- PBMC_F.1_cox_p3_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_p3_models2 <- PBMC_F.1_cox_p3_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_p3_models2 <- PBMC_F.1_cox_p3_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times--
PBMC_F.1_cox_p3_it1_mod <- select(PBMC_F.1_cox_p3_models, Well, slope)
PBMC_F.1_cox_p3_it2_mod <- select(PBMC_F.1_cox_p3_models2, Well, slope)

PBMC_F.1_cox_p3_allmod <- full_join(PBMC_F.1_cox_p3_it1_mod, PBMC_F.1_cox_p3_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_F.1_cox_p3_excl ~ NA,
    Well%in%PBMC_F.1_cox_p3_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_F.1_cox_p3_allmod$final_it_act <- -1*1000*60*(PBMC_F.1_cox_p3_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_F.1_cox_p3_allmod$Well)
PBMC_F.1_cox_p3_allmod <- PBMC_F.1_cox_p3_allmod[well_order, ]

#extract activities per well
PBMC_F.1_cox_p3_act <- select(PBMC_F.1_cox_p3_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(p3, t1a, t1b, t2a, t2b,
   PBMC_F.1_cox_p3_excl, PBMC_F.1_cox_p3_it_diff, 
   PBMC_F.1_cox_p3_it1, PBMC_F.1_cox_p3_it2, 
   PBMC_F.1_cox_p3_models, PBMC_F.1_cox_p3_models2,
   PBMC_F.1_cox_p3_it1_mod, PBMC_F.1_cox_p3_it2_mod,
   PBMC_F.1_cox_p3_allmod
)











# ################## NEGATIVE ALL ########################

#changing spectramax data into long format for 1 plate (n1 in this case)
PBMC_F.1_cox_n1 <- PBMC_F.1_cox[91:117, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_F.1_cox_n1 <- PBMC_F.1_cox_n1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_F.1_cox_n1$Column <- as.numeric(PBMC_F.1_cox_n1$Column)
PBMC_F.1_cox_n1$mOD <- as.numeric(PBMC_F.1_cox_n1$mOD)

PBMC_F.1_cox_n1$Time <- period_to_seconds(hms(PBMC_F.1_cox_n1$Time))

PBMC_F.1_cox_n1 <- arrange(PBMC_F.1_cox_n1, Row)
PBMC_F.1_cox_n1 <- arrange(PBMC_F.1_cox_n1, Column)


#plot all traces as interactive
n1 <- ggplot(data = PBMC_F.1_cox_n1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n1)


##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

#PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 598

t2a <- 6
t2b <- 161

PBMC_F.1_cox_n1_excl <- c()
PBMC_F.1_cox_n1_it_diff <- c()

PBMC_F.1_cox_n1_it1 <- PBMC_F.1_cox_n1[PBMC_F.1_cox_n1$Time>=t1a & PBMC_F.1_cox_n1$Time<=t1b,] 
PBMC_F.1_cox_n1_it2 <- PBMC_F.1_cox_n1[PBMC_F.1_cox_n1$Time>=t2a & PBMC_F.1_cox_n1$Time<=t2b,] 

#time1-- 
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_n1_models <- PBMC_F.1_cox_n1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_n1_models <- PBMC_F.1_cox_n1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2--
# Perform linear regression for each group using dplyr and purrr
PBMC_F.1_cox_n1_models2 <- PBMC_F.1_cox_n1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
PBMC_F.1_cox_n1_models2 <- PBMC_F.1_cox_n1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times--
PBMC_F.1_cox_n1_it1_mod <- select(PBMC_F.1_cox_n1_models, Well, slope)
PBMC_F.1_cox_n1_it2_mod <- select(PBMC_F.1_cox_n1_models2, Well, slope)

PBMC_F.1_cox_n1_allmod <- full_join(PBMC_F.1_cox_n1_it1_mod, PBMC_F.1_cox_n1_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_F.1_cox_n1_excl ~ NA,
    Well%in%PBMC_F.1_cox_n1_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_F.1_cox_n1_allmod$final_it_act <- -1*1000*60*(PBMC_F.1_cox_n1_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_F.1_cox_n1_allmod$Well)
PBMC_F.1_cox_n1_allmod <- PBMC_F.1_cox_n1_allmod[well_order, ]

#extract activities per well
PBMC_F.1_cox_n1_act <- select(PBMC_F.1_cox_n1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(n1, t1a, t1b, t2a, t2b,
   PBMC_F.1_cox_n1_excl, PBMC_F.1_cox_n1_it_diff, 
   PBMC_F.1_cox_n1_it1, PBMC_F.1_cox_n1_it2, 
   PBMC_F.1_cox_n1_models, PBMC_F.1_cox_n1_models2,
   PBMC_F.1_cox_n1_it1_mod, PBMC_F.1_cox_n1_it2_mod,
   PBMC_F.1_cox_n1_allmod
)

#######################################################



# join data together and calculate activities ----

cutoff_cv <- 0.35


molar_ext_coef <- 29.5
hom_buff_per_millcell <- 100
vol_sample_assay <-  20
vol_assay_mix <- 200
vol_well_total <- vol_sample_assay + vol_assay_mix


# nonspec
nonspec_raw_mean = mean(PBMC_F.1_cox_n1_act$Activity,na.rm = TRUE)
nonspec_raw_sd = sd(PBMC_F.1_cox_n1_act$Activity,na.rm = TRUE)
nonspec_raw_cv = nonspec_raw_sd / nonspec_raw_mean

if(abs(nonspec_raw_cv) > cutoff_cv){ 
  PBMC_F.1_cox_n1_act <- PBMC_F.1_cox_n1_act %>%
    mutate(activity_cleaned = ifelse(abs(Activity - nonspec_raw_mean) > nonspec_raw_sd, NA, Activity))
} else {
  PBMC_F.1_cox_n1_act <- PBMC_F.1_cox_n1_act %>%
    mutate(activity_cleaned = Activity)
}

nonspec_cleaned_mean = mean(PBMC_F.1_cox_n1_act$activity_cleaned,na.rm = TRUE)
nonspec_cleaned_sd = sd(PBMC_F.1_cox_n1_act$activity_cleaned,na.rm = TRUE)
nonspec_cleaned_cv = nonspec_cleaned_sd / nonspec_cleaned_mean


#pos & specific
PBMC_F.1_cox_all <- full_join(PBMC_F.1_cox_p1_act, PBMC_F.1_cox_p2_act, by = "Well") %>%
  rename(p1 = Activity.x, p2 = Activity.y) %>%
  full_join(., PBMC_F.1_cox_p3_act, by= "Well") %>%
  rename(p3 = Activity) %>%
  full_join(., PBMC_F.1_cox_n1_act, by= "Well") %>%
  rename(n1 = Activity)%>%
  column_to_rownames(var="Well") %>%
  mutate(total_raw_mean = rowMeans(select(., starts_with("p")), na.rm = TRUE)) %>%
  mutate(total_raw_sd = rowSds(as.matrix(select(., starts_with("p"))), na.rm = TRUE)) %>%
  mutate(total_raw_cv = total_raw_sd/total_raw_mean) %>%
  mutate(spec_raw = total_raw_mean - nonspec_raw_mean) %>%
  mutate(p1_cleaned = ifelse(total_raw_cv > cutoff_cv & abs(p1-total_raw_mean) > abs(p2-total_raw_mean) & abs(p1-total_raw_mean) > abs(p3-total_raw_mean) & !is.na(p2) & !is.na(p3),NA, p1)) %>%
  mutate(p2_cleaned = ifelse(total_raw_cv > cutoff_cv & abs(p2-total_raw_mean) > abs(p1-total_raw_mean) & abs(p2-total_raw_mean) > abs(p3-total_raw_mean) & !is.na(p1) & !is.na(p3),NA, p2)) %>%
  mutate(p3_cleaned = ifelse(total_raw_cv > cutoff_cv & abs(p3-total_raw_mean) > abs(p1-total_raw_mean) & abs(p3-total_raw_mean) > abs(p2-total_raw_mean) & !is.na(p1) & !is.na(p2),NA, p3)) %>%
  mutate(total_cleaned_mean = rowMeans(select(., starts_with("p") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(total_cleaned_sd = rowSds(as.matrix(select(., starts_with("p") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(total_cleaned_cv = total_cleaned_sd/total_cleaned_mean) %>%
  mutate(PBMC_F.1_cox_spec_cleaned = ifelse(total_cleaned_mean - nonspec_cleaned_mean>=0, total_cleaned_mean - nonspec_cleaned_mean, 0)) %>%
  mutate(PBMC_F.1_cox_spec_cleaned_cell = (((PBMC_F.1_cox_spec_cleaned/1000)*hom_buff_per_millcell)/(0.552*molar_ext_coef*vol_sample_assay*(vol_sample_assay/vol_well_total)))*1000)%>%
  rownames_to_column(var = "Well")

# write_xlsx(PBMC_F.1_cox_all, "interim_data/PBMC_F.1_cox_all.xlsx")
pre_av_tot_cv <- mean(abs(PBMC_F.1_cox_all$total_raw_cv), na.rm = T)
pre_av_non_cv <- mean(abs(PBMC_F.1_cox_all$nonspec_raw_cv), na.rm = T)
post_av_tot_cv <- mean(abs(PBMC_F.1_cox_all$total_cleaned_cv), na.rm = T)
post_av_non_cv <- mean(abs(PBMC_F.1_cox_all$nonspec_cleaned_cv), na.rm = T)
excl <- as.numeric(sum(is.na(PBMC_F.1_cox_all$p1_cleaned)) + sum(is.na(PBMC_F.1_cox_all$p2_cleaned)) + sum(is.na(PBMC_F.1_cox_all$p3_cleaned)) +
                     sum(is.na(PBMC_F.1_cox_all$activity_cleaned)))

32/(4*96)

PBMC_F.1_cox_byid <- PBMC_F.1_cox_all %>%
  select(Well, PBMC_F.1_cox_spec_cleaned, PBMC_F.1_cox_spec_cleaned_cell) %>%
  full_join(.,PBMC_F.1, by="Well") %>%
  select(Well, sampleName, PBMC_F.1_cox_spec_cleaned, PBMC_F.1_cox_spec_cleaned_cell)
write_xlsx(PBMC_F.1_cox_byid, "interim_data/manual/PBMC_F.1_cox.xlsx")


rm(PBMC_F.1_cox_p1, PBMC_F.1_cox_p2, PBMC_F.1_cox_p3, PBMC_F.1_cox_n1, PBMC_F.1_cox_n2, PBMC_F.1_cox_n3,
   PBMC_F.1_cox_p1_act, PBMC_F.1_cox_p2_act, PBMC_F.1_cox_p3_act, PBMC_F.1_cox_n1_act, PBMC_F.1_cox_n2_act, PBMC_F.1_cox_n3_act
)    
# 
# #visualize distributions -----
# 
# 
# 
# pos_mean <- PBMC_F.1_cox_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_mean), color = "green")+
#   geom_density(aes(x=total_cleaned_mean), color = "blue")
# neg_mean <- PBMC_F.1_cox_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_mean), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_mean), color = "orange")
# 
# 
# pos_cv <- PBMC_F.1_cox_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_cv), color = "green")+
#   geom_density(aes(x=total_cleaned_cv), color = "blue")
# print(pos_cv)
# neg_cv <- PBMC_F.1_cox_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_cv), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_cv), color = "orange")
# print(neg_cv)
# 
# 
# fd <- ggplot(PBMC_F.1_cox_byid, mapping = aes(x = PBMC_F.1_cox_spec_cleaned)) +
#   geom_histogram(alpha = 0.3, fill = "dodgerblue", show.legend = F, bins = 50) +
#   geom_density(colour = "blue")
# print (fd)
# 
# PBMC_F.1_cox_byid %>%
#   ggplot()+
#   geom_density(aes(x=PBMC_F.1_cox_spec_cleaned))
# 
# plotly::ggplotly(fd)
# 
# rm(PBMC_F.1_cox, PBMC_F.1_cox_notemp, PBMC_F.1_cox_all, PBMC_F.1_ci_byid)


