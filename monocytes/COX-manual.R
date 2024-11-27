library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)



#read file, get rid of temperature column ----

mono_cox <- read.delim ("raw_data/spectramax_data/monocytes/COX-monocytes.xls", header = TRUE, sep = "\t", fileEncoding = "UTF-16LE", skip = 2)[,-2]
mono <- readxl::read_xlsx("plate_maps/MiSBIE_mono.xlsx", sheet = 1)



################### POSITIVE 1 ########################

#changing spectramax data into long format for 1 plate (p1 in this case)
mono_cox_p1 <- mono_cox[1:53, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
mono_cox_p1 <- mono_cox_p1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
mono_cox_p1$Column <- as.numeric(mono_cox_p1$Column)
mono_cox_p1$mOD <- as.numeric(mono_cox_p1$mOD)

mono_cox_p1$Time <- period_to_seconds(hms(mono_cox_p1$Time))

mono_cox_p1 <- arrange(mono_cox_p1, Row)
mono_cox_p1 <- arrange(mono_cox_p1, Column)


#plot all traces as interactive
p1 <- ggplot(data = mono_cox_p1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row)) +
  geom_point(aes(color = Row))

plotly::ggplotly(p1)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 1196

t2a <- 0
t2b <- 368

mono_cox_p1_excl <- c()
mono_cox_p1_it_diff <- c("A1","A12","D6","E6","H1","H12","H7","H8")

mono_cox_p1_it1 <- mono_cox_p1[mono_cox_p1$Time>=t1a & mono_cox_p1$Time<=t1b,] 
mono_cox_p1_it2 <- mono_cox_p1[mono_cox_p1$Time>=t2a & mono_cox_p1$Time<=t2b,] 

#time1--
# Perform linear regression for each group using dplyr and purrr
mono_cox_p1_models <- mono_cox_p1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
mono_cox_p1_models <- mono_cox_p1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
mono_cox_p1_models2 <- mono_cox_p1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
mono_cox_p1_models2 <- mono_cox_p1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times- 
mono_cox_p1_it1_mod <- select(mono_cox_p1_models, Well, slope)
mono_cox_p1_it2_mod <- select(mono_cox_p1_models2, Well, slope)


mono_cox_p1_allmod <- full_join(mono_cox_p1_it1_mod, mono_cox_p1_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% mono_cox_p1_excl ~ NA,
    Well%in%mono_cox_p1_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
mono_cox_p1_allmod$final_it_act <- -1*1000*60*(mono_cox_p1_allmod$final_it_act)

#order by well
well_order <- mixedorder(mono_cox_p1_allmod$Well)
mono_cox_p1_allmod <- mono_cox_p1_allmod[well_order, ]

#extract activities per well
mono_cox_p1_act <- select(mono_cox_p1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)

rm(p1, t1a, t1b, t2a, t2b,
   mono_cox_p1_excl, mono_cox_p1_it_diff, 
   mono_cox_p1_it1, mono_cox_p1_it2, 
   mono_cox_p1_models, mono_cox_p1_models2,
   mono_cox_p1_it1_mod, mono_cox_p1_it2_mod,
   mono_cox_p1_allmod
)




# ################## POSITIVE 2 ########################

#changing spectramax data into long format for 1 plate (p2 in this case)
mono_cox_p2 <- mono_cox[57:109, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
mono_cox_p2 <- mono_cox_p2[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
mono_cox_p2$Column <- as.numeric(mono_cox_p2$Column)
mono_cox_p2$mOD <- as.numeric(mono_cox_p2$mOD)

mono_cox_p2$Time <- period_to_seconds(hms(mono_cox_p2$Time))

mono_cox_p2 <- arrange(mono_cox_p2, Row)
mono_cox_p2 <- arrange(mono_cox_p2, Column)


#plot all traces as interactive
p2 <- ggplot(data = mono_cox_p2, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p2)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 1196

t2a <- 0
t2b <- 368

mono_cox_p2_excl <- c()
mono_cox_p2_it_diff <- c("A1","A12","D6","E6","H1","H12","H7","H8")


mono_cox_p2_it1 <- mono_cox_p2[mono_cox_p2$Time>=t1a & mono_cox_p2$Time<=t1b,] 
mono_cox_p2_it2 <- mono_cox_p2[mono_cox_p2$Time>=t2a & mono_cox_p2$Time<=t2b,] 

#time1--
# Perform linear regression for each group using dplyr and purrr
mono_cox_p2_models <- mono_cox_p2_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
mono_cox_p2_models <- mono_cox_p2_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
mono_cox_p2_models2 <- mono_cox_p2_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
mono_cox_p2_models2 <- mono_cox_p2_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times- 
mono_cox_p2_it1_mod <- select(mono_cox_p2_models, Well, slope)
mono_cox_p2_it2_mod <- select(mono_cox_p2_models2, Well, slope)

mono_cox_p2_allmod <- full_join(mono_cox_p2_it1_mod, mono_cox_p2_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% mono_cox_p2_excl ~ NA,
    Well%in%mono_cox_p2_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
mono_cox_p2_allmod$final_it_act <- -1*1000*60*(mono_cox_p2_allmod$final_it_act)

#order by well
well_order <- mixedorder(mono_cox_p2_allmod$Well)
mono_cox_p2_allmod <- mono_cox_p2_allmod[well_order, ]

#extract activities per well
mono_cox_p2_act <- select(mono_cox_p2_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(p2, t1a, t1b, t2a, t2b,
   mono_cox_p2_excl, mono_cox_p2_it_diff, 
   mono_cox_p2_it1, mono_cox_p2_it2, 
   mono_cox_p2_models, mono_cox_p2_models2,
   mono_cox_p2_it1_mod, mono_cox_p2_it2_mod,
   mono_cox_p2_allmod
)



# ################## POSITIVE 3 ########################

#changing spectramax data into long format for 1 plate (p3 in this case)
mono_cox_p3 <- mono_cox[113:165, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
mono_cox_p3 <- mono_cox_p3[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
mono_cox_p3$Column <- as.numeric(mono_cox_p3$Column)
mono_cox_p3$mOD <- as.numeric(mono_cox_p3$mOD)

mono_cox_p3$Time <- period_to_seconds(hms(mono_cox_p3$Time))

mono_cox_p3 <- arrange(mono_cox_p3, Row)
mono_cox_p3 <- arrange(mono_cox_p3, Column)


#plot all traces as interactive
p3 <- ggplot(data = mono_cox_p3, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p3)


##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 1196

t2a <- 0
t2b <- 368

t3a <- 1
t3b <- 1

mono_cox_p3_excl <- c()
mono_cox_p3_it_diff <- c("A1","A12","D6","E6","H1","H12","H7","H8")
mono_cox_p3_it_diff2 <- c()

mono_cox_p3_it1 <- mono_cox_p3[mono_cox_p3$Time>=t1a & mono_cox_p3$Time<=t1b,] 
mono_cox_p3_it2 <- mono_cox_p3[mono_cox_p3$Time>=t2a & mono_cox_p3$Time<=t2b,]
mono_cox_p3_it3 <- mono_cox_p3[mono_cox_p3$Time>=t3a & mono_cox_p3$Time<=t3b,] 


#time1-- 
# Perform linear regression for each group using dplyr and purrr
mono_cox_p3_models <- mono_cox_p3_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
mono_cox_p3_models <- mono_cox_p3_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
mono_cox_p3_models2 <- mono_cox_p3_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
mono_cox_p3_models2 <- mono_cox_p3_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time3-
# Perform linear regression for each group using dplyr and purrr
mono_cox_p3_models3 <- mono_cox_p3_it3 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
mono_cox_p3_models3 <- mono_cox_p3_models3 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times------------------
mono_cox_p3_it1_mod <- select(mono_cox_p3_models, Well, slope)
mono_cox_p3_it2_mod <- select(mono_cox_p3_models2, Well, slope)
mono_cox_p3_it3_mod <- select(mono_cox_p3_models3, Well, slope)


mono_cox_p3_allmod <- full_join(mono_cox_p3_it1_mod, mono_cox_p3_it2_mod, by = "Well")%>%
  rename(act_it1 = slope.x, act_it2 = slope.y) %>%
  full_join(., mono_cox_p3_it3_mod, by = "Well") %>%
  rename(act_it3 = slope) %>%
  mutate(final_it_act = case_when(
    Well %in% mono_cox_p3_excl ~ NA,
    Well%in%mono_cox_p3_it_diff ~ act_it2,
    Well%in%mono_cox_p3_it_diff2 ~ act_it3,
    TRUE ~ act_it1
  ))


# change slope into correct units
mono_cox_p3_allmod$final_it_act <- -1*1000*60*(mono_cox_p3_allmod$final_it_act)

#order by well
well_order <- mixedorder(mono_cox_p3_allmod$Well)
mono_cox_p3_allmod <- mono_cox_p3_allmod[well_order, ]

#extract activities per well
mono_cox_p3_act <- select(mono_cox_p3_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(t1a, t1b, t2a, t2b, t3a, t3b,
   mono_cox_p3_excl, mono_cox_p3_it_diff, mono_cox_p3_it_diff2,
   mono_cox_p3_it1, mono_cox_p3_it2, mono_cox_p3_it3,
   mono_cox_p3_models, mono_cox_p3_models2, mono_cox_p3_models3,
   mono_cox_p3_it1_mod, mono_cox_p3_it2_mod, mono_cox_p3_it3_mod,
   mono_cox_p3_allmod
)











# ################## NEGATIVE ALL ########################

#changing spectramax data into long format for 1 plate (n1 in this case)
mono_cox_n1 <- mono_cox[169:221, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
mono_cox_n1 <- mono_cox_n1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
mono_cox_n1$Column <- as.numeric(mono_cox_n1$Column)
mono_cox_n1$mOD <- as.numeric(mono_cox_n1$mOD)

mono_cox_n1$Time <- period_to_seconds(hms(mono_cox_n1$Time))

mono_cox_n1 <- arrange(mono_cox_n1, Row)
mono_cox_n1 <- arrange(mono_cox_n1, Column)


#plot all traces as interactive
n1 <- ggplot(data = mono_cox_n1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n1)


##REDUCTION##------------------

# to get reduction: per well, make linear regression between specoxfied time points and then get slope

#subsetting data for correct int times

#PICK INTEGRATION TIMES------------------
t1a <- 184 
t1b <- 1173

t2a <- 6
t2b <- 161

mono_cox_n1_excl <- c("C6","C4","C5","D4","B4")
mono_cox_n1_it_diff <- c()

mono_cox_n1_it1 <- mono_cox_n1[mono_cox_n1$Time>=t1a & mono_cox_n1$Time<=t1b,] 
mono_cox_n1_it2 <- mono_cox_n1[mono_cox_n1$Time>=t2a & mono_cox_n1$Time<=t2b,] 

#time1-- 
# Perform linear regression for each group using dplyr and purrr
mono_cox_n1_models <- mono_cox_n1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
mono_cox_n1_models <- mono_cox_n1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2--
# Perform linear regression for each group using dplyr and purrr
mono_cox_n1_models2 <- mono_cox_n1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficoxents for each group
mono_cox_n1_models2 <- mono_cox_n1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times--
mono_cox_n1_it1_mod <- select(mono_cox_n1_models, Well, slope)
mono_cox_n1_it2_mod <- select(mono_cox_n1_models2, Well, slope)

mono_cox_n1_allmod <- full_join(mono_cox_n1_it1_mod, mono_cox_n1_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% mono_cox_n1_excl ~ NA,
    Well%in%mono_cox_n1_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
mono_cox_n1_allmod$final_it_act <- -1*1000*60*(mono_cox_n1_allmod$final_it_act)

#order by well
well_order <- mixedorder(mono_cox_n1_allmod$Well)
mono_cox_n1_allmod <- mono_cox_n1_allmod[well_order, ]

#extract activities per well
mono_cox_n1_act <- select(mono_cox_n1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(n1, t1a, t1b, t2a, t2b,
   mono_cox_n1_excl, mono_cox_n1_it_diff, 
   mono_cox_n1_it1, mono_cox_n1_it2, 
   mono_cox_n1_models, mono_cox_n1_models2,
   mono_cox_n1_it1_mod, mono_cox_n1_it2_mod,
   mono_cox_n1_allmod
)

#######################################################



# join data together and calculate activities ----

cutoff_cv <- 0.3


molar_ext_coef <- 29.5
hom_buff_per_millcell <- 100
vol_sample_assay <-  20
vol_assay_mix <- 200
vol_well_total <- vol_sample_assay + vol_assay_mix


# nonspec
nonspec_raw_mean = mean(mono_cox_n1_act$Activity,na.rm = TRUE)
nonspec_raw_sd = sd(mono_cox_n1_act$Activity,na.rm = TRUE)
nonspec_raw_cv = nonspec_raw_sd / nonspec_raw_mean

if(abs(nonspec_raw_cv) > cutoff_cv){ 
  mono_cox_n1_act <- mono_cox_n1_act %>%
    mutate(n1_cleaned = ifelse(abs(Activity - nonspec_raw_mean) > nonspec_raw_sd, NA, Activity))
} else {
  mono_cox_n1_act <- mono_cox_n1_act %>%
    mutate(n1_cleaned = Activity)
}

nonspec_cleaned_mean = mean(mono_cox_n1_act$n1_cleaned,na.rm = TRUE)
nonspec_cleaned_sd = sd(mono_cox_n1_act$n1_cleaned,na.rm = TRUE)
nonspec_cleaned_cv = nonspec_cleaned_sd / nonspec_cleaned_mean


#pos & specific
mono_cox_all <- full_join(mono_cox_p1_act, mono_cox_p2_act, by = "Well") %>%
  rename(p1 = Activity.x, p2 = Activity.y) %>%
  full_join(., mono_cox_p3_act, by= "Well") %>%
  rename(p3 = Activity) %>%
  full_join(., mono_cox_n1_act, by= "Well") %>%
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
  mutate(mono_cox_spec_cleaned = ifelse(total_cleaned_mean - nonspec_cleaned_mean>=0, total_cleaned_mean - nonspec_cleaned_mean, 0)) %>%
  mutate(mono_cox_spec_cleaned_cell = (((mono_cox_spec_cleaned/1000)*hom_buff_per_millcell)/(0.552*molar_ext_coef*vol_sample_assay*(vol_sample_assay/vol_well_total)))*1000)%>%
  rownames_to_column(var = "Well")%>%
  filter(!Well %in% c("B12","C12","D12","E12","F12", "G12",
                      "A10","B10","C10","D10","E10","F10","G10","H12",
                      "A11","B11","C11","D11","E11","F11","G11","H11",
                      "D9","E9","F9","G9","H9",
                      "D1","E1","H2","A3"))

# write_xlsx(mono_cox_all, "interim_data/mono_cox_all.xlsx")
pre_av_tot_cv <- mean(abs(mono_cox_all$total_raw_cv), na.rm = T)
pre_av_non_cv <- mean(abs(mono_cox_all$nonspec_raw_cv), na.rm = T)
post_av_tot_cv <- mean(abs(mono_cox_all$total_cleaned_cv), na.rm = T)
post_av_non_cv <- mean(abs(mono_cox_all$nonspec_cleaned_cv), na.rm = T)
excl <- as.numeric(sum(is.na(mono_cox_all$p1_cleaned)) + sum(is.na(mono_cox_all$p2_cleaned)) + sum(is.na(mono_cox_all$p3_cleaned)) +
                     sum(is.na(mono_cox_all$n1_cleaned)) + sum(is.na(mono_cox_all$n2_cleaned)) + sum(is.na(mono_cox_all$n3_cleaned)))

16/(4*65)

mono_cox_byid <- mono_cox_all %>%
  select(Well, mono_cox_spec_cleaned, mono_cox_spec_cleaned_cell) %>%
  full_join(.,mono, by="Well") %>%
  select(Well, sampleName, mono_cox_spec_cleaned, mono_cox_spec_cleaned_cell)
write_xlsx(mono_cox_byid, "interim_data/manual/mono_cox.xlsx")


rm(mono_cox_p1, mono_cox_p2, mono_cox_p3, mono_cox_n1, mono_cox_n2, mono_cox_n3,
   mono_cox_p1_act, mono_cox_p2_act, mono_cox_p3_act, mono_cox_n1_act, mono_cox_n2_act, mono_cox_n3_act
)    
# 
# #visualize distributions -----
# 
# 
# 
# pos_mean <- mono_cox_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_mean), color = "green")+
#   geom_density(aes(x=total_cleaned_mean), color = "blue")
# neg_mean <- mono_cox_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_mean), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_mean), color = "orange")
# 
# 
# pos_cv <- mono_cox_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_cv), color = "green")+
#   geom_density(aes(x=total_cleaned_cv), color = "blue")
# print(pos_cv)
# neg_cv <- mono_cox_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_cv), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_cv), color = "orange")
# print(neg_cv)
# 
# 
# fd <- ggplot(mono_cox_byid, mapping = aes(x = mono_cox_spec_cleaned)) +
#   geom_histogram(alpha = 0.3, fill = "dodgerblue", show.legend = F, bins = 50) +
#   geom_density(colour = "blue")
# print (fd)
# 
# mono_cox_byid %>%
#   ggplot()+
#   geom_density(aes(x=mono_cox_spec_cleaned))
# 
# plotly::ggplotly(fd)
# 
# rm(mono_cox, mono_cox_notemp, mono_cox_all, mono_ci_byid)


