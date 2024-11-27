library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)


#read file, get rid of temperature column ----

PBMC_S.1_ci <- read.delim ("raw_data/spectramax_data/PBMC_S.1/C1-PBMC_S.1.xls", header = TRUE, sep = "\t", fileEncoding = "UTF-16LE", skip = 2)[,-2]
PBMC_S.1 <- readxl::read_xlsx("plate_maps/MiSBIE_PBMC_S.1.xlsx", sheet = 1)

################### POSITIVE 1 ########################

#changing spectramax data into long format for 1 plate (p1 in this case)
PBMC_S.1_ci_p1 <- PBMC_S.1_ci[0:27, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_S.1_ci_p1 <- PBMC_S.1_ci_p1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_S.1_ci_p1$Column <- as.numeric(PBMC_S.1_ci_p1$Column)
PBMC_S.1_ci_p1$mOD <- as.numeric(PBMC_S.1_ci_p1$mOD)

PBMC_S.1_ci_p1$Time <- period_to_seconds(hms(PBMC_S.1_ci_p1$Time))

PBMC_S.1_ci_p1 <- arrange(PBMC_S.1_ci_p1, Row)
PBMC_S.1_ci_p1 <- arrange(PBMC_S.1_ci_p1, Column)


#plot all traces as interactive
p1 <- ggplot(data = PBMC_S.1_ci_p1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row)) +
  geom_point(aes(color = Row))

plotly::ggplotly(p1)

##REDUCTION##

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES----
t1a <- 161  
t1b <- 598

t2a <- 391
t2b <- 575

PBMC_S.1_ci_p1_excl <- c("B12","C12","D12","E12","F12","G12")
PBMC_S.1_ci_p1_it_diff <- c()

PBMC_S.1_ci_p1_it1 <- PBMC_S.1_ci_p1[PBMC_S.1_ci_p1$Time>=t1a & PBMC_S.1_ci_p1$Time<=t1b,] 
PBMC_S.1_ci_p1_it2 <- PBMC_S.1_ci_p1[PBMC_S.1_ci_p1$Time>=t2a & PBMC_S.1_ci_p1$Time<=t2b,] 

#time1
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_p1_models <- PBMC_S.1_ci_p1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_p1_models <- PBMC_S.1_ci_p1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_p1_models2 <- PBMC_S.1_ci_p1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_p1_models2 <- PBMC_S.1_ci_p1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times
PBMC_S.1_ci_p1_it1_mod <- select(PBMC_S.1_ci_p1_models, Well, slope)
PBMC_S.1_ci_p1_it2_mod <- select(PBMC_S.1_ci_p1_models2, Well, slope)


PBMC_S.1_ci_p1_allmod <- full_join(PBMC_S.1_ci_p1_it1_mod, PBMC_S.1_ci_p1_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_S.1_ci_p1_excl ~ NA,
    Well%in%PBMC_S.1_ci_p1_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_S.1_ci_p1_allmod$final_it_act <- -1*1000*60*(PBMC_S.1_ci_p1_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_S.1_ci_p1_allmod$Well)
PBMC_S.1_ci_p1_allmod <- PBMC_S.1_ci_p1_allmod[well_order, ]

#extract activities per well
PBMC_S.1_ci_p1_act <- select(PBMC_S.1_ci_p1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(p1, t1a, t1b, t2a, t2b,
   PBMC_S.1_ci_p1_excl, PBMC_S.1_ci_p1_it_diff, 
   PBMC_S.1_ci_p1_it1, PBMC_S.1_ci_p1_it2, 
   PBMC_S.1_ci_p1_models, PBMC_S.1_ci_p1_models2,
   PBMC_S.1_ci_p1_it1_mod, PBMC_S.1_ci_p1_it2_mod,
   PBMC_S.1_ci_p1_allmod
)

#######################################################

################### POSITIVE 2 ########################

#changing spectramax data into long format for 1 plate (p2 in this case)
PBMC_S.1_ci_p2 <- PBMC_S.1_ci[31:57, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_S.1_ci_p2 <- PBMC_S.1_ci_p2[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_S.1_ci_p2$Column <- as.numeric(PBMC_S.1_ci_p2$Column)
PBMC_S.1_ci_p2$mOD <- as.numeric(PBMC_S.1_ci_p2$mOD)

PBMC_S.1_ci_p2$Time <- period_to_seconds(hms(PBMC_S.1_ci_p2$Time))

PBMC_S.1_ci_p2 <- arrange(PBMC_S.1_ci_p2, Row)
PBMC_S.1_ci_p2 <- arrange(PBMC_S.1_ci_p2, Column)


#plot all traces as interactive
p2 <- ggplot(data = PBMC_S.1_ci_p2, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p2)


##REDUCTION##

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES ----
t1a <- 161
t1b <- 598

t2a <- 414
t2b <- 575

PBMC_S.1_ci_p2_excl <-  c("B12","C12","D12","E12","F12","G12")
PBMC_S.1_ci_p2_it_diff <- c()

PBMC_S.1_ci_p2_it1 <- PBMC_S.1_ci_p2[PBMC_S.1_ci_p2$Time>=t1a & PBMC_S.1_ci_p2$Time<=t1b,] 
PBMC_S.1_ci_p2_it2 <- PBMC_S.1_ci_p2[PBMC_S.1_ci_p2$Time>=t2a & PBMC_S.1_ci_p2$Time<=t2b,] 

#time1
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_p2_models <- PBMC_S.1_ci_p2_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_p2_models <- PBMC_S.1_ci_p2_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_p2_models2 <- PBMC_S.1_ci_p2_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_p2_models2 <- PBMC_S.1_ci_p2_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times
PBMC_S.1_ci_p2_it1_mod <- select(PBMC_S.1_ci_p2_models, Well, slope)
PBMC_S.1_ci_p2_it2_mod <- select(PBMC_S.1_ci_p2_models2, Well, slope)

PBMC_S.1_ci_p2_allmod <- full_join(PBMC_S.1_ci_p2_it1_mod, PBMC_S.1_ci_p2_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_S.1_ci_p2_excl ~ NA,
    Well%in%PBMC_S.1_ci_p2_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_S.1_ci_p2_allmod$final_it_act <- -1*1000*60*(PBMC_S.1_ci_p2_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_S.1_ci_p2_allmod$Well)
PBMC_S.1_ci_p2_allmod <- PBMC_S.1_ci_p2_allmod[well_order, ]

#extract activities per well
PBMC_S.1_ci_p2_act <- select(PBMC_S.1_ci_p2_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(p2, t1a, t1b, t2a, t2b,
   PBMC_S.1_ci_p2_excl, PBMC_S.1_ci_p2_it_diff, 
   PBMC_S.1_ci_p2_it1, PBMC_S.1_ci_p2_it2, 
   PBMC_S.1_ci_p2_models, PBMC_S.1_ci_p2_models2,
   PBMC_S.1_ci_p2_it1_mod, PBMC_S.1_ci_p2_it2_mod,
   PBMC_S.1_ci_p2_allmod
)

#######################################################

################### POSITIVE 3 ########################

#changing spectramax data into long format for 1 plate (p3 in this case)
PBMC_S.1_ci_p3 <- PBMC_S.1_ci[61:87, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_S.1_ci_p3 <- PBMC_S.1_ci_p3[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_S.1_ci_p3$Column <- as.numeric(PBMC_S.1_ci_p3$Column)
PBMC_S.1_ci_p3$mOD <- as.numeric(PBMC_S.1_ci_p3$mOD)

PBMC_S.1_ci_p3$Time <- period_to_seconds(hms(PBMC_S.1_ci_p3$Time))

PBMC_S.1_ci_p3 <- arrange(PBMC_S.1_ci_p3, Row)
PBMC_S.1_ci_p3 <- arrange(PBMC_S.1_ci_p3, Column)


#plot all traces as interactive
p3 <- ggplot(data = PBMC_S.1_ci_p3, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p3)

##REDUCTION##-

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES----
t1a <- 299  
t1b <- 598

t2a <- 391
t2b <- 575

PBMC_S.1_ci_p3_excl <-  c("B12","C12","D12","E12","F12","G12","B11")
PBMC_S.1_ci_p3_it_diff <- c()

PBMC_S.1_ci_p3_it1 <- PBMC_S.1_ci_p3[PBMC_S.1_ci_p3$Time>=t1a & PBMC_S.1_ci_p3$Time<=t1b,] 
PBMC_S.1_ci_p3_it2 <- PBMC_S.1_ci_p3[PBMC_S.1_ci_p3$Time>=t2a & PBMC_S.1_ci_p3$Time<=t2b,] 

#time1
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_p3_models <- PBMC_S.1_ci_p3_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_p3_models <- PBMC_S.1_ci_p3_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_p3_models2 <- PBMC_S.1_ci_p3_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_p3_models2 <- PBMC_S.1_ci_p3_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times
PBMC_S.1_ci_p3_it1_mod <- select(PBMC_S.1_ci_p3_models, Well, slope)
PBMC_S.1_ci_p3_it2_mod <- select(PBMC_S.1_ci_p3_models2, Well, slope)

PBMC_S.1_ci_p3_allmod <- full_join(PBMC_S.1_ci_p3_it1_mod, PBMC_S.1_ci_p3_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_S.1_ci_p3_excl ~ NA,
    Well%in%PBMC_S.1_ci_p3_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_S.1_ci_p3_allmod$final_it_act <- -1*1000*60*(PBMC_S.1_ci_p3_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_S.1_ci_p3_allmod$Well)
PBMC_S.1_ci_p3_allmod <- PBMC_S.1_ci_p3_allmod[well_order, ]

#extract activities per well
PBMC_S.1_ci_p3_act <- select(PBMC_S.1_ci_p3_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(p3, t1a, t1b, t2a, t2b,
   PBMC_S.1_ci_p3_excl, PBMC_S.1_ci_p3_it_diff, 
   PBMC_S.1_ci_p3_it1, PBMC_S.1_ci_p3_it2, 
   PBMC_S.1_ci_p3_models, PBMC_S.1_ci_p3_models2,
   PBMC_S.1_ci_p3_it1_mod, PBMC_S.1_ci_p3_it2_mod,
   PBMC_S.1_ci_p3_allmod
)

#######################################################








################### NEGATIVE 1 ########################

#changing spectramax data into long format for 1 plate (n1 in this case)
PBMC_S.1_ci_n1 <- PBMC_S.1_ci[91:117, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_S.1_ci_n1 <- PBMC_S.1_ci_n1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_S.1_ci_n1$Column <- as.numeric(PBMC_S.1_ci_n1$Column)
PBMC_S.1_ci_n1$mOD <- as.numeric(PBMC_S.1_ci_n1$mOD)

PBMC_S.1_ci_n1$Time <- period_to_seconds(hms(PBMC_S.1_ci_n1$Time))

PBMC_S.1_ci_n1 <- arrange(PBMC_S.1_ci_n1, Row)
PBMC_S.1_ci_n1 <- arrange(PBMC_S.1_ci_n1, Column)


#plot all traces as interactive
n1 <- ggplot(data = PBMC_S.1_ci_n1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n1)

##REDUCTION##

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES----
t1a <- 161
t1b <- 598

t2a <- 322
t2b <- 483

PBMC_S.1_ci_n1_excl <-  c("B12","C12","D12","E12","F12","G12")
PBMC_S.1_ci_n1_it_diff <- c()

PBMC_S.1_ci_n1_it1 <- PBMC_S.1_ci_n1[PBMC_S.1_ci_n1$Time>=t1a & PBMC_S.1_ci_n1$Time<=t1b,] 
PBMC_S.1_ci_n1_it2 <- PBMC_S.1_ci_n1[PBMC_S.1_ci_n1$Time>=t2a & PBMC_S.1_ci_n1$Time<=t2b,] 

#time1
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_n1_models <- PBMC_S.1_ci_n1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_n1_models <- PBMC_S.1_ci_n1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_n1_models2 <- PBMC_S.1_ci_n1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_n1_models2 <- PBMC_S.1_ci_n1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times
PBMC_S.1_ci_n1_it1_mod <- select(PBMC_S.1_ci_n1_models, Well, slope)
PBMC_S.1_ci_n1_it2_mod <- select(PBMC_S.1_ci_n1_models2, Well, slope)

PBMC_S.1_ci_n1_allmod <- full_join(PBMC_S.1_ci_n1_it1_mod, PBMC_S.1_ci_n1_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_S.1_ci_n1_excl ~ NA,
    Well%in%PBMC_S.1_ci_n1_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_S.1_ci_n1_allmod$final_it_act <- -1*1000*60*(PBMC_S.1_ci_n1_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_S.1_ci_n1_allmod$Well)
PBMC_S.1_ci_n1_allmod <- PBMC_S.1_ci_n1_allmod[well_order, ]

#extract activities per well
PBMC_S.1_ci_n1_act <- select(PBMC_S.1_ci_n1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(n1, t1a, t1b, t2a, t2b,
   PBMC_S.1_ci_n1_excl, PBMC_S.1_ci_n1_it_diff, 
   PBMC_S.1_ci_n1_it1, PBMC_S.1_ci_n1_it2, 
   PBMC_S.1_ci_n1_models, PBMC_S.1_ci_n1_models2,
   PBMC_S.1_ci_n1_it1_mod, PBMC_S.1_ci_n1_it2_mod,
   PBMC_S.1_ci_n1_allmod
)

#######################################################

################### NEGATIVE 2 ########################

#changing spectramax data into long format for 1 plate (n2 in this case)
PBMC_S.1_ci_n2 <- PBMC_S.1_ci[121:147, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_S.1_ci_n2 <- PBMC_S.1_ci_n2[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_S.1_ci_n2$Column <- as.numeric(PBMC_S.1_ci_n2$Column)
PBMC_S.1_ci_n2$mOD <- as.numeric(PBMC_S.1_ci_n2$mOD)

PBMC_S.1_ci_n2$Time <- period_to_seconds(hms(PBMC_S.1_ci_n2$Time))

PBMC_S.1_ci_n2 <- arrange(PBMC_S.1_ci_n2, Row)
PBMC_S.1_ci_n2 <- arrange(PBMC_S.1_ci_n2, Column)


#plot all traces as interactive
n2 <- ggplot(data = PBMC_S.1_ci_n2, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n2)

##REDUCTION

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES----
t1a <- 161
t1b <- 598

t2a <- 1
t2b <- 1

PBMC_S.1_ci_n2_excl <-  c("B12","C12","D12","E12","F12","G12")
PBMC_S.1_ci_n2_it_diff <- c()

PBMC_S.1_ci_n2_it1 <- PBMC_S.1_ci_n2[PBMC_S.1_ci_n2$Time>=t1a & PBMC_S.1_ci_n2$Time<=t1b,] 
PBMC_S.1_ci_n2_it2 <- PBMC_S.1_ci_n2[PBMC_S.1_ci_n2$Time>=t2a & PBMC_S.1_ci_n2$Time<=t2b,] 

#time1
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_n2_models <- PBMC_S.1_ci_n2_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_n2_models <- PBMC_S.1_ci_n2_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_n2_models2 <- PBMC_S.1_ci_n2_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_n2_models2 <- PBMC_S.1_ci_n2_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times
PBMC_S.1_ci_n2_it1_mod <- select(PBMC_S.1_ci_n2_models, Well, slope)
PBMC_S.1_ci_n2_it2_mod <- select(PBMC_S.1_ci_n2_models2, Well, slope)


PBMC_S.1_ci_n2_allmod <- full_join(PBMC_S.1_ci_n2_it1_mod, PBMC_S.1_ci_n2_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_S.1_ci_n2_excl ~ NA,
    Well%in%PBMC_S.1_ci_n2_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
PBMC_S.1_ci_n2_allmod$final_it_act <- -1*1000*60*(PBMC_S.1_ci_n2_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_S.1_ci_n2_allmod$Well)
PBMC_S.1_ci_n2_allmod <- PBMC_S.1_ci_n2_allmod[well_order, ]

#extract activities per well
PBMC_S.1_ci_n2_act <- select(PBMC_S.1_ci_n2_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(n2, t1a, t1b, t2a, t2b,
   PBMC_S.1_ci_n2_excl, PBMC_S.1_ci_n2_it_diff, 
   PBMC_S.1_ci_n2_it1, PBMC_S.1_ci_n2_it2, 
   PBMC_S.1_ci_n2_models, PBMC_S.1_ci_n2_models2,
   PBMC_S.1_ci_n2_it1_mod, PBMC_S.1_ci_n2_it2_mod,
   PBMC_S.1_ci_n2_allmod
)

#######################################################


################### NEGATIVE 3 ########################

#changing spectramax data into long format for 1 plate (n3 in this case)
PBMC_S.1_ci_n3 <- PBMC_S.1_ci[151:177, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
PBMC_S.1_ci_n3 <- PBMC_S.1_ci_n3[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
PBMC_S.1_ci_n3$Column <- as.numeric(PBMC_S.1_ci_n3$Column)
PBMC_S.1_ci_n3$mOD <- as.numeric(PBMC_S.1_ci_n3$mOD)

PBMC_S.1_ci_n3$Time <- period_to_seconds(hms(PBMC_S.1_ci_n3$Time))

PBMC_S.1_ci_n3 <- arrange(PBMC_S.1_ci_n3, Row)
PBMC_S.1_ci_n3 <- arrange(PBMC_S.1_ci_n3, Column)


#plot all traces as interactive
n3 <- ggplot(data = PBMC_S.1_ci_n3, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n3)

##REDUCTION

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES----
t1a <- 299
t1b <- 598

t2a <- 1
t2b <- 1

PBMC_S.1_ci_n3_excl <-  c("B12","C12","D12","E12","F12","G12")
PBMC_S.1_ci_n3_it_diff <- c()

PBMC_S.1_ci_n3_it1 <- PBMC_S.1_ci_n3[PBMC_S.1_ci_n3$Time>=t1a & PBMC_S.1_ci_n3$Time<=t1b,] 
PBMC_S.1_ci_n3_it2 <- PBMC_S.1_ci_n3[PBMC_S.1_ci_n3$Time>=t2a & PBMC_S.1_ci_n3$Time<=t2b,] 

#time1
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_n3_models <- PBMC_S.1_ci_n3_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_n3_models <- PBMC_S.1_ci_n3_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2
# Perform linear regression for each group using dplyr and purrr
PBMC_S.1_ci_n3_models2 <- PBMC_S.1_ci_n3_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
PBMC_S.1_ci_n3_models2 <- PBMC_S.1_ci_n3_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times
PBMC_S.1_ci_n3_it1_mod <- select(PBMC_S.1_ci_n3_models, Well, slope)
PBMC_S.1_ci_n3_it2_mod <- select(PBMC_S.1_ci_n3_models2, Well, slope)

PBMC_S.1_ci_n3_allmod <- full_join(PBMC_S.1_ci_n3_it1_mod, PBMC_S.1_ci_n3_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% PBMC_S.1_ci_n3_excl ~ NA,
    Well%in%PBMC_S.1_ci_n3_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))


# change slope into correct units
PBMC_S.1_ci_n3_allmod$final_it_act <- -1*1000*60*(PBMC_S.1_ci_n3_allmod$final_it_act)

#order by well
well_order <- mixedorder(PBMC_S.1_ci_n3_allmod$Well)
PBMC_S.1_ci_n3_allmod <- PBMC_S.1_ci_n3_allmod[well_order, ]

#extract activities per well
PBMC_S.1_ci_n3_act <- select(PBMC_S.1_ci_n3_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)

rm(n3, t1a, t1b, t2a, t2b,
   PBMC_S.1_ci_n3_excl, PBMC_S.1_ci_n3_it_diff, 
   PBMC_S.1_ci_n3_it1, PBMC_S.1_ci_n3_it2, 
   PBMC_S.1_ci_n3_models, PBMC_S.1_ci_n3_models2,
   PBMC_S.1_ci_n3_it1_mod, PBMC_S.1_ci_n3_it2_mod,
   PBMC_S.1_ci_n3_allmod
)

#######################################################


# join data together and calculate activities

cutoff_cv <- 0.15


molar_ext_coef <- 16.3
hom_buff_per_millcell <- 100
vol_sample_assay <-  15
vol_assay_mix <- 200
vol_well_total <- vol_sample_assay + vol_assay_mix

PBMC_S.1_ci_all <- full_join(PBMC_S.1_ci_p1_act, PBMC_S.1_ci_p2_act, by = "Well") %>%
  rename(p1 = Activity.x, p2 = Activity.y) %>%
  full_join(., PBMC_S.1_ci_p3_act, by= "Well") %>%
  rename(p3 = Activity) %>%
  full_join(., PBMC_S.1_ci_n1_act, by= "Well") %>%
  rename(n1 = Activity) %>%
  full_join(., PBMC_S.1_ci_n2_act, by= "Well") %>%
  rename(n2 = Activity) %>%
  full_join(., PBMC_S.1_ci_n3_act, by= "Well") %>%
  rename(n3 = Activity) %>%
  column_to_rownames(var="Well") %>%
  mutate(total_raw_mean = rowMeans(select(., starts_with("p")), na.rm = TRUE)) %>%
  mutate(total_raw_sd = rowSds(as.matrix(select(., starts_with("p"))), na.rm = TRUE)) %>%
  mutate(total_raw_cv = total_raw_sd/total_raw_mean) %>%
  mutate(nonspec_raw_mean = rowMeans(select(., starts_with("n")), na.rm = TRUE)) %>%
  mutate(nonspec_raw_sd = rowSds(as.matrix(select(., starts_with("n"))), na.rm = TRUE)) %>%
  mutate(nonspec_raw_cv = nonspec_raw_sd/nonspec_raw_mean) %>%
  mutate(spec_raw = total_raw_mean - nonspec_raw_mean) %>%
  mutate(p1_cleaned = ifelse(total_raw_cv > cutoff_cv & abs(p1-total_raw_mean) > abs(p2-total_raw_mean) & abs(p1-total_raw_mean) > abs(p3-total_raw_mean) & !is.na(p2) & !is.na(p3),NA, p1)) %>%
  mutate(p2_cleaned = ifelse(total_raw_cv > cutoff_cv & abs(p2-total_raw_mean) > abs(p1-total_raw_mean) & abs(p2-total_raw_mean) > abs(p3-total_raw_mean) & !is.na(p1) & !is.na(p3),NA, p2)) %>%
  mutate(p3_cleaned = ifelse(total_raw_cv > cutoff_cv & abs(p3-total_raw_mean) > abs(p1-total_raw_mean) & abs(p3-total_raw_mean) > abs(p2-total_raw_mean) & !is.na(p1) & !is.na(p2),NA, p3)) %>%
  mutate(n1_cleaned = ifelse(nonspec_raw_cv > cutoff_cv & abs(n1-nonspec_raw_mean) > abs(n2-nonspec_raw_mean) & abs(n1-nonspec_raw_mean) > abs(n3-nonspec_raw_mean) & !is.na(n2) & !is.na(n3),NA, n1)) %>%
  mutate(n2_cleaned = ifelse(nonspec_raw_cv > cutoff_cv & abs(n2-nonspec_raw_mean) > abs(n1-nonspec_raw_mean) & abs(n2-nonspec_raw_mean) > abs(n3-nonspec_raw_mean) & !is.na(n1) & !is.na(n3),NA, n2)) %>%
  mutate(n3_cleaned = ifelse(nonspec_raw_cv > cutoff_cv & abs(n3-nonspec_raw_mean) > abs(n1-nonspec_raw_mean) & abs(n3-nonspec_raw_mean) > abs(n2-nonspec_raw_mean) & !is.na(n1) & !is.na(n2),NA, n3)) %>%
  mutate(total_cleaned_mean = rowMeans(select(., starts_with("p") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(total_cleaned_sd = rowSds(as.matrix(select(., starts_with("p") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(total_cleaned_cv = total_cleaned_sd/total_cleaned_mean) %>%
  mutate(nonspec_cleaned_mean = rowMeans(select(., starts_with("n") & ends_with("cleaned")), na.rm = TRUE)) %>%
  mutate(nonspec_cleaned_sd = rowSds(as.matrix(select(., starts_with("n") & ends_with("cleaned"))), na.rm = TRUE)) %>%
  mutate(nonspec_cleaned_cv = nonspec_cleaned_sd/nonspec_cleaned_mean) %>%
  mutate(PBMC_S.1_ci_spec_cleaned = ifelse(total_cleaned_mean - nonspec_cleaned_mean>=0, total_cleaned_mean - nonspec_cleaned_mean, 0)) %>%
  mutate(PBMC_S.1_ci_spec_cleaned_cell = (((PBMC_S.1_ci_spec_cleaned/1000)*hom_buff_per_millcell)/(0.552*molar_ext_coef*vol_sample_assay*(vol_sample_assay/vol_well_total)))*1000)%>%
  rownames_to_column(var = "Well") %>%
  filter(!Well %in% c("B12","C12","D12","E12","F12", "G12"))
# write_xlsx(PBMC_S.1_ci_all, "interim_data/PBMC_S.1_ci_all.xlsx")
pre_av_tot_cv <- mean(abs(PBMC_S.1_ci_all$total_raw_cv), na.rm = T)
pre_av_non_cv <- mean(abs(PBMC_S.1_ci_all$nonspec_raw_cv), na.rm = T)
post_av_tot_cv <- mean(abs(PBMC_S.1_ci_all$total_cleaned_cv), na.rm = T)
post_av_non_cv <- mean(abs(PBMC_S.1_ci_all$nonspec_cleaned_cv), na.rm = T)
excl <- as.numeric(sum(is.na(PBMC_S.1_ci_all$p1_cleaned)) + sum(is.na(PBMC_S.1_ci_all$p2_cleaned)) + sum(is.na(PBMC_S.1_ci_all$p3_cleaned)) +
                     sum(is.na(PBMC_S.1_ci_all$n1_cleaned)) + sum(is.na(PBMC_S.1_ci_all$n2_cleaned)) + sum(is.na(PBMC_S.1_ci_all$n3_cleaned)))

25/(6*90)


PBMC_S.1_ci_byid <- PBMC_S.1_ci_all %>%
  select(Well, PBMC_S.1_ci_spec_cleaned, PBMC_S.1_ci_spec_cleaned_cell) %>%
  full_join(.,PBMC_S.1, by="Well") %>%
  select(Well, sampleName, PBMC_S.1_ci_spec_cleaned, PBMC_S.1_ci_spec_cleaned_cell)
write_xlsx(PBMC_S.1_ci_byid, "interim_data/manual/PBMC_S.1_ci.xlsx")

rm(PBMC_S.1_ci_p1, PBMC_S.1_ci_p2, PBMC_S.1_ci_p3, PBMC_S.1_ci_n1, PBMC_S.1_ci_n2, PBMC_S.1_ci_n3,
   PBMC_S.1_ci_p1_act, PBMC_S.1_ci_p2_act, PBMC_S.1_ci_p3_act, PBMC_S.1_ci_n1_act, PBMC_S.1_ci_n2_act, PBMC_S.1_ci_n3_act
)    


#visualize distributions -----

# 
# 
# pos_mean <- PBMC_S.1_ci_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_mean), color = "green")+
#   geom_density(aes(x=total_cleaned_mean), color = "blue")
# neg_mean <- PBMC_S.1_ci_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_mean), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_mean), color = "orange")
# 
# 
# pos_cv <- PBMC_S.1_ci_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_cv), color = "green")+
#   geom_density(aes(x=total_cleaned_cv), color = "blue")
# print(pos_cv)
# neg_cv <- PBMC_S.1_ci_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_cv), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_cv), color = "orange")
# print(neg_cv)
# 
# 
# fd <- ggplot(PBMC_S.1_ci_byid, mapping = aes(x = PBMC_S.1_ci_spec_cleaned)) +
#   geom_histogram(alpha = 0.3, fill = "dodgerblue", show.legend = F, bins = 50) +
#   geom_density(colour = "blue")
# print (fd)
# 
# PBMC_S.1_ci_byid %>%
#   ggplot()+
#   geom_density(aes(x=PBMC_S.1_ci_spec_cleaned))
# 
# plotly::ggplotly(fd)
# 
# rm(PBMC_S.1_ci, PBMC_S.1_ci_notemp, PBMC_S.1_ci_all, PBMC_S.1_ci_byid)

