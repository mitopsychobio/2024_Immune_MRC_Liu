library(tidyverse)
library(janitor)
library(plotly)
library(gtools)
library(matrixStats)
library(writexl)

### not this has 3 int times for NEG 3 ONLY

#read file, get rid of temperature column ----

plt_sdh <- read.delim ("raw_data/spectramax_data/platelets/SDH-platelets.xls", header = TRUE, sep = "\t", fileEncoding = "UTF-16LE", skip = 2)[,-2]
plt <- readxl::read_xlsx("plate_maps/MiSBIE_plt.xlsx", sheet = 1)




################### POSITIVE 1 ########################

#changing spectramax data into long format for 1 plate (p1 in this case)
plt_sdh_p1 <- plt_sdh[0:31, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
plt_sdh_p1 <- plt_sdh_p1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
plt_sdh_p1$Column <- as.numeric(plt_sdh_p1$Column)
plt_sdh_p1$mOD <- as.numeric(plt_sdh_p1$mOD)

plt_sdh_p1$Time <- period_to_seconds(hms(plt_sdh_p1$Time))

plt_sdh_p1 <- arrange(plt_sdh_p1, Row)
plt_sdh_p1 <- arrange(plt_sdh_p1, Column)


#plot all traces as interactive
p1 <- ggplot(data = plt_sdh_p1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row)) +
  geom_point(aes(color = Row))

plotly::ggplotly(p1)

rm(p1)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 420
t1b <- 870

t2a <- 150
t2b <- 540

t3a <- 270
t3b <- 330

plt_sdh_p1_excl <- c()
plt_sdh_p1_it_diff <- c()
plt_sdh_p1_it_diff2 <-c()


plt_sdh_p1_it1 <- plt_sdh_p1[plt_sdh_p1$Time>=t1a & plt_sdh_p1$Time<=t1b,] 
plt_sdh_p1_it2 <- plt_sdh_p1[plt_sdh_p1$Time>=t2a & plt_sdh_p1$Time<=t2b,]
plt_sdh_p1_it3 <- plt_sdh_p1[plt_sdh_p1$Time>=t3a & plt_sdh_p1$Time<=t3b,] 


#time1-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p1_models <- plt_sdh_p1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p1_models <- plt_sdh_p1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p1_models2 <- plt_sdh_p1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p1_models2 <- plt_sdh_p1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time3-
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p1_models3 <- plt_sdh_p1_it3 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p1_models3 <- plt_sdh_p1_models3 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times------------------
plt_sdh_p1_it1_mod <- select(plt_sdh_p1_models, Well, slope)
plt_sdh_p1_it2_mod <- select(plt_sdh_p1_models2, Well, slope)
plt_sdh_p1_it3_mod <- select(plt_sdh_p1_models3, Well, slope)


plt_sdh_p1_allmod <- full_join(plt_sdh_p1_it1_mod, plt_sdh_p1_it2_mod, by = "Well")%>%
  rename(act_it1 = slope.x, act_it2 = slope.y) %>%
  full_join(., plt_sdh_p1_it3_mod, by = "Well") %>%
  rename(act_it3 = slope) %>%
  mutate(final_it_act = case_when(
    Well %in% plt_sdh_p1_excl ~ NA,
    Well%in%plt_sdh_p1_it_diff ~ act_it2,
    Well%in%plt_sdh_p1_it_diff2 ~ act_it3,
    TRUE ~ act_it1
  ))


# change slope into correct units
plt_sdh_p1_allmod$final_it_act <- -1*1000*60*(plt_sdh_p1_allmod$final_it_act)

#order by well
well_order <- mixedorder(plt_sdh_p1_allmod$Well)
plt_sdh_p1_allmod <- plt_sdh_p1_allmod[well_order, ]

#extract activities per well
plt_sdh_p1_act <- select(plt_sdh_p1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(t1a, t1b, t2a, t2b, t3a, t3b,
   plt_sdh_p1_excl, plt_sdh_p1_it_diff, plt_sdh_p1_it_diff2,
   plt_sdh_p1_it1, plt_sdh_p1_it2, plt_sdh_p1_it3,
   plt_sdh_p1_models, plt_sdh_p1_models2, plt_sdh_p1_models3,
   plt_sdh_p1_it1_mod, plt_sdh_p1_it2_mod, plt_sdh_p1_it3_mod,
   plt_sdh_p1_allmod
)




# ################## POSITIVE 2 ########################

#changing spectramax data into long format for 1 plate (p2 in this case)
plt_sdh_p2 <- plt_sdh[35:65, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
plt_sdh_p2 <- plt_sdh_p2[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
plt_sdh_p2$Column <- as.numeric(plt_sdh_p2$Column)
plt_sdh_p2$mOD <- as.numeric(plt_sdh_p2$mOD)

plt_sdh_p2$Time <- period_to_seconds(hms(plt_sdh_p2$Time))

plt_sdh_p2 <- arrange(plt_sdh_p2, Row)
plt_sdh_p2 <- arrange(plt_sdh_p2, Column)


#plot all traces as interactive
p2 <- ggplot(data = plt_sdh_p2, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p2)

rm(p2)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 480
t1b <- 900

t2a <- 0
t2b <- 120

t3a <- 150
t3b <- 330

plt_sdh_p2_excl <- c()
plt_sdh_p2_it_diff <- c()
plt_sdh_p2_it_diff2 <-c()


plt_sdh_p2_it1 <- plt_sdh_p2[plt_sdh_p2$Time>=t1a & plt_sdh_p2$Time<=t1b,] 
plt_sdh_p2_it2 <- plt_sdh_p2[plt_sdh_p2$Time>=t2a & plt_sdh_p2$Time<=t2b,]
plt_sdh_p2_it3 <- plt_sdh_p2[plt_sdh_p2$Time>=t3a & plt_sdh_p2$Time<=t3b,] 


#time1-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p2_models <- plt_sdh_p2_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p2_models <- plt_sdh_p2_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p2_models2 <- plt_sdh_p2_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p2_models2 <- plt_sdh_p2_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time3-
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p2_models3 <- plt_sdh_p2_it3 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p2_models3 <- plt_sdh_p2_models3 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times------------------
plt_sdh_p2_it1_mod <- select(plt_sdh_p2_models, Well, slope)
plt_sdh_p2_it2_mod <- select(plt_sdh_p2_models2, Well, slope)
plt_sdh_p2_it3_mod <- select(plt_sdh_p2_models3, Well, slope)


plt_sdh_p2_allmod <- full_join(plt_sdh_p2_it1_mod, plt_sdh_p2_it2_mod, by = "Well")%>%
  rename(act_it1 = slope.x, act_it2 = slope.y) %>%
  full_join(., plt_sdh_p2_it3_mod, by = "Well") %>%
  rename(act_it3 = slope) %>%
  mutate(final_it_act = case_when(
    Well %in% plt_sdh_p2_excl ~ NA,
    Well%in%plt_sdh_p2_it_diff ~ act_it2,
    Well%in%plt_sdh_p2_it_diff2 ~ act_it3,
    TRUE ~ act_it1
  ))


# change slope into correct units
plt_sdh_p2_allmod$final_it_act <- -1*1000*60*(plt_sdh_p2_allmod$final_it_act)

#order by well
well_order <- mixedorder(plt_sdh_p2_allmod$Well)
plt_sdh_p2_allmod <- plt_sdh_p2_allmod[well_order, ]

#extract activities per well
plt_sdh_p2_act <- select(plt_sdh_p2_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(t1a, t1b, t2a, t2b, t3a, t3b,
   plt_sdh_p2_excl, plt_sdh_p2_it_diff, plt_sdh_p2_it_diff2,
   plt_sdh_p2_it1, plt_sdh_p2_it2, plt_sdh_p2_it3,
   plt_sdh_p2_models, plt_sdh_p2_models2, plt_sdh_p2_models3,
   plt_sdh_p2_it1_mod, plt_sdh_p2_it2_mod, plt_sdh_p2_it3_mod,
   plt_sdh_p2_allmod
)





# ################## POSITIVE 3 ########################

#changing spectramax data into long format for 1 plate (p3 in this case)
plt_sdh_p3 <- plt_sdh[69:99, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
plt_sdh_p3 <- plt_sdh_p3[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
plt_sdh_p3$Column <- as.numeric(plt_sdh_p3$Column)
plt_sdh_p3$mOD <- as.numeric(plt_sdh_p3$mOD)

plt_sdh_p3$Time <- period_to_seconds(hms(plt_sdh_p3$Time))

plt_sdh_p3 <- arrange(plt_sdh_p3, Row)
plt_sdh_p3 <- arrange(plt_sdh_p3, Column)


#plot all traces as interactive
p3 <- ggplot(data = plt_sdh_p3, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(p3)

rm(p3)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 480
t1b <- 900

t2a <- 0
t2b <- 120

t3a <- 180
t3b <- 330

plt_sdh_p3_excl <- c("A7","B7","C7","D7","E7","F7","G7","G6")
plt_sdh_p3_it_diff <- c()
plt_sdh_p3_it_diff2 <-c()


plt_sdh_p3_it1 <- plt_sdh_p3[plt_sdh_p3$Time>=t1a & plt_sdh_p3$Time<=t1b,] 
plt_sdh_p3_it2 <- plt_sdh_p3[plt_sdh_p3$Time>=t2a & plt_sdh_p3$Time<=t2b,]
plt_sdh_p3_it3 <- plt_sdh_p3[plt_sdh_p3$Time>=t3a & plt_sdh_p3$Time<=t3b,] 


#time1-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p3_models <- plt_sdh_p3_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p3_models <- plt_sdh_p3_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p3_models2 <- plt_sdh_p3_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p3_models2 <- plt_sdh_p3_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time3-
# Perform linear regression for each group using dplyr and purrr
plt_sdh_p3_models3 <- plt_sdh_p3_it3 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_p3_models3 <- plt_sdh_p3_models3 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times------------------
plt_sdh_p3_it1_mod <- select(plt_sdh_p3_models, Well, slope)
plt_sdh_p3_it2_mod <- select(plt_sdh_p3_models2, Well, slope)
plt_sdh_p3_it3_mod <- select(plt_sdh_p3_models3, Well, slope)


plt_sdh_p3_allmod <- full_join(plt_sdh_p3_it1_mod, plt_sdh_p3_it2_mod, by = "Well")%>%
  rename(act_it1 = slope.x, act_it2 = slope.y) %>%
  full_join(., plt_sdh_p3_it3_mod, by = "Well") %>%
  rename(act_it3 = slope) %>%
  mutate(final_it_act = case_when(
    Well %in% plt_sdh_p3_excl ~ NA,
    Well%in%plt_sdh_p3_it_diff ~ act_it2,
    Well%in%plt_sdh_p3_it_diff2 ~ act_it3,
    TRUE ~ act_it1
  ))


# change slope into correct units
plt_sdh_p3_allmod$final_it_act <- -1*1000*60*(plt_sdh_p3_allmod$final_it_act)

#order by well
well_order <- mixedorder(plt_sdh_p3_allmod$Well)
plt_sdh_p3_allmod <- plt_sdh_p3_allmod[well_order, ]

#extract activities per well
plt_sdh_p3_act <- select(plt_sdh_p3_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(t1a, t1b, t2a, t2b, t3a, t3b,
   plt_sdh_p3_excl, plt_sdh_p3_it_diff, plt_sdh_p3_it_diff2,
   plt_sdh_p3_it1, plt_sdh_p3_it2, plt_sdh_p3_it3,
   plt_sdh_p3_models, plt_sdh_p3_models2, plt_sdh_p3_models3,
   plt_sdh_p3_it1_mod, plt_sdh_p3_it2_mod, plt_sdh_p3_it3_mod,
   plt_sdh_p3_allmod
)




#######################################################








# ################## NEGATIVE 1 ########################

#changing spectramax data into long format for 1 plate (n1 in this case)
plt_sdh_n1 <- plt_sdh[103:133, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
plt_sdh_n1 <- plt_sdh_n1[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
plt_sdh_n1$Column <- as.numeric(plt_sdh_n1$Column)
plt_sdh_n1$mOD <- as.numeric(plt_sdh_n1$mOD)

plt_sdh_n1$Time <- period_to_seconds(hms(plt_sdh_n1$Time))

plt_sdh_n1 <- arrange(plt_sdh_n1, Row)
plt_sdh_n1 <- arrange(plt_sdh_n1, Column)


#plot all traces as interactive
n1 <- ggplot(data = plt_sdh_n1, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n1)


rm(n1)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 300
t1b <- 900

t2a <- 30
t2b <- 150

plt_sdh_n1_excl <- c()
plt_sdh_n1_it_diff <- c()

plt_sdh_n1_it1 <- plt_sdh_n1[plt_sdh_n1$Time>=t1a & plt_sdh_n1$Time<=t1b,] 
plt_sdh_n1_it2 <- plt_sdh_n1[plt_sdh_n1$Time>=t2a & plt_sdh_n1$Time<=t2b,] 

#time1--
# Perform linear regression for each group using dplyr and purrr
plt_sdh_n1_models <- plt_sdh_n1_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_n1_models <- plt_sdh_n1_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_n1_models2 <- plt_sdh_n1_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_n1_models2 <- plt_sdh_n1_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times- 
plt_sdh_n1_it1_mod <- select(plt_sdh_n1_models, Well, slope)
plt_sdh_n1_it2_mod <- select(plt_sdh_n1_models2, Well, slope)

plt_sdh_n1_allmod <- full_join(plt_sdh_n1_it1_mod, plt_sdh_n1_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% plt_sdh_n1_excl ~ NA,
    Well%in%plt_sdh_n1_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
plt_sdh_n1_allmod$final_it_act <- -1*1000*60*(plt_sdh_n1_allmod$final_it_act)

#order by well
well_order <- mixedorder(plt_sdh_n1_allmod$Well)
plt_sdh_n1_allmod <- plt_sdh_n1_allmod[well_order, ]

#extract activities per well
plt_sdh_n1_act <- select(plt_sdh_n1_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(t1a, t1b, t2a, t2b,
   plt_sdh_n1_excl, plt_sdh_n1_it_diff, 
   plt_sdh_n1_it1, plt_sdh_n1_it2, 
   plt_sdh_n1_models, plt_sdh_n1_models2,
   plt_sdh_n1_it1_mod, plt_sdh_n1_it2_mod,
   plt_sdh_n1_allmod
)



# ################## NEGATIVE 2 ########################

#changing spectramax data into long format for 1 plate (n2 in this case)
plt_sdh_n2 <- plt_sdh[137:167, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
plt_sdh_n2 <- plt_sdh_n2[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
plt_sdh_n2$Column <- as.numeric(plt_sdh_n2$Column)
plt_sdh_n2$mOD <- as.numeric(plt_sdh_n2$mOD)

plt_sdh_n2$Time <- period_to_seconds(hms(plt_sdh_n2$Time))

plt_sdh_n2 <- arrange(plt_sdh_n2, Row)
plt_sdh_n2 <- arrange(plt_sdh_n2, Column)


#plot all traces as interactive
n2 <- ggplot(data = plt_sdh_n2, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n2)

rm(n3)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 0
t1b <- 900

t2a <- 30
t2b <- 150

plt_sdh_n2_excl <- c()
plt_sdh_n2_it_diff <- c()

plt_sdh_n2_it1 <- plt_sdh_n2[plt_sdh_n2$Time>=t1a & plt_sdh_n2$Time<=t1b,] 
plt_sdh_n2_it2 <- plt_sdh_n2[plt_sdh_n2$Time>=t2a & plt_sdh_n2$Time<=t2b,] 

#time1-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_n2_models <- plt_sdh_n2_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_n2_models <- plt_sdh_n2_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2--
# Perform linear regression for each group using dplyr and purrr
plt_sdh_n2_models2 <- plt_sdh_n2_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_n2_models2 <- plt_sdh_n2_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times-- 
plt_sdh_n2_it1_mod <- select(plt_sdh_n2_models, Well, slope)
plt_sdh_n2_it2_mod <- select(plt_sdh_n2_models2, Well, slope)

plt_sdh_n2_allmod <- full_join(plt_sdh_n2_it1_mod, plt_sdh_n2_it2_mod, by = "Well") %>%
  rename(act_it1 = slope.x, act_it2 = slope.y)%>%
  mutate(final_it_act = case_when(
    Well %in% plt_sdh_n2_excl ~ NA,
    Well%in%plt_sdh_n2_it_diff ~ act_it2,
    TRUE ~ act_it1
  ))

# change slope into correct units
plt_sdh_n2_allmod$final_it_act <- -1*1000*60*(plt_sdh_n2_allmod$final_it_act)

#order by well
well_order <- mixedorder(plt_sdh_n2_allmod$Well)
plt_sdh_n2_allmod <- plt_sdh_n2_allmod[well_order, ]

#extract activities per well
plt_sdh_n2_act <- select(plt_sdh_n2_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(t1a, t1b, t2a, t2b,
   plt_sdh_n2_excl, plt_sdh_n2_it_diff, 
   plt_sdh_n2_it1, plt_sdh_n2_it2, 
   plt_sdh_n2_models, plt_sdh_n2_models2,
   plt_sdh_n2_it1_mod, plt_sdh_n2_it2_mod,
   plt_sdh_n2_allmod
)




# ################## NEGATIVE 3 ########################

#changing spectramax data into long format for 1 plate (n3 in this case)
plt_sdh_n3 <- plt_sdh[171:201, 1:97]%>%
  pivot_longer(cols = -"Time", names_to = "Well", values_to="mOD")

#adding row/column to make sorting easier, cleaning time column so in sec
plt_sdh_n3 <- plt_sdh_n3[,c("Well","Time","mOD")]%>%
  separate(Well,
           into = c("Row", "Column"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           remove = FALSE
  )
plt_sdh_n3$Column <- as.numeric(plt_sdh_n3$Column)
plt_sdh_n3$mOD <- as.numeric(plt_sdh_n3$mOD)

plt_sdh_n3$Time <- period_to_seconds(hms(plt_sdh_n3$Time))

plt_sdh_n3 <- arrange(plt_sdh_n3, Row)
plt_sdh_n3 <- arrange(plt_sdh_n3, Column)


#plot all traces as interactive
n3 <- ggplot(data = plt_sdh_n3, mapping = aes(x=Time, y=mOD, line=Well)) +
  geom_line(aes(color = Row))+
  geom_point(aes(color = Row))

plotly::ggplotly(n3)
rm(n3)

##REDUCTION##------------------

# to get reduction: per well, make linear regression between specified time points and then get slope

#subsetting data for correct int times

##PICK INTEGRATION TIMES------------------
t1a <- 300
t1b <- 780

t2a <- 30
t2b <- 150

t3a <- 30
t3b <- 510

plt_sdh_n3_excl <- c()
plt_sdh_n3_it_diff <- c()
plt_sdh_n3_it_diff2 <- c()


plt_sdh_n3_it1 <- plt_sdh_n3[plt_sdh_n3$Time>=t1a & plt_sdh_n3$Time<=t1b,] 
plt_sdh_n3_it2 <- plt_sdh_n3[plt_sdh_n3$Time>=t2a & plt_sdh_n3$Time<=t2b,]
plt_sdh_n3_it3 <- plt_sdh_n3[plt_sdh_n3$Time>=t3a & plt_sdh_n3$Time<=t3b,] 


#time1-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_n3_models <- plt_sdh_n3_it1 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_n3_models <- plt_sdh_n3_models %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time2-- 
# Perform linear regression for each group using dplyr and purrr
plt_sdh_n3_models2 <- plt_sdh_n3_it2 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_n3_models2 <- plt_sdh_n3_models2 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#time3-
# Perform linear regression for each group using dplyr and purrr
plt_sdh_n3_models3 <- plt_sdh_n3_it3 %>%
  group_by(Well)%>%
  nest() %>%
  mutate(model = map(data, ~lm(mOD ~ Time, data = .))) 

# Access coefficients for each group
plt_sdh_n3_models3 <- plt_sdh_n3_models3 %>%
  mutate(coef = map(model, coef)) %>%
  mutate(slope = map_dbl(coef, ~ .x[2]))

#all int times------------------
plt_sdh_n3_it1_mod <- select(plt_sdh_n3_models, Well, slope)
plt_sdh_n3_it2_mod <- select(plt_sdh_n3_models2, Well, slope)
plt_sdh_n3_it3_mod <- select(plt_sdh_n3_models3, Well, slope)


plt_sdh_n3_allmod <- full_join(plt_sdh_n3_it1_mod, plt_sdh_n3_it2_mod, by = "Well")%>%
  rename(act_it1 = slope.x, act_it2 = slope.y) %>%
  full_join(., plt_sdh_n3_it3_mod, by = "Well") %>%
  rename(act_it3 = slope) %>%
  mutate(final_it_act = case_when(
    Well %in% plt_sdh_n3_excl ~ NA,
    Well%in%plt_sdh_n3_it_diff ~ act_it2,
    Well%in%plt_sdh_n3_it_diff2 ~ act_it3,
    TRUE ~ act_it1
  ))


# change slope into correct units
plt_sdh_n3_allmod$final_it_act <- -1*1000*60*(plt_sdh_n3_allmod$final_it_act)

#order by well
well_order <- mixedorder(plt_sdh_n3_allmod$Well)
plt_sdh_n3_allmod <- plt_sdh_n3_allmod[well_order, ]

#extract activities per well
plt_sdh_n3_act <- select(plt_sdh_n3_allmod, Well, final_it_act) %>%
  rename(Activity = final_it_act)


rm(t1a, t1b, t2a, t2b, t3a, t3b,
   plt_sdh_n3_excl, plt_sdh_n3_it_diff, plt_sdh_n3_it_diff2,
   plt_sdh_n3_it1, plt_sdh_n3_it2, plt_sdh_n3_it3,
   plt_sdh_n3_models, plt_sdh_n3_models2, plt_sdh_n3_models3,
   plt_sdh_n3_it1_mod, plt_sdh_n3_it2_mod, plt_sdh_n3_it3_mod,
   plt_sdh_n3_allmod
)


#######################################################


# join data together and calculate activities ----

cutoff_cv <- 0.15


molar_ext_coef <- 16.3
hom_buff_per_millcell <- 100
vol_sample_assay <-  20
vol_assay_mix <- 200
vol_well_total <- vol_sample_assay + vol_assay_mix


plt_sdh_all <- full_join(plt_sdh_p1_act, plt_sdh_p2_act, by = "Well") %>%
  rename(p1 = Activity.x, p2 = Activity.y) %>%
  full_join(., plt_sdh_p3_act, by= "Well") %>%
  rename(p3 = Activity) %>%
  full_join(., plt_sdh_n1_act, by= "Well") %>%
  rename(n1 = Activity) %>%
  full_join(., plt_sdh_n2_act, by= "Well") %>%
  rename(n2 = Activity) %>%
  full_join(., plt_sdh_n3_act, by= "Well") %>%
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
  mutate(plt_sdh_spec_cleaned = ifelse(total_cleaned_mean - nonspec_cleaned_mean>=0, total_cleaned_mean - nonspec_cleaned_mean, 0)) %>%
  mutate(plt_sdh_spec_cleaned_cell = (((plt_sdh_spec_cleaned/1000)*hom_buff_per_millcell)/(0.552*molar_ext_coef*vol_sample_assay*(vol_sample_assay/vol_well_total)))*1000)%>%
  rownames_to_column(var = "Well") %>%
  filter(!Well %in% c("B12","C12","D12","E12","F12", "G12",
                      "A10","B10","C10","D10","E10","F10","G10","H12",
                      "A11","B11","C11","D11","E11","F11","G11","H11",
                      "A9","B9","C9","D9","E9","F9","G9","H9",
                      "A8","B8","C8","D8","E8","F8","G8","H8",
                      "H7"))

# write_xlsx(plt_sdh_all, "interim_data/plt_sdh_all.xlsx")
pre_av_tot_cv <- mean(abs(plt_sdh_all$total_raw_cv), na.rm = T)
pre_av_non_cv <- mean(abs(plt_sdh_all$nonspec_raw_cv), na.rm = T)
post_av_tot_cv <- mean(abs(plt_sdh_all$total_cleaned_cv), na.rm = T)
post_av_non_cv <- mean(abs(plt_sdh_all$nonspec_cleaned_cv), na.rm = T)
excl <- as.numeric(sum(is.na(plt_sdh_all$p1_cleaned)) + sum(is.na(plt_sdh_all$p2_cleaned)) + sum(is.na(plt_sdh_all$p3_cleaned)) +
                     sum(is.na(plt_sdh_all$n1_cleaned)) + sum(is.na(plt_sdh_all$n2_cleaned)) + sum(is.na(plt_sdh_all$n3_cleaned)))

57/(6*57)



plt_sdh_byid <- plt_sdh_all %>%
  select(Well, plt_sdh_spec_cleaned, plt_sdh_spec_cleaned_cell) %>%
  full_join(.,plt, by="Well") %>%
  select(Well, sampleName, plt_sdh_spec_cleaned, plt_sdh_spec_cleaned_cell)
write_xlsx(plt_sdh_byid, "interim_data/manual/plt_sdh.xlsx")

rm(plt_sdh_p1, plt_sdh_p2, plt_sdh_p3, plt_sdh_n1, plt_sdh_n2, plt_sdh_n3,
   plt_sdh_p1_act, plt_sdh_p2_act, plt_sdh_p3_act, plt_sdh_n1_act, plt_sdh_n2_act, plt_sdh_n3_act
)   

# #visualize distributions -----
# 
# 
# 
# pos_mean <- plt_sdh_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_mean), color = "green")+
#   geom_density(aes(x=total_cleaned_mean), color = "blue")
# neg_mean <- plt_sdh_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_mean), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_mean), color = "orange")
# 
# 
# pos_cv <- plt_sdh_all %>%
#   ggplot() +
#   geom_density(aes(x=total_raw_cv), color = "green")+
#   geom_density(aes(x=total_cleaned_cv), color = "blue")
# print(pos_cv)
# neg_cv <- plt_sdh_all %>%
#   ggplot() + 
#   geom_density(aes(x=nonspec_raw_cv), color = "yellow")+
#   geom_density(aes(x=nonspec_cleaned_cv), color = "orange")
# print(neg_cv)
# 
# 
# fd <- ggplot(plt_sdh_byid, mapping = aes(x = plt_sdh_spec_cleaned)) +
#   geom_histogram(alpha = 0.3, fill = "dodgerblue", show.legend = F, bins = 50) +
#   geom_density(colour = "blue")
# print (fd)
# 
# plt_sdh_byid %>%
#   ggplot()+
#   geom_density(aes(x=plt_sdh_spec_cleaned))
# 
# plotly::ggplotly(fd)
# 
# rm(plt_sdh, plt_sdh_notemp, plt_sdh_all, plt_sdh_byid)
