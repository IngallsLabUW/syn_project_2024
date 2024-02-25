#Claudia Luthy
#Feb 2024
#Syn salinity project

#section I. to clean and start re-arranging skyline data

#plan:select the data I do want and reformat, also add new column for salinity?
#don't forget abt the standards vs internal standards (IS?) bit, could use select, mutate, and ifelse? Just live in deplyr land?
#backup plan would be grep

#need the right packages and to import the data
library(tidyverse)
library(dplyr)
library(stringr)

transition_list <- read.csv("data_raw/HILIC_QE_POS_SynSalinityExperiment_SkylineReport.csv")
View(transition_list)

#want to trim the data down
trimmed <- transition_list %>%
  select(Replicate.Name= 'Replicate.Name', Compound.Name = 'Precursor.Ion.Name', Retention.Time='Retention.Time', Area='Area', Height = 'Height',
         Background = 'Background', Molecule.List = 'Molecule.List') %>%
  #so far this worked! now area as number
  mutate(Area= as.numeric(Area)) %>%
  mutate(Height= as.numeric(Height)) %>%
  mutate(Background= as.numeric(Background)) %>%
  mutate(Retention.Time= as.numeric(Retention.Time)) %>%
  #now want cmpd_type as S or IS for standard versus internal standard
  mutate(Molecule.List=ifelse(str_detect(Molecule.List, "HILIC_Pos_InternalStandards"), "IS", "S")) %>%
  #next will be to do sample vs blank vs pooled vs standard
  mutate(Precursor.Ion.Name=str_extract(Replicate.Name, "Poo|Blk|Smp|Std")) %>%
  #want to make salinity column, can do this with also str_extract
  mutate(Sal=str_extract(Replicate.Name, '25ppt|30ppt|35ppt|40ppt'), Sal=str_remove(Sal, "ppt")) %>%
  mutate(Sal=as.numeric(Sal)) %>%
  arrange(Compound.Name)
view(trimmed)

#Save trimmed skyline data
#save(trimmed, file="~/Desktop/syn_project_2024/trimmed_HILIC_QE_POS_SynSalinityExperiment_SkylineReport.Rdata")
write_csv(trimmed, file="~/Desktop/syn_project_2024/trimmed_HILIC_QE_POS_SynSalinityExperiment_SkylineReport.csv")

#notes for myself later:
#okay remember the error for the numerical area/rt/sal is just saying that N/As are present in the column now
# ask Susan if we even want the sample type as standard versus internal standard bit. I remember she had that but do we need it?
#you had first tried grep and filter functions for salinity, I think this is actually the easiest now. Ask Susan if it looks fine.

