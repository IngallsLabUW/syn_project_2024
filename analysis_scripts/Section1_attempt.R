#Claudia Luthy
#Feb 2024
#Syn salinity project

#section I. to clean and start re-arranging skyline data

#goals:have the replicate name as basic name?, ion name as compound, retention time + area stay + molecule.list,
#sample type for samp vs std vs blank sorta deal, and maybe salinity?
#and then remove molecule list name, height, background, and mass error ppm.... I think

#plan:select the data I do want and reformat, also add new column for salinity?
#don't forget abt the standards vs internal standards (IS?) bit, could use select, mutate, and ifelse? Just live in deplyr land?
#backup plan would be grep

#need the right packages and to import the data
library(tidyverse)
library(dplyr)
library(stringr)

transition_list <- read.csv("data_raw/HILIC_QE_POS_SynSalinityExperiment_SkylineReport.csv")
view(transition_list)

#want to trim the data down
trimmed <- transition_list %>%
  select(samp_name= 'Replicate.Name', cmpd= 'Precursor.Ion.Name', rt='Retention.Time', area='Area', cmpd_type = 'Molecule.List') %>%
  #so far this worked! now area as number
  mutate(area= as.numeric(area)) %>%
  mutate(rt= as.numeric(rt)) %>%
  #now want cmpd_type as S or IS for standard versus internal standard
  mutate(cmpd_type=ifelse(str_detect(cmpd_type, "HILIC_Pos_InternalStandards"), "IS", "S")) %>%
  #next will be to do sample vs blank vs pooled vs standard
  mutate(samp_type=str_extract(samp_name, "Poo|Blk|Smp|Std")) %>%
  #want to make salinity column, can do this with also str_extract
  mutate(sal=str_extract(samp_name, '25ppt|30ppt|35ppt|40ppt'), sal=str_remove(sal, "ppt")) %>%
  mutate(sal=as.numeric(sal))
  arrange(cmpd)
view(trimmed)



#notes for myself later:
#okay remember the error for the numerical area/rt/sal is just saying that N/As are present in the column now
# ask Susan if we even want the sample type as standard versus internal standard bit. I remember she had that but do we need it?
#you had first tried grep and filter functions for salinity, I think this is actually the easiest now. Ask Susan if it looks fine.







#failed attempts:
#    sal <- filter(str_extract(samp_name, '25ppt'))
#    mutate(sal=ifelse(samp_name= "/25ppt+")) %>%
#    mutate(sal=ifelse(str_detect(samp_name, "25ppt|30ppt|35ppt|40ppt"))) %>%
#    mutate(sal=str_extract(samp_name, '25'="25")) %>%
#    grep1(sal=samp_name)
  # mutate(sal=str_extract(samp_name, '25ppt|30ppt|35ppt|40ppt')str_remove(sal, "ppt")) %>%
  # mutate(sal=str_extract(samp_name, '25ppt|30ppt|35ppt|40ppt', str_remove(sal, "ppt"))) %>%
  # sal=str_remove(sal, "ppt") %>%


