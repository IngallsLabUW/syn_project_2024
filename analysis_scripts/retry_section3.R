# Section 3: B-MIS
# Modified by Claudia Luthy
# February 2024 for Syn salinity project

#upload packages and data! using trimmed (result of section 1) here bc
#QC (section 2) didn't provide anything we need to worry about
library(tidyverse)
library(dplyr)
library(stringr)

clean_dat_1 <- read_csv("~/Desktop/syn_project_2024/trimmed_HILIC_QE_POS_SynSalinityExperiment_SkylineReport.csv")

clean_data_1 <- clean_dat_1 %>%
  rename(samp= 'Replicate.Name', cmpd_name= 'Compound.Name', area= 'Area',
         cmpd_type= 'Molecule.List', samp_type= 'Precursor.Ion.Name', sal= 'Sal') %>%
  select(samp, cmpd_name, area, cmpd_type, samp_type, sal) %>%
  arrange(cmpd_name)

#Removing redundant mixes (lines 33-46 of Section2_CL.R)
mix_data_1 <- if (TRUE %in% grepl("Mix", clean_data_1$samp)) {
  Ingalls.Standards <- read_csv("~/Desktop/syn_project_2024/data_raw/Ingalls_Lab_Standards.csv") %>%
    filter(Compound_Name %in% clean_data_1$cmpd_name) %>%
    select(cmpd_name = Compound_Name, HILIC_Mix, conc = Concentration_uM) %>%
    unique()
  
  clean_data_std_1 <- clean_data_1 %>%
    filter(str_detect(samp, "Std")) %>%
    left_join(Ingalls.Standards, relationship = "many-to-many") %>%
    filter(str_detect(samp, as.character(HILIC_Mix)) | str_detect(samp, regex("H2OinMatrix", ignore_case = TRUE))) %>%
    select(-HILIC_Mix) %>%
    arrange(cmpd_name)
}

#Oscar's combination of clean data and mix data
tmp_1 <- clean_data_1 %>% filter(!samp_type=="Std") #This removes all rows of standard sample type
tmp_1$conc <- 0
clean_data_std_1 <- rbind(tmp_1,mix_data_1) %>% arrange(cmpd_name) %>% as.data.frame

#create IS subset
IS_data_1 <- clean_data_std_1 %>%
  filter(cmpd_type== "IS")

keep_IS_1 <- unique(IS_data_1$cmpd_name)
keep_IS_1 <- sapply(stringr::str_split(keep_IS_1, ","), "[", 1)

#so this is creating a table that calculates each best matched standard with the least variance
#for every single compound
select_data_1 <- clean_data_std_1 %>%
  select(samp, cmpd_name, area, cmpd_type, samp_type, sal, conc) %>%
  filter(cmpd_type=="S") %>%
  filter(samp_type=="Poo") %>%
  left_join(IS_data_1, by="samp", suffix=c("", "_IS"), relationship = "many-to-many") %>%
  select(samp, cmpd_name, cmpd_name_IS, area, area_IS, sal) %>%
  group_by(cmpd_name, cmpd_name_IS) %>%
  mutate(norm_area=area/area_IS*mean(area_IS)) %>%
  summarise(cv_IS=sd(norm_area)/mean(norm_area)) %>%
  arrange(cmpd_name, cv_IS) %>%
  slice(1)

#this makes a table with the bmised areas for each sample by matching it with the closest standard
#from the select_data calcs above
match_data_1 <- clean_data_std_1 %>%
  filter(cmpd_type=="S") %>%
  select(samp, cmpd_name, area, sal, conc) %>%
  left_join(select_data_1) %>%
  select(-cv_IS) %>%
  left_join(clean_data_std_1, by=c("samp", cmpd_name_IS="cmpd_name"), suffix=c("", "_IS")) %>%
  select(samp, cmpd_name, area, area_IS, sal, conc) %>%
  group_by(cmpd_name) %>%
  mutate(bmis_area=(area/area_IS)*mean(area_IS[1:112], na.rm = TRUE)) %>%
  select(samp, cmpd_name, bmis_area, sal, conc)

match_data_1 <- match_data_1 %>% filter(!cmpd_name %in% keep_IS_1)


#Manually calculate the area of the 13 compounds which did not pair correctly with their internal standard
#Use the clean_data_std1 to divide peak area of labeld compound by peak area of unlabeled

#Make a list of the 13 unique compounds that definitely did not pick correctly


keep_IS_1 <- c(keep_IS_1,"L-Alanine","L-Histidine", "L-Valine", "L-Proline")


tmp2 <- clean_data_std_1 %>% filter(grepl(paste(keep_IS_1,collapse = "|"),cmpd_name), samp_type=="Smp")
tmp2 <- tmp2 %>% filter(!cmpd_name=="L-Methionine S-oxide")
tmp2$cmpd_name_2 <- stringr::str_remove_all(tmp2$cmpd_name,"DL-")
tmp2$labeled <- sapply(stringr::str_split(tmp2$cmpd_name_2, ","), "[", 2)
tmp2.labeled <- tmp2 %>% filter(!is.na(labeled))
tmp2.unlabeled <- tmp2 %>% filter(is.na(labeled))
tmp2.unlabeled$bmis <- tmp2.unlabeled$area / tmp2.labeled$area

tmp2_fixed <- tmp2.unlabeled %>% select(samp, cmpd_name, bmis, sal, conc) %>%
  rename(bmis_area= 'bmis')


#create subset with standards for the 13 IS's (since they were taken out but we want them back now)
special_std <- mix_data_1 %>%
  filter(cmpd_name %in% keep_IS_1) %>%
  select(samp, cmpd_name, area, sal, conc) %>%
  rename(bmis_area= 'area') 
  #mix_data_1 %>% filter(cmpd_name %in% keep_IS_1) %>% select(cmpd_name) %>% distinct()
  


#need to somehow filter mix data to just have the compounds listed in tmp2 as cmpd_name
  #merge(mix_data_1, tmp2_fixed, by = "cmpd_name")
fixed <- rbind(special_std, tmp2_fixed)


  
  
  

# special <- clean_data_1 %>%
#   filter(samp_type== "Std")
# 
# special <- filter(unique(tmp2$cmpd_name))
# special <- sapply(stringr::str_split(special, ","), "[", 1)
# 
# tmp3 <- clean_data_std_1 %>% filter(grepl(paste(special,collapse = "|"),cmpd_name), samp_type=="Std")



corrected_IS_1 <- rbind(fixed, match_data_1) %>% as.data.frame()# %>%
  # left_join(tmp3)
# corrected_IS_2 <- rbind(corrected_IS_1, tmp3, by = 'cmpd_name', relationship = "many-to-many") 


