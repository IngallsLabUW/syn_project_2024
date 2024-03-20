# Section 3: B-MIS
# Modified by Claudia Luthy
# Spring 2024 for Syn salinity project

#inputs: trimmed file from section 1
#outputs: clean_data_std which has raw areas and everything,
#also match_data which has BMISed areas

#upload packages and data! using trimmed (result of section 1) here bc
#QC (section 2) didn't provide anything we need to worry about
library(tidyverse)
library(dplyr)
library(stringr)

clean_dat <- read_csv("~/Desktop/syn_project_2024/trimmed_HILIC_QE_POS_SynSalinityExperiment_SkylineReport.csv")

clean_data <- clean_dat %>%
  rename(samp= 'Replicate.Name', cmpd_name= 'Compound.Name', area= 'Area',
         cmpd_type= 'Molecule.List', samp_type= 'Precursor.Ion.Name', sal= 'Sal') %>%
  select(samp, cmpd_name, area, cmpd_type, samp_type, sal) %>%
  arrange(cmpd_name)

#sanity check for a single compound to see if the chosen IS matches
#so can switch out Homarine for any compound
clean_data %>%
  filter(cmpd_name=="Glycine betaine" | cmpd_type=="IS") %>%
  filter(samp_type=="Poo") %>%
  group_by(cmpd_name) %>%
  ggplot() +
  geom_col(aes(x=samp, y=area)) +
  facet_wrap(~cmpd_name, ncol = 2, scales = "free_y")

#Removing redundant mixes (lines 33-46 of Section2_CL.R)
mix_data <- if (TRUE %in% grepl("Mix", clean_data$samp)) {
  Ingalls.Standards <- read_csv("~/Desktop/syn_project_2024/data_raw/Ingalls_Lab_Standards.csv") %>%
    filter(Compound_Name %in% clean_data$cmpd_name) %>%
    select(cmpd_name = Compound_Name, HILIC_Mix, conc = Concentration_uM) %>%
    unique()

  clean_data_std <- clean_data %>%
    filter(str_detect(samp, "Std")) %>%
    left_join(Ingalls.Standards, relationship = "many-to-many") %>%
    filter(str_detect(samp, as.character(HILIC_Mix)) | str_detect(samp, regex("H2OinMatrix", ignore_case = TRUE))) %>%
    select(-HILIC_Mix) %>%
    arrange(cmpd_name)
}

#Oscar's combination of clean data and mix data
tmp <- clean_data %>% filter(!samp_type=="Std") #This removes all rows of standard sample type
#semi_join(Ingalls.Standards, mix_data_1, by = "cmpd_name", copy = FALSE)
#clean_data_std_1 <- rbind(tmp_1,mix_data_1) %>% arrange(cmpd_name) %>% as.data.frame

#tmp_1 <- left_join(tmp_1,mix_data_1%>%select(samp,cmpd_name,conc))
tmp$conc <- 0

clean_data_std <- rbind(tmp,mix_data) %>% arrange(cmpd_name) %>% as.data.frame

write_csv(clean_data_std, file="~/Desktop/syn_project_2024/intermediates/clean_data_std.csv")

#create IS subset
IS_data <- clean_data_std %>%
  filter(cmpd_type== "IS")

keep_IS <- unique(IS_data$cmpd_name)
keep_IS <- sapply(stringr::str_split(keep_IS, ","), "[", 1)

#so this is creating a table that calculates each best matched standard with the least variance
#for every single compound
select_data <- clean_data_std %>%
  select(samp, cmpd_name, area, cmpd_type, samp_type, sal, conc) %>%
  filter(cmpd_type=="S") %>%
  filter(samp_type=="Poo") %>%
  left_join(IS_data, by="samp", suffix=c("", "_IS"), relationship = "many-to-many") %>%
  select(samp, cmpd_name, cmpd_name_IS, area, area_IS, sal) %>%
  group_by(cmpd_name, cmpd_name_IS) %>%
  mutate(norm_area=area/area_IS*mean(area_IS)) %>%
  summarise(cv_IS=sd(norm_area)/mean(norm_area)) %>%
  arrange(cmpd_name, cv_IS) %>%
  slice(1)

#this makes a table with the bmised areas for each sample by matching it with the closest standard
#from the select_data calcs above
match_data <- clean_data_std %>%
  filter(cmpd_type=="S") %>%
  select(samp, cmpd_name, area, sal, conc) %>%
  left_join(select_data) %>%
  select(-cv_IS) %>%
  left_join(clean_data_std, by=c("samp", cmpd_name_IS="cmpd_name"), suffix=c("", "_IS")) %>%
  select(samp, cmpd_name, area, area_IS, sal, conc) %>%
  group_by(cmpd_name) %>%
  mutate(bmis_area=(area/area_IS)*mean(area_IS[1:112], na.rm = TRUE)) %>%
  select(samp, cmpd_name, bmis_area, sal, conc)

match_data <- match_data %>% filter(!cmpd_name %in% keep_IS)

write_csv(match_data, file="~/Desktop/syn_project_2024/intermediates/BMISed_areas.csv")





#these seem to be the important parts of Regina's code, general idea being:
#1) subset internal standards
#2)calculate mean values and normalize the mean values of the IS
#3)do BMIS!

#part 1 here:
# Create Internal Standard data
# Int.Stds.data <- Data.withIS %>%
#   select(Replicate.Name, Metabolite.Name, Area.with.QC) %>%
#   rename(Mass.Feature = Metabolite.Name)

#part 2 here:
# Calculate mean values for each Internal Standard--------------------------------------------------------
# Int.Stds.means <- Int.Stds.data %>%
#   select(-c("a", "RunType", "c", "d")) %>%
#   group_by(Mass.Feature) %>%
#   summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE))

# # Normalize to each internal Standard--------------------------------------------------------------------
# Data.bound <- rbind(Int.Stds.data %>% select(-c("a", "RunType", "c", "d")), Data.long) %>%
#   arrange(Mass.Feature)

#part 3 here:
# # Find the B-MIS for each MassFeature-------------------------------------------------------------------
#
# # Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
# # then choose which IS reduces the RSD the most (Poo.Picked.IS)
# Poodata1 <- Mydata.new %>%
#   filter(type == "Poo") %>%
#   group_by(SampID, Mass.Feature, MIS) %>%
#   mutate(RSD_ofPoo_IND = sd(Adjusted.Area, na.rm = TRUE) / mean(Adjusted.Area, na.rm = TRUE)) %>%
#   mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
#   group_by(Mass.Feature, MIS) %>%
#   summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
#   mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo))
#
#
# Poodata2 <- Poodata1 %>%
#   left_join(Poodata1 %>% group_by(Mass.Feature) %>%
#               summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))
#
#
# # Get the original RSD, calculate RSD change, decide if MIS is acceptable -------------------------------
# Poodata3 <- left_join(Poodata2, Poodata2 %>%
#                        filter(MIS == "Inj_vol" ) %>%
#                        mutate(Orig_RSD = RSD_ofPoo) %>%
#                        select(-RSD_ofPoo, -MIS)) %>%
#   mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
#   mutate(percent.Change = del_RSD/Orig_RSD) %>%
#   mutate(accept_MIS = (percent.Change > cut.off & Orig_RSD > cut.off2))








#this is the original script from Regina:

# Match QC'd data with Internal Standards list ----------------------------------------------------------
# Data.withIS <- QCd.data %>%
#   filter(Metabolite.Name %in% Internal.Standards$Compound_Name)
#
# Data.NoIS <- QCd.data %>%
#   filter(!Metabolite.Name %in% Internal.Standards$Compound_Name)

# Create Internal Standard data
# Int.Stds.data <- Data.withIS %>%
#   select(Replicate.Name, Metabolite.Name, Area.with.QC) %>%
#   rename(Mass.Feature = Metabolite.Name)

# # Add injection volume ---------------------------------------------------------------------------------
# SampKey <- QCd.data %>%
#   select(Replicate.Name) %>%
#   mutate(Area.with.QC = ifelse(str_detect(Replicate.Name, "Half"), 0.5, 1.0)) %>%
#   mutate(Mass.Feature = "Inj_vol")

# Create Internal standard data to identify problematic compounds/replicates ---------------------------
# Int.Stds.data <- rbind(Int.Stds.data, SampKey) %>%
#   ####
#   separate(Replicate.Name, into = c("a", "RunType", "c", "d"), sep = "_", remove = FALSE)

# # Identify internal standards without an Area, i.e. any NA values.
# IS.Issues <- Int.Stds.data[is.na(Int.Stds.data$Area.with.QC), ]
# write.csv(IS.Issues, paste("data_intermediate/MSDial_InternalStdIssues_", currentDate, ".csv", sep = ""))

# # Visualize raw areas of Internal Standards -------------------------------------------------------------
# IS.Raw.Area.Plot <- ggplot(Int.Stds.data, aes(x = Replicate.Name, y = Area.with.QC, color = RunType)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap( ~Mass.Feature, scales = "free_y") +
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 10),
#         legend.position = "top",
#         strip.text = element_text(size = 10)) +
#   ggtitle("Internal Standard Raw Areas")
#
# currentDate <- Sys.Date()
# plotFileName <- paste("figures/IS.Raw.Areas_", currentDate, ".png", sep = "")
#
# ggsave(file = plotFileName, dpi = 600, width = 8, height = 6, units = "in")
# print(IS.Raw.Area.Plot)

# Edit data so names match, test that names are equal across sample sets---------------------------------
# Data.long  <- Data.NoIS %>%
#   rename(Mass.Feature = Metabolite.Name) %>%
#   select(Replicate.Name, Mass.Feature, Area.with.QC) %>%
#   arrange(Replicate.Name)
#
# test_isdata <- as.data.frame(sort(unique(Int.Stds.data$Replicate.Name)), stringsAsFactors = FALSE)
# test_long <- as.data.frame(sort(unique(Data.long$Replicate.Name)), stringsAsFactors = FALSE)
# test <- identical(test_isdata[[1]], test_long[[1]])
# print(paste("Your replicate names are identical:", test))
#
# if(test == FALSE)
#   stop("Error: Your replicate names are not matched across datasets!")

# Calculate mean values for each Internal Standard--------------------------------------------------------
# Int.Stds.means <- Int.Stds.data %>%
#   select(-c("a", "RunType", "c", "d")) %>%
#   group_by(Mass.Feature) %>%
#   summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE))

# # Normalize to each internal Standard--------------------------------------------------------------------
# Data.bound <- rbind(Int.Stds.data %>% select(-c("a", "RunType", "c", "d")), Data.long) %>%
#   arrange(Mass.Feature)
#
# Split_Dat <- list()
# # MIS stands for "Matched Internal Standard"
# for (i in 1:length(unique(Int.Stds.data$Mass.Feature))) {
#   Split_Dat[[i]] <- Data.bound %>%
#     mutate(MIS = unique(Int.Stds.data$Mass.Feature)[i]) %>%
#     left_join(Int.Stds.data %>%
#                 rename(MIS = Mass.Feature, IS_Area = Area.with.QC) %>%
#                 select(MIS, Replicate.Name, IS_Area), by = c("Replicate.Name", "MIS")) %>%
#     left_join(Int.Stds.means %>%
#                 rename(MIS = Mass.Feature), by = "MIS") %>%
#     mutate(Adjusted.Area = Area.with.QC/IS_Area*Average.Area)
# }
#
# Data.area.norm <- do.call(rbind, Split_Dat) %>%
#   select(-IS_Area, -Average.Area)
#
# # Standardize name structure to: Date_type_ID_replicate_anythingextra -----------------------------------
# Mydata.new <- Data.area.norm %>%
#   separate(Replicate.Name, c("runDate", "type", "SampID", "replicate"), "_") %>%
#   mutate(Run.Cmpd = paste(Data.area.norm$Replicate.Name, Data.area.norm$Mass.Feature))
#
#
# # Find the B-MIS for each MassFeature-------------------------------------------------------------------
#
# # Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
# # then choose which IS reduces the RSD the most (Poo.Picked.IS)
# Poodata1 <- Mydata.new %>%
#   filter(type == "Poo") %>%
#   group_by(SampID, Mass.Feature, MIS) %>%
#   mutate(RSD_ofPoo_IND = sd(Adjusted.Area, na.rm = TRUE) / mean(Adjusted.Area, na.rm = TRUE)) %>%
#   mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
#   group_by(Mass.Feature, MIS) %>%
#   summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
#   mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo))
#
#
# Poodata2 <- Poodata1 %>%
#   left_join(Poodata1 %>% group_by(Mass.Feature) %>%
#               summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))
#
#
# # Get the original RSD, calculate RSD change, decide if MIS is acceptable -------------------------------
# Poodata3 <- left_join(Poodata2, Poodata2 %>%
#                        filter(MIS == "Inj_vol" ) %>%
#                        mutate(Orig_RSD = RSD_ofPoo) %>%
#                        select(-RSD_ofPoo, -MIS)) %>%
#   mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
#   mutate(percent.Change = del_RSD/Orig_RSD) %>%
#   mutate(accept_MIS = (percent.Change > cut.off & Orig_RSD > cut.off2))
#
# # Change the BMIS to "Inj_vol" if the BMIS is not an acceptable -----------------------------------------
#
# # Adds a column that has the BMIS, not just Poo.Picked.IS
# # Changes the FinalBMIS to inject_volume if it's no good
# Fixed.poodata <- Poodata3 %>%
#   filter(MIS == Poo.Picked.IS) %>%
#   mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
#   mutate(FinalRSD = RSD_ofPoo)
#
# New.poodata <- Poodata3 %>%
#   left_join(Fixed.poodata %>% select(Mass.Feature, FinalBMIS)) %>%
#   filter(MIS == FinalBMIS) %>%
#   mutate(FinalRSD = RSD_ofPoo)
#
# Try <- New.poodata %>%
#   filter(FinalBMIS != "Inj_vol")
#
# QuickReport <- print(paste("Percent of Mass Features that picked a BMIS:",
#                            length(Try$Mass.Feature) / length(New.poodata$Mass.Feature), "|",
#                            "RSD improvement cutoff", cut.off, "|",
#                            "RSD minimum cutoff", cut.off2,
#                            sep = " "))
#
# reportFileName = paste("data_intermediate/MSDial_QuickReport", file.pattern, "_", currentDate, ".txt", sep = "")
# cat(QuickReport, file = reportFileName)
#
# # Evaluate and visualize the results of your BMIS cutoff-------------------------------------------------
# IS_toISdat <- Mydata.new %>%
#   filter(Mass.Feature %in% Int.Stds.data$Mass.Feature) %>%
#   select(Mass.Feature, MIS, Adjusted.Area, type) %>%
#   filter(type == "Smp") %>%
#   group_by(Mass.Feature, MIS) %>%
#   summarize(RSD_of_Smp = sd(Adjusted.Area, na.rm = TRUE) / mean(Adjusted.Area, na.rm = TRUE)) %>%
#   left_join(Poodata3 %>% select(Mass.Feature, MIS, RSD_ofPoo, accept_MIS))
#
# injectONlY_toPlot <- IS_toISdat %>%
#   filter(MIS == "Inj_vol")
#
# ISTest_plot <- ggplot() +
#   geom_point(dat = IS_toISdat, shape = 21, size = 2, aes(x = RSD_ofPoo, y = RSD_of_Smp, fill = accept_MIS)) +
#   geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_of_Smp), size = 3) +
#   facet_wrap(~ Mass.Feature) +
#   ggtitle(paste("Results of BMIS Cutoff:", cut.off, "RSD decrease,", cut.off2, "RSD minimum."))
#
# plotFileName <- paste("figures/BMIS_Evaluation_", currentDate, ".png", sep = "")
#
# ggsave(file = plotFileName, dpi = 600, width = 8, height = 6, units = "in")
# print(ISTest_plot)
#
#
# # Return data that is normalized via BMIS----------------------------------------------------------------
# BMIS.normalized.data <- New.poodata %>% select(Mass.Feature, FinalBMIS, Orig_RSD, FinalRSD) %>%
#   left_join(Mydata.new, by = "Mass.Feature") %>%
#   filter(MIS == FinalBMIS)
#
# currentDate <- Sys.Date()
# csvFileName <- paste("data_processed/MSDial_BMIS_Output_", Column.Type, "_", currentDate, ".csv", sep = "")
#
#
# write.csv(BMIS.normalized.data, csvFileName, row.names = FALSE)
#
# rm(list = setdiff(ls()[!ls() %in% c("file.pattern")], lsf.str()))
