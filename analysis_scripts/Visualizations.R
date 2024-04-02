# Visualization
#Claudia Luthy
# March 2024 for Syn salinity project

#inputs are final conc which has everything corrected in nm
#outputs are TBD
require(ggplot2)

abs_data <- read_csv("~/Desktop/syn_project_2024/data_raw/absorbance_measurements_syn.csv")
final_concs <- read_csv("~/Desktop/syn_project_2024/intermediates/final_concentrations.csv")
Ingalls <- read_csv("~/Desktop/syn_project_2024/data_raw/Ingalls_Lab_Standards.csv")


#trimming down the sample name to A, B, blank
final_concs$samp <- sapply(stringr::str_split(final_concs$samp, "_"), "[", 4)

#now normalizing to cell density absorbance measurements
normalized <- left_join(final_concs, abs_data) %>%
  mutate(norm_area = nm/avg_abs) %>%
  filter(!is.na(nm), !is.na(norm_area)) #this then removes rows with concs as NA

#this calculates means and std dev for replicates
concs_summary <- normalized %>% filter(!samp == 'Blk') %>%
  group_by(cmpd_name, sal) %>%
  summarise(avg = mean(norm_area), sd = sd(norm_area))

#this creates a table of blank values
blanks <- normalized %>% filter(samp == 'Blk') %>%
  select(cmpd_name, sal, norm_area) %>%
  rename(blk_area= 'norm_area')

#this merges the two above tables and subtracts the blank measurements
concs_summary <- left_join(concs_summary, blanks) %>%
  mutate(corrected_area = avg - blk_area)



#now for ANOVA time
anova_values <- left_join(normalized, blanks, relationship = "many-to-many") %>%
  mutate(subtracted_area = norm_area - blk_area) %>%
  filter(!stringr::str_detect(samp, "Blk")) %>%
  filter(subtracted_area > 0) #might not need this, calcs might be wrong

#run ANOVA's by compound, trying to see differences between salinities and then do post hocs
anova_results <- lapply(split(anova_values, anova_values$cmpd_name), function(i){
  anova(lm(subtracted_area ~ sal, data = i))
})

anova.table <- data.frame(compound=names(anova_results))

fun1 <- function(lst,n){
  sapply(lst, `[`, n)
}

anova.table$p.value <- fun1(anova_results,5) %>% fun1(.,1) %>% as.numeric
anova.table$f <- fun1(anova_results,4) %>% fun1(.,1) %>% as.numeric

#first visualization just to get an idea of what we are working with
concs_summary_1 <- concs_summary %>% filter(!is.na(corrected_area)) #this filters out everything with NA in the area column
p1 <- ggplot(concs_summary_1, aes(corrected_area, cmpd_name)) +
  geom_point() +
  scale_y_discrete() +
  scale_x_log10()
p1 <- p1 + facet_wrap(~sal, ncol=4)

ggsave(plot=p1, '~/Desktop/syn_concs_no.pdf', width=12, height=12, units='in')




#now going to start trying to do plots by compound type
#this first trims down the entire ingalls stds sheet to just name and type
ingalls_trimmed <- Ingalls %>% 
  select(Compound_Name, Compound_Type) %>%
  rename(cmpd_name = 'Compound_Name', cmpd_type = 'Compound_Type')
 
#now we are merging our final nice simple table with the trimmed stds
cmpd_by_type <- left_join(concs_summary_1, ingalls_trimmed, relationship = "many-to-many")


#betaine land
betaine <- cmpd_by_type %>%
  filter(cmpd_type %in% c("Betaine")) %>%
  filter(corrected_area > 0) %>%
  select(cmpd_name, sal, corrected_area)

betaine_plot <- ggplot(betaine, aes(fill=cmpd_name, y=corrected_area, x=sal))+ 
  geom_bar(position="stack", stat="identity")
ggsave(plot=betaine_plot, '~/Desktop/betaine_1.pdf', width=12, height=12, units='in')


#osmolyte land
osmolyte <- cmpd_by_type %>%
  filter(cmpd_type %in% c("Osmolyte")) %>%
  filter(corrected_area > 0) %>%
  select(cmpd_name, sal, corrected_area)

osmolyte_plot <- ggplot(osmolyte, aes(fill=cmpd_name, y=corrected_area, x=sal))+ 
  geom_bar(position="stack", stat="identity")
ggsave(plot=osmolyte_plot, '~/Desktop/syn_project_2024/visualizations/osmolyte_1.pdf', width=12, height=12, units='in')




#nucleic acid land
nucleic_acid <- cmpd_by_type %>%
  filter(cmpd_type %in% c("Nucleic acid")) %>%
  filter(corrected_area > 0) %>%
  select(cmpd_name, sal, corrected_area)

nucleic_acid_plot <- ggplot(nucleic_acid, aes(fill=cmpd_name, y=corrected_area, x=sal))+ 
  geom_bar(position="stack", stat="identity")
ggsave(plot=nucleic_acid_plot, '~/Desktop/syn_project_2024/visualizations/nucleic_acid_1.pdf', width=12, height=12, units='in')

