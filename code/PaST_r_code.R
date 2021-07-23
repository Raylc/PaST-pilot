##Libraries####
detachAllPackages <- function() {
      
      basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
      
      package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
      
      package.list <- setdiff(package.list,basic.packages)
      
      if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
      
}

detachAllPackages()
library(car)
library(ggplot2)
library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)
library(DescTools)
library(stringr)
library("factoextra")
library("psych") 
library(Hmisc)
library(glmulti)
library(report)
library(ggeffects) 
library(sjPlot)
library(rcompanion)
library(caret)
library(permuco)

#### Import datasets #### ####

experiment_data<-read.csv("/Users/JustinPargeter/Google Drive/Projects/PaST/Pilot experiment/Lithics/Experiment database/past_stats/past_stats/past_cores_result.csv") %>%
      rename(subject_id = Subject,
             core_number = nodule_number) %>%
      mutate(subject_id=as.factor(subject_id),
             core_number=as.factor(core_number),
             Condition=as.factor(Condition))

experiment_core_usedata<-read.csv("/Users/JustinPargeter/Google Drive/Projects/PaST/Pilot experiment/Lithics/Experiment database/past_stats/past_stats/core_use_records.csv") %>%
      mutate(subject=as.factor(subject),
             core_number=as.factor(core_number))

lithic_data<-read.csv("past_pilot_lithics.csv") %>%
      rename(flake_mass = flake_mass,
             flake_max_thickness = max_thickness,
             core_final_mass=core_mass) %>%
      mutate(core_number_nyu=core_number) %>%
      rename(subject_id_nyu=subject_id) 

extra_thickness_flakes<-read.csv("extra_flakes.csv") 
  
pre_exp_cores_data<-read.csv("/Users/JustinPargeter/Google Drive/Projects/PaST/Pilot experiment/Lithics/Experiment database/past_stats/past_stats/pre_exp_cores_data.csv") %>%
      rename(core_number = nodule_number) %>%
      mutate(core_number=as.factor(core_number))

incremental_width<-list.files("/Users/JustinPargeter/Google Drive/Projects/PaST/Pilot experiment/Lithics/Experiment database/past_stats/past_stats/combined_txt_files", pattern="*.txt")
area_size<- read.csv(file = "/Users/JustinPargeter/Google Drive/Projects/PaST/Pilot experiment/Lithics/Experiment database/past_stats/past_stats/area_size_combined_cheng_aditi_megan.csv")

psychometrics<-read.csv(file = "/Users/JustinPargeter/Google Drive/Projects/PaST/Pilot experiment/Lithics/Experiment database/past_stats/past_stats/past_psychometrics.csv") %>%
      rename(subject = Subject) %>%
      mutate(subject=as.factor(subject))

#### create merged dataset #### ####

# txt file opening function
open_read<-function(filename, foldername){
      open_file<-read.delim(paste(foldername,"/",filename,sep=""),header=F,sep='_')
      colnames(open_file)<-c("variable","measurement_point","measurement")
      open_file$individual<-gsub(".txt","",filename)
      return(open_file)
}

# execute the txt file opening function
flake_plan_measurements<-do.call(rbind,lapply(incremental_width,open_read,"/Users/JustinPargeter/Google Drive/Projects/PaST/Pilot experiment/Lithics/Experiment database/past_stats/past_stats/combined_txt_files"))

# reshape the herbi measurement data
cols <- c("subject","core_number","flake_number")

flake_plan_measurements_reshape<-flake_plan_measurements %>%
      mutate(variable=recode(flake_plan_measurements$variable, Width = "width"),
             measurement_point = paste(variable, measurement_point, sep="_"),
             shape_width_mm=measurement*10) %>% #convert cm to mm
      separate(individual, into = c("delete", "subject","core_number","delete_2","flake_number"), sep="_") %>%
      dplyr::select(-c(delete,delete_2,variable,measurement)) %>%
      dcast(... ~ measurement_point,value.var="shape_width_mm") %>%
      mutate_at(cols, funs(factor(.)))

# reformat area size data
area_size_data<-area_size %>%
      mutate(area=area*100,
             max_width=max_width*10,
             max_length=max_length*10) %>% #convert cm2 to mm2
      separate(individual, into = c("delete", "subject","core_number","delete_2","flake_number"), sep="_") %>%
      dplyr::select(-c(delete,delete_2)) %>%
      mutate_at(cols, funs(factor(.)))

#find duplicates if necessary (if dcast is returning 1's instead of measurements)
duplicates<-flake_plan_measurements_reshape[duplicated(flake_plan_measurements_reshape[,1:4]),]

unique_data <- unique(flake_plan_measurements_reshape)

## join the two data sets
merged_shape_data<- full_join(area_size_data,flake_plan_measurements_reshape, by=c("subject","core_number","flake_number")) %>%
      rename(max_flake_length=max_length) 

# remove leading zeros
merged_shape_data$subject<-str_remove(merged_shape_data$subject, "^0+")
merged_shape_data$core_number<-str_remove(merged_shape_data$core_number, "^0+")
merged_shape_data$core_number<-str_remove(merged_shape_data$core_number, "^0+")
merged_shape_data$flake_number<-str_remove(merged_shape_data$flake_number, "^0+")

## Fill in missing expert cores data
#create unique ID column for demos to extract relevant core information
merged_shape_data_demo_identifier<-merged_shape_data %>% #a bunch of data missing for the expert cores 236-246 and 35
      mutate(subject=as.factor(subject)) %>%
      distinct(core_number,.keep_all = TRUE) %>%
      select(c(subject,core_number))

pre_exp_cores_data_demo_ID<-merge(pre_exp_cores_data,merged_shape_data_demo_identifier,
                                  by="core_number",
                                  all=T)

experiment_data_filled<-merge(x = experiment_data, y = pre_exp_cores_data_demo_ID 
                              [ , c("subject","core_number","mass")], 
                              by = "core_number", all=TRUE) %>%
      select(-("start_mass")) %>%
      mutate(subject_id = coalesce(subject, subject_id)) %>%
      rename(start_mass=mass) %>%
      filter(!core_number %in% c(2,22,42,149,168, 169, 174, 204, 207, 223, 228)) %>%
      subset(!core_number == "sample") %>% #no starting mass-flakes shown to subjects prior to the experiment
      select(-("subject"))

# merge with NYU lithic data
complete_past_data<-full_join(experiment_data_filled,lithic_data,by=c("core_number")) %>%
      subset(Used==1 | subject_id == "demo") %>% #only cores actually used in the experiment
      select(c("subject_id", "subject_id_nyu","recorder_name","Condition", "flake_number","core_number","core_number_nyu",
               "Total_Cores_Used","Order","Used",
               "maxlength","maxwidth","maxthickness","edge_angle_degree",
               "start_mass","core_final_mass","Proportion_Mass_Removed",
               "flake_mass","flake_max_thickness","flake_count",
               "Highest_nBack","RPM","BEAST","Grip_Strength")) %>%
      arrange(subject_id,core_number) %>%
      mutate(delta_mass=(start_mass-core_final_mass)/start_mass) %>%
      group_by(subject_id) %>%
      mutate(core_count=n_distinct(core_number),
             core_counts_compared=ifelse(core_count==Total_Cores_Used,1,0),
             subject_id=as.factor(subject_id),
             core_number=as.factor(core_number),
             flake_number=as.factor(flake_number)) %>%
      rename(subject=subject_id) %>%
      ungroup() %>%
      select(-c("recorder_name","subject_id_nyu","core_number_nyu","core_count","core_counts_compared","Proportion_Mass_Removed"))

# Check missing data
new_DF <- complete_past_data[is.na(complete_past_data$subject),]
new_DF <- complete_past_data[is.na(complete_past_data$flake_number),] #core 210 missing flakes from the NYU database

# merge with flake size and shape data
#issue is with the merging column names (subject, core number, flake number)
complete_past_data_plusflakes<-full_join(complete_past_data,merged_shape_data,by=c("subject","core_number","flake_number"))

# Check missing data (all flakes below herbi measurements just not working)
# S4 core 39 flake 2 
# s4 29 14, 15
# S4 core 45 flake 6 herbi measurements just not working
# S9 core 129 flake 7 herbi measurements just not working
# S9 core 129 flake 8
# S9 core 129 flake 18, 19, 21, 25, 27
# S20 212 002, 006, 014, 20
# S21 210 001, 007, 011, 012, 013, 015
new_DF <- complete_past_data_plusflakes[is.na(complete_past_data_plusflakes$delta_mass),] %>%
      distinct(core_number,.keep_all = TRUE)

new_DF <- complete_past_data_plusflakes[!complete.cases(complete_past_data_plusflakes),] %>%
      distinct(core_number,.keep_all = TRUE)

# Check against list of cores knapped during the experiment
test_core_comparison<-merge(x = complete_past_data_plusflakes[ , c("subject","core_number")], 
                            y = experiment_core_usedata[ , c("subject","core_number")], 
                            by = "core_number", all=TRUE) %>%
      rename(merged_subject=subject.x,
             exp_notes_subject=subject.y) %>%
      distinct(core_number,.keep_all = TRUE) %>%
      mutate(subject_id_compared=ifelse(exp_notes_subject==merged_subject,1,0))

#### Flake shape PCA ####

# adjust linear measurements to ensure top is pointed, base is wider
complete_past_data_plusflakes = complete_past_data_plusflakes %>% 
      mutate(index = as.numeric(row.names(complete_past_data_plusflakes)))
complete_past_data_plusflakes_diff <- complete_past_data_plusflakes %>% filter(width_0.1 > width_0.9 &
                                                                                     width_0.2 > width_0.8 &
                                                                                     width_0.3 > width_0.7)
complete_past_data_plusflakes_no_diff = complete_past_data_plusflakes %>% 
      filter(!(index %in% complete_past_data_plusflakes_diff$index)) %>%
      rename_with( tolower)
complete_past_data_plusflakes_updated <- complete_past_data_plusflakes_diff %>% mutate(
      Width_0.1 = width_0.9,
      Width_0.2 = width_0.8,
      Width_0.3 = width_0.7,
      Width_0.9 = width_0.1,
      Width_0.8 = width_0.2,
      Width_0.7 = width_0.3
) %>%
      select(-c(width_0.9,width_0.8,width_0.7,width_0.1,width_0.2,width_0.3)) %>%
      rename_with( tolower)
complete_past_data_plusflakes_updated = union_all(complete_past_data_plusflakes_no_diff,
                                                  complete_past_data_plusflakes_updated)

complete_past_data_plusflakes_updated = complete_past_data_plusflakes_updated %>% select(-index)

# subset without demo

pca_flakes<- complete_past_data_plusflakes_updated[!is.na(complete_past_data_plusflakes_updated$width_0.1), ] %>%
      subset(!subject %in% c("demo")) 

# determine factorability
library(psych)
bartlett.test(pca_flakes[,c(15,24,26:34)])

# create identifiers datasets for rejoining

subject_core_flakes <- pca_flakes[, 1:4] 

# Calculate PCA

pca_flakes_complete<-prcomp(pca_flakes[,c(15,24,26:34)], scale. = T)

# Biplot of variables and factors

fviz_pca_var(pca_flakes_complete,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Eigenvalues

eig.val_flakes <- get_eigenvalue(pca_flakes_complete)

# Results for Variables

res.var_flakes<- get_pca_var(pca_flakes_complete)
data.frame(res.var_flakes$coord) %>% select(c("Dim.1","Dim.2","Dim.3"))

# results for individuals
res.ind_flakes <- get_pca_ind(pca_flakes_complete)
indiv_flake_scores<-data.frame(res.ind_flakes$coord)        # Coordinates

# add data back to dataset 
indiv_subject_scores<-data.frame(res.ind_flakes$coord) %>%
      select(c("Dim.1","Dim.2","Dim.3")) %>%
      dplyr::rename(flake_pca_dim1 = Dim.1,
                    flake_pca_dim2 = Dim.2,
                    flake_pca_dim3= Dim.3)

joined_to_identifiers<-cbind(indiv_subject_scores,subject_core_flakes) %>%
      dplyr::rename(Condition = condition)

complete_past_data_shape_pca<-join(complete_past_data_plusflakes,
                                   joined_to_identifiers,
                                by=c("subject","Condition","core_number","flake_number")) %>%
  mutate(flake_pca_dim1=flake_pca_dim1*-1) #change sign so larger values = larger flakes

# plot subjects' data against absolute size measure to track allometry in variance. 

ggplot(complete_past_data_shape_pca, aes(y = log(flake_mass), x = flake_pca_dim2)) +
      geom_point(alpha = 0.5) +
      #geom_smooth(method='loess') +
      theme_classic() +
      ylab("Flake mass (log)") +
      xlab("Flake shape Factor 2 (thickness relative to width)")+
  theme(text = element_text(size=20))

ggplot(complete_past_data_shape_pca, aes(y = log(flake_mass), x = flake_pca_dim3)) +
      geom_point(alpha = 0.5) +
      #geom_smooth(method='lm') +
      theme_classic()+
      ylab("Flake mass (log)") +
      xlab("Flake shape Factor 3 (basal relative to tip width, thickness)")+
  theme(text = element_text(size=18))

ggplot(complete_past_data_shape_pca, aes(x = Condition, y = flake_pca_dim3)) +
  geom_boxplot()

ggplot(complete_past_data_shape_pca, aes(x = Condition, y = flake_pca_dim2)) +
  geom_boxplot()

## apply to Dietz's flake data

pca_flakes_expert<- complete_past_data_plusflakes_updated[!is.na(complete_past_data_plusflakes_updated$width_0.1), ] %>%
      subset(subject %in% c("demo"))

# create identifiers datasets for rejoining

subject_core_flakes_expert <- pca_flakes_expert[, c(1,3:4)]
subject_core_flakes_expert$FIPS<-rownames(subject_core_flakes_expert)

# PCA

dietz_data<-data.frame(scale(pca_flakes_expert[c(15,24,26:34)], 
                             pca_flakes_complete$center, 
                             pca_flakes_complete$scale) %*% pca_flakes_complete$rotation) %>%
      tibble::rowid_to_column( "FIPS") %>%
      select(c("FIPS","PC1","PC2","PC3")) %>%
      mutate(PC1=PC1*-1) #change sign so larger values = larger flakes

# join back

expert_joined_to_identifiers<-merge(subject_core_flakes_expert,dietz_data, by="FIPS") %>%
      select(-FIPS)

# create expert included dataset for use in group comparisons
complete_past_data_shape_pca_expert<-join(complete_past_data_shape_pca,
                                          expert_joined_to_identifiers,
                                          by=c("subject","core_number","flake_number")) %>%
      unite('flake_pca_dim1',PC1,flake_pca_dim1) %>%
      unite('flake_pca_dim2',PC2,flake_pca_dim2) %>%
      unite('flake_pca_dim3',PC3,flake_pca_dim3) %>%
      mutate(flake_pca_dim1=gsub("NA_", "", flake_pca_dim1),
             flake_pca_dim1=gsub("_NA", "", flake_pca_dim1),
             flake_pca_dim2=gsub("_NA", "", flake_pca_dim2),
             flake_pca_dim2=gsub("NA_", "", flake_pca_dim2),
             flake_pca_dim3=gsub("_NA", "", flake_pca_dim3),
             flake_pca_dim3=gsub("NA_", "", flake_pca_dim3),
             flake_pca_dim1=as.numeric(flake_pca_dim1),
             flake_pca_dim2=as.numeric(flake_pca_dim2),
             flake_pca_dim3=as.numeric(flake_pca_dim3))

#### Skill metric PCA ####

# Calculate flake skill summary variables for skill PC 

large_flake_count<-complete_past_data_shape_pca_expert %>% group_by(core_number) %>% tally()
complete_past_data_shape_pca_expert<-join(complete_past_data_shape_pca_expert,large_flake_count,by="core_number")
names(complete_past_data_shape_pca_expert)[names(complete_past_data_shape_pca_expert) == "n"] <- "total_flakes_above40mm_5g_bycore"

complete_past_data_shape_pca_expert_summary<- complete_past_data_shape_pca_expert %>%
      rename(maxlength_core = maxlength,
             maxwidth_core = maxwidth,
             maxthickness_core=maxthickness,
             core_start_mass=start_mass,
             total_flake_count_over_10mm=flake_count) %>%
      mutate(flaked_mass=core_start_mass-core_final_mass,
             flake_area_mass=area/flake_mass) %>%
      group_by(core_number) %>%
      mutate(total_mass_flakes_percore_greater40mm_5g=sum(flake_mass),
             flake_area_mass_bycore=mean(flake_area_mass)) %>%
      filter(!core_number %in% c("39") & !flake_number %in% c("2")) %>% #flakes with missing area measurements
      filter(!core_number %in% c("45") & !flake_number %in% c("6")) %>%
      mutate(total_mass_allflakes_greater40mm_5g=sum(flake_mass),
             flake_mass_bycore=sum(flake_mass),
             total_flaked_mass_bycore=sum(unique(flaked_mass)),
             mass_flakes_flaked_mass_bycore=(total_mass_allflakes_greater40mm_5g/total_flaked_mass_bycore),
             flake_area_mass_bycore=mean(flake_area_mass),
             mass_removed_bycore=core_start_mass-core_final_mass,
             flake_pca_dim1_bycore=mean(flake_pca_dim1, na.rm=T),
             flake_pca_dim2_bycore=mean(flake_pca_dim2,na.rm=T),
             flake_pca_dim3_bycore=mean(flake_pca_dim3,na.rm=T)) %>%
      ungroup() %>%
      group_by(subject) %>%
      mutate(total_mass_allflakes_greater40mm_5g_bysubject=sum(flake_mass),
             total_flaked_mass_bysubject=sum(unique(flaked_mass)),
             mass_flakes_flaked_mass_bysubject=total_mass_allflakes_greater40mm_5g,
             flake_area_mass_bysubject=mean(flake_area_mass),
             delta_mass_bysubject=mean(delta_mass))

# Calculate PCA

#### check possible errors in s15-clear that core 113 is an extreme outlier ####

ggplot(complete_past_data_shape_pca_expert_summary, aes(x=Condition, y=mass_flakes_flaked_mass_bycore))+
  geom_boxplot()+
  geom_point(data=subset(complete_past_data_shape_pca_expert_summary, subject==15), 
             aes(x=Condition,y=mass_flakes_flaked_mass_bycore), 
             color='red',
             size=3)

ggplot(complete_past_data_shape_pca_expert_summary, aes(x=Condition, y=total_flakes_above40mm_5g_bycore))+
  geom_boxplot()+
  geom_point(data=subset(complete_past_data_shape_pca_expert_summary, subject==15), 
             aes(x=Condition,y=total_flakes_above40mm_5g_bycore), 
             color='red',
             size=3)

ggplot(complete_past_data_shape_pca_expert_summary, aes(x=Condition, y=delta_mass))+
  geom_boxplot()+
  geom_point(data=subset(complete_past_data_shape_pca_expert_summary, subject==15), 
             aes(x=Condition,y=delta_mass), 
             color='red',
             size=3)

ggplot(complete_past_data_shape_pca_expert_summary, aes(x=Condition, y=flake_pca_dim1_bycore))+
  geom_boxplot()+
  geom_point(data=subset(complete_past_data_shape_pca_expert_summary, subject==15), 
             aes(x=Condition,y=flake_pca_dim1_bycore), 
             color='red',
             size=3)

ggplot(complete_past_data_shape_pca_expert_summary, aes(x=Condition, y=flake_pca_dim2_bycore))+
  geom_boxplot()+
  geom_point(data=subset(complete_past_data_shape_pca_expert_summary, subject==15), 
             aes(x=Condition,y=flake_pca_dim2_bycore), 
             color='red',
             size=3)

ggplot(complete_past_data_shape_pca_expert_summary, aes(x=Condition, y=flake_pca_dim3_bycore))+
  geom_boxplot()+
  geom_point(data=subset(complete_past_data_shape_pca_expert_summary, subject==15), 
             aes(x=Condition,y=flake_pca_dim3_bycore), 
             color='red',
             size=3)

ggplot(complete_past_data_shape_pca_expert_summary, aes(x=Condition, y=Total_Cores_Used))+
  geom_boxplot()+
  geom_point(data=subset(complete_past_data_shape_pca_expert_summary, subject==15), 
             aes(x=Condition,y=Total_Cores_Used), 
             color='red',
             size=3)


### Create PCA dataset ####

subject_flake_core_skill_pca <- complete_past_data_shape_pca_expert_summary %>% subset(!subject=="demo") %>% 
      ungroup() %>%
      select(core_number) %>%
      distinct(core_number,.keep_all=T)

pca_skill_data_bycore<-complete_past_data_shape_pca_expert_summary %>% 
      subset(!subject=="demo" & !core_number==113) %>% #removing core 113 
      ungroup() %>%
      select(c("core_number","mass_flakes_flaked_mass_bycore","total_flakes_above40mm_5g_bycore",
                "delta_mass","flake_pca_dim1_bycore","flake_pca_dim2_bycore","flake_pca_dim3_bycore",
                  "Total_Cores_Used"))  %>%
      dplyr::rename(Mass_flakes_flaked_mass = mass_flakes_flaked_mass_bycore,
                    large_flake_count=total_flakes_above40mm_5g_bycore,
                    Delta_mass = delta_mass,
                    Shape_PC2_relative_thickness = flake_pca_dim2_bycore,
                    Shape_PC1_size = flake_pca_dim1_bycore,
                    Shape_PC3_thick_basal_tip_width = flake_pca_dim3_bycore,
                    Total_cores_used = Total_Cores_Used) %>%
      distinct(core_number,.keep_all=T) %>%
      select(-core_number)

pca_skill_metrics<-prcomp(na.omit(pca_skill_data_bycore), scale. = T)

# determine factorability
library(psych)
bartlett.test(pca_skill_data_bycore)

# Graph of variables
library(factoextra)
fviz_pca_var(pca_skill_metrics,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

#screeplot
screeplot(pca_skill_metrics,type='line')

# Eigenvalues
eig.val <- get_eigenvalue(pca_skill_metrics)

# Results for Variables
res.var_skill_metrics<- get_pca_var(pca_skill_metrics)
data.frame(res.var_skill_metrics$coord)         # Coordinates

# results for individual cores
res.ind_skill_metrics <- get_pca_ind(pca_skill_metrics)
indiv_skill_scores<-data.frame(res.ind_skill_metrics$coord)  

## add data back to dataset 
indiv_core_skill_scores<-indiv_skill_scores %>%
      tibble::rowid_to_column( "FIPS") %>%
      select(c("FIPS","Dim.1","Dim.2")) %>%
      dplyr::rename(skill_pca_dim1_quantity_flaking = Dim.1,
                    skill_pca_dim2_quality_flaking = Dim.2)

subject_flake_core_skill_pca$FIPS<-rownames(subject_flake_core_skill_pca)

skill_joined_to_identifiers<-join(indiv_core_skill_scores,subject_flake_core_skill_pca, by="FIPS") %>%
      select(-FIPS)

complete_past_data_shape_skill_pca_combined<-join(complete_past_data_shape_pca_expert_summary,
                                skill_joined_to_identifiers,
                                by=c("core_number")) %>%
      mutate(skill_pca_dim2_quality_flaking=skill_pca_dim2_quality_flaking*-1) %>% #transform sign so higher = better
      distinct(core_number,.keep_all=T)

#Are PC 1 and 2 related?

subject_skillPC_data<-complete_past_data_shape_skill_pca_combined %>%
  group_by(subject) %>%
  mutate(skill_pca_dim1_quantity_flaking=mean(skill_pca_dim1_quantity_flaking),
         skill_pca_dim2_quality_flaking=mean(skill_pca_dim2_quality_flaking)) %>%
  distinct(subject,.keep_all=T) %>%
  select(c(subject, Condition, core_number, flake_number, skill_pca_dim1_quantity_flaking, skill_pca_dim2_quality_flaking)) %>%
  na.omit(skill_pca_dim1_quantity_flaking)

ggplot(subset(subject_skillPC_data), 
       aes(x=skill_pca_dim1_quantity_flaking, y=skill_pca_dim2_quality_flaking,color=Condition)) +
      geom_point() +
      geom_smooth(method='lm', se=T)+ 
      scale_color_manual(name = "Training condition", labels = c("Untrained", "Trained"), values = c("blue", "red"))+
      theme_classic() +
      ylab("Factor 1 'Quality flaking' (higher = higher quality flaking)") +
      xlab("Factor 2 'Quantity flaking' (higher = higher quantity flaking)")+
      theme(text = element_text(size=16))

summary(lm(skill_pca_dim1_quantity_flaking~skill_pca_dim2_quality_flaking, data=subset(subject_skillPC_data,Condition==1 &skill_pca_dim2_quality_flaking>-1)))
summary(lm(skill_pca_dim1_quantity_flaking~skill_pca_dim2_quality_flaking, data=subset(subject_skillPC_data,Condition==0 &skill_pca_dim2_quality_flaking>-1)))

#### Group level comparisons ####
## Two shape PCs, Core comparisons, Skill PC comparisons

complete_past_data_shape_skill_pca_combined_bycore<-complete_past_data_shape_skill_pca_combined %>%
      mutate(condition_threeway = case_when(subject == 'demo' ~ "expert",
                                            Condition == '0' ~ "untrained",
                                            Condition == '1' ~ "trained",
                                          TRUE ~ "NA"),
             core_shape=maxlength_core*maxwidth_core*maxthickness_core)

complete_past_data_shape_skill_pca_combined_byflake<-complete_past_data_shape_pca_expert %>%
      mutate(condition_threeway = case_when(subject == 'demo' ~ "expert",
                                            Condition == '0' ~ "untrained",
                                            Condition == '1' ~ "trained",
                                            TRUE ~ "NA"))

# Skill PC comparisons between novice groups
ggplot(filter(model_data), aes(x=Condition, y=skill_pca_dim1_quantity_flaking)) +
  geom_boxplot()+
  geom_jitter()

leveneTest(skill_pca_dim1_quantity_flaking ~ Condition, data = subject_skillPC_data)
aovperm(skill_pca_dim1_quantity_flaking ~ Condition,
        data = subject_skillPC_data)

ggplot(filter(subject_skillPC_data,!subject=="demo"), 
       aes(x=Condition, y=skill_pca_dim2_quality_flaking)) +
  geom_boxplot()+
  geom_jitter()

leveneTest(skill_pca_dim2_quality_flaking ~ Condition, data = filter(complete_past_data_shape_skill_pca_combined_bycore,!subject=="demo" ))
aovperm(skill_pca_dim2_quality_flaking ~ Condition,
        data = subject_skillPC_data)

# Flake shape
ggplot(complete_past_data_shape_skill_pca_combined_byflake, aes(x=condition_threeway, y=flake_pca_dim1))+
  geom_boxplot()

leveneTest(flake_pca_dim1 ~ condition_threeway, data = complete_past_data_shape_skill_pca_combined_byflake) #check assumption of equal variance
aov_1<-aov(flake_pca_dim1~condition_threeway, data=complete_past_data_shape_skill_pca_combined_byflake)
aov_1 %>% report() %>% as.report_table(summary=TRUE)

ggplot(complete_past_data_shape_skill_pca_combined_byflake, aes(x=condition_threeway, y=flake_pca_dim2))+
  geom_boxplot()
leveneTest(flake_pca_dim2 ~ condition_threeway, data = complete_past_data_shape_skill_pca_combined_byflake)
aov_2<-aov(flake_pca_dim2~condition_threeway, data=complete_past_data_shape_skill_pca_combined_byflake)
aov_2 %>% report() %>% as.report_table(summary=TRUE)

ggplot(complete_past_data_shape_skill_pca_combined_byflake, aes(x=condition_threeway, y=flake_pca_dim3))+
  geom_boxplot()
leveneTest(flake_pca_dim3 ~ condition_threeway, data = complete_past_data_shape_skill_pca_combined_byflake)
aovperm(flake_pca_dim3 ~ condition_threeway, #implement permutation tests for un-equal variances
        data = complete_past_data_shape_skill_pca_combined_byflake)

# Performance skill, compare by attribute to allow for comparisons between expert and novice, and between PC scores for novice groups

ggplot(complete_past_data_shape_skill_pca_combined_bycore, aes(x=condition_threeway, y=total_flakes_above40mm_5g_bycore)) +
  geom_boxplot()+
  geom_jitter()+
  theme_classic() +
  xlab("") +
  ylab("Large (>40mm and 5g) flake count")+
  theme(text = element_text(size=20))

leveneTest(total_flakes_above40mm_5g_bycore ~ condition_threeway, data = complete_past_data_shape_skill_pca_combined_bycore)
aovperm(total_flakes_above40mm_5g_bycore ~ condition_threeway,
        data = complete_past_data_shape_skill_pca_combined_bycore)

aovperm(total_flakes_above40mm_5g_bycore ~ Condition,
        data = complete_past_data_shape_skill_pca_combined_bycore)

epsilonSquared(x = complete_past_data_shape_skill_pca_combined_bycore$total_flakes_above40mm_5g_bycore,
               g = complete_past_data_shape_skill_pca_combined_bycore$condition_threeway)

ggplot(complete_past_data_shape_skill_pca_combined_bycore, aes(x=condition_threeway, y=mass_flakes_flaked_mass_bycore)) +
      geom_boxplot() +
      geom_jitter()
leveneTest(mass_flakes_flaked_mass_bycore ~ condition_threeway, data = complete_past_data_shape_skill_pca_combined_bycore) #test for homogeneity of variances
aov_3<-aov(mass_flakes_flaked_mass_bycore~condition_threeway, data=complete_past_data_shape_skill_pca_combined_bycore)
aov_3 %>% report() %>% as.report_table(summary=TRUE)

ggplot(complete_past_data_shape_skill_pca_combined_bycore, aes(x=condition_threeway, y=delta_mass)) +
      geom_boxplot()+
      geom_jitter()+
  theme_classic() +
  xlab("") +
  ylab("Core delta mass")+
  theme(text = element_text(size=20))

leveneTest(delta_mass ~ condition_threeway, data = complete_past_data_shape_skill_pca_combined_bycore)
aovperm(delta_mass ~ condition_threeway,
        data = complete_past_data_shape_skill_pca_combined_bycore)

aovperm(delta_mass ~ Condition,
        data = complete_past_data_shape_skill_pca_combined_bycore)

epsilonSquared(x = complete_past_data_shape_skill_pca_combined_bycore$delta_mass,
               g = complete_past_data_shape_skill_pca_combined_bycore$condition_threeway)


# Core final mass comparisons
ggplot(complete_past_data_shape_skill_pca_combined_bycore, aes(x=condition_threeway, y=log(core_final_mass))) +
      geom_boxplot()+
  geom_jitter()+
  theme_classic() +
  xlab("") +
  ylab("Core final mass (g)")+
  theme(text = element_text(size=20))

leveneTest(core_final_mass ~ condition_threeway, data = filter(complete_past_data_shape_skill_pca_combined_bycore,!subject=="demo"))
aovperm(core_final_mass ~ condition_threeway,
        data = complete_past_data_shape_skill_pca_combined_bycore)

aovperm(core_final_mass ~ Condition,
        data = complete_past_data_shape_skill_pca_combined_bycore)

epsilonSquared(x = complete_past_data_shape_skill_pca_combined_bycore$core_final_mass,
               g = complete_past_data_shape_skill_pca_combined_bycore$condition_threeway)

#### learning rates ####

# creat grouping variable for learning rates

learning_rate_data<-complete_past_data_shape_skill_pca_combined_bycore %>%
      filter(!subject=="demo") %>%
      distinct(core_number,.keep_all=T) %>%
      group_by(subject,Order) %>%
      mutate(adjusted_cores_used=(Order/Total_Cores_Used)*100) %>%
      mutate(core_order_stage = case_when(adjusted_cores_used > 0 & adjusted_cores_used < 20 ~ "0_20",
                                          adjusted_cores_used > 20 & adjusted_cores_used < 40 ~ "20_40",
                                          adjusted_cores_used > 40 & adjusted_cores_used < 60 ~ "40_60",
                                          adjusted_cores_used > 60 & adjusted_cores_used < 80 ~ "60_80",
                                          adjusted_cores_used > 80 & adjusted_cores_used < 100 ~ "80-100",
                                          TRUE ~ "80-100"),
             Condition=dplyr::recode(Condition, '0' = "Untrained", '1'="Trained")) %>%
      ungroup() %>%
      select(c(subject, Condition, core_number, adjusted_cores_used,core_order_stage,
               flake_pca_dim1,flake_pca_dim2,flake_pca_dim3,
               skill_pca_dim1_quantity_flaking,
               skill_pca_dim2_quality_flaking,
               core_start_mass))

# skill PCS, core size

aov(flake_pca_dim1~Condition*core_order_stage, data=learning_rate_data) %>% report() %>% as.report_table(summary=TRUE)

aov(flake_pca_dim2~Condition*core_order_stage, data=learning_rate_data) %>% report() %>% as.report_table(summary=TRUE)

aov(flake_pca_dim3~Condition*core_order_stage, data=learning_rate_data) %>% report() %>% as.report_table(summary=TRUE)
ggplot(subset(learning_rate_data,flake_pca_dim3 <3.5 ),aes(x=core_order_stage,y=flake_pca_dim3))+
  geom_boxplot()+
  theme_classic() +
  xlab("Relative core order category (%)") +
  ylab("Flake factor 3 (basal relative to tip width)")+
  theme(text = element_text(size=20))

aov(skill_pca_dim1_quantity_flaking~Condition*core_order_stage, data=learning_rate_data) %>% report() %>% as.report_table(summary=TRUE)

aov(skill_pca_dim2_quality_flaking~Condition*core_order_stage, data=learning_rate_data) %>% report() %>% as.report_table(summary=TRUE)

aov(log(core_start_mass)~Condition*core_order_stage, data=learning_rate_data) %>% report() %>% as.report_table(summary=TRUE)

ggplot(learning_rate_data,aes(x=core_order_stage,y=log(core_start_mass), color=Condition))+
      geom_boxplot()+
      theme_classic() +
      xlab("Relative core order category (%)") +
      ylab("Nodule starting mass (log)")+
      theme(text = element_text(size=20))

ggplot(learning_rate_data,aes(x=adjusted_cores_used,y=skill_pca_dim1_quantity_flaking, color=Condition))+
  geom_point()+
  facet_wrap(~subject)+
  xlab("Relative core order %") +
  ylab("Skill PC 1 Quantity flaking")+
  geom_smooth(method="loess",se=F)

ggplot(learning_rate_data,aes(x=adjusted_cores_used,y=skill_pca_dim2_quality_flaking, color=Condition))+
  geom_point()+
  facet_wrap(~subject)+
  xlab("Relative core order %") +
  ylab("Skill PC 2 Quality flaking")+
  geom_smooth(method="loess",se=F)

#### Performance and individual variation models ####

## Add psychometric data

complete_past_data_psycho<-full_join(complete_past_data_shape_skill_pca_combined,psychometrics,by=c("subject")) %>%
      subset(!subject %in% c("6","9", "demo")) #missing psych measures, maybe also 21 as outlier

## Generate model data

model_data<-complete_past_data_psycho %>%
      drop_na(skill_pca_dim1_quantity_flaking) %>%
      mutate(core_shape=maxlength_core*maxwidth_core*maxthickness_core) %>%
      select(c(subject,Condition,core_number,skill_pca_dim1_quantity_flaking,
               skill_pca_dim2_quality_flaking,core_shape,
               Highest_nBack,RPM,BEAST,
               Grip_Strength, fitts_mtavg_ms,
               core_start_mass)) %>%
      group_by(subject) %>%
      mutate(skill_pca_dim1_quantity_flaking=mean(skill_pca_dim1_quantity_flaking),
             skill_pca_dim2_quality_flaking=mean(skill_pca_dim2_quality_flaking)) %>%
      ungroup() %>%
      mutate(Highest_nBack=as.numeric(scale(Highest_nBack,center = TRUE)),
             RPM=as.numeric(scale(RPM,center = TRUE)),
             BEAST=as.numeric(scale(BEAST,center = TRUE)),
             core_start_mass=as.numeric(scale(core_start_mass,center = TRUE)),
             Grip_Strength=as.numeric(scale(Grip_Strength,center = TRUE)),
             fitts_mtavg_ms=as.numeric(scale(fitts_mtavg_ms,center = TRUE))) %>%
      distinct(subject, .keep_all=T) 
  
### Skill PC 1 ####

#generate complete model

PC1_complete<- lm(skill_pca_dim1_quantity_flaking ~  Condition + Highest_nBack +
                            RPM + BEAST  + 
                            fitts_mtavg_ms + 
                            Grip_Strength,
                      data = model_data)

#the multimodal selection

PC1_multi <- glmulti(PC1_complete, # use the model with built as a starting point
                   level = 2, #1=just look at main effects, 2=look at main and interraction effects
                   method='g',
                   crit="aic",minsize = 0, maxsize = 7, minK = 0, maxK = -1) 

print(PC1_multi) #all models

summary(PC1_multi@objects[[1]]) #best model

# add mass to final model, but not in model selection
PC1_best_model<-lm(skill_pca_dim1_quantity_flaking~1+
                     Condition + 
                     Highest_nBack + 
                     BEAST+
                     Grip_Strength + 
                     fitts_mtavg_ms:RPM +
                     Condition:Highest_nBack +
                     Condition:BEAST +
                     core_start_mass,
                   data = model_data)

summary(PC1_best_model)

PC1_best_model %>%
      report() %>%
      as.report_table(summary=TRUE)

plot(PC1_best_model)

library(olsrr)
ols_test_normality(PC1_best_model) #compare residuals against normal distribution
library(lmtest)
bptest(PC1_best_model) #test for heteroskedasticity

#### check for overfitting ####

set.seed(3456)
trainIndex_PC1 <- createDataPartition(model_data$Condition, p = .7, 
                                  list = FALSE, 
                                  times = 1)

PC1Train <- model_data[ trainIndex_PC1,]
PC1Test  <- model_data[-trainIndex_PC1,]

PC1_best_model_train<-lm(skill_pca_dim1_quantity_flaking~1+
                           Condition + 
                           Highest_nBack + 
                           BEAST+
                           Grip_Strength + 
                           Condition:Highest_nBack +
                           Condition:BEAST +
                           core_start_mass,
                   data = PC1Train)

plot(PC1_best_model_train)

summary(PC1_best_model_train)

PC1.mod.test.pred <- predict(PC1_best_model_train, newdata=PC1Test)

PC1.mod.test.R2 <- cor(PC1.mod.test.pred, PC1Test$skill_pca_dim1_quantity_flaking)^2
PC1.mod.test.R2

Overfitting.PC1.mod <- 0.83 - PC1.mod.test.R2
Overfitting.PC1.mod

##### plot PC 1 model ####

pred.mm_condition <- ggpredict(PC1_best_model, terms = c("Condition"))  # this gives overall predictions for the model
levels(pred.mm_condition$group) <- c("Untrained", "Trained")
plot(pred.mm_condition,facet=F, add.data = T,show.title =F, alpha=0.07) + 
      labs(x = "Training condition", 
           y = "Quantity flaking (higher = more quantity flaking)")+
      theme_classic()+
      theme(text = element_text(size=20))

pred.mm_nback <- ggpredict(PC1_best_model, terms = c("Highest_nBack"))  # this gives overall predictions for the model
plot(pred.mm_nback,facet=F, add.data = T,show.title =F, alpha=0.07) + 
      labs(x = "Highest n-back level", 
           y = "Quantity flaking (higher = more quantity flaking)")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_beast<- ggpredict(PC1_best_model, terms = c("BEAST"))  # this gives overall predictions for the model
plot(pred.mm_beast,facet=F, add.data = T,show.title =F, alpha=0.07) + 
      labs(x = "BEAST", 
           y = "Quantity flaking (higher = more quantity flaking)")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_grip<- ggpredict(PC1_best_model, terms = c("Grip_Strength"))  # this gives overall predictions for the model
plot(pred.mm_grip,facet=F, add.data = T,show.title =F, alpha=0.07) + 
      labs(x = "Grip strength", 
           y = "Quantity flaking (higher = more quantity flaking)")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_cond_nback<- ggpredict(PC1_best_model, terms = c("Highest_nBack","Condition"))  # this gives overall predictions for the model
levels(pred.mm_cond_nback$group) <- c("Untrained", "Trained")
plot(pred.mm_cond_nback,facet=F, add.data = F,show.title =F, alpha=0.07) + 
  labs(x = "Highest n-back level", 
       y = "Quantity flaking (higher = more quantity flaking)")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_beast_condition<- ggpredict(PC1_best_model, terms = c("BEAST","Condition"))  # this gives overall predictions for the model
levels(pred.mm_beast_condition$group) <- c("Untrained", "Trained")
plot(pred.mm_beast_condition,facet=F, add.data = F,show.title =F, alpha=0.07) + 
  labs(x = "BEAST score", 
       y = "Quantity flaking (higher = more quantity flaking)")+
  theme_classic()+
  theme(text = element_text(size=20))

##### Skill PC 2 ####

#generate complete model

PC2_complete<- lm(skill_pca_dim2_quality_flaking ~  Condition + Highest_nBack +
                        RPM + BEAST + 
                        fitts_mtavg_ms +
                        Grip_Strength,
                  data = model_data)

#the multimodal selection

PC2_multi <- glmulti(PC2_complete, # use the model with built as a starting point
                   level = 2, #1=just look at main effects, 2=look at main and interraction effects
                   method='g',
                   crit="aic",minsize = 0, maxsize = 7, minK = 0, maxK = -1) 

print(PC2_multi) #all models

summary(PC2_multi@objects[[1]]) #best model

# add mass to final model, but not in model selection
PC2_best_model<-lm(skill_pca_dim2_quality_flaking ~
                         1+Highest_nBack+
                     fitts_mtavg_ms+
                     Grip_Strength+
                     fitts_mtavg_ms:BEAST+
                     Grip_Strength:BEAST+
                     Grip_Strength:fitts_mtavg_ms+
                     Condition:Grip_Strength+
                         core_start_mass,
                   data = model_data)

plot(PC2_best_model)
ols_test_normality(PC2_best_model) #compare residuals against normal distribution
bptest(PC2_best_model) #test for heteroskedasticity

summary(PC2_best_model)

PC2_best_model %>%
      report() %>%
      as.report_table(summary=TRUE)

#### check for overfitting ####

set.seed(3456)
trainIndex <- createDataPartition(model_data$Condition, p = .7, 
                                  list = FALSE, 
                                  times = 1)

PC2Train <- model_data[ trainIndex,]
PC2Test  <- model_data[-trainIndex,]

PC2_best_model_train<-lm(skill_pca_dim2_quality_flaking ~
                           1+Highest_nBack+
                           fitts_mtavg_ms+
                           Grip_Strength+
                           fitts_mtavg_ms:BEAST+
                           Grip_Strength:BEAST+
                           Grip_Strength:fitts_mtavg_ms+
                           Condition:Grip_Strength+
                           core_start_mass,
                   data = PC2Train)

summary(PC2_best_model_train)

perfect.mod.test.pred <- predict(PC2_best_model_train, newdata=PC2Test)
summary(PC2_best_model_train)

perfect.mod.test.R2 <- cor(perfect.mod.test.pred, PC2Test$skill_pca_dim2_quality_flaking)^2
perfect.mod.test.R2

Overfitting.perfect.mod <- 0.75 - perfect.mod.test.R2
Overfitting.perfect.mod #but should report the adjust

#### plot PC 2 model ####

pred.mm_fitts<- ggpredict(PC2_best_model, terms = c("fitts_mtavg_ms"))  # this gives overall predictions for the model
plot(pred.mm_fitts,facet=F, add.data = T,show.title =F, alpha=0.07) + 
  labs(x = "Fitt's score", 
       y = "Quality flaking (higher = increased quality flaking")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_grip<- ggpredict(PC2_best_model, terms = c("Grip_Strength"))  # this gives overall predictions for the model
plot(pred.mm_grip,facet=F, add.data = T,show.title =F, alpha=0.07) + 
  labs(x = "Grip strength", 
       y = "Quality flaking (higher = increased quality flaking")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_nback<- ggpredict(PC2_best_model, terms = c("Highest_nBack"))  # this gives overall predictions for the model
plot(pred.mm_nback,facet=F, add.data = T,show.title =F, alpha=0.07) + 
      labs(x = "Highest n-back", 
           y = "Quality flaking (higher = increased quality flaking")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_mass<- ggpredict(PC2_best_model, terms = c("core_start_mass"))  # this gives overall predictions for the model
plot(pred.mm_mass,facet=F, add.data =T,show.title =F, alpha=0.07) + 
      labs(x = "Nodule starting mass", 
           y = "Quality flaking (higher = increased quality flaking")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_grip_beast<- ggpredict(PC2_best_model, terms = c("Grip_Strength","BEAST"))  # this gives overall predictions for the model
plot(pred.mm_grip_beast,facet=F, add.data =F,show.title =F, alpha=0.07) + 
  labs(x = "Grip strength", 
       y = "Quality flaking (higher = increased quality flaking")+
  theme_classic()+
  theme(text = element_text(size=20))

pred.mm_grip_cond<- ggpredict(PC2_best_model, terms = c("Condition","Grip_Strength"))  # this gives overall predictions for the model
levels(pred.mm_grip_cond$x) <- c("Untrained", "Trained")
plot(pred.mm_grip_cond,facet=F, add.data =F,show.title =F, alpha=0.07) + 
  labs(x = "Training condition", 
       y = "Quality flaking (higher = increased quality flaking")+
  theme_classic()+
  theme(text = element_text(size=20))


