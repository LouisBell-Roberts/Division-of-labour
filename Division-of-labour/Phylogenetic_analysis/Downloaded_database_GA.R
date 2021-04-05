#Online database manipulation - Global ants
#Louis Bell-Roberts
#02/04/2021

library(tidyverse)

#Read in datafiles
data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Downloaded online databases/GlobalAnts.csv", header = T)
primary_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)

#Create a new column which is the combination of genus and species name = "animal"
data1 <- data %>% unite("animal", Genus:Species, remove = F)

#Remove species which don't have HW data
##Convert data to numeric
data2 <- data1
data2$Head.width.across.eyes..mm. <- as.numeric(as.character(data1$Head.width.across.eyes..mm.))
##Filter out species which have NA
data3 <- filter(data2, Head.width.across.eyes..mm. >= 0)
View(data3)

#Calculate a mean for each species by splitting up the HW column by species and claculating a mean
mean_splitHW <- sapply(split(data3$Head.width.across.eyes..mm.,data3$animal), mean)
sd_splitHW <- sapply(split(data3$Head.width.across.eyes..mm.,data3$animal), sd)

#Convert to a dataframe
mean_splitHW_df<-as.data.frame(mean_splitHW)
sd_splitHW_df<-as.data.frame(sd_splitHW)

#Give dataframe an animal column name
mean_splitHW_df.1<-rownames_to_column(mean_splitHW_df, var = "animal")
sd_splitHW_df.1<-rownames_to_column(sd_splitHW_df, var = "animal")

#Check how much overlap with my existing database
common_sp <- intersect(mean_splitHW_df.1$animal,primary_data$animal)
common_sp1 <- intersect(sd_splitHW_df.1$animal,primary_data$animal)

#Filter for these rows from my primary dataset
filter_primary<-primary_data[match(common_sp, primary_data$animal), ]





#Create combined dataset - creates dataset with 270 species even though only 269 species match???
joined_df <- left_join(primary_data, sd_splitHW_df.1)
View(select(joined_df, animal, sd_splitHW, eff.mating.freq.MEAN.harmonic))

#Filter to select only species with data on body size, mating frequency , polymorphism, colony size.
joined_filtered <- filter(joined_df, sd_splitHW >= 0, eff.mating.freq.MEAN.harmonic >=0)
View(joined_filtered)

joined_filtered1 <- filter(joined_df, sd_splitHW >= 0, colony.size >=0)

#Create linear models - don;t show any real trends
lm <- lm(sd_splitHW ~ eff.mating.freq.MEAN.harmonic, data = joined_filtered)
summary(lm)
plot(sd_splitHW ~ eff.mating.freq.MEAN.harmonic, joined_filtered)
abline(lm)                      


lm2 <- lm(sd_splitHW ~ log(colony.size), data = joined_filtered1) 
summary(lm2)                      
plot(sd_splitHW ~ log(colony.size), joined_filtered)
abline(lm2)










