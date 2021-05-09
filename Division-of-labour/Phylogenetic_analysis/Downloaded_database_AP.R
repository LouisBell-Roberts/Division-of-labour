#Online database manipulation - Ant profiler
#Louis Bell-Roberts
#02/04/2021


library(tidyverse)

#Read in datafiles
data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Downloaded online databases/AntProfiler.csv", header = T)
primary_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)

#Remove " " from the "Name" column
data1 <- data
data1$Name <- str_replace(data$Name, " ", "_")
View(data1)

#Filter data so I only get rows with Min and Max values
data2 <- filter(data1, MinBodySize >=0, MaxBodySize >= 0)

#Create new column which is a measure of worker variation = Max-Min/Min
data3 <- mutate(data2, Variation = (MaxBodySize-MinBodySize)/MinBodySize)
View(select(data3, Name, MinBodySize, MaxBodySize, Variation, WorkerPolymorphism))

#Check how much overlap with my existing database - 269
common_sp <- intersect(data3$Name,primary_data$animal)

#Filter for these rows from my primary dataset
filter_primary<-primary_data[match(common_sp, primary_data$animal), ]

#Create combined dataset - creates dataset with 270 species even though only 269 species match???
joined_df <- left_join(primary_data, data3, 
                       by = c("animal" = "Name"))
View(select(joined_df, animal, MinBodySize, MaxBodySize, Variation, WorkerPolymorphism))

#Filter to select only species with data on body size, mating frequency , polymorphism, colony size.
joined_filtered <- filter(joined_df, Variation >= 0.0, eff.mating.freq.MEAN.harmonic >=0)
View(joined_filtered)

joined_filtered1 <- filter(joined_df, Variation >= 0, colony.size >=0)
joined_filtered2 <- filter(joined_df, Variation >= 0.0, W.policing.clean >=0)
joined_filtered3 <- filter(joined_df, Variation >= 0.0, WorkerPolymorphism >=0)


View(select(joined_filtered, animal, MinBodySize, MaxBodySize, Variation, WorkerPolymorphism))

#Create linear models
lm <- lm(Variation ~ eff.mating.freq.MEAN.harmonic, data = joined_filtered)
summary(lm)
plot(Variation ~ eff.mating.freq.MEAN.harmonic, joined_filtered)
abline(lm)                      


lm2 <- lm(Variation ~ log(colony.size), data = joined_filtered1) 
summary(lm2)                      
plot(Variation ~ log(colony.size), joined_filtered)
abline(lm2)

lm3 <- lm(Variation ~ W.policing.clean, data = joined_filtered2) 
summary(lm3)                      
plot(Variation ~ W.policing.clean, joined_filtered2)
abline(lm3)

lm4 <- lm(Variation ~ WorkerPolymorphism, data = joined_filtered3) 
summary(lm4)                      
plot(Variation ~ WorkerPolymorphism, joined_filtered3)
abline(lm4)
                      