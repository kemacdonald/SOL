source("/Users/kmacdonald/Documents/programming/rscripts/RScripts_v_3.6/libraries_v_3.6.R")

## This script analyzes and plots VLP data for SOL Project 

## First, we read in the data
iChart <- readiChart("")

#### Set Responses

## Option 1 -> Changes all Aways and Centers to Ds for computing reaction time later

iChart$Response[iChart$Response == "A"] <- "D"
iChart$Response[iChart$Response == "C"] <- "D"

## Option 2 ->  If at F0 you are on T, D, or A, then your response becomes an A, and that trial is not included in the accuracy or RT computation. 
## This method only analyzes trials on which the child was on Center at F0.
iChart$Response[iChart$Response == "T"] <- "A"
iChart$Response[iChart$Response == "D"] <- "A"
iChart$Response[iChart$Response == "C"] <- "D"
	
## Check responses 
sum(iChart$Response=="C")
sum(iChart$Response=="D")
sum(iChart$Response=="T")
sum(iChart$Response=="A")

## Change all trials to the same condition -> 'van' for vanilla.

iChart$Condition <- "van"

#### Explore data

## At 1800 ms after the start of the sign, what is the disribution of C, T, and D?

sum(iChart$"1800" == "0.5")
sum(iChart$"1800" == "1")
sum(iChart$"1800" == "0")

## Checks the distribution of C, D, T, and A at F0

sum(iChart$"0" == "0.5")
sum(iChart$"0" == "1")
sum(iChart$"0" == "0")
sum(iChart$"0" == ".")
sum(iChart$"0" == "-")

######### Compute accuracy #########

## define onset and reject prescreen

iChart <- defineOnset(iChart[iChart$Prescreen.Notes == "",], critonset=0, includeAways=TRUE)

# compute gaps
iChart <- computeStatistics(iChart, startWindow=0, endWindow=4000)

# reject trials with extreme RT and gaps
iChart <- filteriChart(iChart, minRT=300, maxRT=3500, maxfirstgap=15, maxlonggap=15)

## Compute accuracy

accuracy <- poolData(meanAccuracy(iChart, startWindowAcc=300, endWindowAcc=3000), RejectFirstGap=FALSE,RejectLongestGap=FALSE, RejectRT=FALSE, color=TRUE, dependent="Accuracy", group="", facet="", dodge="", xlab="", ylab= "Proportion\n  Looking\n  to target", paired=TRUE, miny = 0.2, maxy = 0.80, size=13, legend.direction="horizontal", legend.position="bottom", breaks=c(0.25, 0.50, 0.75))
                     
## change months to a numeric to do median split
iChart$Months <- as.numeric(iChart$Months) 

## median split by age (months)
iChart$Condition <- ifelse(iChart$Months >= median(iChart$Months, na.rm=TRUE), "H", "L") 

## compute accuracy for each age group (high and low)

accuracy.by.condition <- poolData(meanAccuracy(iChart, startWindowAcc=300, endWindowAcc=3000), RejectFirstGap=FALSE, 
                     RejectLongestGap=FALSE, RejectRT=FALSE, color=TRUE, dependent="Accuracy", group="", 						 			 facet="Condition", dodge="", xlab="", 
                     ylab= "Proportion\n  Looking\n  to target", paired=TRUE, miny = 0.2, maxy = 0.80, size=13, 
                     legend.direction="horizontal", legend.position="bottom", breaks=c(0.25, 0.50, 0.75))

## Plot by condition 

createPlots(iChart, startWindow=300, endWindow=3000, RejectLongestGap=FALSE, 
            RejectFirstGap=FALSE, RejectRT=FALSE, color=TRUE, smooth=400, 
            targetEnd=800, carrier="", targets=c(""), 
            group="",  plotStats="PP", miny = 0.4, maxy=0.85, size=15, 
            legend.direction = "horizontal", legend.position=c(0.7, 0.95), 
            breaks=c(0.25, 0.50, 0.75), x.target=0.33)


####### RT ###########

## What to do with the 0.5 after the shift? 

## Read in iChart
iChart <- readiChart("/Users/kmacdonald/Documents/Projects/SOL/SOL_Data/SOL_Child_Data_n23/SOL.n23.ichart.txt")

iChart[iChart$"0" == 0.5, "0"] <- "0" ## At F0, changes all 0.5s to 0

## for loop: 
col_i <- match("0",names(iChart)):match("3500",names(iChart))
for(i in col_i) {
  iChart[iChart[,i]== 0.5, i] <- 0
}

## Set all trials to condition 'van' for vanilla
iChart$Condition <- "van"

## Read Months variable as numeric
iChart$Months <- as.numeric(iChart$Months)

## Split data set by median age
iChart$Condition <- ifelse(iChart$Months >= median(iChart$Months, na.rm=TRUE), "H", "L")

## Define critical onset, in this case we set it to 0 
iChart <- defineOnset(iChart[iChart$Prescreen.Notes == "",], critonset=0, includeAways=FALSE)

## Compute statistics, setting the start and end windows 
iChart <- computeStatistics(iChart, startWindow=0, endWindow=4000)

## Filter iChart, setting your filtering criteria 
iChart <- filteriChart(iChart, minRT=300, maxRT=4000, maxfirstgap=15, maxlonggap=15)

## Compute RT for each subjet
RT <- poolData(iChart[iChart$Response == "D",], RejectFirstGap=FALSE, RejectLongestGap=FALSE, RejectRT=TRUE, color=FALSE, dependent="RT", 
               group="", facet="", dodge="Response", xlab="", ylab="mean RT (ms)", paired=TRUE, 
               miny = 400, maxy=1300, size=13, legend.direction = "horizontal", legend.position="bottom", breaks=c(400, 800, 1200))

cor.test(as.numeric(as.character(RT$Months)), RT$Accuracy)

## Create RT plots (OC Plots)
createPlots(iChart, startWindow=0, endWindow=3500, RejectLongestGap=FALSE, RejectFirstGap=FALSE, RejectRT=FALSE, 
            color=TRUE, smooth=200, targetEnd=800, carrier="Where's the", targets=c("b a l l", "d o f a"), 
            group="",  plotStats="OC_D", miny = 0, maxy=0.85, size=15, legend.direction = "horizontal", 
            legend.position=c(0.7, 0.95), breaks=c(0.25, 0.50, 0.75), x.target=0.33)


##### Create aggregated table of acc, rt, cdi, and age by subject ######

ages = aggregate(Months~Sub.Num, data=iChart, mean)

RT <- merge(RT, ages, by="Sub.Num")
colnames(RT)[5] <- "RT"


all <- merge(RT, accuracy, by="Sub.Num")

colnames(all)[5] <- "van_d"
colnames(all)[11] <- "acc"

###### Analysis #######

## correlation between reaction time and age

cor.test(RT$RT, RT$Months.y)
plot(RT$RT, RT$Months.y)

## correlation between accuracy and age

cor.test(all$van, all$Months)
plot(all$van, all$Months)

#### Subset data: only include subjects over 16 months of age and remove kid with CI (30030)

all_over16 <- all[which(all$Months >= 16 & all$Sub.Num != 30030),]

## correlation between reaction time and age with subsetted data 
plot(all_over18$van_D, all_over18$Months)
cor.test(all_over18$van_D, all_over18$Months)

## accuracy analysis with subsetted data: acc by age
plot(all_over18$van, all_over18$Months)
cor.test(all_over18$van, all_over18$Months)


#### Age analyses ######

## T-test of reaction time difference by age median split ##

t.test(all$H_H_D, all$L_L_D)

## T-test of accuracy difference by age median split ##

t.test(all$H_H, all$L_L)

#### CDI Analyses ######

read.csv("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /sol_cdi_data.csv") -> CDI
colnames(CDI)[1] <- "Sub.Num"

CDI <- merge(all, CDI, by="Sub.Num")

write.csv(CDI, file="SOL.n23.alldata.csv")

#### Subset data: only include subjects under 30 months of age and remove kid with CI (30030). We do this because there are ceiling effects on the CDI for kids over 30 months.

CDI_under30 <- CDI[which(CDI$Months <= 30 & all$Sub.Num != 30030),]
cor.test(CDI_under30$signs_produced, CDI_under30$vanilla)
cor.test(CDI_under30$signs_produced, CDI_under30$Months)

plot(CDI_under30$signs_produced, CDI_under30$vanilla)
plot(CDI_under30$signs_produced, CDI_under30$Months)


#### Subset data: only include subjects under 32 months of age and remove kid with CI (30030).  We do this because there are ceiling effects on the CDI for kids over 30 months.

CDI_under32 <- CDI[which(CDI$Months <= 32 & all$Sub.Num != 30030),]
cor.test(CDI_under32$signs_produced, CDI_under32$vanilla)
cor.test(CDI_under32$signs_produced, CDI_under32$Months)

plot(CDI_under32$signs_produced, CDI_under32$vanilla)
plot(CDI_under32$signs_produced, CDI_under32$Months)


## correlations
cor.test(CDI$phrases, CDI$vanilla)
cor.test(CDI$signs_produced, CDI$vanilla)
cor.test(CDI$signs_produced, CDI$Months)
cor.test(CDI$signs_produced, CDI$van_D)

## plots 
plot(CDI$phrases, CDI$vanilla)
plot(CDI$vanilla, CDI$signs_produced)
plot(CDI$Months, CDI$signs_produced)
plot(CDI$signs_produced, CDI$van_D)


# REMOVE UNKNOWN
read.table("/Users/ricardobion/Downloads/sol_knows.txt", header=TRUE) -> known

known$birdy <- 1

melt(known, id=c("Sub.Num", "response")) ->known
known$variable <- as.character(known$variable)
known$variable[known$variable == "cat"] <- "kitty"
known$variable[known$variable == "bear"] <- "teddy"

names(known)[3] <- "Target.Image"

iChart$Target.Image <- sub("[1-4][RL].pct", "", iChart$Target.Image)
iChart$Target.Image <- sub("[1-4].pct", "", iChart$Target.Image)
merge(iChart, known, by = c("Sub.Num", "Target.Image"), all.x=TRUE, all.y=FALSE) -> iChart

iChart <- iChart[iChart$value == 1 & !is.na(iChart$value),]


########## ADULT DATA ################

iChart <- readiChart("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /SOL_Adult_Data_n11/SOL.adult.iChart.n11.txt")

## Remove subject 30002 -> not performing above chance and very slow reaction times

iChart <- iChart[which(iChart$Sub.Num!= 30002),]

#### Set Responses

	## Option 1 -> Changes all Aways and Centers to Ds for computing reaction time later

	iChart$Response[iChart$Response == "A"] <- "D"
	iChart$Response[iChart$Response == "C"] <- "D"

	## Option 2 ->  If at F0 you are on T, D, or A, then your response becomes an A, and that trial is not included in the accuracy or RT computation. This method only analyzes trials on which the child was on Center at F0.
	iChart$Response[iChart$Response == "T"] <- "A"
	iChart$Response[iChart$Response == "D"] <- "A"
	iChart$Response[iChart$Response == "C"] <- "D"
	
## Check responses 
sum(iChart$Response=="C")
sum(iChart$Response=="D")
sum(iChart$Response=="T")
sum(iChart$Response=="A")

## Change all trials to the same condition -> 'van' for vanilla.

iChart$Condition <- "van"

#### Explore data

	## At 1800 ms after the start of the sign, what is the disribution of C, T, and D?

	sum(iChart$"1800" == "0.5")
	sum(iChart$"1800" == "1")
	sum(iChart$"1800" == "0")

	## Checks the distribution of C, D, T, and A at F0

	sum(iChart$"0" == "0.5")
	sum(iChart$"0" == "1")
	sum(iChart$"0" == "0")
	sum(iChart$"0" == ".")
	sum(iChart$"0" == "-")

######### Compute accuracy #########

## Get regex to find C-T and C-D and plot these as different funcitons. STILL NEED TO DO THIS ##


## define onset and reject prescreen

iChart <- defineOnset(iChart[iChart$Prescreen.Notes == "",], critonset=0, includeAways=TRUE)

# compute gaps
iChart <- computeStatistics(iChart, startWindow=0, endWindow=3500)

# reject trials with extreme RT and gaps
iChart <- filteriChart(iChart, minRT=300, maxRT=2000, maxfirstgap=15, maxlonggap=15)

## Compute accuracy

accuracy <- poolData(meanAccuracy(iChart, startWindowAcc=300, endWindowAcc=2000), RejectFirstGap=FALSE,RejectLongestGap=FALSE, RejectRT=FALSE, color=TRUE, dependent="Accuracy", group="", facet="", dodge="", xlab="", ylab= "Proportion\n  Looking\n  to target", paired=TRUE, miny = 0.2, maxy = 0.75, size=13, legend.direction="horizontal", legend.position="bottom", breaks=c(0.25, 0.50, 0.75))
                     
## change months to a numeric to do median split
iChart$Months <- as.numeric(iChart$Months) 

## PP Plot

createPlots(iChart, startWindow=300, endWindow=3000, RejectLongestGap=FALSE, 
            RejectFirstGap=FALSE, RejectRT=FALSE, color=TRUE, smooth=400, 
            targetEnd=800, carrier="", targets=c(""), 
            group="",  plotStats="PP", miny = 0.4, maxy=1.0, size=15, 
            legend.direction = "horizontal", legend.position=c(0.7, 0.95), 
            breaks=c(0.25, 0.50, 0.75), x.target=0.33)


####### RT ###########

## What to do with the 0.5 after the shift? 

## Read in iChart
iChart <- readiChart("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /SOL_Adult_Data_n11/SOL.adult.iChart.n11.txt")

iChart[iChart$"0" == 0.5, "0"] <- "0" ## At F0, changes all 0.5s to 0

## for loop: 
col_i <- match("0",names(iChart)):match("3500",names(iChart))
for(i in col_i) {
  iChart[iChart[,i]== 0.5, i] <- 0
}

## Set all trials to condition 'van' for vanilla
iChart$Condition <- "van"

## Read Months variable as numeric
iChart$Months <- as.numeric(iChart$Months)

## Define critical onset, in this case we set it to 0 
iChart <- defineOnset(iChart[iChart$Prescreen.Notes == "",], critonset=0, includeAways=FALSE)

## Compute statistics, setting the start and end windows 
iChart <- computeStatistics(iChart, startWindow=0, endWindow=3500)

## Filter iChart, setting your filtering criteria 
iChart <- filteriChart(iChart, minRT=300, maxRT=4000, maxfirstgap=15, maxlonggap=15)

## Compute RT for each subjet
RT <- poolData(iChart[iChart$Response == "D",], RejectFirstGap=FALSE, RejectLongestGap=FALSE, RejectRT=TRUE, color=FALSE, dependent="RT", 
               group="", facet="", dodge="Response", xlab="", ylab="mean RT (ms)", paired=TRUE, 
               miny = 400, maxy=1300, size=13, legend.direction = "horizontal", legend.position="bottom", breaks=c(400, 800, 1200))

## Create RT plots (OC Plots)
createPlots(iChart, startWindow=0, endWindow=3500, RejectLongestGap=FALSE, RejectFirstGap=FALSE, RejectRT=FALSE, 
            color=TRUE, smooth=200, targetEnd=800, carrier="Where's the", targets=c("b a l l", "d o f a"), 
            group="",  plotStats="OC_D", miny = 0, maxy=0.85, size=15, legend.direction = "horizontal", 
            legend.position=c(0.7, 0.95), breaks=c(0.25, 0.50, 0.75), x.target=0.33)
            


##### Create aggregated table of acc, rt, cdi, and age by subject ######

write.csv(RT, file="SOL.n10.adult.csv")

## Comparing ages ANOVA ####

d1 <- read.csv("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /SOL.n33.alldata.csv"

rs.rt <- anova(lm(van_D ~ group, d1))
rs.rt

## Direct comparison of RT: children vs. adults ##

d1$group <- as.character(d1$group)

	d1$group[d1$group == "H"] <- "C"
	d1$group[d1$group == "L"] <- "C"

t.test(van_D~group, data=d1)

######### Adult Prime Data ############

iChart <- readiChart("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /SOL_Adult_Data_n11/SOL.adult.iChart.n11.txt")
	
## Remove subject 30002 -> not performing above chance and very slow reaction times

iChart <- iChart[which(iChart$Sub.Num!= 30002),]



####### Plot adult and child data together ########

## First, we read in the data
iChart <- readiChart("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /SOL_All_Data_n33/SOL.ALL.iChart.n33.txt")

#### Set Responses

	## Option 1 -> Changes all Aways and Centers to Ds for computing reaction time later

	iChart$Response[iChart$Response == "A"] <- "D"
	iChart$Response[iChart$Response == "C"] <- "D"

	## Option 2 ->  If at F0 you are on T, D, or A, then your response becomes an A, and that trial is not included in the accuracy or RT computation. This method only analyzes trials on which the child was on Center at F0.
	iChart$Response[iChart$Response == "T"] <- "A"
	iChart$Response[iChart$Response == "D"] <- "A"
	iChart$Response[iChart$Response == "C"] <- "D"
	
## Check responses 
sum(iChart$Response=="C")
sum(iChart$Response=="D")
sum(iChart$Response=="T")
sum(iChart$Response=="A")

## Change all trials to the same condition -> 'van' for vanilla.

iChart$Condition <- "van"

#### Explore data

	## At 1800 ms after the start of the sign, what is the disribution of C, T, and D?

	sum(iChart$"1800" == "0.5")
	sum(iChart$"1800" == "1")
	sum(iChart$"1800" == "0")

	## Checks the distribution of C, D, T, and A at F0

	sum(iChart$"0" == "0.5")
	sum(iChart$"0" == "1")
	sum(iChart$"0" == "0")
	sum(iChart$"0" == ".")
	sum(iChart$"0" == "-")

######### Compute accuracy #########

## Get regex to find C-T and C-D and plot these as different funcitons. STILL NEED TO DO THIS ##


## define onset and reject prescreen

iChart <- defineOnset(iChart[iChart$Prescreen.Notes == "",], critonset=0, includeAways=TRUE)

# compute gaps
iChart <- computeStatistics(iChart, startWindow=0, endWindow=3500)

# reject trials with extreme RT and gaps
iChart <- filteriChart(iChart, minRT=300, maxRT=3500, maxfirstgap=15, maxlonggap=15)

## Compute accuracy

accuracy <- poolData(meanAccuracy(iChart, startWindowAcc=300, endWindowAcc=3000), RejectFirstGap=FALSE,RejectLongestGap=FALSE, RejectRT=FALSE, color=TRUE, dependent="Accuracy", group="", facet="", dodge="", xlab="", ylab= "Proportion\n  Looking\n  to target", paired=TRUE, miny = 0.2, maxy = 0.80, size=13, legend.direction="horizontal", legend.position="bottom", breaks=c(0.25, 0.50, 0.75))
                     
## change months to a numeric to split our conditions
iChart$Months <- as.numeric(iChart$Months) 

## split sample into three groups based on age: young (< 24months, > 24months, adults)

iChart$Condition <- ifelse(iChart$Months < 24, "< 24 Months", 
               ifelse(iChart$Months >= 24 & iChart$Months < 50, "> 24 Months ", "Adults"))


## compute accuracy for each age group

accuracy.by.condition <- poolData(meanAccuracy(iChart, startWindowAcc=300, endWindowAcc=3000), RejectFirstGap=FALSE, 
                     RejectLongestGap=FALSE, RejectRT=FALSE, color=TRUE, dependent="Accuracy", group="", 						 			 facet="Condition", dodge="", xlab="", 
                     ylab= "Proportion\n  Looking\n  to target", paired=TRUE, miny = 0.2, maxy = 0.80, size=13, 
                     legend.direction="horizontal", legend.position="bottom", breaks=c(0.25, 0.50, 0.75))

## Plot by condition 

createPlots(iChart, startWindow=300, endWindow=3000, RejectLongestGap=FALSE, 
            RejectFirstGap=FALSE, RejectRT=FALSE, color=TRUE, smooth=400, 
            targetEnd=800, carrier="", targets=c(""), 
            group="",  plotStats="PP", miny = 0.4, maxy=0.95, size=15, 
            legend.direction = "vertical", legend.position=c(0.85, 0.9), 
            breaks=c(0.25, 0.50, 0.75), x.target=0.33)


####### RT ###########

## What to do with the 0.5 after the shift? 

## Read in iChart
iChart <- readiChart("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /SOL_Child_Data_n23/SOL.n23.ichart.txt")

iChart[iChart$"0" == 0.5, "0"] <- "0" ## At F0, changes all 0.5s to 0

## for loop: 
col_i <- match("0",names(iChart)):match("3500",names(iChart))
for(i in col_i) {
  iChart[iChart[,i]== 0.5, i] <- 0
}

## Set all trials to condition 'van' for vanilla
iChart$Condition <- "van"

## Read Months variable as numeric
iChart$Months <- as.numeric(iChart$Months)

## Split data set by median age
iChart$Condition <- ifelse(iChart$Months >= median(iChart$Months, na.rm=TRUE), "H", "L")

## Define critical onset, in this case we set it to 0 
iChart <- defineOnset(iChart[iChart$Prescreen.Notes == "",], critonset=0, includeAways=FALSE)

## Compute statistics, setting the start and end windows 
iChart <- computeStatistics(iChart, startWindow=0, endWindow=4000)

## Filter iChart, setting your filtering criteria 
iChart <- filteriChart(iChart, minRT=300, maxRT=4000, maxfirstgap=15, maxlonggap=15)

## Compute RT for each subjet
RT <- poolData(iChart[iChart$Response == "D",], RejectFirstGap=FALSE, RejectLongestGap=FALSE, RejectRT=TRUE, color=FALSE, dependent="RT", 
               group="", facet="", dodge="Response", xlab="", ylab="mean RT (ms)", paired=TRUE, 
               miny = 400, maxy=1300, size=13, legend.direction = "horizontal", legend.position="bottom", breaks=c(400, 800, 1200))

cor.test(as.numeric(as.character(RT$Months)), RT$Accuracy)

## Create RT plots (OC Plots)
createPlots(iChart, startWindow=0, endWindow=3500, RejectLongestGap=FALSE, RejectFirstGap=FALSE, RejectRT=FALSE, 
            color=TRUE, smooth=200, targetEnd=800, carrier="Where's the", targets=c("b a l l", "d o f a"), 
            group="",  plotStats="OC_D", miny = 0, maxy=0.85, size=15, legend.direction = "horizontal", 
            legend.position=c(0.7, 0.95), breaks=c(0.25, 0.50, 0.75), x.target=0.33)






