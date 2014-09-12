source("/Users/kylemacdonaldadmin/RScripts_v_3.6/libraries_v_3.6.R")
iChart <- readiChart("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data/SOL.n15.ichart.txt")

iChart2 <- readiChart("/Volumes/Babylab/General/Studies/SOL/SOL Data/SOL Online Data/SOL-INC Adult Prime Noun Onset/SOL-INC-Adult.iChart.n8.prime.nounonset.txt")

iChart[iChart$Sub.Num != "30008",] -> iChart ## removes 85-month-old

sum(iChart2$Response.Center=="C")
sum(iChart2$Response=="D")
sum(iChart2$Response=="T")
sum(iChart2$Response=="A")

## Changes all As and Cs to Ds

iChart2$Response[iChart$Response == "A"] <- "D"
iChart2$Response.Center[iChart2$Response.Center == "C"] <- "D"

## This line changes all trials to the same condition -> 'van' for vanilla.

iChart2$Condition <- "van"
str(iChart2$Condition)

## At 1800 ms after the start of the sign, what is the disribution of C, T, and D?

sum(iChart2$"1800" == "0.5")
sum(iChart2$"1800" == "1")
sum(iChart2$"1800" == "0")

## This checks the distribution of C, D, T, and A at F0

sum(iChart2$"0" == "0.5")
sum(iChart2$"0" == "1")
sum(iChart2$"0" == "0")

######### Compute accuracy #########

## Get regex to find C-T and C-D and plot these as different funcitons.

## If at F0 you are on the T or D, then your response becomes an A, and that trial is not included in the accuracy or RT computation.

accuracy <- poolData(meanAccuracy(iChart2, startWindowAcc=300, endWindowAcc=3300), RejectFirstGap=FALSE, 
                     RejectLongestGap=FALSE, RejectRT=FALSE, color=TRUE, dependent="Accuracy", group="", facet="", dodge="", xlab="", 
                     ylab= "Proportion\n  Looking\n  to target", paired=TRUE, miny = 0.2, maxy = 0.80, size=13, 
                     legend.direction="horizontal", legend.position="bottom", breaks=c(0.25, 0.50, 0.75))

iChart$Months <- as.numeric(iChart$Months) ## change months to a numeric

iChart$Condition <- ifelse(iChart$Months >= median(iChart$Months, na.rm=TRUE), "H", "L") ## divide groups into high and low based on median age

accuracy.by.condition <- poolData(meanAccuracy(iChart, startWindowAcc=300, endWindowAcc=3300), RejectFirstGap=FALSE, 
                     RejectLongestGap=FALSE, RejectRT=FALSE, color=TRUE, dependent="Accuracy", group="", facet="Condition", dodge="", xlab="", 
                     ylab= "Proportion\n  Looking\n  to target", paired=TRUE, miny = 0.2, maxy = 0.80, size=13, 
                     legend.direction="horizontal", legend.position="bottom", breaks=c(0.25, 0.50, 0.75))

createPlots(iChart, startWindow=300, endWindow=3500, RejectLongestGap=FALSE, 
            RejectFirstGap=FALSE, RejectRT=FALSE, color=TRUE, smooth=400, 
            targetEnd=800, carrier="", targets=c(""), 
            group="",  plotStats="PP", miny = 0.4, maxy=0.85, size=15, 
            legend.direction = "horizontal", legend.position=c(0.7, 0.95), 
            breaks=c(0.25, 0.50, 0.75), x.target=0.33)


####### RT ###########

## What to do with the 0.5 after the shift? 

## Read in iChart
iChart <- readiChart("/Volumes/Babylab/General/Studies/SOL/SOL Data/SOL Online Data/SOL-INC Adult Prime Noun Onset/SOL-INC-Adult.iChart.n8.prime.nounonset.txt")
iChart[iChart$Sub.Num != "30008",] -> iChart ## Remove 85-month-old

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
iChart <- filteriChart(iChart, minRT=300, maxRT=3500, maxfirstgap=15, maxlonggap=15)

## Compute RT for each subjet
RT <- poolData(iChart, RejectFirstGap=FALSE, RejectLongestGap=FALSE, RejectRT=TRUE, color=FALSE, dependent="RT", 
               group="", facet="", dodge="Response", xlab="", ylab="mean RT (ms)", paired=TRUE, 
               miny = 400, maxy=1300, size=13, legend.direction = "horizontal", legend.position="bottom", breaks=c(400, 800, 1200))

## Create RT plots (OC Plots)
createPlots(iChart, startWindow=0, endWindow=3500, RejectLongestGap=FALSE, RejectFirstGap=FALSE, RejectRT=FALSE, 
            color=TRUE, smooth=200, targetEnd=800, carrier="Where's the", targets=c("b a l l", "d o f a"), 
            group="",  plotStats="OC_D", miny = 0, maxy=0.85, size=15, legend.direction = "horizontal", 
            legend.position=c(0.7, 0.95), breaks=c(0.25, 0.50, 0.75), x.target=0.33)

## Create aggregated table of acc and rt for each sub

ages = aggregate(Months~Sub.Num, data=iChart, mean)

RT <- merge(RT, ages, by="Sub.Num")

all <- merge(RT, accuracy.by.condition, by="Sub.Num")

## correlation between accuracy and age

cor.test(RT$van_D[all$Sub.Num != "30008"], RT$Months[all$Sub.Num != "30008"])
plot(RT$van_D[all$Sub.Num != "30008"], RT$Months[all$Sub.Num != "30008"])

## T-test of reaction time difference by age median split ##

t.test(all$H_D, all$L_D)

## T-test of accuracy difference by age median split ##

t.test(all$H_H, all$L_L)

## Add cdi data

CDI <- merge(RT, CDI, by="Sub.Num")


cor.test(all$van[all$Sub.Num != "30008"], all$Months[all$Sub.Num != "30008"])
plot(all$van[all$Sub.Num != "30008"], all$Months[all$Sub.Num != "30008"])


read.csv("/Users/ricardobion/Downloads/sol_cdi.csv") -> CDI



cor.test(CDI$phrases[CDI$Sub.Num != "30008"], CDI$van[CDI$Sub.Num != "30008"])

plot(CDI$produced[CDI$Sub.Num != "30008"], CDI$Months[CDI$Sub.Num != "30008"])

plot(all$van[all$Sub.Num != "30008"], all$Months[all$Sub.Num != "30008"])


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
