# this script removes all current objects and loads the necessary libraries 

# create function for flagging these trial types
trial_type_fun <- function(trial) { 
    start_col <- which(names(trial)=="0")
    end_col <- which(names(trial)=="2700")
    trial_type <- "no-shift"
    for (col in start_col:end_col) {
        curr_val <- trial[col]
        next_val <- trial[col+1]
        
        if(curr_val == "." & next_val %in% c("0", "1", "0.5", ".5")) {
            if(next_val == "1") {
                trial_type <- "C-T"
            } else if(next_val == "0.5" | next_val == ".5") {
                trial_type <- "C-C"
            } else {
                trial_type <- "C-D"    
            }
            break
        }
    }
    return(trial_type)
}

# removes all current objects
#rm(list=ls()) 

print("loading library")

options(scipen=10)

temp <- installed.packages()

if(!("reshape" %in% temp)){
	install.packages("reshape", repos = "http://cran.r-project.org")
	}
	
if(!("ggplot2" %in% temp)){
	install.packages("ggplot2", repos = "http://cran.r-project.org")
	}	
	
if(!("lme4" %in% temp)){
	install.packages("lme4", repos = "http://cran.r-project.org")
	}

# load libraries
library(reshape)
library(ggplot2)
library(grid)
library(lme4)

msToFrames <- function(value) { return(value*3/100) }
framesToMs <- function(value) { return(value*100/3) }

minimumStare <-  function(iChart) {
	endFrame<- msToFrames(startWindow)
	startFrame <- which(names(iChart)=="0")
	responses <- iChart[, startFrame + endFrame]
	backfilledResponses <- rep(responses, endFrame)
	backfilledMatrix <- matrix(backfilledResponses, nrow = length(responses), ncol = endFrame)
	same <- iChart[, (startFrame+1):(startFrame+ endFrame)] == 	backfilledMatrix
	sumSame <- rowSums(same)
	iChart <- iChart[sumSame == endFrame,]
	return(iChart)
	}
	

convertiChart <- function() {
	
	}

backFill <- function(iChart, startWindow, framesToBackFill) {
	
	startWindow <- which(names(iChart)==startWindow)
	iChart$Response <- ifelse(iChart[,startWindow] == "1", "T", ifelse(iChart[,startWindow] == "0", "D", "A"))
		
	framesToBeBackfilled <- ((iChart[,startWindow]=='-'|iChart[,startWindow]=='.')&(iChart[,startWindow+framesToBackFill]=='0'|iChart[,startWindow+framesToBackFill]=='1'))

	responses <- iChart[framesToBeBackfilled, startWindow+framesToBackFill]
	backfilledResponses <- rep(responses, framesToBackFill)
	backfilledMatrix <- matrix(backfilledResponses, nrow = length(responses), ncol = framesToBackFill)
	iChart[framesToBeBackfilled, startWindow:(startWindow+framesToBackFill-1)] <- backfilledMatrix
	iChart$BackFilled <- framesToBeBackfilled
	iChart$Response[framesToBeBackfilled] <- ifelse(iChart[framesToBeBackfilled, startWindow+framesToBackFill] == "1", "T", "D")
	return(iChart)
	}
	
	


poolRT <- function(iChart, GoodFirstGap, GoodLongestGap, GoodRT) {
	
	# pool RT accross subject
	filteredRT <- iChart
	if(GoodFirstGap) filteredRT <- filteredRT[filteredRT$GoodFirstGap | is.na(filteredRT$GoodFirstGap),]
	if(GoodLongestGap) filteredRT <- filteredRT[filteredRT$GoodLongestGap | is.na(filteredRT$GoodLongestGap),]
	if(GoodRT) filteredRT <- filteredRT[filteredRT$GoodRT,]
	filteredRT <- filteredRT[filteredRT$Response == "D" | filteredRT$Response == "T",]

	sumRT <- aggregate(list(Sum_RT=as.numeric(filteredRT$RT)), list(Sub.Num = filteredRT$Sub.Num, Response = filteredRT$Response, Condition = filteredRT$Condition), function(value){sum(!is.na(value))})
	pooledRT <- aggregate(list(RT=as.numeric(filteredRT$RT)), list(Sub.Num = filteredRT$Sub.Num, Response = filteredRT$Response, Condition = filteredRT$Condition), mean, na.rm=T)
	sdRT <- aggregate(list(SE_RT=as.numeric(filteredRT$RT)), list(Sub.Num = filteredRT$Sub.Num, Response = filteredRT$Response, Condition = filteredRT$Condition), se)
	
	tableRT <- cbind(pooledRT, sumRT["Sum_RT"], sdRT["SE_RT"])
	
		longestFirstGap <- max(filteredRT$firstGap, na.rm=T)
		longestLongestGap <- max(filteredRT$longestGap, na.rm=T)
		longestRT_D <- max(filteredRT$RT, na.rm=T)
		shortestRT_D <- min(filteredRT$RT, na.rm=T)
		
		
	npar <- length(unique(iChart$Sub.Num))
	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_meanRT_by_subs_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_n_", npar, ".txt", sep="")
	
	write.table(tableRT, save_as, sep="\t", row.names=F)

	# save paired
	tableRTPaired <- cast(pooledRT, Sub.Num~Response+Condition, mean)
	
	
	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_meanRT_by_subs_paired_", iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_n_", npar, ".txt", sep="")
		
	
	write.table(tableRTPaired, save_as, sep="\t", row.names=F)


### create graphs


meanRT <- aggregate(list(mean= tableRT$RT), list(Response = tableRT$Response, Condition = tableRT$Condition), mean, na.rm=T)
seRT <- aggregate(list(se= tableRT$RT), list(Response = tableRT$Response, Condition = tableRT$Condition), se)

plotRT <- cbind(meanRT, seRT["se"])

	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_meanRT_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_n_", npar, ".txt", sep="")
		
	
	write.table(plotRT, save_as, sep="\t", row.names=F)



dodge <- position_dodge(width=0.9)
barRT <- ggplot(plotRT, aes(Condition, mean, fill=Response)) + scale_fill_grey() + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=dodge, width=0.2) + labs(x = "Condition", y = "mean RT (ms)", fill = "Trial Type:") + labs(title="Mean RT by Condition") + theme(panel.background = element_blank()) + theme(axis.line = element_line())

		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_meanRT_barchart_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_lg_", longestLongestGap, "_n_", npar, ".pdf", sep="")

ggsave(save_as, barRT)


	return(tableRTPaired)
	}
		
poolData <- function(iChart, RejectFirstGap =FALSE, RejectLongestGap =FALSE, RejectRT =FALSE, color = FALSE, dependent = "Accuracy", group = "", facet = "", dodge="", xlab="", ylab="", paired=TRUE, miny = -1, maxy=-1, size=13, legend.direction = "vertical", legend.position="right", breaks=-1) {
	
	GoodFirstGap <- RejectFirstGap
	GoodLongestGap <- RejectLongestGap
	GoodRT <- RejectRT
	dodg <- dodge
	
	
		iChart <- iChart[iChart$Response == "D" | iChart$Response == "T",]
	
	dep <- which(names(iChart)==dependent)
	iChart$Accuracy <- iChart[,dep]
	
	
	groupName <- group
	facetName <- facet
	dodgName <- dodg
	
	if(group=="") {group <- FALSE} else {group <- as.factor(iChart[, group])}

	if(facet=="") {facet <- FALSE} else {facet <- as.factor(iChart[, facet])}		
	if(dodg == "") {dodg <- FALSE} else {dodg <- as.factor(iChart[, dodg])}
	

	# pool Accuracy across subjects
	
	filteredAccuracy <- iChart
		
	if(group[1]!=FALSE) filteredAccuracy <- cbind(filteredAccuracy, group)
	if(facet[1]!=FALSE) filteredAccuracy <- cbind(filteredAccuracy, facet)
	if(dodg[1]!=FALSE) filteredAccuracy <- cbind(filteredAccuracy, dodg)

		
	if(GoodFirstGap) filteredAccuracy <- filteredAccuracy[filteredAccuracy$GoodFirstGap | is.na(filteredAccuracy$GoodFirstGap),]
	if(GoodLongestGap) filteredAccuracy <- filteredAccuracy[filteredAccuracy$GoodLongestGap | is.na(filteredAccuracy$GoodLongestGap),]
	if(GoodRT) filteredAccuracy <- filteredAccuracy[filteredAccuracy$GoodRT,]

valuesBySub <- list(Sub.Num = filteredAccuracy$Sub.Num, Condition = filteredAccuracy$Condition)
if(group[1]!=FALSE) valuesBySub <- c(valuesBySub, list(group = filteredAccuracy$group))
if(facet[1]!=FALSE) valuesBySub <- c(valuesBySub, list(facet = filteredAccuracy$facet))
if(dodg[1]!=FALSE) valuesBySub <- c(valuesBySub, list(dodg = filteredAccuracy$dodg))

	sumAccuracy <- aggregate(list(n= filteredAccuracy$Accuracy), valuesBySub, function(value){sum(!is.na(value))})
	
	pooledAccuracy <- aggregate(list(Accuracy=as.numeric(filteredAccuracy$Accuracy)), valuesBySub, mean, na.rm=T)	
	sdAccuracy <- aggregate(list(sd=as.numeric(filteredAccuracy$Accuracy)), valuesBySub, sd)
	
	
	tableAccuracy <- cbind(pooledAccuracy, sumAccuracy["n"], sdAccuracy["sd"])
	
		if ("firstGap" %in% names(filteredAccuracy)) {longestFirstGap <- max(filteredAccuracy$firstGap, na.rm=T) } else {longestFirstGap <-NA }
		
		if ("longestGap" %in% names(filteredAccuracy)) {longestLongestGap <- max(filteredAccuracy$longestGap, na.rm=T) } else { longestLongestGap <- NA }
		
		longestRT_D <- max(filteredAccuracy$RT, na.rm=T)
		shortestRT_D <- min(filteredAccuracy$RT, na.rm=T)
		npar <- length(unique(iChart$Sub.Num))

	save_as_ta <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_mean_", dependent, "_by_subs_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_lg_", longestLongestGap, "_n_", npar, ".txt", sep="")

	
	# save paired

averaging <- "Sub.Num~Condition"
if(group[1]!=FALSE) averaging <- paste(averaging, "+group")
if(facet[1]!=FALSE) averaging <- paste(averaging, "+facet")
if(dodg[1]!=FALSE) averaging <- paste(averaging, "+dodg")


	tableAccuracyPaired_2 <- cast(pooledAccuracy,as.formula(averaging), mean,  value = 'Accuracy')
	
	tableAccuracyPairedse <- cast(sdAccuracy,as.formula(averaging), mean,  value = 'sd')
			
	tableAccuracyPairedNs <- cast(sumAccuracy,as.formula(averaging), mean, value = 'n')
	
	
	
	names(tableAccuracyPairedNs)[2:ncol(tableAccuracyPairedNs)] <- paste(names(tableAccuracyPairedNs)[2:ncol(tableAccuracyPairedNs)], "n", sep="_")
	
	names(tableAccuracyPairedse)[2:ncol(tableAccuracyPairedse)] <- paste(names(tableAccuracyPairedse)[2:ncol(tableAccuracyPairedse)], "sd", sep="_")	
	
	tableAccuracyPaired <- merge(tableAccuracyPaired_2, tableAccuracyPairedse, by="Sub.Num")
	
	tableAccuracyPaired <- merge(tableAccuracyPaired, tableAccuracyPairedNs, by="Sub.Num")	
	
			
		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_mean_", dependent, "_by_subs_paired_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_lg_", longestLongestGap, "_n_", npar, ".txt", sep="")
			
	### create graphs

values <- list(Condition = tableAccuracy$Condition)
if(group[1]!=FALSE) values <- c(values, list(group = tableAccuracy$group))
if(facet[1]!=FALSE) values <- c(values, list(facet = tableAccuracy$facet))
if(dodg[1]!=FALSE) values <- c(values, list(dodg = tableAccuracy$dodg))


meanAcc <- aggregate(list(mean= tableAccuracy$Accuracy), values, mean, na.rm=T)
seAcc <- aggregate(list(se= tableAccuracy$Accuracy), values, se)	
plotAcc <- cbind(meanAcc, seAcc["se"])

		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_mean_", dependent, "_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_lg_", longestLongestGap, "_n_", npar, ".txt", sep="")
		



write.table(plotAcc, save_as, sep="\t", row.names=F)

dodge <- position_dodge(width=0.9)

if(dodg[1]==FALSE){
	
		if(group[1]==FALSE) {
barAcc <- ggplot(plotAcc, aes(Condition, mean, fill=Condition)) + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=dodge, width=0.2) + labs(x = xlab, y = ylab, fill = "Condition:") + labs(title="") + theme(panel.background = element_blank()) + theme(axis.line = element_line())
	} else {
barAcc <- ggplot(plotAcc, aes(as.factor(group), mean, fill=Condition)) + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=dodge, width=0.2) + labs(x = xlab, y = ylab, fill = "Condition:") + labs(title="") + theme(panel.background = element_blank()) + theme(axis.line = element_line())
		}	
	
	
	} else {
		
		if(group[1]==FALSE) {
barAcc <- ggplot(plotAcc, aes(Condition, mean, fill=dodg)) + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=dodge, width=0.2) + labs(x = xlab, y = ylab, fill = "Condition:") + labs(title="") + theme(panel.background = element_blank()) + theme(axis.line = element_line())
	} else {
barAcc <- ggplot(plotAcc, aes(as.factor(group), mean, fill=dodg)) + geom_bar(position=dodge, stat="identity") + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=dodge, width=0.2) + labs(x = xlab, y = ylab, fill = "Condition:") + labs(title="") + theme(panel.background = element_blank()) + theme(axis.line = element_line())
		}
		
		}




barAcc <- barAcc + geom_abline(intercept = 0.5, slope = 0, size=0.1, linetype = 2)  + theme(axis.title.x = element_text(vjust = -0.5, size=size))  + theme(axis.title.y = element_text(hjust=-0.5, vjust=0.5, size=size)) + theme(axis.line = element_line()) + theme(strip.text.y = element_text(angle = 30))+ theme(strip.text.x = element_text(angle = 0, size=size))+theme(axis.text.x=element_text(size=size, vjust=1))+theme(axis.text.y=element_text(size=size-2, hjust=1, colour="black"))

barAcc <- barAcc + theme(legend.position="top")+ scale_fill_discrete("")

barAcc <- barAcc + theme(legend.direction = legend.direction, legend.position=legend.position) + theme(legend.text= element_text(size = size))

if(maxy != -1) barAcc <- barAcc + coord_cartesian(ylim = c(miny, maxy))

if(breaks != -1) barAcc <- barAcc + scale_y_continuous(breaks = breaks)

if(!color) barAcc <- barAcc + scale_fill_grey()
if(facet[1]!=FALSE) barAcc <- barAcc + facet_grid(. ~ facet)

		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_mean_", dependent, "_barchart_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_lg_", longestLongestGap, "_n_", npar, ".pdf", sep="")
		
ggsave(save_as, barAcc)

	

indexGroup <- which(names(tableAccuracy)=="group")
if("group" %in% names(tableAccuracy)) {names(tableAccuracy)[indexGroup] <- groupName} 

indexFacet <- which(names(tableAccuracy)=="facet")
if("facet" %in% names(tableAccuracy)) {names(tableAccuracy)[indexFacet] <- facetName} 

indexDodg <- which(names(tableAccuracy)=="dodg")
if("dodg" %in% names(tableAccuracy)) {names(tableAccuracy)[indexDodg] <- dodgName} 


tableAccuracyPaired <- tableAccuracyPaired[,c(1,order(names(tableAccuracyPaired)[2:ncol(tableAccuracyPaired)])+1)]

if(paired==TRUE) {
		write.table(tableAccuracyPaired, save_as_ta, sep="\t", row.names=F)

	return(tableAccuracyPaired_2)
	} else {
		
	write.table(tableAccuracy, save_as_ta, sep="\t", row.names=F)

	return(tableAccuracy)
		
		}

	}

plotData <- function(variable, iChart_quantile, label) {
	
	variable <- variable[!is.na(variable)]

	h = hist(variable, plot = FALSE)
	d = density(variable)
	
	xlimit = range(h$breaks, d$x)
	ylimit = range(0, h$density, d$y)

	hist(variable, freq=FALSE, xlim=xlimit, ylim=ylimit, main="", xlab=label, ylab="", yaxt = "n", col="lightgrey", border="darkgrey")
	lines(d$x, d$y)
	abline(v = quantile(variable, iChart_quantile), lty = 2)
	
	}

getSummaryStats <- function(iChart, iChart_quantile) {
	
	
	npar <- length(unique(iChart$Sub.Num))
	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_descriptives_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_n_", npar, sep="")
	
	
	pdf(paste(save_as, ".pdf", sep = ""))
	
	par(mfrow = c(3, 2))
	
	goodRTs <- (iChart$Response=="D" | iChart$Response=="T") & !is.na(iChart$RT)
	Dtrials <- (iChart$Response=="D") & !is.na(iChart$RT)
	Ttrials <- (iChart$Response=="T") & !is.na(iChart$RT)

	
	plotData(iChart$RT[Dtrials], iChart_quantile, "RT (ms) D trials")
	plotData(iChart$RT[Ttrials], iChart_quantile, "RT (ms) T trials")
	
	plotData(iChart$firstGap[goodRTs], iChart_quantile, "Duration of first Gap (frames)")
	plotData(iChart$durationLongestLook[goodRTs], iChart_quantile, "Duration of longest fixation (frames)")
	
	goodLongestGap <- ((iChart$Response=="D" | iChart$Response=="T") & !is.na(iChart$longestGap))
	plotData(iChart$longestGap[goodLongestGap], iChart_quantile, "Duration of Longest Gap (frames)")
	
	goodLongestGapPosition <- ((iChart$Response=="D" | iChart$Response=="T") & !is.na(iChart$longestGapPosition))
	plotData(iChart$longestGapPosition[goodLongestGapPosition], 0.5, "Start of Longest Gap (ms)")

	sink(paste(save_as, ".txt", sep = ""))
	
	quantileRT <- quantile(iChart$RT[goodRTs], iChart_quantile, na.rm=T)
	print("")
	print(paste(iChart_quantile, "quantile of RT for D and T trials =", round(quantileRT,0)))
	print("Summary RT all Trials")
	print(summary(iChart$RT[goodRTs]))
	print("trials used")
	print(sum(!is.na(iChart$RT[goodRTs])))
	
	quantileRT_D <- quantile(iChart$RT[Dtrials], iChart_quantile, na.rm=T)
	print("")
	print(paste(iChart_quantile, "quantile of RT for D trials =", round(quantileRT_D,0)))
	print("Summary RT D Trials")
	print(summary(iChart$RT[Dtrials]))
	print("trials used")
	print(sum(!is.na(iChart$RT[Dtrials])))

	
	quantileRT_T <- quantile(iChart$RT[Ttrials], iChart_quantile, na.rm=T)
	print("")
	print(paste(iChart_quantile, "quantile of RT for T trials =", round(quantileRT_T, 0)))
	print("Summary RT T Trials")
	print(summary(iChart$RT[Ttrials]))
	print("trials used")
	print(sum(!is.na(iChart$RT[Ttrials])))
	
	quantileFirstGap <- quantile(iChart$firstGap[goodRTs], iChart_quantile, na.rm=T)
	print("")
	print(paste(iChart_quantile, "quantile of First GAPS =", round(quantileFirstGap, 0)))
	print("Summary first GAP")
	print(summary(iChart$firstGap[goodRTs]))
	print("trials used")
	print(sum(!is.na(iChart$firstGap[goodRTs])))

	quantileLongestFixation <- quantile(iChart$durationLongestLook[goodRTs], iChart_quantile, na.rm=T)
	print("")
	print(paste(iChart_quantile, "quantile of longest fixation =", round(quantileLongestFixation, 0)))
	print("Summary longest Fixation")
	print(summary(iChart$durationLongestLook[goodRTs]))
	print("trials used")
	print(sum(!is.na(iChart$durationLongestLook[goodRTs])))

	print("")
	quantileLongestGap_no_0 <- quantile(iChart$longestGap[goodLongestGap & iChart$longestGap != 0], iChart_quantile, na.rm=T)
	print("")
	print(paste(quantileLongestGap_no_0, "quantile of longest GAPS (without 0) =", round(quantileLongestGap_no_0, 0)))
	print("Summary longest GAP (without 0)")
	print(summary(iChart$longestGap[goodLongestGap & iChart$longestGap != 0]))
	print("trials used")
	print(sum(!is.na(iChart$longestGap[goodLongestGap & iChart$longestGap != 0])))

	print("")
	quantileLongestGap <- quantile(iChart$longestGap[goodLongestGap], iChart_quantile, na.rm=T)
	print("")
	print(paste(iChart_quantile, "quantile of longest GAPS =", round(quantileLongestGap, 0)))
	print("Summary longest GAP")
	print(summary(iChart$longestGap[goodLongestGap]))
	print("trials used")
	print(sum(!is.na(iChart$longestGap[goodLongestGap])))
		
	quantileLongestGapPosition <- quantile(iChart$longestGapPosition[goodLongestGapPosition], iChart_quantile, na.rm=T)
	print("")
	print(paste(iChart_quantile, "quantile of start of longest GAPS =", round(quantileLongestGapPosition, 0)))
	print("Summary start of longest GAP")
	print(summary(iChart$longestGapPosition[goodLongestGapPosition]))
	print("trials used")
	print(sum(!is.na(iChart$longestGapPosition[goodLongestGapPosition])))
	
	
	# check how many trials per condition each participant has 
	trials_per_condition <- aggregate(list(Trials = iChart$Tr.Num), list(Months = iChart$Months, Sub.Num = iChart$Sub.Num, Condition = iChart$Condition), length)
	
	
		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_descriptives_trials_per_condition_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_n_", npar, ".txt", sep="")
		
	write.table(trials_per_condition[order(trials_per_condition$Trials),], save_as, sep="\t", row.names=F)

	# check how many trials each participant has 
	trials_per_subject <- aggregate(list(Trials = iChart$Tr.Num), list(Months = iChart$Months, Sub.Num = iChart$Sub.Num), length)
	
			save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_descriptives_trials_by_subs_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_n_", npar, ".txt", sep="")
			
	write.table(trials_per_subject[order(trials_per_subject$Trials),], save_as, sep="\t", row.names=F)
	
	sink()
	graphics.off()
	
	}

readCDI <- function() {
	
		CDI <- read.table(file.choose(), header=T, fill=T, as.is=T, sep="\t")
		return(CDI)
	}
	
defineOnset <- function(iChart, critonset, includeAways=FALSE) {
	
		F0 <- critonset
		F0 <- which(names(iChart)==F0)
	iChart$Response <- ifelse(iChart[,F0] == "1", "T", ifelse(iChart[,F0] == "0", "D", "A"))
	
	
	if(includeAways) 	iChart$Response <- ifelse(iChart[,F0] == "1", "T", ifelse(iChart[,F0] == "0", "D", "T"))

	return(iChart)
	
	}

## includeAways == FALSE -> only include trials child was looking at center at F0
## includeOffCenter == TRUE -> include trials child was looking at center, target, or distractor at F0 
defineOnsetSOL <- function(iChart, critonset, includeOffCenter=FALSE) {
      
      F0 <- critonset
      F0 <- which(names(iChart)==F0)
      iChart$Response <- ifelse(iChart[,F0] == ".5" | iChart[,F0] == 0.5, "D", 
                                ifelse(iChart[,F0] == "", "A", "A")) 
      
      if(includeOffCenter) iChart$Response <- ifelse(iChart[,F0] == ".5" | iChart[,F0] == "0.5", "D", 
                                                 ifelse(iChart[,F0] == "1" | 
                                                              iChart[,F0] == "0", "D", "A")) 
      
      return(iChart)
      
}

re.readiChart <- function(iChartFile = file.choose()) {
	
	iChart <- read.table(iChartFile, header=T, fill=T, as.is=T, sep="\t", check.names=FALSE)
	return(iChart)
	
	}

readiChart <- function(iChartFile=file.choose()) {
			
	#temp <- gregexpr("/", iChartFile)
	
	#last_separator <- temp[[1]][length(temp[[1]])]

	StudyName <- substr(basename(iChartFile), 1, nchar(iChartFile)-4)

	Directory <- paste(dirname(iChartFile),"/" , sep="")	

	iChart <- read.table(iChartFile, header=T, fill=T, as.is=T, sep="\t")
				
	iChart$RTOld <- iChart$RT
	
	startWindow <- match("F0", colnames(iChart))
	if(is.na(startWindow)) startWindow <- match("Word.Onset.Frame", colnames(iChart))

	# rename frame columns
	colnames(iChart)[startWindow:length(iChart)] <- round(framesToMs(round(0:(length(iChart)-startWindow))))
	

	# the original file had one header per subject
	# this line removes all rows that are the same as the header
	# and selects only the columns of interest
	iChart <- iChart[iChart[,"Months"]!="Months",]
	

	iChart$Block <- rep(NA, nrow(iChart))
	iChart$Offset <- rep(0, nrow(iChart))
	iChart$StartWindowRT <- rep(NA, nrow(iChart))
	iChart$EndWindowRT <- rep(NA, nrow(iChart))
	iChart$StartWindowAcc <- rep(NA, nrow(iChart))
	iChart$EndWindowAcc <- rep(NA, nrow(iChart))
	iChart$StudyName <- rep(StudyName, nrow(iChart))
	iChart$Directory <- rep(Directory, nrow(iChart))

	iChart <- iChart[iChart$Condition != "",]

	return(iChart)
	
	}


computeStatistics <- function(iChart, startWindow, endWindow) {


	iChart$StartWindowRT <- rep(startWindow, nrow(iChart))
	iChart$EndWindowRT <- rep(endWindow, nrow(iChart))
	iChart$StartWindowAcc <- rep(startWindow, nrow(iChart))
	iChart$EndWindowAcc <- rep(endWindow, nrow(iChart))
	
	startWindow <- which(names(iChart)==startWindow)
	
	if (endWindow == FALSE) endWindow <- length(iChart)
	else endWindow <- which(names(iChart)==endWindow)
	
	iChart$longestLook <- rep(NA, nrow(iChart))
	iChart$longestGap <- rep(NA, nrow(iChart))
	iChart$Gap <- rep(NA, nrow(iChart))
	iChart$longestGapPosition <- rep(NA, nrow(iChart))
	iChart$firstGap <- rep(NA, nrow(iChart))
	iChart$RT <- rep(NA, nrow(iChart))
	iChart$durationLongestLook <- rep(NA, nrow(iChart))
	iChart$StartlongestLook_T <- rep(NA, nrow(iChart))
	iChart$NumberOfShifts <- rep(NA, nrow(iChart))
	iChart$FirstFixation <- rep(NA, nrow(iChart))
	
	endRows <- nrow(iChart)
	print("### Trials left ###")

	
	for (row in 1:endRows) {
		
		if((endRows-row)%%50 == 0)print(endRows-row)
		consecutiveAways <- 0
		containsAway <- FALSE
		currentLook <- 0
		longestLook <- 0
		firstGap <- 0
		longestGap <- 0
		longestGap_col <- NA
		Gap <- 0
		longestLook_T <- 0
		StartlongestLook_T <- NA
		NumberOfShifts <- 0
		FirstFixation <- 0
		isFirst <- TRUE
		isShift <- TRUE
	
	    start_elem <- iChart[row, startWindow]
	    	
	    rt_notdetected <- TRUE
		for (col in startWindow:endWindow) {
			curr_elem <- iChart[row, col]
			
			if (is.na(curr_elem) || curr_elem == ""){
			col = endWindow
			} else if (curr_elem == "." || curr_elem == "-") {
			
				currentLook <- 0
				
				if(isShift) {
					NumberOfShifts <- NumberOfShifts + 1
					isShift <- FALSE
				}

				if(FirstFixation > 0) {
					isFirst <- FALSE
				}
				
				
				if(consecutiveAways==0) rt <- col
				consecutiveAways <- consecutiveAways + 1
				Gap <- Gap + 1
			
				if(!containsAway && (curr_elem == "." || curr_elem == ".")) containsAway <- TRUE
				
				if (containsAway && (consecutiveAways > longestGap)) {
					longestGap_col <- col
					longestGap <- consecutiveAways
				}
			
			
			} else if (curr_elem== "1" || curr_elem == "0") {
			
				isShift <- TRUE
				

				#define longest look
				currentLook <- currentLook + 1
				
					
				if (currentLook > longestLook) {
						longestLook <- currentLook
						longestLook_elem <- curr_elem
					}
    		if(rt_notdetected){
    			if(curr_elem != start_elem && iChart$Response[row] != "R") {
    					iChart$RT[row] <- as.numeric(colnames(iChart)[rt])
    					iChart$firstGap[row] <- consecutiveAways
    					rt_notdetected <- FALSE	
    				}
				
				}
				
				if(!rt_notdetected & isFirst){
					FirstFixation <- currentLook
				}
				
				if(!rt_notdetected && (currentLook > longestLook_T) && curr_elem == "1") {
						longestLook_T <- currentLook
						StartlongestLook_T <- rt
					}
				
				consecutiveAways <- 0
				containsAway <- FALSE

			} # if it is an away or off trial
		} # col
		   iChart$durationLongestLook[row] <- longestLook
		   iChart$longestLook[row] <- longestLook_elem
			iChart$longestGap[row] <- longestGap
			iChart$Gap[row] <- Gap
			
			iChart$NumberOfShifts[row] <- NumberOfShifts
			
			iChart$FirstFixation[row] <- FirstFixation
			
			
			if(!is.na(StartlongestLook_T)) iChart$StartlongestLook_T[row] <- as.numeric(colnames(iChart)[StartlongestLook_T])
			if(!is.na(longestGap_col))iChart$longestGapPosition[row] <- as.numeric(colnames(iChart)[longestGap_col-longestGap+1])
	   }
	

	iChart$Accuracy <- rowMeans(data.matrix(iChart[,startWindow:endWindow]), na.rm=T)
		
		longestFirstGap <- max(iChart$firstGap[iChart$Response == "D" | iChart$Response == "T"], na.rm=T)
		longestLongestGap <- max(iChart$longestGap[iChart$Response == "D" | iChart$Response == "T"], na.rm=T)
		longestRT_D <- max(iChart$RT[iChart$Response == "D"], na.rm=T)
		shortestRT_D <- min(iChart$RT[iChart$Response == "D"], na.rm=T)
		
		
	npar <- length(unique(iChart$Sub.Num))
	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_originaliChart", "_RT_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_Acc_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_n_", npar, ".txt", sep="")
	
	write.table(iChart, save_as, sep="\t", row.names=F)
	
	return(iChart)
	
	
	}
	
	
	meanAccuracy <- function(iChart, startWindowAcc, endWindowAcc) {
		
		startWindow <- startWindowAcc
		endWindow <- endWindowAcc
		
		iChart$StartWindowAcc <- rep(startWindow, nrow(iChart))
		iChart$EndWindowAcc <- rep(endWindow, nrow(iChart))
	
		startWindow <- which(names(iChart)==startWindow)
		endWindow <- which(names(iChart)==endWindow)
	
			
			iChart$Accuracy <- rowMeans(data.matrix(iChart[,startWindow:endWindow]), na.rm=T)
	
			longestFirstGap <- max(iChart$firstGap[iChart$Response == "D" | iChart$Response == "T"], na.rm=T)
		longestLongestGap <- max(iChart$longestGap[iChart$Response == "D" | iChart$Response == "T"], na.rm=T)
		longestRT_D <- max(iChart$RT[iChart$Response == "D"], na.rm=T)
		shortestRT_D <- min(iChart$RT[iChart$Response == "D"], na.rm=T)
		
		
	npar <- length(unique(iChart$Sub.Num))
	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_filterediChart", "_RT_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_Acc_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_n_", npar, ".txt", sep="")
	
	write.table(iChart, save_as, sep="\t", row.names=F)
	return(iChart)
		
		}
	
	
	
	

filteriChart <- function(iChart, minRT, maxRT, maxfirstgap, maxlonggap) {
	
	high_percentileRT <- maxRT
	low_percentileRT <- minRT
	percentileFirstGap <- maxfirstgap
	percentileLongestGap <- maxlonggap

	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_filtering_Criteria.txt", sep="")
	
	sink(paste(save_as, sep=""))
	
	
	goodTrials <- iChart$Response == "D" | iChart$Response == "T"
	T_Trials <- iChart$Response == "T"
	D_Trials <- iChart$Response == "D"

	# first Gap
	if(percentileFirstGap < 1) {
		longestFirstGap <- quantile(iChart$firstGap[goodTrials], percentileFirstGap, na.rm=T)
	print("")
		print(paste(percentileFirstGap, "quantile of First Gaps =", longestFirstGap))
		} else {
		longestFirstGap <- percentileFirstGap
			}

	print("First GAP")
	print(paste("Longest First Gap = ", longestFirstGap))
	print(summary(iChart$firstGap[goodTrials]))

	# longest Gap
	if(percentileLongestGap < 1) {
		longestLongestGap <- quantile(iChart$longestGap[goodTrials], percentileLongestGap, na.rm=T)
	print("")
		print(paste(percentileLongestGap, "quantile of Longest Gaps =", longestLongestGap))
		} else {
		longestLongestGap <- percentileLongestGap
			}
	
	print("Longest GAP")
	print(paste("Longest Gap = ", longestLongestGap))
	print(summary(iChart$longestGap[goodTrials]))

	# RTs
	
	if(low_percentileRT < 1) {
		shortestRT_D <- quantile(iChart$RT[D_Trials], low_percentileRT, na.rm=T)
		print(paste(low_percentileRT, "quantile of RTs (D) =", round(shortestRT_D)))

		shortestRT_T <- quantile(iChart$RT[T_Trials], low_percentileRT, na.rm=T)
		print(paste(low_percentileRT, "quantile of RTs (T) =", round(shortestRT_T)))

		} else {
			shortestRT_D <- low_percentileRT
			shortestRT_T <- low_percentileRT
			
			}
	
	if(high_percentileRT < 1) {
		longestRT_D <- quantile(iChart$RT[D_Trials], high_percentileRT, na.rm=T)
		print("")
		print(paste(high_percentileRT, "quantile of RTs (D) =", round(longestRT_D)))
		
		longestRT_T <- quantile(iChart$RT[T_Trials], high_percentileRT, na.rm=T)
		print("")
		print(paste(high_percentileRT, "quantile of RTs (T) =", round(longestRT_T)))
		
		} else {
			
			longestRT_D <- high_percentileRT
			longestRT_T <- high_percentileRT
			
			}

	print("RT D Trials (filtering based on these statistics)")
	print(paste("Longest RT D trials: ", longestRT_D))
	print(paste("Shortest RT D trials: ", shortestRT_D))
	print(summary(iChart$RT[D_Trials]))
	
	print("RT T Trials")
	print(paste("Longest RT T trials: ", longestRT_T))
	print(paste("Shortest RT T trials: ", shortestRT_T))
	print(summary(iChart$RT[T_Trials]))	

	iChart$GoodFirstGap <- iChart$firstGap <= longestFirstGap
	
	print("Trials Rejected")
	print("")
	print(paste("First Gap = ", sum(iChart$GoodFirstGap,na.rm=T), " out of ", sum(iChart$GoodFirstGap,na.rm=T)+sum(!iChart$GoodFirstGap,na.rm=T), sep=""))
	
	
	iChart$GoodRT <- (iChart$RT <= longestRT_D)&(iChart$RT >= shortestRT_D)
	
	print("")
	print(paste("RT = ", sum(iChart$GoodRT,na.rm=T), " out of ", sum(iChart$GoodRT,na.rm=T)+sum(!iChart$GoodRT,na.rm=T),sep=""))
	
	
	iChart$GoodLongestGap <- iChart$longestGap <= longestLongestGap
	
	print("")
	print(paste("Longest Gap = ", sum(iChart$GoodLongestGap,na.rm=T), " out of ", sum(iChart$GoodLongestGap,na.rm=T)+sum(!iChart$GoodLongestGap, na.rm=T), sep=""))
	
	sink()

		
	npar <- length(unique(iChart$Sub.Num))
	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_filterediChart", "_RT_",iChart[1, "StartWindowRT"], "_", iChart[1, "EndWindowRT"], "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_Acc_", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], "_n_", npar, ".txt", sep="")

	
	write.table(iChart, save_as, sep="\t", row.names=F)

	
	return(iChart)
	
	}	
	
	
	createPlots <- function(iChart, startWindow, endWindow, RejectLongestGap, RejectFirstGap, RejectRT, color = FALSE, smooth=33, targetEnd = 800, carrier = "", targets = "", group="", plotStats = "OC", miny = -1, maxy = -1, legendPosition = FALSE, size=13, legend.direction = "vertical", legend.position="right", breaks=-1, x.target=0.30) {
				
				GoodLongestGap <- RejectLongestGap
				GoodFirstGap <- RejectFirstGap
				GoodRT <- RejectRT
				endTarget <- targetEnd
				
		if(group == "") {group <- FALSE} else {group <- iChart[,group]}
		
		startW <- startWindow
				

		if(group[1]!=FALSE) iChart <- cbind(iChart, group)
		
		if(GoodLongestGap) iChart <- iChart[iChart$GoodLongestGap | is.na(iChart$GoodLongestGap),]
		if(GoodFirstGap) iChart <- iChart[iChart$GoodFirstGap | is.na(iChart$GoodFirstGap),]
		if(GoodRT) iChart <- iChart[iChart$GoodRT,]

		
		sw <- startWindow
		ew <- endWindow
		
		onsetBy <- startWindow
		startWindow <- which(names(iChart)==startWindow)
		endWindow <- which(names(iChart)==endWindow)
	
		
		iChart <- iChart[iChart$Response == "D" | iChart$Response == "T", ]
		
		if(plotStats == "OC_D") iChart <- iChart[iChart$Response == "D", ]
		if(plotStats == "OC_T") iChart <- iChart[iChart$Response == "T", ]
		
		
		longestFirstGap <- max(iChart$firstGap, na.rm=T)
		longestLongestGap <- max(iChart$longestGap, na.rm=T)
		longestRT_D <- max(iChart$RT, na.rm=T)
		shortestRT_D <- min(iChart$RT, na.rm=T)
		
		
	npar <- length(unique(iChart$Sub.Num))	
		
				
		# Profile
		
		value <- list(Sub.Num = iChart$Sub.Num, Condition = iChart$Condition)
		if(group[1]!=FALSE) value <- c(value, list(groupping=iChart$group))
		
		if(plotStats=="OC") value <- c(value, list(Response=iChart$Response))
		
		
		pooledPlotAccuracy <- aggregate(data.matrix(iChart[,startWindow:endWindow]), value, mean, na.rm=T)
		
		startPlot <- which(names(pooledPlotAccuracy)==sw)
		endPlot <- which(names(pooledPlotAccuracy)==ew)
		
		if(plotStats=="OC") pooledPlotAccuracy[pooledPlotAccuracy$Response == "T",startPlot:endPlot] <- 1 - pooledPlotAccuracy[pooledPlotAccuracy$Response == "T",startPlot:endPlot]
		
		if(plotStats=="OC_T") pooledPlotAccuracy[,startPlot:endPlot] <- 1 - pooledPlotAccuracy[,startPlot:endPlot]

		
		
		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_", plotStats,"_graphValues_by_subs_",sw, "_", ew, "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_n_", npar, ".txt", sep="")
		
		if(plotStats %in% c("PP","OC")) {

		write.table(pooledPlotAccuracy, save_as, sep="\t", row.names=F)
		
		}
		
		value <- list(Condition = pooledPlotAccuracy$Condition)
		if(group[1]!=FALSE) value <- c(value, list(groupping= pooledPlotAccuracy$group))
		if(plotStats=="OC") value <- c(value, list(Response= pooledPlotAccuracy$Response))

		
		startPlot <- which(names(pooledPlotAccuracy)==sw)
		endPlot <- which(names(pooledPlotAccuracy)==ew)
		
		plotProfilemean <- aggregate(data.matrix(pooledPlotAccuracy[,startPlot:endPlot]), value, mean, na.rm=T)
		plotProfilemean <- cbind(statistic = rep("mean", length(plotProfilemean$Condition)), plotProfilemean)

		plotProfilese <- aggregate(data.matrix(pooledPlotAccuracy[,startPlot:endPlot]), value, function(value) { sd(value, na.rm=TRUE)/sqrt(sum(!is.na(value)))})

		plotProfilese <- cbind(statistic = rep("se", length(plotProfilese$Condition)), plotProfilese)

		plotProfile <- rbind(plotProfilemean, plotProfilese)


		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_", plotStats, "_graphValues_",sw, "_", ew, "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_n_", npar, ".txt", sep="")
		
		if(plotStats %in% c("PP","OC")) {

		write.table(plotProfile, save_as, sep="\t", row.names=F)
		
		}

idvalue <- c("Condition")
if(group[1]!=FALSE) idvalue <- c("Condition", "groupping")
if(plotStats == "OC") idvalue <- c(idvalue, "Response")

		tablePlotProfile = melt(plotProfilemean[,2:endPlot], id=idvalue)
		tablePlotProfilese = melt(plotProfilese[,2:endPlot], id=idvalue)
		tablePlotProfile$ucl <- tablePlotProfile$value + tablePlotProfilese$value
		tablePlotProfile$lcl <- tablePlotProfile$value - tablePlotProfilese$value
		
		tablePlotProfile$variable <- as.numeric(as.character(tablePlotProfile$variable))
		
		
		
		if (smooth > 33) {
			tablePlotProfile$lcl[tablePlotProfile$variable %% smooth != 0] <- NA
			tablePlotProfile$ucl[tablePlotProfile$variable %% smooth != 0] <- NA

			}

if(plotStats == "OC") {
				
				
			plotMeansAccuracy <- qplot(variable, value, data = tablePlotProfile, shape= Response, colour=Condition, group=interaction(Condition, Response),geom="blank") + geom_errorbar(aes(ymin = lcl, ymax = ucl, fill=Condition), data=tablePlotProfile, stat="identity", width=100) + xlab("Time (ms) from noun onset") + labs(colour = "", shape = "", fill = "") + ylab("Proportion\n  Looking\n  to target") + theme(panel.background = element_blank()) + theme(axis.title.x = element_text(vjust = -0.5, size=13))  + theme(axis.title.y = element_text(hjust=-0.5, vjust=0.5, size=13)) + theme(panel.grid.minor = element_blank()) + theme(axis.line = element_line()) + theme(plot.title = element_text(vjust = 1.3, size=18)) + theme(plot.margin = unit(c(1, 1, 1, 3), "lines")) + theme(legend.key = element_blank()) + geom_line(size = 0.8) + theme(strip.text.y = element_text(angle = -90, size=13)) + theme(legend.text = element_text(size=12))+geom_point(aes(x=variable, y=value, colour=Condition, shape=Response), size=4, data= tablePlotProfile[!is.na(tablePlotProfile$ucl),]) + geom_vline(xintercept = startW, linetype = "dashed", size=0.3)+ geom_vline(xintercept = endTarget, linetype = "dashed", size=0.3)			
	
	
	} else {
		
				if (length(unique(iChart$Condition)) > 6) {
			
			plotMeansAccuracy <- qplot(variable, value, data = tablePlotProfile, colour=Condition, group=Condition,geom="line") + geom_errorbar(aes(ymin = lcl, ymax = ucl, fill=Condition), data=tablePlotProfile, stat="identity", width=100) + xlab("Time (ms) from noun onset") + labs(colour = "", shape = "", fill = "") + ylab("Proportion\n  Looking\n  to target") + theme(panel.background = element_blank()) + theme(axis.title.x = element_text(vjust = -0.5, size=13))  + theme(axis.title.y = element_text(hjust=-0.5, vjust=0.5, size=13)) + theme(panel.grid.minor = element_blank()) + theme(axis.line = element_line()) + theme(plot.title = element_text(vjust = 1.3, size=18)) + theme(plot.margin = unit(c(1, 1, 1, 3), "lines")) + theme(legend.key = element_blank()) + geom_line(size = 0.8)  + theme(strip.text.y = element_text(angle = -90, size=13)) + theme(legend.text = element_text(size=12))+geom_point(aes(x=variable, y=value, colour=Condition), size=4, data= tablePlotProfile[!is.na(tablePlotProfile$ucl),]) + geom_vline(xintercept = startW, linetype = "dashed", size=0.3)+ geom_vline(xintercept = endTarget, linetype = "dashed", size=0.3)		
				
			} else {				
				
			plotMeansAccuracy <- qplot(variable, value, data = tablePlotProfile, colour=Condition, group=Condition, geom="line") + geom_errorbar(aes(ymin = lcl, ymax = ucl, fill=Condition), data=tablePlotProfile, stat="identity", width=100) + xlab("Time (ms) from picture onset") + labs(colour = "", shape = "", fill = "") + ylab("Proportion looking to target") + theme(panel.background = element_blank()) + theme(axis.title.x = element_text(vjust = -0.5, size=13))  + theme(axis.title.y = element_text(hjust=-0.5, vjust=0.5, size=13)) + theme(panel.grid.minor = element_blank()) + theme(axis.line = element_line()) + theme(plot.title = element_text(vjust = 1.3, size=18)) + theme(plot.margin = unit(c(1, 1, 1, 3), "lines")) + theme(legend.key = element_blank()) + geom_line(size = 0.8) + theme(strip.text.y = element_text(angle = -90, size=13)) + theme(legend.text = element_text(size=12))+geom_point(aes(x=variable, y=value, colour=Condition, shape=Condition), size=4, data= tablePlotProfile[!is.na(tablePlotProfile$ucl),]) + geom_vline(xintercept = startW, linetype = "dashed", size=0.3)+ geom_vline(xintercept = endTarget, linetype = "dashed", size=0.3)
							}

		
		
		}
				
							
				

if(group[1]!=FALSE) plotMeansAccuracy <- plotMeansAccuracy + facet_grid(groupping~.) + theme(panel.background = element_rect()) + theme(axis.line = element_blank()) + theme(axis.title.y = element_text(hjust=-0.5, vjust=0.8, size=13))

if(legendPosition == TRUE) plotMeansAccuracy <-  plotMeansAccuracy + theme(legend.position=c(0.86,0.90))


if(plotStats == "PP") plotMeansAccuracy <-  plotMeansAccuracy +  geom_hline(yintercept = 0.5, linetype = "dashed", size=0.3)

if(breaks != -1) plotMeansAccuracy <- plotMeansAccuracy + scale_y_continuous(breaks = breaks)

if(maxy != -1) plotMeansAccuracy <-  plotMeansAccuracy + coord_cartesian(ylim = c(miny, maxy))


if(!color) {plotMeansAccuracy <- plotMeansAccuracy + scale_colour_grey()} else { plotMeansAccuracy <- plotMeansAccuracy + scale_colour_hue() 
}

plotMeansAccuracy <- plotMeansAccuracy + theme(axis.title.x = element_text(vjust = -0.5, size=size))  + theme(axis.title.y = element_text(hjust=-0.5, vjust=0.5, size=size)) + theme(axis.line = element_line()) + theme(strip.text.y = element_text(angle = 30))+ theme(strip.text.x = element_text(angle = 0, size=size))+theme(axis.text.x=element_text(size=size-2, vjust=1))+theme(axis.text.y=element_text(size=size-2, hjust=1, colour="black"))

plotMeansAccuracy <- plotMeansAccuracy + theme(legend.direction = legend.direction, legend.position=legend.position) + theme(legend.text= element_text(size = size))

plotMeansAccuracy <- plotMeansAccuracy + theme(legend.text= element_text(size = size))

if(group!="")plotMeansAccuracy <- plotMeansAccuracy +  theme(strip.text.y = element_text(size = size, angle = 90))+theme(strip.text.x = element_text(angle = 0,size=size))+theme(axis.title.y = element_text(size=size, angle = 90, vjust=-0.2, hjust= - 0.1))



			save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_", plotStats, "_plot_",sw, "_", ew, "_minRT_",  shortestRT_D, "_maxRT_", longestRT_D, "_lg_", longestLongestGap, "_fg_", longestFirstGap, "_n_", npar, ".pdf", sep="")
			
			
pdf(file=save_as, width=10, height=5)
print(plotMeansAccuracy)

x <- 0.07
y <- 0.93
grid.text(carrier, x=x, y=y, gp=gpar(fontsize=21, col="black", fontface="italic"))


x <- x.target
y <- 0.99
for(i_word in 1:length(targets)) {
	
	y <- y - 0.06
	grid.text(targets[i_word], x=x, y=y, gp=gpar(fontsize=21, col="black"))
	
		if(i_word > 1){grid.text(targets[i_word], x=x, y=y, gp=gpar(fontsize=21, col="grey"))} else {grid.text(targets[i_word], x=x, y=y, gp=gpar(fontsize=21, col="black"))}

	
}	

dev.off()

		
		#return(tablePlotProfile)
		
		}



plotPercentageAways <- function(iChart, startWindow, endWindow, save_as) {
	
	pdf(paste(save_as, ".pdf", sep = ""))

	temp <- iChart[,which(names(iChart)==startWindow):which(names(iChart)==endWindow)]
	aways <- ifelse(temp == "-", 1, ifelse(!is.na(temp), 0, NA))
	means <- colMeans(aways, na.rm=T)

	write.table(means, paste(save_as, ".txt", sep=""), sep="\t", col.names=F)
	
	plot(as.numeric(names(means)), means, xlab = "Time from noun onset (ms)", ylab = "Percentage of trials away")

	graphics.off()
	return(means)


	}


renameCondition <- function(iChart, oldCondition, newCondition) {
	
	for (i in 1:length(oldCondition)) {
			iChart$Condition[iChart$Condition == oldCondition[i]] <- newCondition[i]
		}
	return(iChart)
	}
	
	renameItem <- function(iChart, oldItem, newItem) {
	
	for (i in 1:length(oldItem)) {
			iChart$Target.Image[iChart$Target.Image == oldItem[i]] <- newItem[i]
		}
	return(iChart)
	}

removeParticipants <- function(iChart, participants) {
		iChart[!iChart$Sub.Num %in% participants,]
	}

pairedMeans <- function(iChart, variable, save_as) {
	if(variable == "RT") {
		
		end <- which(names(iChart) == "RT")
		pooled <- cast(iChart[,1:end], Sub.Num+Sex~Condition+Response+Months, mean, na.rm=TRUE)

		} else {
		
		end <- which(names(iChart) == variable)
		pooled <- cast(iChart[,1:end], Sub.Num+Sex~Condition+Months, mean, na.rm=TRUE)
			}
	
		write.table(pooled, paste(save_as, ".txt", sep=""), sep="\t", row.names=F)

	return(pooled)
	
	}

addCDI <- function(iChart, save_as){
	
	CDI <- read.table(file.choose(), header=T, fill=T, as.is=T, sep="\t")

	pooled <- merge(iChart, CDI, by.x = "Sub.Num", all = TRUE)
	
	write.table(pooled, paste(save_as, ".txt", sep=""), sep="\t", row.names=F)
	
	return(pooled)
	}
	
	
	se <- function(x){sd(x, na.rm=T)/sqrt(sum(!is.na(x)))}
	
	
	
	combineData <- function(iChart, accuracy, RT, CDI=FALSE) {
		
		RTandAcc <- 	merge(accuracy, RT, by = "Sub.Num", all = TRUE)
		
		save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_meanAcc_and_RT_by_subs_paired", iChart[1, "StartWindowAcc"], "_", iChart[1, "EndWindowAcc"], ".txt", sep="")


	if(!(CDI == FALSE)) {
		
	RTandAcc <- merge(RTandAcc, CDI, by.x = "Sub.Num", by.y = names(CDI)[1], all = TRUE)	
	
		}
		
		write.table(RTandAcc, save_as, sep="\t", row.names=F)


		return(RTandAcc)
		
		}
	
	

useKnown <- function(iChart, minUnderstands, minSays) {
	
	if("GoodUnderstands" %in% names(iChart)) iChart <- iChart[,-which(names(iChart)%in%c("GoodSays", "GoodUnderstands", "says", "understands"))]
	
	file_raw <- file.choose()
	temp <- gregexpr("/", file_raw)
	
	last_separator <- temp[[1]][length(temp[[1]])]

	Directory <- substr(file_raw, 1, last_separator)
	file <- substr(file_raw, last_separator+1, nchar(file_raw))
	

	knows <- read.table(file_raw, sep="\t", as.is=T, header=T, fill=T)

names_knows <- names(knows)[3:length(names(knows))]
extra_knows <- setdiff(names_knows, unique(iChart$Target.Image))
extra_iChart <- setdiff(unique(iChart$Target.Image), names_knows)

print(paste("The following words are listed in", file, "but do not occur as target words in the working file:"))
print(extra_knows)

print(paste("The following words are listed as target words in the working file but do not occur in ", file, ":", sep=""))
print(extra_iChart)

subs <- unique(iChart$Sub.Num)
subs_understands<- unique(knows$Sub.Num[knows$response=="understands"])
subs_says <- unique(knows$Sub.Num[knows$response=="says"])

	extra_understands <- setdiff(subs_understands, subs)
	extra_says <- setdiff(subs_says, subs)

	print(paste("The following Sub.Num are listed in the understands columns of ", file,  "but are not listed in the working file:"))
	print(extra_understands)
	
	print(paste("The following Sub.Num are listed in the says columns of", file,  "but are not listed in the working file:"))
	print(extra_says)
		
	extra_i_understands <- setdiff(subs, subs_understands)
	extra_i_says <- setdiff(subs, subs_says)
	
	print(paste("The following Sub.Num are listed in the working file but do not occur in the understands column of ", file, ":", sep=""))
	
	print(extra_i_understands)
	
	print(paste("The following Sub.Num are listed in the working file but do not occur in the says column of ", file, ":", sep=""))
	
	print(extra_i_says)


	temp_1 <- melt(knows, id=c("Sub.Num", "response"))
	temp_2 <- cast(temp_1, Sub.Num+variable~response)

	subjects <- unique(iChart$Sub.Num)

	iChart <- merge(iChart, temp_2, by.x = c("Sub.Num", "Target.Image"), by.y = c("Sub.Num", "variable"), all.x=TRUE, all.y=FALSE)
	

	iChart <- iChart[!is.na(iChart$Condition),]
	
	iChart$GoodSays <- iChart$says >= minSays
	iChart$GoodUnderstands <- iChart$understands >= minUnderstands

	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_iChart_known_words", ".txt", sep="")
	
	write.table(iChart, save_as, sep="\t", row.names=F)
	
	iChart <- iChart[iChart$GoodSays & iChart$GoodUnderstands & !is.na(iChart$GoodSays) & !is.na(iChart$GoodUnderstands),]
		
	return(iChart)
	
	}
	
	
	# rrapply(my.data, FUN=cor.test2, STATS=c("estimate", "p.value", "ci1","ci2")) should do what you want.


t.tests <- function(x, extra = "", FUN=t.test, STATS=c("statistic", "p.value"), ...){
   # use FUN=cor.test, STATS=c("estimate", "p.value") for matrix ofcorrelations and p-values
   
   x <- x[,2:ncol(x)]
   
   K <- ncol(x)
   RESU <- list()
   my.fun <- FUN
       fuf<-my.fun(x[,2],x[,1], ...)
       if(!is.list(fuf)) {
             STATS <- "A"
             my.fun <- function(x,y,...) list(A=FUN(x,y,...))
             }
   neimz <- colnames(x)
   for(i in 1:length(STATS)) {
         fof <- STATS[i]
         RESU[[fof]] <- matrix(NA, ncol=K, nrow=K)
         colnames(RESU[[fof]]) <- neimz
         rownames(RESU[[fof]]) <- neimz
         names(RESU)[i] <- fof
         }

   for(i in 2:K) for(j in 1:(i-1)) {
        foo<-my.fun(x[,i],x[,j], ...)
          for(h in 1:length(STATS)) {
               fof <- STATS[h]
               RESU[[fof]][i,j]<- foo[[fof]]
               }
   }
   
   save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_", extra,  "_t.test_matrix.txt", sep="")
   
   if (file.exists(save_as)) { file.remove(save_as) }
   

     	
  if(is.list(fuf)) {
  	
  	out<-capture.output(RESU)
	cat(out,file=save_as,sep="\n",append=TRUE)
	return(RESU)		

  	
  	} else {
  		
  	out<-capture.output(RESU$A)
	cat(out,file=save_as,sep="\n",append=TRUE)  
	return(RESU$A)		
  		}	
	
}


correlations <- function(x, extra = "", FUN=cor.test, STATS=c("estimate", "p.value"), ...){
   # use FUN=cor.test, STATS=c("estimate", "p.value") for matrix ofcorrelations and p-values
   
   x <- x[,2:ncol(x)]
   
   K <- ncol(x)
   RESU <- list()
   my.fun <- FUN
       fuf<-my.fun(x[,2],x[,1], ...)
       if(!is.list(fuf)) {
             STATS <- "A"
             my.fun <- function(x,y,...) list(A=FUN(x,y,...))
             }
   neimz <- colnames(x)
   for(i in 1:length(STATS)) {
         fof <- STATS[i]
         RESU[[fof]] <- matrix(NA, ncol=K, nrow=K)
         colnames(RESU[[fof]]) <- neimz
         rownames(RESU[[fof]]) <- neimz
         names(RESU)[i] <- fof
         }

   for(i in 2:K) for(j in 1:(i-1)) {
        foo<-my.fun(x[,i],x[,j], ...)
          for(h in 1:length(STATS)) {
               fof <- STATS[h]
               RESU[[fof]][i,j]<- foo[[fof]]
               }
   }
   
   save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_", extra, "_correlation_matrix.txt", sep="")
   
   if (file.exists(save_as)) { file.remove(save_as) }
     	
  if(is.list(fuf)) {
  	
  	out<-capture.output(RESU)
	cat(out,file=save_as,sep="\n",append=TRUE)
  	return(RESU)		

  	} else {
  		
  	out<-capture.output(RESU$A)
	cat(out,file=save_as,sep="\n",append=TRUE)
	RESU$A

		
  		}	
	
}



chance <- function(accuracy) {
	
	columns <- ncol(accuracy)

p <- 0
t <- 0
cond <- 0
mean_acc <- 0

for (i in 2:columns) {
	
	test <- t.test(accuracy[,i]-0.5)
	cond[i-1] <- names(accuracy)[i]
	t[i-1] <- test$statistic
	p[i-1] <- test$p.value
	mean_acc[i-1] <- test$estimate+0.5
	
	}
	
	chance <- data.frame(condition = cond, mean = mean_acc, t = t, p = p)
	
	save_as <- paste(iChart[1, "Directory"], iChart[1, "StudyName"], "_accuracy_against_chance.txt", sep="")
	
	write.table(chance, save_as, sep="\t", row.names=F)
	
	return(chance)
	
	
	
	}
	
	
	deleteTrials <- function(iChart, GoodFirstGap, GoodLongestGap, GoodRT) {
	
	# pool RT accross subject
	filteredRT <- iChart
	if(GoodFirstGap) filteredRT <- filteredRT[filteredRT$GoodFirstGap | is.na(filteredRT$GoodFirstGap),]
	if(GoodLongestGap) filteredRT <- filteredRT[filteredRT$GoodLongestGap | is.na(filteredRT$GoodLongestGap),]
	if(GoodRT) filteredRT <- filteredRT[filteredRT$GoodRT,]
	filteredRT <- filteredRT[filteredRT$Response == "D" | filteredRT$Response == "T",]
	
	return(filteredRT)
}

check <- function (iChart) {
	
	print(paste(nrow(iChart),"rows -", "(", nrow(iChart[iChart$Response=="D",]), "Ds )", "(", nrow(iChart[iChart$Response=="T",]), "Ts )","(", nrow(iChart[iChart$Response=="A",]), "As )"))
	
	print(paste(length(unique(iChart$Sub.Num)), "Subjects -", length(unique(iChart$Condition)), "Conditions"))	
	print(xtabs(~Sub.Num+Condition, data=iChart))
	
	}
	
	
	combineFiles <- function(file1= file.choose(), file2 = file.choose(), header1 = "", header2 = "", save_as = "combinedFile.txt") {
	
	
	if (is.character(file1)) {
		
			file1 <- read.table(file1, header=T, fill=T, as.is=T, sep="\t")

		
		}
		
	if (is.character(file2)) {
			file2 <- read.table(file2, header=T, fill=T, as.is=T, sep="\t")
		
		}
		
		names(file1)[2:ncol(file1)] <- paste(names(file1)[2:ncol(file1)], header1, sep="_")
		
		names(file2)[2:ncol(file2)] <- paste(names(file2)[2:ncol(file2)], header2, sep="_")


		combined <- merge(file1, file2, by = "Sub.Num", all = TRUE)


		write.table(combined, save_as, sep="\t", row.names=F)

		return(combined)
	
	}
	
	
	
	
	

reliability <- function(iChart.x, iChart.y, rejectionAccuracy, rejectionRT, agreement, RejectFirstGap, RejectLongestGap) {


if (RejectFirstGap) {
	
iChart.x <- iChart.x[iChart.x$GoodFirstGap| is.na(iChart.x$GoodFirstGap ),]
iChart.y <- iChart.y[iChart.y$GoodFirstGap| is.na(iChart.y$GoodFirstGap ),]

}


if (RejectLongestGap) {
	
iChart.x <- iChart.x[iChart.x$GoodLongestGap| is.na(iChart.x$GoodLongestGap ),]
iChart.y <- iChart.y[iChart.y$GoodLongestGap| is.na(iChart.y$GoodLongestGap ),]


}


######CREATE HEADER WITH DATE AND INFORMATION

file.remove(paste(iChart.y$Directory[1],"Reliability_Report.txt"))
sink(paste(iChart.y$Directory[1],"Reliability_Report.txt"))

Sys.time()

cat("\n\nx = ",iChart.x$StudyName[1], "\n")
cat("y = ",iChart.y$StudyName[1], "\n\n\n")

cat("Rejection Accuracy = ", round(rejectionAccuracy, digits = 2), "\n")
cat("Rejection RT = ", rejectionRT, "\n\n")


cat("StartWindowRT = ", iChart.x$StartWindowRT[1], "\nEndWindowRT = ", iChart.x$EndWindowRT[1], "\nStartWindowAccuracy = ",  iChart.x$StartWindowAcc[1], "\nEndWindowAccuracy = ", iChart.x$EndWindowAcc[1], "\n\n")


####### CLEAN AND MERGE TWO ICHARTS

# recompute accuracy on all trials
iChart.x.temp <- iChart.x[,names(iChart.x) %in% paste(c(as.numeric(iChart.x$StartWindowAcc[1]):as.numeric(iChart.x$EndWindowAcc[1])))]
iChart.y.temp <- iChart.y[,names(iChart.y) %in% paste(c(as.numeric(iChart.y$StartWindowAcc[1]):as.numeric(iChart.y$EndWindowAcc[1])))]

iChart.x$Accuracy <- rowSums(iChart.x.temp=="1")/(rowSums(iChart.x.temp=="1")+rowSums(iChart.x.temp=="0"))
iChart.y$Accuracy <- rowSums(iChart.y.temp=="1")/(rowSums(iChart.y.temp=="1")+rowSums(iChart.y.temp=="0"))


iChart.x[,c("Sub.Num", "Tr.Num", "Response", "Accuracy", "RT", "Condition", "GoodRT", "GoodLongestGap", "GoodFirstGap")] -> iChart.xR
iChart.y[,c("Sub.Num", "Tr.Num", "Response", "Accuracy", "RT", "Condition", "GoodRT", "GoodLongestGap", "GoodFirstGap")] -> iChart.yR


iChart.xR$Accuracy <- round(iChart.xR$Accuracy, digits = 2)
iChart.yR$Accuracy <- round(iChart.yR$Accuracy, digits = 2)

merge(iChart.xR, iChart.yR, by=c("Sub.Num", "Tr.Num", "Condition")) -> iChartR

# remove trials in which both coders marked respose as A
iChartR <- iChartR[!(iChartR$Response.x == "A" & iChartR$Response.y == "A"),]



######## CHECK RESPONSES
correctResponses <- sum(iChartR$Response.y == iChartR$Response.x)/nrow(iChartR)

cat("\nReliability for Responses = ", round(correctResponses, digits = 2))

table(iChartR$Response.x, iChartR$Response.y)

cat("\n\n")	

if(correctResponses != 1) print(iChartR[iChartR$Response.y != iChartR$Response.x, c("Sub.Num", "Tr.Num", "Condition","Response.x", "Response.y")])

######## CHECK NA TRIALS

sumNA <- sum(is.na(iChartR$Accuracy.x) == is.na(iChartR$Accuracy.y))
totalNA <- nrow(iChartR)
correctNA <- sumNA/totalNA

cat("\n\nACCURACY\n\nNA Accuracy trials, including aways")
table(is.na(iChartR$Accuracy.x), is.na(iChartR$Accuracy.y))

if(correctNA != 1) {
	
	cat("\n\n")
	
	print(iChartR[is.na(iChartR$Accuracy.x) != is.na(iChartR$Accuracy.y), c("Sub.Num","Tr.Num","Accuracy.x", "Accuracy.y")])

}

######## CHECK ACCURACY


names(iChartR) %in% paste(c(as.numeric(iChart.y$StartWindowAcc[1]):as.numeric(iChart.y$EndWindowAcc[1])))


sumAccuracy <- sum(ifelse(abs(iChartR$Accuracy.x-iChartR$Accuracy.y) <= rejectionAccuracy, TRUE, FALSE), na.rm=TRUE)
totalAccuracy <- sum(!is.na(iChartR$Accuracy.x) & !is.na(iChartR$Accuracy.y))

correctAccuracy <- sumAccuracy/totalAccuracy

cat("\n\n\nReliability for Accuracy = ", round(correctAccuracy, digits = 2), "\n")
cat("Correct trials = ", sumAccuracy, "\n")
cat("Out of  = ", totalAccuracy, "\n")

totalAccuracyNA <- nrow(iChartR) - sum(is.na(iChartR$Accuracy.x) & is.na(iChartR$Accuracy.y))
correctAccuracyNA <- sumAccuracy/totalAccuracyNA 

cat("\n\n\nReliability for Accuracy including NAs = ", round(correctAccuracyNA, digits = 2), "\n")
cat("Correct trials = ", sumAccuracy, "\n")
cat("Out of  = ", totalAccuracyNA , "\n")
cat("\n\n")

if(correctAccuracy != 1) print(iChartR[!ifelse(abs(iChartR$Accuracy.x-iChartR$Accuracy.y) <= rejectionAccuracy, TRUE, FALSE)
	& !is.na(iChartR$Accuracy.x) & !is.na(iChartR$Accuracy.y), c("Sub.Num", "Tr.Num", "Condition", "Accuracy.x", "Accuracy.y")])


######## CHECK RT FOR D TRIALS ###############

iChartD <- iChartR[iChartR$Response.x == "D" & iChartR$Response.y =="D",]

RejectRT <- TRUE	# Reject trials with min/max values outside of RT window (defined above, inclusive)
# RejectRT <- FALSE # DO NOT reject trials with min/max values outside of RT window (defined above, inclusive)

iChartD <- iChartD[!(!iChartD$GoodRT.x & !iChartD$GoodRT.y) | is.na(iChartD$GoodRT.x) | is.na(iChartD$GoodRT.y),]


######## CHECK NA TRIALS

sumNA <- sum(is.na(iChartD$RT.x) == is.na(iChartD$RT.y))
totalNA <- nrow(iChartD)
correctNA <- sumNA/totalNA

cat("\n\nRT (D) \n\nNA RT trials (D)")
table(is.na(iChartD$RT.x), is.na(iChartD$RT.y))

if(correctNA != 1) {
	
	cat("\n\n")
	print(iChartD[is.na(iChartD$RT.x) != is.na(iChartD$RT.y), c("Sub.Num", "Tr.Num", "RT.x", "RT.y")])
	
	
	}

#### COMPUTE RT


sumRT <- sum(ifelse(abs(iChartD$RT.x-iChartD$RT.y) <= rejectionRT, TRUE, FALSE), na.rm=TRUE)
totalRT <- sum(!is.na(iChartD$RT.x))
correctRT <- sumRT/totalRT

cat("\n\nReliability for RT (D) = ", round(correctRT, digits = 2), "\n")
cat("Correct trials = ", sumRT, "\n")
cat("Out of  = ", totalRT, "\n\n\n")

if(correctRT != 1) print(iChartD[!ifelse(abs(iChartD$RT.x-iChartD$RT.y) <= rejectionRT, TRUE, FALSE) & !is.na(iChartD$RT.x) & !is.na(iChartD$RT.y), c("Sub.Num", "Tr.Num", "Condition", "RT.x", "RT.y")])


######## CHECK RT FOR T TRIALS ###############

iChartT <- iChartR[iChartR$Response.x == "T" & iChartR$Response.y =="T",]


RejectRT <- TRUE	# Reject trials with min/max values outside of RT window (defined above, inclusive)
# RejectRT <- FALSE # DO NOT reject trials with min/max values outside of RT window (defined above, inclusive)


iChartT <- iChartT[!(!iChartT$GoodRT.x & !iChartT$GoodRT.y) | is.na(iChartT$GoodRT.x) | is.na(iChartT$GoodRT.y),]


######## CHECK NA TRIALS

sumNA <- sum(is.na(iChartT$RT.x) == is.na(iChartT$RT.y))
totalNA <- nrow(iChartT)
correctNA <- sumNA/totalNA

cat("\n\nRT (T) \n\nNA RT trials (T)")
table(is.na(iChartT$RT.x), is.na(iChartT$RT.y))

if(correctNA != 1) {
	
	cat("\n\n")
	
	print(iChartT[is.na(iChartT$RT.x) != is.na(iChartT$RT.y), c("Sub.Num", "Tr.Num", "RT.x", "RT.y")])
	
	}

#### COMPUTE RT
sumRT <- sum(ifelse(abs(iChartT$RT.x-iChartT$RT.y) <= rejectionRT, TRUE, FALSE), na.rm=TRUE)
totalRT <- sum(!is.na(iChartT$RT.x))
correctRT <- sumRT/totalRT

cat("\n\nReliability for RT (T) = ", round(correctRT, digits = 2), "\n")
cat("Correct trials = ", sumRT, "\n")
cat("Out of  = ", totalRT, "\n\n\n")

if(correctRT != 1) print(iChartT[!ifelse(abs(iChartT$RT.x-iChartT$RT.y) <= rejectionRT, TRUE, FALSE) & !is.na(iChartT$RT.x) & !is.na(iChartT$RT.y), c("Sub.Num", "Tr.Num", "Condition", "RT.x", "RT.y")])


#### COMPUTE FRAME AGREEMENT



iChart.x[,names(iChart.x) %in% c(as.numeric(iChart.x$StartWindowAcc[1]):as.numeric(iChart.x$EndWindowAcc[1]), "Sub.Num", "Tr.Num", "Response", "Accuracy", "RT", "Condition")] -> iChart.xR


iChart.y[,names(iChart.y) %in% c(as.numeric(iChart.y$StartWindowAcc[1]):as.numeric(iChart.y$EndWindowAcc[1]), "Sub.Num", "Tr.Num", "Response", "Accuracy", "RT", "Condition")] -> iChart.yR


merge(iChart.xR, iChart.yR, by=c("Sub.Num", "Tr.Num", "Condition")) -> iChartR

iChartR <- iChartR[(iChartR$Response.x != "A") & iChartR$Response.y != "A",]

iChartR[iChartR == "."] <- "-"

iChart.xR <- iChartR[,names(iChartR) %in% paste(c(as.numeric(iChart.y$StartWindowAcc[1]):as.numeric(iChart.y$EndWindowAcc[1])), ".x", sep="")]

iChart.yR <- iChartR[,names(iChartR) %in% paste(c(as.numeric(iChart.y$StartWindowAcc[1]):as.numeric(iChart.y$EndWindowAcc[1])), ".y", sep="")]

sumRows <- rowSums((iChart.xR == "0" & iChart.yR == "0") | (iChart.xR == "1" & iChart.yR == "1") | (iChart.xR == "-" & iChart.yR == "-"))

corrects <- round(sumRows/(ncol(iChart.yR) - rowSums(iChart.xR == "" & iChart.yR == "")), digits=2)
iChartR$FrameAgreement <- corrects

cat("\n\nMean frame agreement = ", round(mean(corrects), digits = 2), "\n\n")

iChartR[corrects <= agreement,]

print(iChartR[corrects <= agreement, c("Sub.Num", "Tr.Num", "Condition", "Accuracy.x", "Accuracy.y", "FrameAgreement")])




sink()

}
	
	
	
	
	
	
	
compare.iChart <- function(iChart.x, iChart.y) {
	
merge(iChart.x[,c("Sub.Num", "Tr.Num", "Response")], iChart.y[,c("Sub.Num", "Tr.Num", "Response")], by=c("Sub.Num", "Tr.Num"), all.x=TRUE, all.y=TRUE) -> comparing

cat("\n\nThe following trials were not coded on iChart.x, but were coded on iChart.y\n")

print(comparing[is.na(comparing$Response.x), c("Sub.Num","Tr.Num")])

cat("\n\nThe following trials were not coded on iChart.y, but were coded on iChart.x\n")

print(comparing[is.na(comparing$Response.y), c("Sub.Num","Tr.Num")])


}

