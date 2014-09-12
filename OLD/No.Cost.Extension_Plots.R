#### Plots for No Cost Extension Report #####

library(ggplot2)
source("http://www.stanford.edu/~ricardoh/fortify.R")

d1 <- read.csv("/Users/kylemacdonaldadmin/Documents/Projects/SOL/SOL_Data /trio_face_plot.csv")

### Plot of proportion shifting to the target image ###

qplot(data=d1, x=Group, y=MeanCT, main="Mean Proportion of First Shifts to Target", fill=Group, label=MeanCT) +
  scale_fill_manual(breqplot(data=subst, x=Group, y=MeanCT, main="Mean Proportion of First Shifts to Target", fill=Group, label=MeanCT) +
  scale_fill_manual(breaks = 1:2, values=c("gray68", "gray30")) +
  scale_x_discrete(labels=c("Non-signers", "Older Signers")) +
  geom_bar(width=.4) +
  coord_cartesian(ylim= c(0, 1)) +
  geom_hline(yintercept=0.5, linetype="dotted") +
  ylab(NULL) + 
  xlab("Group") +
  theme_bw(base_size = 16) + 
  theme(axis.title.x=element_text(vjust=-0.5)) +
  theme(title=element_text(vjust=1.5)) + 
  theme(axis.text.x=element_text(size=14)) + 
  geom_text(vjust=-2)aks = 1:3, values=c("gray68", "gray46", "gray30")) +
  scale_x_discrete(labels=c("Non-signers", "Younger Signers", "Older Signers")) +
  geom_bar(width=.4) +
  coord_cartesian(ylim= c(0, 1)) +
  geom_hline(yintercept=0.5, linetype="dotted") +
  ylab(NULL) + 
  xlab("Group") +
  theme_bw(base_size = 16) + 
  theme(axis.title.x=element_text(vjust=-0.5)) +
  theme(title=element_text(vjust=1.5)) + 
  theme(axis.text.x=element_text(size=14)) + 
  geom_text(vjust=-2)

## Subset data to make different graphs

subst <- subset(d1, d1$Group != "Signers")

 
#### Correlation Plots ####

d2 <- read.csv("/Users/kmacdonald/Documents/Projects/SOL/SOL_Data/SOL.n33.alldata.csv")

#### Accuracy by Age ####

cor.test(d2$vanilla, d2$Months)

age_acc_plot <- ggplot(data=d2, aes(x=Months, y=vanilla)) +
                        geom_point(size=4) + 
                        xlab("Age (months)") +
                        ylab("Mean Accuracy") +
                        ylim(0.45, 0.75) +
                        theme(axis.title.x=element_text(vjust=-0.5)) +
                        theme(axis.title.y=element_text(vjust=0.1)) +
                        stat_smooth(method="lm", se=FALSE) +
                        ggtitle("Positive Relationship Between Age and Accuracy") +
                        geom_text(x=40, y=0.54, label="r(21) = .60", size = 8, face="bold")

d2_18_30 <- d2[which(d2$Months <= 30 & d2$Months >= 18),]

d2_all_good_subs <- subset(d2, d2$X <= 23m)

write.table(d2_all_good_subs, file="SOL_n23", na="NA", sep="\t", row.names=F)

#### Accuracy by CDI ####

## First we have to subset the data to remove kids older than 30 months

d2_under30 <- d2[which(d2$Months <= 30 & d2$Sub.Num != 30030),]

## Correlations: Acc_CDI, RT_CDI ##

cor.test(d2_under30$vanilla, d2_under30$signs_produced)

cdi_acc_plot <- ggplot(data=d2_under30, aes(x=signs_produced, y=vanilla)) +
                        geom_point(size=4) + 
                        xlab("Signs Produced") +
                        ylab("Mean Accuracy") +
                        ylim(0.45, 0.75) +
                        theme(axis.title.x=element_text(vjust=-0.5)) +
                        theme(axis.title.y=element_text(vjust=0.3)) +
                        stat_smooth(method="lm", se=FALSE) +
                        ggtitle("Positive Relationship Between Vocabulary and Accuracy") +
                        geom_text(x=65, y=0.54, label="r(10) = .72", size = 8, face="bold")
                

## Now plot correlations together

multiplot(age_acc_plot, cdi_acc_plot, cols = 2)

## RT by age and cdi ##

cor.test(d2$van_D, d2$Months)
cor.test(d2$van_D, d2$signs_produced)

## Remove outlier ##

d2_no_out <- d2[which(d2$Sub.Num != 30011 & d2$Sub.Num != 30010),]

cor.test(d2_no_out$van_D, d2_no_out$Months)
plot(d2_no_out$van_D, d2_no_out$Months)

cor.test(d2_no_out$van_D, d2_no_out$signs_produced)

#### Profile Plots of the VLP data ####


