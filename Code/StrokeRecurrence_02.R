###############################################################################
#
# Project: 	ACO Projects
# Script:	StrokeRecurrence_02.R
# Version:	
# Created:	Oct 14, 2014
# Updated:	Oct 16, 2014
# Author: 	ry5t
###############################################################################

library("pec")
library("timereg")
library("xtable")
library("plyr")
library("reshape")
library("survival")
library("rms")
library("randomForestSRC")
library("party")
library("prodlim")
library("ggplot2")
library("caret")
data("cost")

##########################
# Set multi-core options #
##########################
options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)

##################################################
# Data set is from the Scandinavian Stroke Study #
##################################################
featureTable1 <- head(subset(cost, select=c(age, sex, hypTen, ihd, prevStroke, othDisease, alcohol, diabetes)))
featureTable2 <- head(subset(cost, select=c(smoke, atrialFib, strokeScore, cholest, time, status)))

table1 <- xtable(featureTable1)
table1

table2 <- xtable(featureTable2)
table2

######################################################
# First stage is to  compare a few algorithm options #
# including Cox Regression, Conditional Forest, and  #
# Random Survival Forest                             #
#													 #
# The initial model includes all the available 		 #
# features in an additive formulation				 #
######################################################

fitformFull <- Surv(time,status) ~ age + sex + hypTen + ihd + prevStroke + othDisease + alcohol + diabetes + smoke + atrialFib + hemor + strokeScore + cholest

##################
# Fit the models #
##################
set.seed(13)
coxModel <- selectCox(fitformFull, data = cost, rule = "aic")
randomForestModel <- rfsrc(fitformFull, data = cost, forest = TRUE, ntree = 10)
conditionalForestModel <- pecCforest(fitformFull, data = cost, controls = cforest_classical(ntree = 10))

##################
# Test the fits  #
##################
modelFitCompAlg <- pec(list("CoxModel" = coxModel, "RandomForestModel" = randomForestModel, "ConditionalForestModel" = conditionalForestModel), data = cost, formula = Surv(time, status) ~ 1,
					splitMethod = "Boot632plus", B = 10, M = 350, keep.index = TRUE, keep.matrix = TRUE)
table3 <- xtable(as.table(crps(modelFitCompAlgo)))
table3

#################################################
# Collate the data representing the probability #
# curves to use in plotting with ggplot2        #
#################################################
curves <- with(modelFitCompAlg,as.data.frame(cbind(as.integer(time), Boot632plusErr$NullModel, Boot632plusErr$CoxModel, Boot632plusErr$RandomForestModel, Boot632plusErr$ConditionalForestModel)))
names(curves) <- c("Days", "Null", "CoxModel", "RandomForestModel", "ConditionalForestModel" )
curves <- melt(curves, id = "Days")
names(curves) <- c("Days", "Model", "Error")
x11()
ggplot(data = curves, aes(x = Days, y = Error)) + 
		geom_smooth(aes(color = Model), size = 2) +		
		theme(plot.title = element_text(lineheight=3, vjust = 1, face="bold.italic", size = 16), axis.title = element_text(size = 12, face = "italic")) + 
		ggtitle("Cross Validation Accuracy") + 
		xlab("Time (Days)") +
		ylab("Error (Integrated Brier Score)") +
		xlim(0, 4000)


table4 <- as.data.frame(sort(vimp(randomForestModel)$importance, decreasing = TRUE))
colnames(table4) <- c("VIMP")
xtable(table4, digits = 4)

fitformSelected <- Surv(time,status) ~ age + sex + prevStroke + othDisease + alcohol + diabetes + atrialFib + strokeScore + cholest

randomForestModelSelected <- rfsrc(fitformSelected, data = cost, forest = TRUE, ntree = 10)
modelFitCompFeat <- pec(list("FullModel" = randomForestModel, "ReducedModel" = randomForestModelSelected), data = cost, formula = Surv(time, status) ~ 1,
					splitMethod = "Boot632plus", B = 10, M = 350, keep.index = TRUE, keep.matrix = TRUE)
table5 <- xtable(as.table(crps(modelFitCompFeat)))
table5
			
#################################################
# Collate the data representing the probability #
# curves to use in plotting with ggplot2        #
#################################################
curves <- with(modelFitCompFeat,as.data.frame(cbind(as.integer(time), Boot632plusErr$NullModel, Boot632plusErr$FullModel, Boot632plusErr$ReducedModel)))
names(curves) <- c("Days", "Null", "FullModel", "ReducedModel")
curves <- melt(curves, id = "Days")
names(curves) <- c("Days", "Model", "Error")
x11()
ggplot(data = curves, aes(x = Days, y = Error)) + 
		geom_smooth(aes(color = Model), size = 2) +		
		theme(plot.title = element_text(lineheight=3, vjust = 1, face="bold.italic", size = 16), axis.title = element_text(size = 12, face = "italic")) + 
		ggtitle("Cross Validation Accuracy") + 
		xlab("Time (Days)") +
		ylab("Error (Integrated Brier Score)") +
		xlim(0, 4000)


table6 <- xtable(as.table(sort(vimp(randomForestModel)$importance, decreasing = TRUE)))
table6
#####################################################
# At this point, I have selected an algorithm and   #
# a set of features. Now I will see if I can refine #
# the selected features through transformations     #
#####################################################

# Histogram overlaid with kernel density curve
x11()
ggplot(cost, aes(x = strokeScore)) + 
		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
				binwidth=.5,
				colour="black", fill="white") +
		geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

x11()
ggplot(cost, aes(x = age)) + 
		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
				binwidth=.5,
				colour="black", fill="white") +
		geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

x11()
ggplot(cost, aes(x = cholest)) + 
		geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
				binwidth=.5,
				colour="black", fill="white") +
		geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

