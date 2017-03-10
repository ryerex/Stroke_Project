###############################################################################
#
# Project: 	ACO Projects
# Script:	StrokeRecurrence.R
# Version:	
# Created:	Oct 14, 2014
# Updated:	Feb 24, 2015
# Author: 	ry5t
###############################################################################

library('pec')
library('timereg')
library('xtable')
library('plyr')
library('reshape')
library('survival')
library('rms')
library('randomForestSRC')
library('party')
library('prodlim')
library('ggplot2')
library('caret')
library('data.table')


options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)


#----------------
# Partition Data
#----------------
set.seed(123)	# set the seed to make partition reproductible
data('cost')	# Scandinavian Stroke Study Data Set
cost <- as.data.table(cost)
smp_size <- floor(0.90 * nrow(cost))
train_ind <- sample(seq_len(nrow(cost)), size = smp_size)
train <- cost[train_ind]
test <- cost[-train_ind]
#-------------------------------------------------------------------------------------------------
# First stage is to  compare a few algorithm options including Cox Regression, Conditional Forest, 
# and Random Survival Forest using all available features in an additive formulation
#-------------------------------------------------------------------------------------------------

fitformFull <- Surv(time,status) ~ age + sex + hypTen + ihd + prevStroke + othDisease + alcohol + diabetes + 
								   smoke + atrialFib + hemor + strokeScore + cholest

set.seed(13)
coxModel <- selectCox(fitformFull, data = train, rule = 'aic')
randomForestModel <- rfsrc(fitformFull, data = train, forest = TRUE, ntree = 10)
conditionalForestModel <- pecCforest(fitformFull, data = train, controls = cforest_classical(ntree = 100))

modelFitCompAlg <- pec(list('CoxModel' = coxModel, 'RandomForestModel' = randomForestModel, 'ConditionalForestModel' = conditionalForestModel), 
						data = train, formula = Surv(time, status) ~ 1, splitMethod = 'Boot632plus', B = 100, M = 350, 
						keep.index = TRUE, keep.matrix = TRUE)

curves <- with(modelFitCompAlg,as.data.frame(cbind(as.integer(time), Boot632plusErr$Reference, 
						Boot632plusErr$CoxModel, Boot632plusErr$RandomForestModel, Boot632plusErr$ConditionalForestModel)))

names(curves) <- c('Days', 'Null', 'CoxModel', 'RandomForestModel', 'ConditionalForestModel' )
curves <- melt(curves, id = 'Days')
names(curves) <- c('Days', 'Model', 'Error')
x11()
ggplot(data = curves, aes(x = Days, y = Error)) + 
		geom_smooth(aes(color = Model), size = 2) +		
		theme(plot.title = element_text(lineheight=3, vjust = 1, face='bold.italic', size = 16), 
				axis.title = element_text(size = 12, face = 'italic')) + 
		ggtitle('Cross Validation Accuracy') + 
		xlab('Time (Days)') +
		ylab('Error (Integrated Brier Score)') +
		xlim(0, 4000)


table4 <- as.data.frame(sort(vimp(randomForestModel)$importance, decreasing = TRUE))
colnames(table4) <- c('VIMP')

#------------------------------------------------------------------
# Fit a reduced model to the best algorithm, Random Survival Forest
#------------------------------------------------------------------
fitformSelected <- Surv(time,status) ~ age + sex + prevStroke + othDisease + alcohol + diabetes + atrialFib + strokeScore + cholest

randomForestModelSelected <- rfsrc(fitformSelected, data = train, forest = TRUE, ntree = 100)
modelFitCompFeat <- pec(list('FullModel' = randomForestModel, 'ReducedModel' = randomForestModelSelected), data = train, formula = Surv(time, status) ~ 1,
					splitMethod = 'Boot632plus', B = 100, M = 350, keep.index = TRUE, keep.matrix = TRUE)
table5 <- as.table(crps(modelFitCompFeat))
curves <- with(modelFitCompFeat,as.data.frame(cbind(as.integer(time), Boot632plusErr$Reference, Boot632plusErr$FullModel, Boot632plusErr$ReducedModel)))
names(curves) <- c('Days', 'Null', 'FullModel', 'ReducedModel')
curves <- melt(curves, id = 'Days')
names(curves) <- c('Days', 'Model', 'Error')
x11()
ggplot(data = curves, aes(x = Days, y = Error)) + 
		geom_smooth(aes(color = Model), size = 2) +		
		theme(plot.title = element_text(lineheight=3, vjust = 1, face='bold.italic', size = 16), axis.title = element_text(size = 12, face = 'italic')) + 
		ggtitle('Cross Validation Accuracy') + 
		xlab('Time (Days)') +
		ylab('Error (Integrated Brier Score)') +
		xlim(0, 4000)

table6 <- xtable(as.table(sort(vimp(randomForestModel)$importance, decreasing = TRUE)))
cindex(randomForestModelSelected, fitformSelected, data = train)

x11()
calPlot(randomForestModelSelected, fitformSelected, time = 30,  data = train)

test$prsf <- predictSurvProb(randomForestModelSelected, newdata = test, times = 10 * 365.25)