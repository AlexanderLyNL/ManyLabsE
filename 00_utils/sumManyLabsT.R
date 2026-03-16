load(file=paste0(myWd, "critcher1", "Result.RData"))

# Scenario 1-----
round(mean(res1$logMetaE), 3)
round(sd(res1$logMetaE), 3)

round(mean(res1$logMetaE >= log(1/alpha))*100, 3)
# round(mean(res1$stopDecision==1)*100, 3)

# res1$logMetaEFut[which(is.infinite(res1$logMetaEFut))] <- -1e27

round(mean(res1$logMetaEFut), 3)
round(sd(res1$logMetaEFut), 3)

# round(mean(res1$stopDecision==-1)*100, 3)

round(mean(res1$logMetaEFut <= log(betaFutility))*100, 3)

round(mean(res1$eValues >= 1/alpha)*100, 3)
round(mean(res1$eValuesFut <= betaFutility)*100, 3)


round(mean(res1$nStudies), 3)
round(sd(res1$nStudies), 3)


res1$nStudiesAlternativeWorstCase
res1$nStudiesFutilityWorstCase

# res1$nSamplesAlternativeWorstCase
# res1$nSamplesFutilityWorstCase


# round(mean(res1$stopDecision==-1)*100, 3)

round(mean(res1$totalStoppingTimes), 3)
round(sd(res1$totalStoppingTimes), 3)

# Scenario 2-----

logMetaE<- rowSums(log(res2$eValues))
round(mean(logMetaE), 3)
round(sd(logMetaE), 3)


round(mean(logMetaE >= log(1/alpha))*100, 3)
# round(sd(logMetaE >= log(1/alpha)), 3)

logMetaEFut <- rowSums(log(res2$eValuesFut))
round(mean(logMetaEFut), 3)
round(sd(logMetaEFut), 3)

round(mean(logMetaEFut <= log(betaFutility))*100, 3)
# round(sd(logMetaEFut <= log(betaFutility))*100, 3)

round(mean(res2$alternativeProportion)*100, 3)
round(sd(res2$alternativeProportion)*100, 3)

round(mean(res2$futilityProportion)*100, 3)
round(sd(res2$futilityProportion)*100, 3)


round(mean(res2$totalStoppingTimes), 3)
round(sd(res2$totalStoppingTimes), 3)

#Scenario 3 ------
round(mean(res3$logMetaE), 3)
round(sd(res3$logMetaE), 3)


# round(mean(rowMeans(res3$nStopDecision==1)*100), 3)


round(mean(res3$logMetaE >= log(1/alpha))*100, 3)
# round(sd(res3$logMetaE >= log(1/alpha)), 3)


round(mean(res3$logMetaEFut), 3)
round(sd(res3$logMetaEFut), 3)


round(mean(res3$logMetaEFut <= log(betaFutility))*100, 3)
# round(sd(res3$logMetaEFut <= log(betaFutility))*100, 3)

round(mean(res3$alternativeProportion)*100, 3)
round(sd(res3$alternativeProportion)*100, 3)

round(mean(res3$futilityProportion), 3)
round(sd(res3$futilityProportion), 3)

round(mean(res3$totalStoppingTimes), 3)
round(sd(res3$totalStoppingTimes), 3)


# PLUS ---------------------
load(file=paste0(myWd, "critcher1", "ResultPlus.RData"))
# Scenario 1-----
round(mean(res1Plus$logMetaE), 3)
round(sd(res1Plus$logMetaE), 3)

round(mean(res1Plus$logMetaE >= log(1/alpha))*100, 3)

round(mean(res1Plus$logMetaEFut), 3)
round(sd(res1Plus$logMetaEFut), 3)

round(mean(res1Plus$logMetaEFut <= log(betaFutility))*100, 3)
# round(mean(res1Plus$stopDecision==-1)*100, 3)


round(mean(res1Plus$eValues >= 1/alpha)*100, 3)
round(mean(res1Plus$eValuesFut <= betaFutility)*100, 3)


round(mean(res1Plus$nStudies), 3)
round(sd(res1Plus$nStudies), 3)


res1Plus$nStudiesAlternativeWorstCase
res1Plus$nStudiesFutilityWorstCase

res1Plus$nSamplesAlternativeWorstCase
res1Plus$nSamplesFutilityWorstCase


# round(mean(res1Plus$stopDecision==-1)*100, 3)

round(mean(res1Plus$totalStoppingTimes), 3)
round(sd(res1Plus$totalStoppingTimes), 3)

# Scenario 2-----

logMetaE<- rowSums(log(res2Plus$eValues))
round(mean(logMetaE), 3)
round(sd(logMetaE), 3)


round(mean(logMetaE >= log(1/alpha)), 3)
# round(sd(logMetaE >= log(1/alpha)), 3)

logMetaEFut <- rowSums(log(res2Plus$eValuesFut))
round(mean(logMetaEFut), 3)
round(sd(logMetaEFut), 3)

round(mean(logMetaEFut <= log(betaFutility))*100, 3)
# round(sd(logMetaEFut <= log(betaFutility))*100, 3)

round(mean(res2Plus$alternativeProportion)*100, 3)
round(sd(res2Plus$alternativeProportion)*100, 3)

round(mean(res2Plus$futilityProportion)*100, 3)
round(sd(res2Plus$futilityProportion)*100, 3)


round(mean(res2Plus$totalStoppingTimes), 3)
round(sd(res2Plus$totalStoppingTimes), 3)

#Scenario 3 ------
round(mean(res3Plus$logMetaE), 3)
round(sd(res3Plus$logMetaE), 3)

# round(mean(res3Plus$nStopDecision==1)*100, 3)

round(mean(res3Plus$logMetaE >= log(1/alpha))*100, 3)
# round(sd(res3Plus$logMetaE >= 1/alpha), 3)

round(mean(res3Plus$logMetaEFut), 3)
round(sd(res3Plus$logMetaEFut), 3)

round(mean(res3Plus$logMetaEFut <= log(betaFutility))*100, 3)
# round(sd(res3Plus$logMetaEFut <= betaFutility)*100, 3)

round(mean(res3Plus$alternativeProportion)*100, 3)
round(sd(res3Plus$alternativeProportion)*100, 3)

round(mean(res3Plus$futilityProportion), 3)
round(sd(res3Plus$futilityProportion), 3)

round(mean(res3Plus$totalStoppingTimes), 3)
round(sd(res3Plus$totalStoppingTimes), 3)
