# Setup paths ----
sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/01_tables/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))

# Data load ----


# fileNeem <- "Tversky_1_clean_tables.csv"

# authorNeem <- "Rottenstreich"
# authorNeem <- "Tversky"
# authorNeem <- "Hauser1"
authorNeem <- "Hauser2"

fileNeem <- switch(authorNeem,
                   Rottenstreich="Rottenstreich_1_clean_tables.csv",
                   Tversky="Tversky_1_clean_tables.csv",
                   Hauser1="Hauser_1_clean_tables.csv",
                   Hauser2="Hauser_4_clean_tables.csv")

# fileNeem <- "Hauser_4_clean_tables.csv"
tableFile <- paste0(myWd, "data/", fileNeem)
dat <- read.csv(file=tableFile)

names(dat) <- c("source", "ya", "yb", "na", "nb")

# save(dat, file="hauser2.RData")

dim(dat)
head(dat)

debugonce(scenario12x2Help)


res1TwoSided <- manyLabsMetaScenarios(1, nSim=1e3L, alternative="twoSided")

res1DeltaMinFactor$resultTable
res1$resultTable
res1TwoSided$resultTable


res1TwoSided$individualResultList$tversky$logMetaE
res1TwoSided$individualResultList$tversky$logMetaEMinEffi

res1TwoSided$individualResultList$rottenstreich$logMetaE
res1TwoSided$individualResultList$rottenstreich$logMetaEMinEffi

res1DeltaMinFactor <- res1

res1$resultTable
res1DeltaMinFactor$resultTable
res1TwoSided$resultTable
res1$resultTableFull

# save(res1, res1DeltaMinFactor, res1TwoSided, file="2x2TestAll.RData")

res1B <- res1
res1DeltaMinFactorB <- res1DeltaMinFactor
res1TwoSidedB <- res1TwoSided

# remove(res1, res1DeltaMinFactor, res1TwoSided)

load(file="2x2TestAll.RData")
# allSources <- checkEnoughDataInTable(dat)$allSources
names(dat)
names(dat) <- c("source", "ya", "yb", "na", "nb")
allSources <- unique(dat$source)
dat <- dat[dat$source %in% allSources, ]

length(allSources)


head(dat)

# Here -------
alpha <- 0.05
betaMinEffi <- alpha
alternative <- switch(authorNeem,
                      Rottenstreich="less",
                      Tversky="greater",
                      Hauser1="twoSided",
                      Hauser2="greater")
# esMin <- 1.08
# deltaMin <- switch(authorNeem,
#                 Rottenstreich=0.74,
#                 Tversky=1.08,
#                 Hauser1=2.5,
#                 Hauser2=0.34)

esMin <- switch(authorNeem,
                Rottenstreich=0.74*pi/sqrt(3),
                Tversky=log(4.96),
                Hauser1=2.5*pi/sqrt(3),
                Hauser2=0.34*pi/sqrt(3))
minEffiParam <- esMin
logOddsRatio <- esMin
#
# esMin <- log(4.96)*sqrt(3)/pi
# esMin <- 4.96




designObj <- list(esMin=esMin, minEffiTestResult=list(parameter=minEffiParam),
                  alternative=alternative,
                  testName="2x2")

# aap <- scenario1Table(dat, allSources, designObj)

aap <- metaScenario1(dat, allSources, designObj)


metaE <- cumsum(log(aap$eValues))
metaEMinEffi <- cumsum(log(aap$eValuesMinEffi))

metaE >= log(1/alpha)
metaEMinEffi <= log(betaMinEffi)

plot(log(aap$eValuesMinEffi), type="l", col="darkgoldenrod")
lines(1:length(aap$eValues), log(aap$eValues), col="blue")

round(aap$eValues-aap$eValuesMinEffi, 2)
# which(round(aap$eValues-1, 6)==0)
# which(round(aap$eValuesMinEffi-1, 6)==0)

yLim <- c(min(c(metaE, metaEMinEffi)), max(c(metaE, metaEMinEffi)))

plot(metaEMinEffi, type="l", col="darkgoldenrod", ylim=yLim)
lines(1:length(aap$eValues), metaE, col="blue")
abline(h=log(1/alpha))
abline(h=log(betaMinEffi))

log(aap$eValuesMinEffi)[which(aap$eValues >= 1/alpha)]
log(aap$eValues)[which(aap$eValuesMinEffi <= betaMinEffi)]



# Aggregate analysis -------
## freq ----
sum(dat$na) + sum(dat$nb)

sumStatAll <- matrix(ncol=2, nrow=2)

sumStatAll[1, 1] <- sum(dat$ya)
sumStatAll[1, 2] <- sum(dat$yb)
sumStatAll[2, 1] <- sum(dat$na)-sum(dat$ya)
sumStatAll[2, 2] <- sum(dat$nb)-sum(dat$yb)

freqTest <- fisher.test(x=sumStatAll, alternative="two.sided")
freqTest <- fisher.test(x=sumStatAll, alternative=alternative)

# Percentage
sum(dat$ya)/sum(dat$na)

# Percentage
sum(dat$yb)/sum(dat$nb)
freqTest$p.value

# Cohen's d
round(log(freqTest$estimate)*sqrt(3)/pi, 2)
round(log(freqTest$conf.int)*sqrt(3)/pi, 2)


# OR
round(freqTest$estimate, 2)
round(freqTest$conf.int, 2)



## eValue ----
s10 <- saviTwoPropConditionalStat(
  ya = sum(dat$ya), na = sum(dat$na), nb = sum(dat$nb),
  n1 = sum(dat$ya)+sum(dat$yb),
  logOddsRatio = designObj$esMin,
  alternative = alternative)

r1f <- saviMinEffiTwoPropConditionalStat(
  ya = sum(dat$ya), na = sum(dat$na), nb = sum(dat$nb),
  n1 = sum(dat$ya)+sum(dat$yb),
  logOddsRatio = designObj[["minEffiTestResult"]][["parameter"]],
  alternative = alternative)

s10$eValue
r1f$eValue


# BiasedUrn check-------

na <- 17
nb <- 34
ya <- 8
n1 <- 41
yb <- n1-ya

n0 <- na+nb-n1

deltaS <- 0.3


nTotal <- na+nb

yMin <- max(0, n1-nb)
yMax <- min(n1, na)

numer <- lfactorial(na)-lfactorial(ya)-lfactorial(na-ya)+lfactorial(nb)-lfactorial(n1-ya)-lfactorial(nb-n1+ya)+deltaS*ya

k <- seq(yMin, yMax)

denom <- sum(exp(lfactorial(na)-lfactorial(k)-lfactorial(na-k)+lfactorial(nb)-lfactorial(n1-k)-lfactorial(nb-n1+k)+deltaS*k))

exp(numer)/denom


saviZTestStat

brie <- saviTwoPropConditionalStat(ya, na, nb, ya+yb, logOddsRatio=1)


sum(dat$na)+sum(dat$nb)

saviTwoPropConditionalStat(sum(dat$ya), sum(dat$na), sum(dat$nb), sum(dat$ya)+sum(dat$yb),
                           logOddsRatio=0.40)

saviMinEffiTwoPropConditionalStat(sum(dat$ya), sum(dat$na), sum(dat$nb), sum(dat$ya)+sum(dat$yb),
                           logOddsRatio=0.40)


#
