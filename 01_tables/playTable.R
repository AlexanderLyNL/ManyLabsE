# Setup paths ----
sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/01_tables/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))

# Data load ----

fileNeem <- "Tversky_1_clean_tables.csv"

dat <- read.csv(file=paste0(myWd, "data/", fileNeem))

dim(dat)
head(dat)
dat$sites

# allSources <- checkEnoughDataInTable(dat)$allSources
allSources <- unique(dat$sites)
dat <- dat[dat$sites %in% allSources, ]

length(allSources)


# Here -------
alpha <- 0.05
betaFutility <- alpha
alternative <- "greater"
esMin <- 1.08
futParam <- esMin
logOddsRatio <- esMin

designObj <- list(esMin=esMin, futilityResult=list(parameter=futParam), alternative=alternative)

aap <- scenario1Table(dat, allSources, designObj)

cumsum(log(aap$eValues)) >= 1/alpha
cumsum(log(aap$eValuesFut)) <= betaFutility

plot(log(aap$eValuesFut), type="l")
lines(1:length(aap$eValues), log(aap$eValues), col="blue")

# plot(round(aap$eValues-aap$eValuesFut, 2), type="l"

# Total sample size -------
sum(dat$na) + sum(dat$nb)
s10 <- saviTwoPropConditionalStat(
  ya = sum(dat$ya), na = sum(dat$na), nb = sum(dat$nb),
  n1 = sum(dat$ya)+sum(dat$yb),
  logOddsRatio = designObj$esMin,
  alternative = alternative)

r1f <- saviFutilityTwoPropConditionalStat(
  ya = sum(dat$ya), na = sum(dat$na), nb = sum(dat$nb),
  n1 = sum(dat$ya)+sum(dat$yb),
  logOddsRatio = designObj$futilityResult$parameter,
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

saviFutilityTwoPropConditionalStat(sum(dat$ya), sum(dat$na), sum(dat$nb), sum(dat$ya)+sum(dat$yb),
                           logOddsRatio=0.40)


