# Table 1 for the "Multiplicative errors-in-variables beta regression" paper.
rm(list = ls())

require(betareg)
require(numDeriv)

source("npl.mle.R")
source("rc.mle.R")

# The National Agency of Supplementary Health data

ans <- readRDS("ans.Rds")
head(ans)

# Multiplicative measurement error mean:

mu.e <- 0.65 # 0.80, 0.95

#  Naïve method

summary(betareg(y~z+w|1,data=ans))

# Pseudo-likelihood method

npl.mle(y~z+w,data=ans,mu.e=mu.e)


# Regression calibration method

rc.mle(y~z+w,data=ans,mu.e=mu.e,B=300)

