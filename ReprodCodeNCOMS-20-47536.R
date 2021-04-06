# ============================= CODE METADATA ================================ #
# DATE CREATED: 2020-12-10
# DESCRIPTION: Demo code for reproducibility of manuscript NCOMS-20-47536, 
#              Nature Communications
# ============================== START CODE ================================== #
# ======================== #
# ==== GENERAL SETUP: ====
# ======================== #
# Working directory:
setwd("Directory where all files are saved")

# Install BaSTA.ZIMS package:
install.packages("BaSTA.ZIMS_1.0.2.tar.gz",
                 repos = NULL, type = "source")

# Install package snowfall for parallel computing:
install.packages("snowfall")

Start <- Sys.time()
# Load libraries:
library(BaSTA.ZIMS)
library(snowfall)

# Load survival test dataset:
dat <- read.table("testData.txt", sep = "\t", header = TRUE)

# ==================== #
# ==== FUNCTIONS: ====
# ==================== #
# Mortality or hazards rate:
mux <- function(x, th) {
  exp(th["a0"] - th["a1"] * x) + th["c"] + exp(th["b0"] + th["b1"] * x)
}

# Cumulative hazards:
Ux <- function(x, th) {
  exp(th["a0"]) / th["a1"] * (1 - exp(- th["a1"] * x)) + th["c"] * x +
    exp(th["b0"]) / th["b1"] * (exp(th["b1"] * x) - 1)
}

# Cumulative survival:
Sx <- function(x, th) exp(- Ux(x, th))

# Life expectancy:
ex <- function(x, th, dx) sum(Sx(x, th) * dx) / Sx(x[1], th)

# Life table entropy:
Hx <- function(x, th, dx) {
  Sx1 <- Sx(x, th)
  Sx1 <- Sx1[Sx1 > 0]
  -sum(Sx1 * log(Sx1) * dx) / ex(x, th, dx)
}

# Lifespan equality (as in Colchero et al 2016 PNAS):
epx <- function(x, th, dx) - log(Hx(x, th, dx))

# Partial derivative of survival with respect to theta_j:
Sth <- function(x, th, idth = "a0") {
  if (idth == "a0") {
    sth <- Sx(x, th) * exp(th["a0"]) / th["a1"] * (exp(-th["a1"] * x) - 1)
  } else if (idth == "a1") {
    sth <- Sx(x, th) * (exp(th["a0"]) / th["a1"]) * 
      (1 / th["a1"] - exp(-th["a1"] * x) * (x + 1 / th["a1"]))
  } else if (idth == "c") {
    sth <- - Sx(x, th) * x
  } else if (idth == "b0") {
    sth <- Sx(x, th) * exp(th["b0"]) / th["b1"] * (1 - exp(th["b1"] * x))
  } else if (idth == "b1") {
    sth <- Sx(x, th) * exp(th["b0"]) / (th["b1"]^2) * 
      (exp(th['b1'] * x) * (1 / th["b1"] - x) - 1 / th["b1"])
  } else {
    stop("Wrong 'idth'. Values should be, 'a0', 'a1', 'c', 'b0', or 'b1'.\n", 
         call. = FALSE)
  }
  return(sth)
}

# Sensitivity or elasticities of e_0 or epsilon_0 to theta_j:
CalcSensElastPaceShape <- function(x, th, dx, idth = "a0", lifeEqual = FALSE,
                                   sensLifeExp = NULL, elast = FALSE) {
  # Calculate sensitivity of e0 to th_j:
  if (is.null(sensLifeExp)) {
    sensLifeExp <- sum(Sth(x, th, idth) * dx)
  }
  # Calculate sensitivity of epsilon to th_j:
  if (lifeEqual) {
    HxInv <- 1 / Hx(x, th, dx)
    sensPS <- (sensLifeExp * (1 + HxInv) -
                 HxInv * sum(Sth(x, th, idth) * Ux(x, th) * dx)) / 
      ex(x, th, dx)
    # Calculate elasticities:
    if (elast) {
      sensPS <- th[idth] / epx(x, th, dx) * sensPS
    } 
  } else {
    sensPS <- sensLifeExp
    # Calculate elasticities:
    if (elast) {
      sensPS <- th[idth] / ex(x, th, dx) * sensPS
    }
  }
  return(sensPS)
}

# Gradient vector for sensitivities of e_0 or epsilon_0:
CalcGradientPS <- function(x, th, dx, ...) {
  grad <- rep(NA, 5)
  parname <- c("a0", "a1", "c", "b0", "b1")
  names(grad) <- parname
  for (thi in 1:5) {
    grad[thi] <- CalcSensElastPaceShape(x, th, dx, idth = parname[thi], ...)
  }
  return(grad)
}

# Sensitivities of theta_j to e_0 or epsilon_0:
CalcSensTh <- function(x, th, dx, idth = "a0", lifeEqual = FALSE, 
                       sensLifeExp = NULL, propSens = FALSE, absVal = FALSE) {
  sensTh <- 1 / CalcSensElastPaceShape(x = x, th = th, dx = dx, idth = idth, 
                                       lifeEqual = lifeEqual,
                                       sensLifeExp = sensLifeExp)
  if (propSens) sensTh <- 1 / th[idth] * sensTh
  if (absVal) sensTh <- abs(sensTh)
  return(sensTh)
}

# Gradient vector for sensitivities of th_j:
CalcGradientTh <- function(x, th, dx, ...) {
  grad <- rep(NA, 2)
  parname <- c("a0", "a1", "c", "b0", "b1")
  psname <- c("e0", "eps0")
  grad <- matrix(NA, 2, 5, dimnames = list(psname, parname))
  for (thi in 1:5) {
    idth <- parname[thi]
    for (psi in 1:2) {
      lifeEqual <- ifelse(psi == 1, FALSE, TRUE)
      if (psi == 1) {
        sensLifeExp <- NULL
      } else {
        sensLifeExp <- grad[1]
      }
      grad[psi, thi] <- CalcSensTh(x, th, dx, idth = idth, 
                                   lifeEqual = lifeEqual, 
                                   sensLifeExp = sensLifeExp, ...)
    }
    
  }
  return(grad)
}

# ============================ #
# ==== SURVIVAL ANALYSIS: ====
# ============================ #
# Run BaSTA on test data for Siler model:
out <- bastaZIMS(dat, model = "GO", shape = "bathtub", ncpus = 4, nsim = 4,
                 niter = 21001, burnin = 1000, parallel = TRUE)

# Check goodness of fit:
plot(out, plot.type = "gof")

# Extract Siler mortality parameters:
theta <- out$coefficients[, 1]

# ================================================= #
# ==== SENSITIVITY OF LIFE EXPECT. AND LIFESPAN EQ. 
#      TO MORTALITY PARAMETERS: ====
# ================================================= #
# ------------------------------------------- #
# Life expectancy and lifespan equality: ----
# ------------------------------------------- #
# Vector of age intervals:
dx <- 0.001
xv <- seq(0, 200, dx)

# Calculate life expectancy:
e0 <- ex(x = xv, dx = dx, th = theta)

# Calculate lifespan equality:
eps0 <- epx(x = xv, dx = dx, th = theta)

# ------------------------------- #
# Examples of sensitivities: ----
# ------------------------------- #
# Sensitivity of life expectancy to a1:
sensE0a1 <- CalcSensElastPaceShape(x = xv, dx = dx, th = theta, idth = "a1")

# Sensitivity of lifespan equality to a1:
sensEps0a1 <- CalcSensElastPaceShape(x = xv, dx = dx, th = theta, idth = "a1",
                                     lifeEqual = TRUE, sensLifeExp = sensE0a1)

# Sensitivity of life expectancy to b1:
sensE0b1 <- CalcSensElastPaceShape(x = xv, dx = dx, th = theta, idth = "b1")

# Sensitivity of lifespan equality to b1:
sensEps0b1 <- CalcSensElastPaceShape(x = xv, dx = dx, th = theta, idth = "b1",
                                     lifeEqual = TRUE, sensLifeExp = sensE0b1)

# Print to the console the results:
cat(sprintf("Sensitivity of life expectancy to:\na1 = %s\nb1 = %s\n", 
            round(abs(sensE0a1), 2), round(abs(sensE0b1), 2)))

cat(sprintf("Sensitivity of lifespan equality to:\na1 = %s\nb1 = %s\n", 
            round(abs(sensEps0a1), 2), round(abs(sensEps0b1), 2)))

# -------------------------- #
# Full gradient vector: ----
# -------------------------- #
# Gradient vector (columns for life expectancy and lifespan equality,
#                  rows for each Siler parameter):
gradVec01 <- cbind(e0 = CalcGradientPS(x = xv, dx = dx, th = theta),
                   eps0 = CalcGradientPS(x = xv, dx = dx, th = theta,
                                         lifeEqual = TRUE))

# Plot vectors:
xlim <- c(0, 40) 
ylim <- eps0 + c(-1, 1)
parCols <- c(a0 = "#A6CEE3", a1 = "#1F78B4", c = "#FF7F00", b0 = "#B2DF8A", 
             b1 = "#33A02C")

par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
plot(xlim, ylim, col = NA, xlab = "Life expectancy", 
     ylab = "lifespan equality")
for (thi in 1:5) {
  arrows(x0 = e0, y0 = eps0, x1 = e0 + gradVec01[thi, 1],
         y1 = eps0 + gradVec01[thi, 2], col = parCols[thi], lwd = 2,
         length = 0.05, angle = 10)
}
legend('bottomright', legend = names(theta), col = parCols, lwd = 2, pch = NA,
       bty = "n")

# ====================================================== #
# ==== PROPORTIONAL SENSITIVITY OF MORTALITY PARAMETERS  
#      TO LIFE EXPECT. AND LIFESPAN EQ.: ====
# ====================================================== #
# -------------------------------------------- #
# Examples of proportional sensitivities: ----
# -------------------------------------------- #
# Sensitivity of a1 to life expectancy:
sensa1E0 <- CalcSensTh(x = xv, dx = dx, th = theta, idth = "a1", 
                       propSens = TRUE)

# Sensitivity of a1 to lifespan equality:
sensa1Eps0 <- CalcSensTh(x = xv, dx = dx, th = theta, idth = "a1", 
                       lifeEqual = TRUE, sensLifeExp = sensa1E0, 
                       propSens = TRUE)

# Sensitivity of b1 to life expectancy:
sensb1E0 <- CalcSensTh(x = xv, dx = dx, th = theta, idth = "b1", 
                       propSens = TRUE)

# Sensitivity of b1 to lifespan equality:
sensb1Eps0 <- CalcSensTh(x = xv, dx = dx, th = theta, idth = "b1", 
                       lifeEqual = TRUE, sensLifeExp = sensb1E0, 
                       propSens = TRUE)

# Print to the console the results:
cat(sprintf("Sensitivity to life expectancy of:\na1 = %s\nb1 = %s\n", 
            round(abs(sensa1E0), 4), round(abs(sensb1E0), 4)))

cat(sprintf("Sensitivity to lifespan equality of:\na1 = %s\nb1 = %s\n", 
            round(abs(sensa1Eps0), 2), round(abs(sensb1Eps0), 2)))

# -------------------------- #
# Full gradient vector: ----
# -------------------------- #
# Gradient vector of proportional changes in parameters
# (columns for each Siler parameter, rows for life
# expectancy and lifespan equality):
gradVec02 <- CalcGradientTh(x = xv, dx = dx, th = theta, propSens = TRUE)

# Plot of changes in each parameter as a function of changes in 
# life expectancy and lifespan equality:
xlim <- max(abs(gradVec02["e0", ])) * c(-1, 1)
ylim <- max(abs(gradVec02["eps0", ])) * c(-1, 1)
par(mfrow = c(3, 2), mar = c(2, 2, 1, 1))
for (thi in 1:5) {
  plot(xlim, ylim, col = NA, xlab = "", ylab = "")
  arrows(x0 = 0, y0 = 0, x1 = gradVec02["e0", thi], 
         y1 = gradVec02["eps0", thi], length = 0.05, angle = 10,
         col = 'dark red', lwd = 2)
  mtext(names(theta)[thi], side = 3, line = -2)
}

End <- Sys.time()
End - Start
# ========================= END OF CODE ====================================== #

# save(list = c("out", "dat", "xv", "dx", "e0", "eps0", "sensE0a1", "sensE0b1",
#               "sensEps0a1", "sensEps0b1", "sensa1E0", "sensb1E0",
#               "sensa1Eps0", "sensb1Eps0", "gradVec01", "gradVec02"),
#      file = "~/FERNANDO/PROJECTS/1.ACTIVE/PrimatePaceShape/Documents/Manuscript/Submissions/READMEforReprodCode/test.RData")
