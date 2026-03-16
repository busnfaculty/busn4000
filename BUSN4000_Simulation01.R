# ==============================================================================
# SamErrSim.R - Sampling Error Simulation Engine
# BUSN 4000 - University of Georgia
# ==============================================================================
# This script is sourced by the student-facing script. It expects these 
# variables to already exist in the environment:
#   - Niter    : Number of iterations (samples to draw)
#   - N        : Sample size for each iteration
#   - Filename : Name of the CSV file (hosted on GitHub)
#   - ShowHist : (Optional) TRUE to display histograms, FALSE to skip
#   - Lag      : (Optional) Seconds to pause between each sample result (default 0.3)
# ==============================================================================

# --- Setup -------------------------------------------------------------------

Lag <- 0.5          # Seconds between each sample result (0 for instant)

# Set defaults if not specified
if (!exists("ShowHist")) ShowHist <- FALSE
if (!exists("Lag")) Lag <- 0.3

# Load data from GitHub
url <- "https://raw.githubusercontent.com/busnfaculty/busn4000/main/"
PopData <- read.csv(paste0(url, Filename))

PopN <- nrow(PopData)

# Validate inputs
if (N >= PopN) {
  stop(paste0("ERROR: Sample size N (", N, ") must be less than population size (", PopN, ")"))
}

if (N < 10) {
  warning("WARNING: Very small sample size. Results may be unstable.")
}

# --- Population Model --------------------------------------------------------

PopModel <- lm(SALES ~ ADV + BONUS, data = PopData)
PopCoefs <- coef(PopModel)

# --- Run Simulation ----------------------------------------------------------

# Storage for results
results <- matrix(NA, nrow = Niter, ncol = 3)
colnames(results) <- c("Beta0", "Beta1_ADV", "Beta2_BONUS")

# Set seed for reproducibility (students get same results each run)
set.seed(4000)

# Draw samples and fit models
for (i in 1:Niter) {
  samp <- PopData[sample(1:PopN, N, replace = FALSE), ]
  mod <- lm(SALES ~ ADV + BONUS, data = samp)
  results[i, ] <- coef(mod)
}

# --- Calculate Summary Statistics --------------------------------------------

means <- colMeans(results)
sds <- apply(results, 2, sd)

# --- Build Output Table ------------------------------------------------------

cat("\n")
cat("================================================================================\n")
cat("                    SAMPLING ERROR EXPLORATION RESULTS                          \n")
cat("================================================================================\n")
cat("\n")
cat("  Population Size (N):", PopN, "\n")
cat("  Sample Size (n):    ", N, "\n")
cat("  Number of Samples:  ", Niter, "\n")
cat("\n")
cat("================================================================================\n")
cat("\n")

# Header
cat(sprintf("%-16s %6s %12s %3s %12s %3s %12s\n", 
            "", "N", "Beta0", "", "Beta1(ADV)", "", "Beta2(BONUS)"))
cat("================================================================================\n")
cat("\n")

# Population model row
cat(sprintf("%-16s %6d %8s %8.3f %3s %8.3f%-4s %3s %8.3f%-4s\n",
            "Population Model", PopN, "y.hat =", PopCoefs[1], "+", PopCoefs[2], "ADV", "+", PopCoefs[3], "BONUS"))
cat("\n")

# Flush output so header appears immediately
flush.console()
Sys.sleep(Lag)

# Sample rows - print each with a lag
for (i in 1:Niter) {
  cat(sprintf("%-16s %6d %8s %8.3f %3s %8.3f%-4s %3s %8.3f%-4s\n",
              paste0("Sample #", i), N, "y.hat =", results[i,1], "+", results[i,2], "ADV", "+", results[i,3], "BONUS"))
  flush.console()
  if (i < Niter) Sys.sleep(Lag)  # Don't lag after the last one
}

cat("\n")

# Summary rows
cat(sprintf("%-16s %6s %12.3f %16.3f %16.3f\n", "Mean", "NA", means[1], means[2], means[3]))
cat(sprintf("%-16s %6s %12.3f %16.3f %16.3f\n", "Std Dev", "NA", sds[1], sds[2], sds[3]))
cat("\n")
cat("================================================================================\n")

# --- Additional Output: Compare to Theoretical SE ----------------------------

cat("\n")
cat("COMPARISON TO REGRESSION OUTPUT:\n")
cat("--------------------------------\n")
PopSE <- summary(PopModel)$coefficients[, "Std. Error"]
cat(sprintf("  Population Model Reported SEs:  SE(b0) = %.3f   SE(b1) = %.3f   SE(b2) = %.3f\n", 
            PopSE[1], PopSE[2], PopSE[3]))
cat(sprintf("  Simulation Std Devs:            SD(b0) = %.3f   SD(b1) = %.3f   SD(b2) = %.3f\n", 
            sds[1], sds[2], sds[3]))
cat("\n")
cat("  Note: The simulation SDs estimate sampling variability empirically.\n")
cat("        The reported SEs are theoretical estimates based on model assumptions.\n")
cat("        With many iterations, these values should be similar.\n")
cat("\n")

# --- Histograms (Optional) ---------------------------------------------------

if (ShowHist) {
  
  # Save current par settings and restore on exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(1, 3), mar = c(5, 4, 4, 2) + 0.1)
  
  # Beta0 histogram
  hist(results[, 1], 
       main = expression(paste("Distribution of ", hat(beta)[0])),
       xlab = expression(hat(beta)[0]),
       col = "lightblue",
       border = "white",
       breaks = min(15, Niter/2))
  abline(v = PopCoefs[1], col = "red", lwd = 2, lty = 2)
  abline(v = means[1], col = "blue", lwd = 2)
  legend("topright", 
         legend = c("Population", "Sample Mean"),
         col = c("red", "blue"), 
         lwd = 2, 
         lty = c(2, 1),
         cex = 0.8)
  
  # Beta1 histogram
  hist(results[, 2], 
       main = expression(paste("Distribution of ", hat(beta)[1], " (ADV)")),
       xlab = expression(hat(beta)[1]),
       col = "lightgreen",
       border = "white",
       breaks = min(15, Niter/2))
  abline(v = PopCoefs[2], col = "red", lwd = 2, lty = 2)
  abline(v = means[2], col = "blue", lwd = 2)
  legend("topright", 
         legend = c("Population", "Sample Mean"),
         col = c("red", "blue"), 
         lwd = 2, 
         lty = c(2, 1),
         cex = 0.8)
  
  # Beta2 histogram
  hist(results[, 3], 
       main = expression(paste("Distribution of ", hat(beta)[2], " (BONUS)")),
       xlab = expression(hat(beta)[2]),
       col = "lightyellow",
       border = "white",
       breaks = min(15, Niter/2))
  abline(v = PopCoefs[3], col = "red", lwd = 2, lty = 2)
  abline(v = means[3], col = "blue", lwd = 2)
  legend("topright", 
         legend = c("Population", "Sample Mean"),
         col = c("red", "blue"), 
         lwd = 2, 
         lty = c(2, 1),
         cex = 0.8)
  
  cat("Histograms displayed. Red dashed line = Population parameter. Blue solid line = Mean of samples.\n\n")
}

# --- Clean Up Environment ----------------------------------------------------
# Remove intermediate objects but keep results available for student exploration

rm(i, samp, mod, url)

cat("================================================================================\n")
cat("  Objects available for further exploration:\n")
cat("    - PopData    : The full population dataset\n")
cat("    - PopModel   : The population regression model (use summary(PopModel))\n")
cat("    - results    : Matrix of all sample coefficients\n")
cat("================================================================================\n")
