#BUSN 4000 - V.1.0.2 - SPR26

regprint <- function(model, conf_level = 95, digits = NULL,
                     robust = c("none","HC0","HC1","HC2","HC3","HC4","HC4m","HC5"),
                     std_beta = FALSE) {

  # ---- validations ----
  if (!inherits(model, "lm")) {
    obj <- deparse(substitute(model))
    stop(sprintf("Error: %s is not a lm object. regprint() requires you fit a linear model first using the lm() function and save it as an object.", obj), call. = FALSE)
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 || !is.finite(conf_level) ||
      conf_level < 1 || conf_level > 99) {
    stop("Error: 'conf_level' must be a single number between 1 and 99 (e.g., 90, 95).", call. = FALSE)
  }
  if (!is.null(digits)) {
    if (!is.numeric(digits) || length(digits) != 1 || !is.finite(digits) ||
        digits < 0 || digits != as.integer(digits)) {
      stop("Error: 'digits' must be a single nonnegative integer (e.g., 3, 4, 5).", call. = FALSE)
    }
  }
  robust <- match.arg(robust)

  level <- conf_level / 100
  lb_name <- sprintf("%d%% CI LB", conf_level)
  ub_name <- sprintf("%d%% CI UB", conf_level)

  # ---- formatting helpers (with thousands separators) ----
  k_p  <- if (is.null(digits)) 4L else digits
  k_r2 <- if (is.null(digits)) 4L else digits
  k_oth<- if (is.null(digits)) 3L else digits

  fmt_fix <- function(x, k) {
    ifelse(is.na(x), "NA",
           formatC(x, format = "f", digits = k, big.mark = ","))
  }
  fmt4  <- function(x) fmt_fix(x, k_r2)   # Multiple R, R^2, Adj R^2
  fmt3  <- function(x) fmt_fix(x, k_oth)  # Everything else
  fmt_p <- function(p) {
    thr <- 10^(-k_p)
    below <- paste0("<", formatC(thr, format="f", digits=k_p, big.mark=""))
    ifelse(p < thr, below, formatC(p, format="f", digits=k_p, big.mark=""))
  }
  fmt_int <- function(x) {
    ifelse(is.na(x), "NA", formatC(as.integer(x), format = "d", big.mark = ","))
  }

  # ---- core pieces ----
  y_name <- as.character(formula(model))[2]
  sm <- summary(model)

  r2 <- sm$r.squared; adjr <- sm$adj.r.squared
  multR <- sqrt(r2); s_err <- sm$sigma
  n_used <- length(sm$residuals)
  n_miss <- if (is.null(model$na.action)) 0L else length(model$na.action)
  n_read <- n_used + n_miss

  mf <- model.frame(model)
  y  <- model.response(mf)
  e  <- resid(model)
  sse <- sum(e^2); sst <- sum((y - mean(y))^2); ssr <- sst - sse

  df_reg <- unname(sm$fstatistic["numdf"])
  df_res <- unname(sm$fstatistic["dendf"])
  df_tot <- df_reg + df_res
  ms_reg <- ssr / df_reg; ms_res <- sse / df_res
  fstat  <- unname(sm$fstatistic["value"])
  p_f    <- pf(fstat, df_reg, df_res, lower.tail = FALSE)

  # ---- coefficients (classic vs robust) ----
  coefs <- coef(model)
  if (robust == "none") {
    se    <- sm$coefficients[, "Std. Error"]
    tvals <- sm$coefficients[, "t value"]
    pvals <- sm$coefficients[, "Pr(>|t|)"]
    ci    <- confint(model, level = level)
  } else {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Robust SEs require package 'sandwich'. Install with install.packages('sandwich').", call. = FALSE)
    }
    V  <- sandwich::vcovHC(model, type = robust)
    se <- sqrt(diag(V))
    tvals <- coefs / se
    pvals <- 2 * pt(abs(tvals), df = df_res, lower.tail = FALSE)
    tcrit <- qt(1 - (1 - level)/2, df = df_res)
    ci <- cbind(coefs - tcrit * se, coefs + tcrit * se)
    colnames(ci) <- c("low","high")
  }

  ctab <- data.frame(
    "Coefficients"    = coefs,
    "Standard Error"  = se,
    "t-test"          = tvals,
    "p-value"         = pvals,
    check.names = FALSE
  )
  ctab[[lb_name]] <- ci[, 1]
  ctab[[ub_name]] <- ci[, 2]
  rownames(ctab) <- sub("^\\(Intercept\\)$", "Intercept", rownames(ctab))

  # ---- standardized betas ----
  if (isTRUE(std_beta)) {
    X <- model.matrix(model)
    sdy <- stats::sd(y)
    sdx <- apply(X, 2, stats::sd)
    stdb <- rep(NA_real_, length(coefs))
    idx <- which(colnames(X) != "(Intercept)")
    if (length(idx)) stdb[idx] <- coefs[idx] * sdx[idx] / sdy
    names(stdb) <- rownames(ctab)
    ctab$`Std. Beta` <- stdb
  }

  # ---- VIFs (exclude intercept) ----
  X <- model.matrix(model)
  pred_cols <- setdiff(colnames(X), "(Intercept)")
  Xp <- if (length(pred_cols)) X[, pred_cols, drop = FALSE] else NULL

  vif <- rep(NA_real_, nrow(ctab))  # NA for Intercept
  if (!is.null(Xp)) {
    if (ncol(Xp) == 1L) {
      vif[match(colnames(Xp), rownames(ctab))] <- 1
    } else {
      v <- numeric(ncol(Xp))
      for (j in seq_len(ncol(Xp))) {
        r2j <- summary(lm(Xp[, j] ~ Xp[, -j]))$r.squared
        v[j] <- 1 / (1 - r2j)
      }
      vif[match(colnames(Xp), rownames(ctab))] <- v
    }
  }
  ctab$VIF <- vif

  # ---- printing (no warnings, explicit formats) ----
  obj_name <- deparse(substitute(model))
  cat(sprintf("Linear Regression Model: %s\n", obj_name))
  cat(sprintf("Dependent Variable: %s\n\n", y_name))

  cat("Regression Statistics Table\n")
  cat(sprintf("%-22s %8s\n", "Multiple R",        fmt4(multR)))
  cat(sprintf("%-22s %8s\n", "R Square",          fmt4(r2)))
  cat(sprintf("%-22s %8s\n", "Adjusted R Square", fmt4(adjr)))
  cat(sprintf("%-22s %8s\n", "Standard Error",    fmt3(s_err)))
  cat(sprintf("%-22s %8s\n", "N Obs Read",        fmt_int(n_read)))
  cat(sprintf("%-22s %8s\n", "N Obs Missing",     fmt_int(n_miss)))
  cat(sprintf("%-22s %8s\n\n","N Obs Used",       fmt_int(n_used)))

  cat("ANOVA Table\n")
  cat(sprintf("%-12s %6s %14s %14s %14s %10s\n",
              "Source","df","SS","MS","F-Statistic","p-value"))
  cat(sprintf("%-12s %6s %14s %14s %14s %10s\n",
              "Regression", fmt_int(df_reg), fmt3(ssr), fmt3(ms_reg), fmt3(fstat), fmt_p(p_f)))
  cat(sprintf("%-12s %6s %14s %14s %14s %10s\n",
              "Residual",   fmt_int(df_res), fmt3(sse), fmt3(ms_res), "", ""))
  cat(sprintf("%-12s %6s %14s %14s %14s %10s\n\n",
              "Total",      fmt_int(df_tot), fmt3(sst), "", "", ""))

  cat("Coefficients Table\n")
  if (isTRUE(std_beta)) {
    cat(sprintf("%-12s %12s %15s %10s %10s %12s %12s %12s %8s\n",
                "Variables","Coefficients","Standard Error","t-test","p-value",
                paste0(conf_level, "% CI LB"), paste0(conf_level, "% CI UB"),
                "Std. Beta","VIF"))
  } else {
    cat(sprintf("%-12s %12s %15s %10s %10s %12s %12s %8s\n",
                "Variables","Coefficients","Standard Error","t-test","p-value",
                paste0(conf_level, "% CI LB"), paste0(conf_level, "% CI UB"),
                "VIF"))
  }

  for (i in seq_len(nrow(ctab))) {
    if (isTRUE(std_beta)) {
      cat(sprintf("%-12s %12s %15s %10s %10s %12s %12s %12s %8s\n",
                  rownames(ctab)[i],
                  fmt3(ctab$Coefficients[i]),
                  fmt3(ctab$`Standard Error`[i]),
                  fmt3(ctab$`t-test`[i]),
                  fmt_p(ctab$`p-value`[i]),
                  fmt3(ctab[[lb_name]][i]),
                  fmt3(ctab[[ub_name]][i]),
                  ifelse(is.null(ctab$`Std. Beta`[i]) || is.na(ctab$`Std. Beta`[i]), "", fmt3(ctab$`Std. Beta`[i])),
                  ifelse(is.na(ctab$VIF[i]), "", fmt3(ctab$VIF[i]))))
    } else {
      cat(sprintf("%-12s %12s %15s %10s %10s %12s %12s %8s\n",
                  rownames(ctab)[i],
                  fmt3(ctab$Coefficients[i]),
                  fmt3(ctab$`Standard Error`[i]),
                  fmt3(ctab$`t-test`[i]),
                  fmt_p(ctab$`p-value`[i]),
                  fmt3(ctab[[lb_name]][i]),
                  fmt3(ctab[[ub_name]][i]),
                  ifelse(is.na(ctab$VIF[i]), "", fmt3(ctab$VIF[i]))))
    }
  }

  invisible(list(
    coefficients = ctab,
    anova = data.frame(
      Source = c("Regression","Residual","Total"),
      df = c(df_reg, df_res, df_tot),
      SS = c(ssr, sse, sst),
      MS = c(ms_reg, ms_res, NA_real_),
      `F-Statistic` = c(fstat, NA_real_, NA_real_),
      `p-value` = c(p_f, NA_real_, NA_real_),
      check.names = FALSE
    ),
    summary = list(
      Multiple_R = multR, R_Square = r2, Adj_R_Square = adjr,
      Standard_Error = s_err, N_Read = n_read, N_Missing = n_miss, N_Used = n_used
    ),
    robust = robust,
    conf_level = conf_level,
    digits = digits,
    std_beta = std_beta
  ))
}
