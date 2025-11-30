# BUSN 4000 Tool Kit
# V1.01 - SPR 2026


# ================================================================================
# ================================================================================


# Standardized Residuals

sresid <- function(x) {
  if (!inherits(x, "lm")) {
    stop('The sresid() function requires a single lm() object as input. "',
         deparse(substitute(x)), '" is not an lm() object.')
  }
  resid(x) / summary(x)$sigma
}


# ================================================================================
# ================================================================================


logY_predict <- function(model,
                         newdata = NULL,
                         na.action = na.pass) {
  # Check if model is an lm object
  if (!inherits(model, "lm")) {
    stop("The 'model' argument must be an lm object created with lm()")
  }
  
  # linear predictor on log scale: η̂ = Xβ̂
  eta_hat <- stats::predict(model, newdata = newdata, na.action = na.action)
  # residual standard error from the log model
  sigma2 <- summary(model)$sigma^2
  # lognormal mean correction: exp(η̂ + σ̂²/2)
  y_hat <- exp(eta_hat + 0.5 * sigma2)
  y_hat
}


# ================================================================================
# ================================================================================


#regprint - V.1.0.9

regprint <- function(model, CL = 95, digits = NULL,
                     robust = c("none","HC0","HC1","HC2","HC3","HC4","HC4m","HC5"),
                     std_beta = FALSE, diag = FALSE) {

  # ---- validations ----
  is_lm <- inherits(model, "lm") && !inherits(model, "glm")
  is_glm <- inherits(model, "glm")
  
  if (!is_lm && !is_glm) {
    obj <- deparse(substitute(model))
    stop(sprintf("Error: %s is not a lm or glm object. regprint() requires you fit a linear or generalized linear model first using the lm() or glm() function and save it as an object.", obj), call. = FALSE)
  }
  
  # Check if glm is binary logistic regression
  if (is_glm) {
    if (!(family(model)$family == "binomial" && family(model)$link == "logit")) {
      stop("Error: regprint() currently only supports binary logistic regression for glm objects (family = binomial with logit link).", call. = FALSE)
    }
  }
  
  if (!is.numeric(CL) || length(CL) != 1 || !is.finite(CL) ||
      CL < 1 || CL > 99) {
    stop("Error: 'CL' must be a single number between 1 and 99 (e.g., 90, 95).", call. = FALSE)
  }
  if (!is.null(digits)) {
    if (!is.numeric(digits) || length(digits) != 1 || !is.finite(digits) ||
        digits < 0 || digits != as.integer(digits)) {
      stop("Error: 'digits' must be a single nonnegative integer (e.g., 3, 4, 5).", call. = FALSE)
    }
  }
  robust <- match.arg(robust)

  level <- CL / 100
  lb_name <- sprintf("%d%% CI LB", CL)
  ub_name <- sprintf("%d%% CI UB", CL)

  # ---- formatting helpers (thousands separators) ----
  k_p  <- if (is.null(digits)) 4L else digits
  k_r2 <- if (is.null(digits)) 4L else digits
  k_oth<- if (is.null(digits)) 3L else digits

  fmt_fix <- function(x, k) ifelse(is.na(x), "NA", formatC(x, format="f", digits=k, big.mark=","))
  fmt4  <- function(x) fmt_fix(x, k_r2)   # Multiple R, R^2, Adj R^2
  fmt3  <- function(x) fmt_fix(x, k_oth)  # Everything else
  fmt_p <- function(p) {
    thr <- 10^(-k_p)
    below <- paste0("<", formatC(thr, format="f", digits=k_p, big.mark=""))
    ifelse(p < thr, below, formatC(p, format="f", digits=k_p, big.mark=""))
  }
  fmt_int <- function(x) ifelse(is.na(x), "NA", formatC(as.integer(x), format="d", big.mark=","))

  # ---- core pieces ----
  y_name <- as.character(formula(model))[2]
  sm <- summary(model)

  # Get number of observations used
  if (is_glm) {
    n_used <- length(model$y)
  } else {
    n_used <- length(sm$residuals)
  }
  n_miss <- if (is.null(model$na.action)) 0L else length(model$na.action)
  n_read <- n_used + n_miss

  mf <- model.frame(model)
  y  <- model.response(mf)
  e  <- resid(model)
  
  if (is_lm) {
    r2 <- sm$r.squared; adjr <- sm$adj.r.squared
    multR <- sqrt(r2); s_err <- sm$sigma
    sse <- sum(e^2); sst <- sum((y - mean(y))^2); ssr <- sst - sse

    df_reg <- unname(sm$fstatistic["numdf"])
    df_res <- unname(sm$fstatistic["dendf"])
    df_tot <- df_reg + df_res
    ms_reg <- ssr / df_reg; ms_res <- sse / df_res
    fstat  <- unname(sm$fstatistic["value"])
    p_f    <- pf(fstat, df_reg, df_res, lower.tail = FALSE)
  } else if (is_glm) {
    df_reg <- sm$df.null - sm$df.residual
    df_res <- sm$df.residual
  }

  # ---- coefficients (classic vs robust) ----
  coefs <- coef(model)
  
  # Check for NA coefficients (perfect multicollinearity / dummy variable trap)
  if (any(is.na(coefs))) {
    msg <- "One or more of the variables selected in the X Range are a perfect linear combination of others (i.e., the design matrix is not full rank). Please respecify the model."
    stop(paste(strwrap(msg), collapse = "\n"), call. = FALSE)
  }
  
  if (robust == "none") {
    se    <- sm$coefficients[, "Std. Error"]
    if (is_lm) {
      tvals <- sm$coefficients[, "t value"]
      pvals <- sm$coefficients[, "Pr(>|t|)"]
    } else if (is_glm) {
      tvals <- sm$coefficients[, "z value"]
      pvals <- sm$coefficients[, "Pr(>|z|)"]
    }
    ci    <- suppressMessages(confint(model, level = level))
  } else {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      stop("Robust SEs require package 'sandwich'. Install with install.packages('sandwich').", call. = FALSE)
    }
    V  <- sandwich::vcovHC(model, type = robust)
    se <- sqrt(diag(V))
    tvals <- coefs / se
    if (is_lm) {
      pvals <- 2 * pt(abs(tvals), df = df_res, lower.tail = FALSE)
      tcrit <- qt(1 - (1 - level)/2, df = df_res)
    } else if (is_glm) {
      pvals <- 2 * pnorm(abs(tvals), lower.tail = FALSE)
      tcrit <- qnorm(1 - (1 - level)/2)
    }
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

  # helper: build format string with ≥10% spacing and a 3-space minimum
  build_fmt <- function(widths) {
    if (length(widths) == 1L) return(paste0("%", widths, "s\n"))
    sp <- pmax(3L, ceiling(0.10 * widths[-1]))
    paste0("%-", widths[1], "s",
           paste(mapply(function(s, w) paste0(strrep(" ", s), "%", w, "s"),
                        sp, widths[-1]), collapse = ""),
           "\n")
  }

  # ---- printing ----
  obj_name <- deparse(substitute(model))
  
  if (is_lm) {
    cat(sprintf("Linear Regression Model: %s\n", obj_name))
    cat(sprintf("Dependent Variable: %s\n\n", y_name))

    # Regression Statistics
    cat("Regression Statistics Table\n")
    rs_names <- c("Multiple R","R Square","Adjusted R Square",
                  "Standard Error","N Obs Read","N Obs Missing","N Obs Used")
    rs_vals  <- c(fmt4(multR), fmt4(r2), fmt4(adjr), fmt3(s_err),
                  fmt_int(n_read), fmt_int(n_miss), fmt_int(n_used))
    rs_w <- c(max(nchar(rs_names)), max(nchar(rs_vals)))
    rs_fmt <- build_fmt(rs_w)
    for (i in seq_along(rs_names)) cat(sprintf(rs_fmt, rs_names[i], rs_vals[i]))
    cat("\n")

    # ANOVA
    cat("ANOVA Table\n")
    a_hdr <- c("Source","df","SS","MS","F-Statistic","p-value")
    a_src <- c("Regression","Residual","Total")
    a_df  <- c(fmt_int(df_reg), fmt_int(df_res), fmt_int(df_tot))
    a_ss  <- c(fmt3(ssr), fmt3(sse), fmt3(sst))
    a_ms  <- c(fmt3(ms_reg), fmt3(ms_res), "")
    a_F   <- c(fmt3(fstat), "", "")
    a_p   <- c(fmt_p(p_f), "", "")
    a_cols <- list(a_src, a_df, a_ss, a_ms, a_F, a_p)
    a_w <- integer(length(a_cols))
    for (i in seq_along(a_cols)) a_w[i] <- max(nchar(a_hdr[i]), max(nchar(a_cols[[i]]), na.rm = TRUE))
    a_fmt <- build_fmt(a_w)
    cat(do.call(sprintf, c(a_fmt, as.list(a_hdr))))
    for (i in seq_along(a_src)) {
      row <- list(a_src[i], a_df[i], a_ss[i], a_ms[i], a_F[i], a_p[i])
      cat(do.call(sprintf, c(a_fmt, row)))
    }
    cat("\n")
  } else if (is_glm) {
    cat(sprintf("Logistic Regression Model: %s\n", obj_name))
    cat(sprintf("Dependent Variable: %s\n\n", y_name))
    
    # Dependent Variable Information
    cat("Dependent Variable Information\n")
    y_table <- table(model$y)
    # Sort so that 1 appears first, then 0
    sort_order <- order(as.numeric(names(y_table)), decreasing = TRUE)
    y_table <- y_table[sort_order]
    y_coding <- as.numeric(names(y_table))
    event_idx <- which.max(y_coding)
    dv_hdr <- c("Variable", "Value", "Coding", "Frequency", "Percent", "")
    dv_var <- c(y_name, rep("", length(y_table) - 1), "Total (N Used)")
    dv_val <- c(names(y_table), "")
    dv_code <- c(names(y_table), "")
    dv_freq <- c(fmt_int(y_table), fmt_int(n_used))
    dv_pct <- c(fmt3(100 * y_table / sum(y_table)), fmt3(100.00))
    dv_event <- c(rep("", length(y_table)), "")
    dv_event[event_idx] <- "(Event)"
    
    dv_cols <- list(dv_var, dv_val, dv_code, dv_freq, dv_pct, dv_event)
    dv_w <- integer(length(dv_cols))
    for (i in seq_along(dv_cols)) dv_w[i] <- max(nchar(dv_hdr[i]), max(nchar(dv_cols[[i]]), na.rm = TRUE))
    dv_fmt <- build_fmt(dv_w)
    cat(do.call(sprintf, c(dv_fmt, as.list(dv_hdr))))
    for (i in seq_along(dv_var)) {
      row <- list(dv_var[i], dv_val[i], dv_code[i], dv_freq[i], dv_pct[i], dv_event[i])
      cat(do.call(sprintf, c(dv_fmt, row)))
    }
    cat("\n")
    
    # Convergence message
    cat(sprintf("**Convergence Criterion Satisfied: Model Converged in %d iterations.**\n\n", 
                model$iter))
    
    # Model Fit Statistics
    cat("Model Fit Statistics\n")
    null_dev <- model$null.deviance
    resid_dev <- model$deviance
    aic_val <- AIC(model)
    
    # Calculate BIC manually: -2*logLik + log(n)*k
    k_full <- length(coefs)
    bic_val <- resid_dev + log(n_used) * k_full
    
    # For intercept-only model: k=1 parameter
    aic_null <- null_dev + 2 * 1
    bic_null <- null_dev + log(n_used) * 1
    
    mf_hdr <- c("Criterion", "Intercept Only", "Intercept and Covariates")
    mf_crit <- c("AIC", "SBC", "-2 log L")
    mf_int <- c(fmt3(aic_null), fmt3(bic_null), fmt3(null_dev))
    mf_full <- c(fmt3(aic_val), fmt3(bic_val), fmt3(resid_dev))
    
    mf_cols <- list(mf_crit, mf_int, mf_full)
    mf_w <- integer(length(mf_cols))
    for (i in seq_along(mf_cols)) mf_w[i] <- max(nchar(mf_hdr[i]), max(nchar(mf_cols[[i]]), na.rm = TRUE))
    mf_fmt <- build_fmt(mf_w)
    cat(do.call(sprintf, c(mf_fmt, as.list(mf_hdr))))
    for (i in seq_along(mf_crit)) {
      row <- list(mf_crit[i], mf_int[i], mf_full[i])
      cat(do.call(sprintf, c(mf_fmt, row)))
    }
    cat("\n")
    
    # Global Test (Likelihood Ratio Test)
    cat("Global Test for All Slopes: BETA = 0\n")
    lr_stat <- null_dev - resid_dev
    lr_df <- df_reg
    lr_p <- pchisq(lr_stat, lr_df, lower.tail = FALSE)
    
    gt_hdr <- c("Test", "Chi-Square", "df", "p-value")
    gt_test <- "Likelihood Ratio"
    gt_chi <- fmt3(lr_stat)
    gt_df <- fmt_int(lr_df)
    gt_p <- fmt_p(lr_p)
    
    gt_cols <- list(gt_test, gt_chi, gt_df, gt_p)
    gt_w <- integer(length(gt_cols))
    for (i in seq_along(gt_cols)) gt_w[i] <- max(nchar(gt_hdr[i]), nchar(gt_cols[[i]]))
    gt_fmt <- build_fmt(gt_w)
    cat(do.call(sprintf, c(gt_fmt, as.list(gt_hdr))))
    cat(do.call(sprintf, c(gt_fmt, gt_cols)))
    cat("\n")
  }

  # Coefficients
  cat("Coefficients Table\n")
  
  if (is_lm) {
    c_hdr <- c("Variables","Coefficients","Standard Error","t-test","p-value",
               paste0(CL, "% CI LB"), paste0(CL, "% CI UB"))
    if (isTRUE(std_beta)) c_hdr <- c(c_hdr, "Std. Beta")
    c_hdr <- c(c_hdr, "VIF")

    c_vars <- rownames(ctab)
    c_coef <- fmt3(ctab$Coefficients)
    c_se   <- fmt3(ctab$`Standard Error`)
    c_t    <- fmt3(ctab$`t-test`)
    c_p    <- fmt_p(ctab$`p-value`)
    c_lb   <- fmt3(ctab[[lb_name]])
    c_ub   <- fmt3(ctab[[ub_name]])
    c_beta <- if (isTRUE(std_beta)) fmt3(ctab$`Std. Beta`) else NULL
    c_vif  <- ifelse(is.na(ctab$VIF), "", fmt3(ctab$VIF))

    c_cols <- list(c_vars, c_coef, c_se, c_t, c_p, c_lb, c_ub)
    if (isTRUE(std_beta)) c_cols <- c(c_cols, list(c_beta))
    c_cols <- c(c_cols, list(c_vif))

    c_w <- integer(length(c_cols))
    for (i in seq_along(c_cols)) c_w[i] <- max(nchar(c_hdr[i]), max(nchar(c_cols[[i]]), na.rm = TRUE))
    c_fmt <- build_fmt(c_w)
    cat(do.call(sprintf, c(c_fmt, as.list(c_hdr))))
    for (i in seq_len(nrow(ctab))) {
      row <- list(c_vars[i], c_coef[i], c_se[i], c_t[i], c_p[i], c_lb[i], c_ub[i])
      if (isTRUE(std_beta)) row <- c(row, c_beta[i])
      row <- c(row, c_vif[i])
      cat(do.call(sprintf, c(c_fmt, row)))
    }
  } else if (is_glm) {
    # GLM Coefficients Table with Odds Ratios
    c_hdr <- c("Parameter", "Estimate", "Standard Error", "Wald Chi-Square", "df", "p-value",
               "Odds Ratio", paste0(CL, "% CI LB OR"), paste0(CL, "% CI UB OR"), "VIF")
    
    c_vars <- rownames(ctab)
    c_coef <- fmt3(ctab$Coefficients)
    c_se   <- fmt3(ctab$`Standard Error`)
    c_wald <- fmt3((ctab$Coefficients / ctab$`Standard Error`)^2)
    c_df   <- rep("1", nrow(ctab))
    c_p    <- fmt_p(ctab$`p-value`)
    
    # Odds Ratios and CIs
    odds_ratio <- exp(ctab$Coefficients)
    or_lb <- exp(ctab[[lb_name]])
    or_ub <- exp(ctab[[ub_name]])
    
    c_or <- fmt3(odds_ratio)
    c_or_lb <- fmt3(or_lb)
    c_or_ub <- fmt3(or_ub)
    c_vif <- ifelse(is.na(ctab$VIF), "", fmt3(ctab$VIF))
    
    # Blank out OR values for intercept
    c_or[1] <- ""
    c_or_lb[1] <- ""
    c_or_ub[1] <- ""

    c_cols <- list(c_vars, c_coef, c_se, c_wald, c_df, c_p, c_or, c_or_lb, c_or_ub, c_vif)

    c_w <- integer(length(c_cols))
    for (i in seq_along(c_cols)) c_w[i] <- max(nchar(c_hdr[i]), max(nchar(c_cols[[i]]), na.rm = TRUE))
    c_fmt <- build_fmt(c_w)
    cat(do.call(sprintf, c(c_fmt, as.list(c_hdr))))
    for (i in seq_len(nrow(ctab))) {
      row <- list(c_vars[i], c_coef[i], c_se[i], c_wald[i], c_df[i], c_p[i], 
                  c_or[i], c_or_lb[i], c_or_ub[i], c_vif[i])
      cat(do.call(sprintf, c(c_fmt, row)))
    }
  }

  # Predictions Table (if diag = TRUE and lm model)
  if (isTRUE(diag) && is_lm) {
    cat("\n")
    cat("Predictions Table\n")
    
    # Get predictions and residuals
    yhat <- fitted(model)
    resids <- resid(model)
    std_resids <- resid(model) / summary(model)$sigma
    
    # Calculate leverage (hat values)
    X <- model.matrix(model)
    H <- X %*% solve(t(X) %*% X) %*% t(X)
    lev <- diag(H)
    
    # Calculate Cook's Distance
    cooks <- cooks.distance(model)
    
    # Calculate DFFITS
    dffits_vals <- dffits(model)
    
    # Get predictor values
    pred_data <- model.frame(model)[, -1, drop = FALSE]
    
    # Build the predictions table
    p_hdr <- c("Observation Number", colnames(pred_data), y_name, 
               "Predicted Values", "Residuals", "Standardized Residuals",
               "Leverage", "Cook's D", "DFFITS")
    
    # Format values
    obs_num <- seq_len(n_used)
    p_obs <- fmt_int(obs_num)
    
    # Format predictor columns
    pred_cols_fmt <- lapply(pred_data, function(col) fmt3(col))
    
    p_y <- fmt3(y)
    p_yhat <- fmt3(yhat)
    p_resids <- fmt3(resids)
    p_std_resids <- fmt3(std_resids)
    p_lev <- fmt3(lev)
    p_cooks <- fmt3(cooks)
    p_dffits <- fmt3(dffits_vals)
    
    # Combine all columns
    p_cols <- c(list(p_obs), pred_cols_fmt, list(p_y, p_yhat, p_resids, 
                p_std_resids, p_lev, p_cooks, p_dffits))
    
    # Calculate widths
    p_w <- integer(length(p_cols))
    for (i in seq_along(p_cols)) {
      p_w[i] <- max(nchar(p_hdr[i]), max(nchar(p_cols[[i]]), na.rm = TRUE))
    }
    
    # Build and print header
    p_fmt <- build_fmt(p_w)
    cat(do.call(sprintf, c(p_fmt, as.list(p_hdr))))
    
    # Print each row
    for (i in seq_len(n_used)) {
      row <- lapply(p_cols, function(col) col[i])
      cat(do.call(sprintf, c(p_fmt, row)))
    }
  }
  
  # Predictions Table (if diag = TRUE and glm model)
  if (isTRUE(diag) && is_glm) {
    cat("\n")
    cat("Predictions Table\n")
    
    # Get predictor values
    pred_data <- model.frame(model)[, -1, drop = FALSE]
    
    # Get observed values, predicted probabilities, and classifications
    obs_y <- model$y
    pred_prob <- fitted(model)
    pred_class <- ifelse(pred_prob >= 0.5, 1, 0)
    
    # Confidence intervals for predicted probabilities
    pred_logit <- predict(model, type = "link", se.fit = TRUE)
    se_logit <- pred_logit$se.fit
    ci_lower_logit <- pred_logit$fit - qnorm(1 - (1 - level)/2) * se_logit
    ci_upper_logit <- pred_logit$fit + qnorm(1 - (1 - level)/2) * se_logit
    ci_lower_prob <- exp(ci_lower_logit) / (1 + exp(ci_lower_logit))
    ci_upper_prob <- exp(ci_upper_logit) / (1 + exp(ci_upper_logit))
    
    # Pearson and Deviance residuals
    pearson_resid <- residuals(model, type = "pearson")
    deviance_resid <- residuals(model, type = "deviance")
    
    # Standardized residuals
    std_pearson <- rstandard(model, type = "pearson")
    std_deviance <- rstandard(model, type = "deviance")
    
    # Calculate leverage (hat values)
    X <- model.matrix(model)
    W <- diag(pred_prob * (1 - pred_prob))
    H <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    lev <- diag(H)
    
    # Calculate Cook's Distance and C-bar (scaled Cook's D)
    cooks <- cooks.distance(model)
    p <- ncol(X)  # number of parameters
    c_bar <- cooks * p
    
    # Calculate DFFITS (Delta Beta)
    dffits_vals <- dffits(model)
    dfbetas_vals <- dfbetas(model)
    delta_deviance <- dfbetas_vals[, 1]  # Just use first column as representative
    
    # Build the predictions table header
    p_hdr <- c("Observation Number", colnames(pred_data), y_name,
               "Predicted Classifications", "Predicted Probabilities",
               paste0(CL, "% CI Mean LB"), paste0(CL, "% CI Mean UB"),
               "Pearson Residuals", "Deviance Residuals",
               "Studentized Pearson Residuals", "Studentized Deviance Residuals",
               "Leverage", "C", "C Bar", "Delta Pearson", "Delta Deviance")
    
    # Format values
    obs_num <- seq_len(n_used)
    p_obs <- fmt_int(obs_num)
    
    # Format predictor columns
    pred_cols_fmt <- lapply(pred_data, function(col) fmt3(col))
    
    p_y <- fmt_int(obs_y)
    p_class <- fmt_int(pred_class)
    p_prob <- fmt3(pred_prob)
    p_ci_lb <- fmt3(ci_lower_prob)
    p_ci_ub <- fmt3(ci_upper_prob)
    p_pearson <- fmt3(pearson_resid)
    p_deviance <- fmt3(deviance_resid)
    p_std_pearson <- fmt3(std_pearson)
    p_std_deviance <- fmt3(std_deviance)
    p_lev <- fmt3(lev)
    p_cooks <- fmt3(cooks)
    p_c_bar <- fmt3(c_bar)
    p_delta_pearson <- fmt3(dffits_vals)
    p_delta_deviance <- fmt3(delta_deviance)
    
    # Combine all columns
    p_cols <- c(list(p_obs), pred_cols_fmt, list(p_y, p_class, p_prob, 
                p_ci_lb, p_ci_ub, p_pearson, p_deviance,
                p_std_pearson, p_std_deviance, p_lev, p_cooks, p_c_bar,
                p_delta_pearson, p_delta_deviance))
    
    # Calculate widths
    p_w <- integer(length(p_cols))
    for (i in seq_along(p_cols)) {
      p_w[i] <- max(nchar(p_hdr[i]), max(nchar(p_cols[[i]]), na.rm = TRUE))
    }
    
    # Build and print header
    p_fmt <- build_fmt(p_w)
    cat(do.call(sprintf, c(p_fmt, as.list(p_hdr))))
    
    # Print each row
    for (i in seq_len(n_used)) {
      row <- lapply(p_cols, function(col) col[i])
      cat(do.call(sprintf, c(p_fmt, row)))
    }
  }

  if (is_lm) {
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
      CL = CL,
      digits = digits,
      std_beta = std_beta
    ))
  } else if (is_glm) {
    invisible(list(
      coefficients = ctab,
      deviance = list(
        null_deviance = model$null.deviance,
        residual_deviance = model$deviance,
        df_null = sm$df.null,
        df_residual = sm$df.residual
      ),
      fit_statistics = list(
        AIC = AIC(model),
        BIC = BIC(model)
      ),
      summary = list(
        N_Read = n_read, N_Missing = n_miss, N_Used = n_used
      ),
      robust = robust,
      CL = CL,
      digits = digits
    ))
  }
}
