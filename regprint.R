regprint <- function(model, conf_level = 95) {
  # friendly error if not an lm
  if (!inherits(model, "lm")) {
    obj <- deparse(substitute(model))
    stop(sprintf(
      "Error: %s is not a lm object. regprint() requires you fit a linear model first using the lm() function and save it as an object.",
      obj
    ), call. = FALSE)
  }
  # validate confidence level (1..99)
  if (!is.numeric(conf_level) || length(conf_level) != 1 || !is.finite(conf_level) ||
      conf_level < 1 || conf_level > 99) {
    stop("Error: 'conf_level' must be a single number between 1 and 99 (e.g., 90, 95).", call. = FALSE)
  }
  level <- conf_level / 100
  lb_name <- sprintf("%d%% CI LB", conf_level)
  ub_name <- sprintf("%d%% CI UB", conf_level)

  # --- rest of the function unchanged ---
  y_name <- as.character(formula(model))[2]
  sm <- summary(model)
  fmt_num <- function(x, k = 4) sprintf(paste0("%.", k, "f"), x)
  fmt_p   <- function(p) ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p))

  r2 <- sm$r.squared; adjr <- sm$adj.r.squared
  multR <- sqrt(r2); s_err <- sm$sigma
  n_used <- length(sm$residuals)
  n_miss <- if (is.null(model$na.action)) 0L else length(model$na.action)
  n_read <- n_used + n_miss

  y <- model.response(model.frame(model))
  e <- resid(model); yhat <- fitted(model) # yhat not used, but kept for clarity
  sse <- sum(e^2); sst <- sum((y - mean(y))^2); ssr <- sst - sse

  df_reg <- unname(sm$fstatistic["numdf"])
  df_res <- unname(sm$fstatistic["dendf"])
  df_tot <- df_reg + df_res
  ms_reg <- ssr / df_reg; ms_res <- sse / df_res
  fstat  <- unname(sm$fstatistic["value"])
  p_f    <- pf(fstat, df_reg, df_res, lower.tail = FALSE)

  ctab <- as.data.frame(sm$coefficients)
  names(ctab) <- c("Coefficients", "Standard Error", "t-test", "p-value")

  ci <- confint(model, level = level)
  ctab[[lb_name]] <- ci[, 1]
  ctab[[ub_name]] <- ci[, 2]

  # ----- VIFs: EXCLUDE INTERCEPT -----
  X <- model.matrix(model)
  pred_cols <- setdiff(colnames(X), "Intercept")   # only predictors
  Xp <- if (length(pred_cols)) X[, pred_cols, drop = FALSE] else NULL

  vif <- rep(NA_real_, nrow(ctab))                   # default NA (incl. intercept)
  if (!is.null(Xp)) {
    if (ncol(Xp) == 1L) {
      # With a single predictor, VIF is 1
      vif[match(colnames(Xp), rownames(ctab))] <- 1
    } else if (ncol(Xp) > 1L) {
      v <- numeric(ncol(Xp))
      for (j in seq_len(ncol(Xp))) {
        r2j <- summary(lm(Xp[, j] ~ Xp[, -j]))$r.squared
        v[j] <- 1 / (1 - r2j)
      }
      vif[match(colnames(Xp), rownames(ctab))] <- v
    }
  }
  ctab$VIF <- vif
  # -----------------------------------

  obj_name <- deparse(substitute(model))
  cat(sprintf("Linear Regression Model: %s\n", obj_name))
  cat(sprintf("Dependent Variable: %s\n\n", y_name))

  cat("Regression Statistics Table\n")
  cat(sprintf("%-22s %8s\n", "Multiple R",        fmt_num(multR, 4)))
  cat(sprintf("%-22s %8s\n", "R Square",          fmt_num(r2, 4)))
  cat(sprintf("%-22s %8s\n", "Adjusted R Square", fmt_num(adjr, 4)))
  cat(sprintf("%-22s %8s\n", "Standard Error",    sprintf("%.3f", s_err)))
  cat(sprintf("%-22s %8d\n", "N Obs Read",        n_read))
  cat(sprintf("%-22s %8d\n", "N Obs Missing",     n_miss))
  cat(sprintf("%-22s %8d\n\n", "N Obs Used",      n_used))

  cat("ANOVA Table\n")
  cat(sprintf("%-12s %6s %10s %10s %12s %10s\n",
              "Source","df","SS","MS","F-Statistic","p-value"))
  cat(sprintf("%-12s %6d %10.3f %10.3f %12.3f %10s\n",
              "Regression", df_reg, ssr, ms_reg, fstat, fmt_p(p_f)))
  cat(sprintf("%-12s %6d %10.3f %10.3f %12s %10s\n",
              "Residual",   df_res, sse, ms_res, "", ""))
  cat(sprintf("%-12s %6d %10.3f %10s %12s %10s\n\n",
              "Total",      df_tot, sst, "", "", ""))

  cat("Coefficients Table\n")
  cat(sprintf("%-12s %12s %15s %10s %10s %12s %12s %8s\n",
              "Variables","Coefficients","Standard Error",
              "t-test","p-value", paste0(conf_level, "% CI LB"),
              paste0(conf_level, "% CI UB"), "VIF"))
  for (i in seq_len(nrow(ctab))) {
    cat(sprintf("%-12s %12s %15s %10s %10s %12s %12s %8s\n",
                rownames(ctab)[i],
                sprintf("%.3f", ctab$Coefficients[i]),
                sprintf("%.3f", ctab$`Standard Error`[i]),
                sprintf("%.3f", ctab$`t-test`[i]),
                fmt_p(ctab$`p-value`[i]),
                sprintf("%.3f", ctab[[lb_name]][i]),
                sprintf("%.3f", ctab[[ub_name]][i]),
                ifelse(is.na(ctab$VIF[i]), "", sprintf("%.3f", ctab$VIF[i]))))
  }
  invisible(NULL)
}
