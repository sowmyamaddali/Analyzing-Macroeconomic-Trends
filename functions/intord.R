

intord <- function(y, freq = 12, year = 1990, period = 1) {
  pacman::p_load(ggplot2, forecast, gridExtra, ggExtra, tidyverse, 
          tsibble)
  #
  # purpose: determine order of integration for variable y.
  # 
  # usage: out = intord(y)
  # This creates time plots, calculates SDs
  #  and performs an ADF unit root test for
  # for level, 1st & 2nd difference of y
  #
  # out contains ADF t-statistics and Mackinnon (2010) critical values.
  #
  # Critical values from:
  # MacKinnon, J.G. (2010) Critical Values for Cointegration Tests.
  # Queen?s Economics Department Working Paper No. 1227, Table 2.
  #
  # Jeff Mills, 2011
  #
  suppressWarnings({
  if (is.ts(y)) {
    data <- as_tsibble(y)
  } else {
    data <- ts(y, freq = freq, start = c(year,period)) %>%
      as_tsibble()
  }
  
  
  data <- data %>%
    mutate(first_diff = difference(value),
           second_diff = difference(value, differences = 2))


  # Plot levels
  p1 <- ggplot(data = data, aes(x = index, y = value)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = mean(data$value), color = "red") +
    ggtitle(paste("SD=", round(sd(data$value), 4))) +
    ylab("Outcome") +
    xlab("Time") + theme_bw() + 
    scale_y_continuous(labels = scales::comma)
  
  # Detrended series
  p2 <- ggplot(data = data, aes(x = index, y = first_diff)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = mean(data$first_diff, na.rm = T), color = "red") +
    ggtitle(paste("SD=", round(sd(data$first_diff, na.rm = T), 4))) +
    ylab(expression(paste(Delta * y))) +
    xlab("Time") + theme_bw() + 
    scale_y_continuous(labels = scales::comma)
  
  # First difference
  p3 <- ggplot(data = data, aes(x = index, y = second_diff)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = mean(data$second_diff, na.rm = T), color = "red") +
    ggtitle(paste("SD=", round(sd(data$second_diff, na.rm = T), 4))) +
    ylab(expression(paste(Delta^2 * y))) +
    xlab("Time")+ theme_bw() + 
    scale_y_continuous(labels = scales::comma)
  
  # ACF of levels
  p4 <- forecast::ggAcf(data$value) + 
    ggtitle("ACF for Level") + theme_bw()
  
  # ACF of first difference
  p5 <- forecast::ggAcf(data$first_diff) +
    ggtitle("ACF for 1st Difference") + theme_bw() 
  
  # ADF test first round
  
  # pmax is selected max lag length, maxp must be at least 2
  pmax <- 12
  maxp <- pmax + 1
  

  y <- data$value
  n <- length(y)
  dy <- diff(y)
  
  z <- embed(dy, maxp)
  
  zz <- embed(y, maxp)
  y1 <- zz[,2]
  xx <- cbind(z[,1], y[maxp:(n-1)], z[,2:maxp])
  nobs <- nrow(xx)
  # DF test (0 lags)
  c <- rep(1, nrow(xx))
  xvars <- cbind(c, y[maxp:(n-1)])
  yvar <- xx[,1]
  ixx <- solve(t(xvars) %*% xvars)
  bh <- ixx %*% t(xvars) %*% yvar
  yh <- xvars %*% bh
  res <- yvar - yh
  rss <- t(res) %*% res
  k <- ncol(xvars)
  s2 <- as.numeric(rss / (nobs - k))
  covb <- s2 * ixx
  seb <- sqrt(diag(covb))
  
  bic <- rep(0, maxp)
  adft <- bic
  adft[1] <- bh[2] / seb[2]
  bic[1] <- log(rss / nobs) + log(nobs) * (k + 1) / nobs
  
  for (i in 3:(maxp + 1)) {
    xvars <- cbind(c, xx[,2:i])
    ixx <- solve(t(xvars) %*% xvars)
    bh <- ixx %*% t(xvars) %*% yvar
    yh <- xvars %*% bh
    res <- yvar - yh
    rss <- t(res) %*% res
    k <- ncol(xvars)
    s2 <- as.numeric(rss / (nobs - k))
    covb <- s2 * ixx
    seb <- sqrt(diag(covb))
    adft[i - 1] <- bh[2] / seb[2]
    bic[i - 1] <- log(rss / nobs) + log(nobs) * (k + 1) / nobs
  }
  
  ind <- which.min(bic)
  # cat("ADF t-value","lags")
  round1 <- c(round(adft[ind], 2))
  
  # ADF test second round
  
  y <- dy
  
  n <- length(y)
  dy <- diff(y)
  
  z <- embed(dy, maxp)
  
  zz <- embed(y, maxp)
  y1 <- zz[,2]
  xx <- cbind(z[,1], y[maxp:(n-1)], z[,2:maxp])
  nobs <- nrow(xx)
  # DF test (0 lags)
  c <- rep(1, nrow(xx))
  xvars <- cbind(c, y[maxp:(n-1)])
  yvar <- xx[,1]
  ixx <- solve(t(xvars) %*% xvars)
  bh <- ixx %*% t(xvars) %*% yvar
  yh <- xvars %*% bh
  res <- yvar - yh
  rss <- t(res) %*% res
  k <- ncol(xvars)
  s2 <- as.numeric(rss / (nobs - k))
  covb <- s2 * ixx
  seb <- sqrt(diag(covb))
  
  bic <- rep(0, maxp)
  adft <- bic
  adft[1] <- bh[2] / seb[2]
  bic[1] <- log(rss / nobs) + log(nobs) * (k + 1) / nobs
  
  for (i in 3:(maxp + 1)) {
    xvars <- cbind(c, xx[,2:i])
    ixx <- solve(t(xvars) %*% xvars)
    bh <- ixx %*% t(xvars) %*% yvar
    yh <- xvars %*% bh
    res <- yvar - yh
    rss <- t(res) %*% res
    k <- ncol(xvars)
    s2 <- as.numeric(rss / (nobs - k))
    covb <- s2 * ixx
    seb <- sqrt(diag(covb))
    adft[i - 1] <- bh[2] / seb[2]
    bic[i - 1] <- log(rss / nobs) + log(nobs) * (k + 1) / nobs
  }
  bic
  ind <- which.min(bic)
  round2 <- c(round(adft[ind], 2))
  rbind(round1, round2)
  
  adf_statistics <- cbind(round1, round2)
  
  # MacKinnon critical values
  c1 <- -3.43035 - 6.5393 / nobs - 16.786 / nobs^2 - 79.433 / nobs^3
  c5 <- -2.86154 - 2.8903 / nobs - 4.234 / nobs^2 - 40.04 / nobs^3
  c10 <- -2.56677 - 1.5384 / nobs - 2.809 / nobs^2
  
  # cat("10%, 5% and 1% critical values")
  test_data <- data.frame(crit_value = round(c(c10, c5, c1), 2)) %>%
    mutate(level = factor(c('10%',
                                   '5%',
                                   '1%'),
                             levels = c("10%", "5%", "1%"))) %>%
    relocate(level, crit_value) %>%
    mutate(ADF_round1_stat = round1,
           ADF_round2_stat = round2) %>%
    mutate(round_1_dec = ifelse(ADF_round1_stat <= crit_value, 'REJ','FTR'),
           round_2_dec = ifelse(ADF_round2_stat <= crit_value, 'REJ','FTR'))
  
  # Window with the test statistic results
  p6 <- ggplot(test_data) + theme_bw() + 
    geom_point(aes(x = level, y = crit_value), size = 4, shape = 'square') + 
    geom_point(aes(x = level, y = ADF_round1_stat, color = round_1_dec),
               size = 4, shape = 'triangle') + 
    geom_point(aes(x = level, y = ADF_round2_stat, color = round_2_dec),
               size = 4, shape = 'circle') +
    scale_color_manual(name = 'Decision',values = c('blue','red')) + 
    theme(legend.position = 'none') + ylab('Statistic') + 
    xlab('Level of Signficance')
    
    
  
  
  # Combine plots into a grid
  grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
  
  
  return(list(test_results = test_data))
  })
}

