quantifyIsoRatio <- function(data, control = "Control", file = "20171226_02_NP02_fit.csv", peak = "13C", ion = "P") {
    peak_ <- peak
    ion_ <- ion
    file_ <- file
    message(file, ",", peak, ",", ion)
    s_control <- data %>% ungroup() %>% filter(group == control & peak == peak_, ion == ion_) %>% filter(ltic > 1)
    s_sample <- data %>% ungroup() %>% filter(file == file_ & peak == peak_, ion == ion_) %>% filter(ltic > 1)
    
    if (nrow(s_control) < 10 | nrow(s_sample) < 10) 
        return(data.frame())
    mm_control <- MASS::rlm(gamma ~ ldI, data = s_control)
    residuals_control <- mm_control$residuals/sd(mm_control$residuals)
    s_control$color <- abs(residuals_control) < 3
    s_control <- s_control %>% filter(color) %>% select(-color)
    
    mm_sample <- MASS::rlm(gamma ~ ldI, data = s_sample)
    residuals_sample <- mm_sample$residuals/sd(mm_sample$residuals)
    s_sample$color <- abs(residuals_sample) < 3
    s_sample <- s_sample %>% filter(color) %>% select(-color)
    
    loess_control <- loess(gamma ~ ldI, data = s_control, weights = ltic)
    s_control$loess <- loess_control$fitted
    loess_sample <- loess(gamma ~ ldI, data = s_sample, weights = ltic)
    s_sample$loess <- loess_sample$fitted
    
    # s_control %>% ggplot(aes(x=ldI)) + geom_point(aes(y=gamma), color='black') + geom_line(aes(y=loess), color='red') + labs(title='Control') + theme_minimal()
    # s_sample %>% ggplot(aes(x=ldI)) + geom_point(aes(y=gamma), color='black') + geom_line(aes(y=loess), color='red') + labs(title='Control') + theme_minimal()
    
    xmin <- max(min(s_control$ldI), min(s_sample$ldI))
    xmax <- min(max(s_control$ldI), max(s_sample$ldI))
    # xrange <- seq(xmin, xmax, length=3000)
    
    xrange <- s_sample$ldI[s_sample$ldI > xmin & s_sample$ldI < xmax]
    control_predicted <- predict(loess_control, data.frame(ldI = xrange), se = T)
    sample_predicted <- predict(loess_sample, data.frame(ldI = xrange), se = T)
    
    res <- data.frame(control = control_predicted$fit, control_se = control_predicted$se.fit, sample = sample_predicted$fit, sample_se = sample_predicted$se.fit, 
        ldI = xrange, diff = (sample_predicted$fit - control_predicted$fit)/control_predicted$fit * 1000, diff_se = sqrt((control_predicted$se.fit^2 + sample_predicted$se.fit^2)/2)/control_predicted$fit * 
            1000)
    # res %>% ggplot(aes(x=control, y=sample, color=ltic)) + geom_point() + theme_minimal()
    
    # res %>% ggplot(aes(x=ldI, y=diff)) + geom_point() + theme_minimal() + geom_line(aes(y=diff-diff_se)) + geom_line(aes(y=diff+diff_se)) res %>%
    # ggplot(aes(x=diff)) + geom_density() + theme_minimal()
    
    return(res)
}

mostprobable <- function(x) {
    xx <- na.exclude(x)
    if (length(xx) > 10) {
        dens <- density(xx)
        return(dens$x[which.max(dens$y)])
    } else {
        return(mean(xx))
    }
}
