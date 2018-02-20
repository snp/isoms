plot_peakfit <- function(spectrum, model) {
    spec_df <- data_frame(mz = spectrum[, 1], i = spectrum[, 2])
    fit_mz = seq(min(spec_df$mz), max(spec_df$mz), length.out = 400)
    fit_df <- data_frame(mz = fit_mz, i = predict(model, newdata = list(x = fit_mz)))
    ggplot(spec_df, aes(x = mz, y = i)) + geom_point() + geom_line(data = fit_df) + theme_minimal()
}
