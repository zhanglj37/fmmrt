# Example Code

```
devtools::install_github('zhanglj37/fmmrt')

# Fit the proposed model
fit <- fmm_rt(responses, rt)

# Simulated data
dat <- read.table(system.file("extdata", "simdata.txt", package = "fmmrt"), header = TRUE)
responses = dat[, grep("y",colnames(dat))]
rt = dat[, grep("rt",colnames(dat))]

# Results Output
fit$epsr
fit$probability

```
