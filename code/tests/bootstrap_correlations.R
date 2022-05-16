library(tidyverse)
library(MASS)

# Set seed for reproducibility
set.seed(1)

# Set target correlation and generate sample with given size
r <- 0.5
dat <- mvrnorm(n=154, mu=c(0,0), Sigma = matrix(c(1, r, r, 1), nrow=2), empirical = TRUE) %>%
  as_tibble

# Bootstrap function
boot_cor <- function(x) {
  cor(sample_frac(x, size=1, replace=T))[1,2]
}

# Test bootstrap function
boot_cor(dat)
bootr <- replicate(n=500, boot_cor(dat)) 
hist(bootr, 100)
median(bootr) # Bootstrap stabilized value

# Let's look at how much this bootstrap stabilized value varies
r <- 0.5
bootn <- 500 # Number of bootstrap samples
distn <- 500 # Number of times to compute bootstrap stabilized value

dat <- mvrnorm(n=154, mu=c(0,0), Sigma = matrix(c(1, r, r, 1), nrow=2), empirical = TRUE) %>%
  as_tibble
prog <- txtProgressBar(max=distn, style=3)
bootdist <- numeric(distn)
for (n in 1:distn) {
  rundist <- replicate(n=bootn, boot_cor(dat)) 
  bootdist[[n]] <- median(rundist)
  setTxtProgressBar(prog, n)
}
close(prog)

# Examine variation in bootstrap stablized value
hist(bootdist, 100)
quantile(bootdist, c(0.025, 0.975))
