# install packages
install.packages("tidyverse")
install.packages("magrittr")
install.packages("forcats")
install.packages("gtools")
install.packages("arrangements")

# import packages
library(tidyverse)
library(magrittr)
library(forcats)
library(gtools)
library(arrangements)

# data format
ss <- data.frame(strain = NA,
                 rep = NA,
                 days_after_bleach = NA,
                 total = NA,
                 alive = NA,
                 proportion = NA)

ss$strain %<>% as.factor()
ss$strain %<>% fct_relevel(levels = c("strain1", "strain2", "strain3")) # relevel in the order you like

# create a dataframe called "result" where you'll save half life values
result <- ss %>%
  select(c("strain", "rep")) %>%
  unique() %>%
  as.data.frame()

# fit data and extract half lives
for (i in 1:nrow(result)) {
  model <- glm(data = filter(ss, strain == result$strain[i] & rep == result$rep[i]),
               formula = proportion ~ days_after_bleach,
               # weights = total,
               family = "quasibinomial")
  
  A = model$coefficients[2]
  B = model$coefficients[1]
  
  fifty_percent_life <- -B/A # time when proportion alive is 50%
  half_life <- -log(2+exp(B), base = exp(A)) # time when proportion alive is half of starting value
  # p < 0.05 for model_sig means there is at least one significant predictor in the model
  model_sig <- pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual, lower.tail = FALSE)
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  
  result$fifty_percent_life[i] <- fifty_percent_life
  result$half_life[i] <- half_life
  result$goodness_of_fit[i] <- goodness_of_fit
  result$model_sig[i] <- model_sig
}
result

# use Bartlett's test to see whether or not can pool variances across strains
bartlett.test(half_life ~ strain, data = result)

# create a dataframe called "stats" where you'll save the stats
stats <- result
stats$strain %<>% factor()
all_strains <- levels(stats$strain)
all_pairs <- combinations(all_strains, 2)

t_test_stats <- data.frame(strain1 = NA,
                           strain2 = NA,
                           t.test.p.unadjusted = NA) 

for (i in 1: dim(all_pairs)[1]) {
  tmp <- t.test(x = filter(stats, strain == all_pairs[i,1])$half_life,
                y = filter(stats, strain == all_pairs[i,2])$half_life,
                var.equal = TRUE) # set to FALSE if Bartlett's p < 0.05
  t_test_stats[i,1] <- all_pairs[i,1]
  t_test_stats[i,2] <- all_pairs[i,2]
  t_test_stats[i,3] <- tmp$p.value
}
t_test_stats %<>% arrange(t.test.p.unadjusted)
t_test_stats$p.signif <- stars.pval(t_test_stats$t.test.p.unadjusted)


