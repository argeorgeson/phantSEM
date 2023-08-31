## code to prepare `Lookup_Table` dataset goes here

usethis::use_data(Lookup_Table, overwrite = TRUE)

# create lookup table for phantom baseline longitudinal model



observedall <- expand.grid(corxm = seq(-.9, .9, .1), corxy = seq(-.9, .9, .1), cormy = seq(-.9, .9, .1))


observedall$cond <- seq(1, nrow(observedall), 1)


test <- data.frame(NA)
for (i in 1:nrow(observedall)) {
  df <- matrix(c(1, observedall[i, 1], observedall[i, 2], observedall[i, 1], 1, observedall[i, 3], observedall[i, 2], observedall[i, 3], 1), nrow = 3)
  test[i, ] <- corpcor::is.positive.definite(df)
}

observedall <- cbind(observedall, test)
observed <- observedall[which(observedall$NA. == TRUE), ]

results <- matrix(nrow = 0, ncol = 1)
for (i in 1:nrow(observed)) {
  df <- matrix(c(1, observed[i, 1], observed[i, 2], observed[i, 1], 1, observed[i, 3], observed[i, 2], observed[i, 3], 1), nrow = 3)
  temp <- cbind(cond = observed[i, 4], SA_tempbias_cov(CorStability = c(-.9, .9), CorCL = c(-.9, .9), corM1Y1 = c(-.9, .9), increments = c(.1, .1, .3), covmat = df, n = 300, extreme = TRUE))
  results <- rbind(results, temp)
}

# remove NA
results_noNA <- results %>% filter(jointsigab == "Yes" | jointsigab == "No")
setwd("G:/My Drive/ASU/Dave lab/temporal bias/Manuscript/Revision/Lookup Table")
write.csv(results, "./lookup_table_all.csv")
write_rds(results, "./lookup_table_all.rds")



write.csv(results_noNA, "./lookup_table_noNA.csv")
write_rds(results_noNA, "./lookup_table_noNA.rds")

# remove empty columns
obs2 <- observed[, 1:4]
write_rds(obs2, "./conditions.rds")


results_noNA <- readRDS("G:/My Drive/ASU/Dave lab/temporal bias/Manuscript/Revision/Lookup Table/lookup_table_noNA_saverds.rds")


#saveRDS(results_smaller,"G:/My Drive/ASU/Dave lab/temporal bias/Manuscript/Revision/Lookup Table/lookup_table_noNA_small.rds",compress="xz")
#saveRDS(results_noNA,"G:/My Drive/ASU/Dave lab/temporal bias/Manuscript/Revision/Lookup Table/lookup_table_noNA_saverds.rds",compress="xz")


results_smaller <- results_noNA[,c(1:5,9,14)]

results_smaller$cond <- as.integer(results_smaller$cond)


#results_smaller <- readRDS("G:/My Drive/ASU/Dave lab/temporal bias/Manuscript/Revision/Lookup Table/lookup_table_noNA_small.rds")

# this puts it in the package
usethis::use_data(results_smaller, obs2, internal = TRUE, overwrite=TRUE, compress = "xz")



# creates data-raw
# usethis::use_data_raw("Lookup_Table")
