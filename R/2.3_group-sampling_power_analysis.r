popsize <- 1000
groupsize <- 10

dyad_replication_calculator <- function(proposed_samplesizes,
                                        # target_samegrp_dyads,
                                        groupsize,
                                        popsize,
                                        n_reps){
  
  # figure out number of groups
  # for every animal, pick a group at random
  # check to be sure no group has more than groupsize animals, if yes, repick for somebody
  # table number of animals in each group
  # replace table values by tabled value choose 2 (this accounts for multiple dyads within the same group)
  # sum those numbers up over all groups
  
  n_groups <- popsize / groupsize
  dyads_within_grps <- rep(NA, n_reps)
  for(i in 1:n_reps){
    group_assignments <- sample(x = 1:n_groups, size = proposed_samplesizes, replace = T)
    counts_within_grps <- table(group_assignments)
    dyads_within_grps[i] <- sum(choose(counts_within_grps, 2))
  }
  
  return(dyads_within_grps)
}

dyad_replication_calculator <- function(proposed_samplesizes,
                                        # target_samegrp_dyads,
                                        groupsize,
                                        popsize,
                                        n_reps){
  
  # figure out number of groups
  # for every animal, pick a group at random
  # check to be sure no group has more than groupsize animals, if yes, repick for somebody
  # table number of animals in each group
  # replace table values by tabled value choose 2 (this accounts for multiple dyads within the same group)
  # sum those numbers up over all groups
  
  simmed_group_sizes <- rpois(n = 100, lambda = groupsize)
  group_number <- seq(1:100)
  group_labels <- rep(group_number, each = simmed_group_sizes)[1:popsize]
  #  cumsum(simmed_group_sizes)
  #  n_groups <- popsize / groupsize
  dyads_within_grps <- rep(NA, n_reps)
  for(i in 1:n_reps){
    group_assignments <- sample(x = group_labels, size = proposed_samplesizes, replace = F)
    counts_within_grps <- table(group_assignments)
    dyads_within_grps[i] <- sum(choose(counts_within_grps, 2))
  }
  
  return(dyads_within_grps)
}



dyad_replication_calculator(proposed_samplesizes = 5,
                            # target_samegrp_dyads = 10,
                            groupsize = 10,
                            popsize = 1000,
                            n_reps = 100)

proposed_ss <- seq(5:200)
n_reps_in <- 10
dyads_within_grps <- dyads_within_grps_gs2 <- dyads_within_grps_gs50 <- proposed_ss_mat <- matrix(NA, nrow = length(proposed_ss), ncol = n_reps_in)
for(i in 1:length(proposed_ss)){
  dyads_within_grps[i, ] <- dyad_replication_calculator(proposed_samplesizes = proposed_ss[i],
                                                        # target_samegrp_dyads = 10,
                                                        groupsize = 10,
                                                        popsize = 1000,
                                                        n_reps = n_reps_in)
  dyads_within_grps_gs2[i, ] <- dyad_replication_calculator(proposed_samplesizes = proposed_ss[i],
                                                            # target_samegrp_dyads = 10,
                                                            groupsize = 2,
                                                            popsize = 1000,
                                                            n_reps = n_reps_in)
  dyads_within_grps_gs50[i, ] <- dyad_replication_calculator(proposed_samplesizes = proposed_ss[i],
                                                             # target_samegrp_dyads = 10,
                                                             groupsize = 50,
                                                             popsize = 1000,
                                                             n_reps = n_reps_in)
  proposed_ss_mat[i, ] <- rep(proposed_ss[i], n_reps_in)
}

dyad_count_vec <- as.vector(dyads_within_grps)
dyad_count_vec_gs2 <- as.vector(dyads_within_grps_gs2)
dyad_count_vec_gs50 <- as.vector(dyads_within_grps_gs50)
ss_vec <- as.vector(proposed_ss_mat)


median_dyads <- apply(dyads_within_grps, 1, median)
median_dyads_gs2 <- apply(dyads_within_grps_gs2, 1, median)
median_dyads_gs50 <- apply(dyads_within_grps_gs50, 1, median)
# ss where median is closest to 10:
closest_to_10 <- which.min(abs(median_dyads - 10))
closest_to_10_gs2 <- which.min(abs(median_dyads_gs2 - 10))
closest_to_10_gs50 <- which.min(abs(median_dyads_gs50 - 10))
#proposed_ss[closest_to_10]

closest_to_20 <- which.min(abs(median_dyads - 20))
closest_to_20_gs2 <- which.min(abs(median_dyads_gs2 - 20))
closest_to_20_gs50 <- which.min(abs(median_dyads_gs50 - 20))
#proposed_ss[closest_to_20]

pdf("CombinatoricsGroupThing_v2.pdf", height = 6, width = 10)
par(mfrow = c(1, 3), mar = c(4, 6, 2, 2))
plot(dyad_count_vec_gs2 ~ ss_vec, xlab = "Animals sampled", ylim = c(1, 1300),
     log = "y", 
     ylab = "Number of tracked dyads\nwithin the same group", las = 1)
abline(h = 10, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_10_gs2])
abline(h = 20, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_20_gs2])
title("Groups of 2, Population of 1000")

plot(dyad_count_vec ~ ss_vec, xlab = "Animals sampled", ylim = c(1, 1300),
     log = "y",
     ylab = "Number of tracked dyads\nwithin the same group", las = 1)
abline(h = 10, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_10])
abline(h = 20, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_20])
title("Groups of 10, Population of 1000")

plot(dyad_count_vec_gs50 ~ ss_vec, xlab = "Animals sampled", ylim = c(1, 1300),
     log = "y",
     ylab = "Number of tracked dyads\nwithin the same group", las = 1)
abline(h = 10, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_10_gs50])
abline(h = 20, lty = 2, col = "grey80")
abline(v = proposed_ss[closest_to_20_gs50])
title("Groups of 50, Population of 1000")
dev.off()