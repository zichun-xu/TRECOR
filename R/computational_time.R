n = 500
p = 25
d = 3
R = 1
times = c()
for(i in 1:100){
  tryCatch({
  load(paste0("~/Research_projects/tree_cov_reg/results/time/n_", n, "_p_", p, "_sparse_0.5_tree_LTNM_R_", R, "_cov_reg_simulation_time_", i))
  times = c(times, as.numeric(time))
  },
  error = function(e){}
  )
}
print(mean(times))
