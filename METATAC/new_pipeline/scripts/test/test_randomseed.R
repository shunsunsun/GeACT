


library(parallel)
# RNGkind("L'Ecuyer-CMRG")

set.seed(0)

# mc.reset.stream()
mclapply(1:10, function(i) {
  # set.seed(0)
  sample(1:10, 3)
})

