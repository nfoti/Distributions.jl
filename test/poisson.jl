
using Distributions
using Base.Test

# Poisson sufficient statistics
x = [1, 2, 5, 7]
ss = suffstats(Poisson, x)

@test ss.sx == 15
@test ss.sw == 4
