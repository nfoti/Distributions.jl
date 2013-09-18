
using Distributions
using Base.Test

bb = BetaBinomial(20, 1., 3.)

@test_approx_eq pdf(bb, 4) gamma(21.)/(gamma(5.)*gamma(17.))*gamma(5.)*gamma(19.)/gamma(24.)*gamma(4.)/(gamma(1.)*gamma(3.))
@test_approx_eq logpdf(bb, 4) lgamma(21.) - lgamma(5.) - lgamma(17.) + lgamma(5.) + lgamma(19.) - lgamma(24.) + lgamma(4.) - lgamma(1.) - lgamma(3.)

@test insupport(BetaBinomial, 5)
@test insupport(BetaBinomial, 0)
@test !insupport(BetaBinomial, -1)

@test_approx_eq mean(bb) 20./4.
@test_approx_eq var(bb) 20.*3.*24./(4.^2*5.)
@test_approx_eq skewness(bb) 44.*2./6.*sqrt(5./(20.*3.*24.))
