
immutable BetaBinomial <: DiscreteUnivariateDistribution
    n::Int64
    a::Float64
    b::Float64
    function BetaBinomial(n::Int, a::Real, b::Real)
        (n > zero(n) && a > zero(a) && b > zero(b)) || error("alpha and beta must be positive")
        new(int64(n), float64(a), float64(b))
    end
end

function pdf(d::BetaBinomial, k::Int)
    n = d.n
    a = d.a
    b = d.b
    nmk = n - k
    ab = a + b
    return gamma(n+1)/(gamma(k+1)*gamma(nmk+1))*gamma(k+a)*gamma(nmk+b)/gamma(n+ab)*gamma(ab)/(gamma(a)*gamma(b))
end

function logpdf(d::BetaBinomial, k::Int)
    n = d.n
    a = d.a
    b = d.b
    nmk = n - k
    ab = a + b
    return lgamma(n+1) - lgamma(k + 1) - gamma(nmk + 1) + lgamma(k + a) + lgamma(nmk + b) - lgamma(n + ab) + lgamma(ab) - lgamma(a) - lgamma(b)
end

insupport(::BetaBinomial, k::Int) = zero(k) <= k && isfinite(k)
insupport(::Type{BetaBinomial}, k::Int) = zero(k) <= k && isfinite(k)

mean(d::BetaBinomial) = d.n*d.a/(d.a + d.b)

var(d::BetaBinomial) = d.n*d.a*d.b*(d.a + d.b + d.n)/((d.a+d.b).^2*(d.a+d.b+1.))

function skewness(d::BetaBinomial)
    (d.a+d.b+2.*d.n)*(d.b-d.a)/(d.a+d.b+2.) * sqrt((1.+d.a+d.b)/(d.n*d.a*d.b*(d.n+d.a+d.b)))
end

# fit_mle can be implemented with an iterative method.  For details, see Minka,
# "Estimating a Dirichlet distribuiton", 2003.
