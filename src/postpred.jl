

# These are kind of uninteresting b/c the data doesn't really show up, but it's
# implemented to keep the interface consistent.
function posterior_pred{T<:Real}(pri::Beta, ::Type{Bernoulli}, x::Vector{T})
    return BetaBernoulli(length(x), pri.alpha, pri.beta)
end

function posterior_pred{T<:Real}(pri::Beta, ::Type{Binomial}, n::Int, x::Vector{T})
    return BetaBernoulli(n, pri.alpha, pri.beta)
end


### Dirichlet: Need DirichletMultinomial distribution

### Gamma-Poisson: Makes NegativeBinomial

### Normal, known sigma

### Normal, known mean

### Normal, unknown mean and variance

### MVN, known Sigma

### MVN, unkown mean and covariance
