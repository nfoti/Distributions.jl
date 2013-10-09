

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

function posterior_pred(pri::Gamma, ss::PoissonStats)
    nb = rate(pri) + ss.sw
    mu = (pri.shape + ss.sx) / nb
    var = mu * (nb+1.) / nb 
    r = mu.^2 / (var - mu)
    return NegativeBinomial(r, r / (r + mu))
end

function posterior_pred{T<:Real}(pri::Gamma, ::Type{Poisson}, x::Vector{T})
    return posterior_pred(pri, suffstats(Poisson, x))
end

function posterior_pred{T<:Real}(pri::Gamma, ::Type{Poisson}, x::Vector{T}, w::Array{Float64})
    return posterior_pred(pri, suffstats(Poisson, x, w))
end

### Normal, known sigma

### Normal, known mean

### Normal, unknown mean and variance

### MVN, known Sigma

### MVN, unkown mean and covariance
