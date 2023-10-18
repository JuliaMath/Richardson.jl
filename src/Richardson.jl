"""
The `Richardson` module provides a function `extrapolate` that
extrapolates a given function `f(x)` to `f(x0)`, evaluating
`f` only  at a geometric sequence of points `> x0`
(or optionally `< x0`).

The key algorithm is Richardson extrapolation using a Neville—Aitken
tableau, which adaptively increases the degree of an extrapolation
polynomial until convergence is achieved to a desired tolerance
(or convergence stalls due to e.g. floating-point errors).  This
allows one to obtain `f(x0)` to high-order accuracy, assuming
that `f(x0+h)` has a Taylor series or some other power
series in `h`.
"""
module Richardson

using LinearAlgebra

export extrapolate, extrapolate!

"""
    extrapolate(f, h; contract=0.125, x0=zero(h), power=1,
                      atol=0, rtol=atol>0 ? 0 : sqrt(ε),
                      maxeval=typemax(Int), breaktol=2)

Extrapolate `f(x)` to `f₀ ≈ f(x0)`, evaluating `f` only at `x > x0` points
(or `x < x0` if `h < 0`) using Richardson extrapolation starting at
`x=x₀+h`.  It returns a tuple `(f₀, err)` of the estimated `f(x0)`
and an error estimate.

The return value of `f` can be any type supporting `±` and `norm`
operations (i.e. a normed vector space).
Similarly, `h` and `x0` can be in any normed vector space,
in which case `extrapolate` performs Richardson extrapolation
of `f(x0+s*h)` to `s=0⁺` (i.e. it takes the limit as `x` goes
to `x0` along the `h` direction).

On each step of Richardson extrapolation, it shrinks `x-x0` by
a factor of `contract`, stopping when the estimated error is
`< max(rtol*norm(f₀), atol)`, when the estimated error
increases by more than `breaktol` (e.g. due to numerical errors in the
computation of `f`), when `f` returns a non-finite value (`NaN` or `Inf`),
 or when `f` has been evaluated `maxeval` times.   Note that
if the function may converge to zero, you may want
specify a nonzero `atol` (which cannot be set by default
because it depends on the scale/units of `f`); alternatively,
in such cases `extrapolate` will halt when it becomes
limited by the floating-point precision.   (Passing `breaktol=Inf`
can be useful to force `extrapolate` to continue shrinking `h` even
if polynomial extrapolation is initially failing to converge,
possibly at the cost of extraneous function evaluations.)

If `x0 = ±∞` (`±Inf`), then `extrapolate` computes the limit of
`f(x)` as `x ⟶ ±∞` using geometrically *increasing* values
of `h` (by factors of `1/contract`).

In general, the starting `h` should be large enough that `f(x0+h)`
can be computed accurately and efficiently (e.g. without
severe cancellation errors), but small enough that `f` does not
oscillate much between `x0` and `h`.  i.e. `h` should be a typical
scale over which the function `f` varies significantly.

Technically, Richardson extrapolation assumes that `f(x0+h)` can
be expanded in a power series in `h^power`, where the default
`power=1` corresponds to an ordinary Taylor series (i.e. assuming
`f` is analytic at `x0`).  If this is not true, you may obtain
slow convergence from `extrapolate`, but you can pass a different
value of `power` (e.g. `power=0.5`) if your `f` has some different
(Puiseux) power-series expansion.   Conversely, if `f` is
an *even* function around `x0`, i.e. `f(x0+h) == f(x0-h)`,
so that its Taylor series contains only *even* powers of `h`,
you can accelerate convergence by passing `power=2`.
"""
function extrapolate(f, h_::Number; contract::Number=oftype(float(real(h_)), 0.125), x0::Number=zero(h_), power::Number=1,
                     atol::Real=0, rtol::Real = atol > zero(atol) ? zero(one(float(real(x0+h_)))) : sqrt(eps(typeof(one(float(real(x0+h_)))))),
                     maxeval::Integer=typemax(Int), breaktol::Real=2)
    if isinf(x0)
        # use a change of variables x = 1/u
        return extrapolate(u -> f(inv(u)), inv(h_); rtol=rtol, atol=atol, maxeval=maxeval, contract = abs(contract) > 1 ? inv(contract) : contract, x0=inv(x0), power=power)
    end
    (rtol ≥ 0 && atol ≥ zero(atol)) || throw(ArgumentError("rtol and atol must be nonnegative"))
    breaktol > 0 || throw(ArgumentError("breaktol must be positive"))
    0 < abs(contract) < 1 || throw(ArgumentError("contract must be in (0,1)"))
    h::typeof(float(x0+h_*contract)) = h_
    invcontract = inv(contract)^power
    neville = [float(f(x0+h))] # the current diagonal of the Neville tableau
    f₀ = neville[1]
    err::typeof(float(norm(f₀))) = Inf
    numeval = 1
    while numeval < maxeval
        numeval += 1
        h *= contract
        push!(neville, f(x0+h))
        c = invcontract
        minerr′ = oftype(err, Inf)
        for i = length(neville)-1:-1:1
            old = neville[i]
            neville[i] = neville[i+1] + (neville[i+1] - neville[i]) / (c - 1)
            err′ = norm(neville[i] - old)
            minerr′ = min(minerr′, err′)
            if err′ < err
                f₀, err = neville[i], err′
            end
            c *= invcontract
        end
        (minerr′ > breaktol*err || !isfinite(minerr′)) && break # stop early if error increases too much
        err ≤ max(rtol*norm(f₀), atol) && break # converged
    end
    return (f₀, err)
end

# support non-numeric h as long as it is in a vector space
extrapolate(f, h; x0=zero(h), kws...) =
    extrapolate(s -> f(x0+s*h), one(norm(h)); kws...)

"""
    extrapolate(fh_itr; power=1,
                        atol=0, rtol=0, maxeval=typemax(Int), breaktol=Inf)

Similar to `extrapolate(f, h)`, performs Richardson extrapolation of a sequence of
values `f(h)` to `h → 0`, but takes an iterable collection `fh_itr` of a
sequence of `(f(h), h)` tuples (in order of decreasing `|h|`).

There is no `contract` keyword argument since the contraction factors are determined
by the sequence of `h` values (which need not contract by the same amount).  The
tolerances `atol` and `rtol` both default to `0` so that by default it examines
*all* of the values in the `fh_itr` collection.   Otherwise, the keyword arguments
have the same meanings as in `extrapolate(f, h)`.
"""
function extrapolate(fh_itr; power::Number=1, atol::Real=0, rtol::Real = 0,
                             breaktol::Real=Inf, maxeval::Integer=typemax(Int))
    (rtol ≥ 0 && atol ≥ zero(atol)) || throw(ArgumentError("rtol and atol must be nonnegative"))
    breaktol > 0 || throw(ArgumentError("breaktol must be positive"))
    Base.IteratorSize(fh_itr) isa Base.IsInfinite && iszero(atol) && iszero(rtol) && breaktol==Inf && maxeval==typemax(Int) &&
        throw(ArgumentError("stopping criteria must be supplied for infinite iterators"))
    itr = iterate(fh_itr)
    itr === nothing && throw(ArgumentError("(f,h) iterator must be non-empty"))
    (f,h), state = itr
    neville = [float(f)] # the current diagonal of the Neville tableau
    f₀ = neville[1]
    hvals = [float(h)]
    if Base.IteratorSize(fh_itr) isa Base.HasLength
        n = length(fh_itr)
        sizehint!(neville, n)
        sizehint!(hvals, n)
    end
    err::typeof(float(norm(f₀))) = Inf
    numeval = 1
    while numeval < maxeval
        numeval += 1
        itr = iterate(fh_itr, state)
        itr === nothing && break
        (f′,h′), state = itr
        push!(neville, f′)
        push!(hvals, h′)
        abs(h) > abs(h′) || throw(ArgumentError("|$h′| ≥ |$h| is not decreasing"))
        h = h′
        minerr′ = oftype(err, Inf)
        for i = length(neville)-1:-1:1
            c = (hvals[i] / h′)^power
            old = neville[i]
            neville[i] = neville[i+1] + (neville[i+1] - neville[i]) / (c - 1)
            err′ = norm(neville[i] - old)
            minerr′ = min(minerr′, err′)
            if err′ < err
                f₀, err = neville[i], err′
            end
        end
        (minerr′ > breaktol*err || !isfinite(minerr′)) && break # stop early if error increases too much
        err ≤ max(rtol*norm(f₀), atol) && break # converged
    end
    return (f₀, err)
end

"""
    extrapolate!(fh::AbstractVector; power=1, atol=0, rtol=0, maxeval=typemax(Int), breaktol=Inf)

Similar to `extrapolate(fh)`, performs Richardson extrapolation on an array `fh`
of `(f(h), h)` tuples (in order of decreasing `|h|`), but overwrites the array
`fh` in-place with intermediate calculations.

(Thus, the array `fh` must be a vector of `Tuple{T,H}` values, where `H<:Number` is
the type of `h` and `T` is the type of the extrapolated `f(0)` **result**.  This `T`
should be a floating-point type, i.e. `fh` should contain `float(f(h))` if the
function you are extrapolating is not already floating-point-valued.)
"""
function extrapolate!(fh::AbstractVector{<:Tuple{<:Any,<:Number}}; power::Number=1, atol::Real=0, rtol::Real = 0,
                      breaktol::Real=Inf, maxeval::Integer=typemax(Int))
    (rtol ≥ 0 && atol ≥ zero(atol)) || throw(ArgumentError("rtol and atol must be nonnegative"))
    breaktol > 0 || throw(ArgumentError("breaktol must be positive"))
    isempty(fh) && throw(ArgumentError("(f,h) array must be non-empty"))
    (f₀,h) = first(fh)
    err::typeof(float(norm(f₀))) = Inf
    numeval = 1
    maxeval = min(maxeval, length(fh))
    while numeval < maxeval
        numeval += 1
        (f′,h′) = fh[numeval + (firstindex(fh)-1)]
        abs(h) > abs(h′) || throw(ArgumentError("|$h′| ≥ |$h| is not decreasing"))
        h = h′
        minerr′ = oftype(err, Inf)
        fᵢ₊₁ = f′
        for i = numeval-1:-1:1
            fᵢ,hᵢ = fh[i + (firstindex(fh)-1)]
            c = (hᵢ / h′)^power
            fᵢ₊₁ += (fᵢ₊₁ - fᵢ) / (c - 1)
            fh[i + (firstindex(fh)-1)] = (fᵢ₊₁,hᵢ)
            err′ = norm(fᵢ₊₁ - fᵢ)
            minerr′ = min(minerr′, err′)
            if err′ < err
                f₀, err = fᵢ₊₁, err′
            end
        end
        (minerr′ > breaktol*err || !isfinite(minerr′)) && break # stop early if error increases too much
        err ≤ max(rtol*norm(f₀), atol) && break # converged
    end
    return (f₀, err)
end

end # module
