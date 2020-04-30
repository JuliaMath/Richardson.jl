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
that `f(x)` is analytic (has a Taylor series) around `x0`.
"""
module Richardson

using LinearAlgebra

export extrapolate

"""
    extrapolate(f, h; contract=0.125, x0=zero(h), power=1,
                      rtol=sqrt(ε), atol=0, maxeval=typemax(Int))

Extrapolate `f(x)` to `f₀ ≈ f(x0)`, evaluating `f` only at `x > x0` points
(or `x < x0` if `h < 0`) using Richardson extrapolation starting at
`x=x₀+h`.  It returns a tuple `(f₀, err)` of the estimated `f(x0)`
and an error estimate.

More generally, `h` and `x0` can be in an arbitrary vector space,
in which case `extrapolate` performs Richardson extrapolation
of `f(x0+s*h)` to `s=0⁺` (i.e. it takes the limit as `x` goes
to `x0` along the `h` direction).

On each step of Richardson extrapolation, it shrinks `x-x0` by
a factor of `contract`, stopping when the estimated error is
`< max(rtol*norm(f₀), atol)`, when the estimated error starts to
increase (e.g. due to numerical errors in the computation of `f`),
or when `f` has been evaluated `maxeval` times.   Note that
if the function may converge to zero, you should probably
specify a nonzero `atol` (which cannot be set by default
because it depends on the scale/units of `f`).

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
so that its Taylor series containsonly *even* powers of `h`,
you can accelerate convergence by passing `power=2`.
"""
function extrapolate(f, h_::Number; contract::Real=0.125, x0::Number=zero(h_), power::Number=1,
                     rtol::Real=sqrt(eps(typeof(float(h_)))), atol::Real=0, maxeval=typemax(Int))
    if isinf(x0)
        # use a change of variables x = 1/u
        return extrapolate(u -> f(inv(u)), inv(h_); rtol=rtol, atol=atol, maxeval=maxeval, contract = contract > 1 ? inv(contract) : contract, x0=inv(x0), power=power)
    end
    (rtol ≥ 0 && atol ≥ 0) || throw(ArgumentError("rtol and atol must be nonnegative"))
    0 < contract < 1 || throw(ArgumentError("contract must be in (0,1)"))
    h = oftype(float(x0+h_), h_)
    invcontract = inv(contract)^power
    neville = [f(x0+h)] # the current diagonal of the Neville tableau
    f₀ = neville[1]
    err = oftype(norm(f₀), Inf)
    numeval = 1
    while numeval < maxeval
        numeval += 1
        h = oftype(h, h * contract)
        push!(neville, f(x0+h))
        c = invcontract
        minerr′ = oftype(err, Inf)
        for i = length(neville)-1:-1:1
            old = neville[i]
            neville[i] = neville[i+1] + (neville[i+1] - neville[i]) / (c - 1)
            err′ = norm(neville[i] - old)
            minerr′ = min(minerr′, err′)
            if err′ < err
                f₀, err = neville[i], oftype(err, err′)
            end
            c *= invcontract
        end
        (minerr′ > 2err || !isfinite(minerr′)) && break # stop early if error increases too much
        err ≤ max(rtol*norm(f₀), atol) && break # converged
    end
    return (f₀, err)
end

# support non-numeric h as long as it is in a vector space
extrapolate(f, h; x0=zero(h), kws...) =
    extrapolate(s -> f(x0+s*h), one(norm(h)); kws...)

end # module
