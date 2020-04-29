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

export extrapolate

"""
    extrapolate(f, h; rtol=sqrt(ε), atol=0, contract=0.125, x0=0)

Extrapolate `f(x)` to `f₀ ≈ f(x0)`, evaluating `f` only at `x > x0` points
(or `x < x0` if `h < 0`) using Richardson's extrapolation starting at
`x=x₀+h`.  It returns a tuple `(f₀, err)` of the estimated `f(x0)`
and an error estimate.

On each step of Richardson's extrapolation, it shrinks `x-x0` by
a factor of `contract`, stopping when the estimated error is
`< max(rtol*f₀, atol)` or when the estimated error starts to
increase (e.g. due to numerical errors in the computation of `f`).

If `x0 = ±∞` (`±Inf`), then `extrapolate` computes the limit of
`f(x)` as `x ⟶ ±∞` using geometrically *increasing* values
of `h` (by factors of `1/contract`).

In general, the starting `h` should be large enough that `f(x0+h)`
can be computed accurately and efficiently (e.g. without
severe cancellation errors), but small enough that `f` does not
oscillate much between `x0` and `h`.  i.e. `h` should be a typical
scale over which the function `f` varies significantly.
"""
function extrapolate(f, h_::Number; rtol::Real=sqrt(eps(typeof(float(h_)))), atol::Real=0, contract::Real=0.125, x0::Number=0)
    if isinf(x0)
        # use a change of variables x = 1/u
        return extrapolate(u -> f(inv(u)), inv(h_); rtol=rtol, atol=atol, contract = contract > 1 ? inv(contract) : contract, x0=inv(x0))
    end
    (rtol ≥ 0 && atol ≥ 0) || throw(ArgumentError("rtol and atol must be nonnegative"))
    0 < contract < 1 || throw(ArgumentError("contract must be in (0,1)"))
    h = float(x0+h_)
    invcontract = inv(contract)
    neville = [f(x0+h)] # the current diagonal of the Neville tableau
    f₀ = neville[1]
    err = oftype(abs(f₀), Inf)
    while true
        h *= contract
        push!(neville, f(x0+h))
        c = invcontract
        minerr′ = oftype(err, Inf)
        for i = length(neville)-1:-1:1
            old = neville[i]
            neville[i] = neville[i+1] + (neville[i+1] - neville[i]) / (c - 1)
            err′ = abs(neville[i] - old)
            minerr′ = min(minerr′, err′)
            if err′ < err
                f₀, err = neville[i], err′
            end
            c *= invcontract
        end
        minerr′ > 2err && break # stop early if error increases too much
        err ≤ max(rtol*f₀, atol) && break # converged
    end
    return (f₀, err)
end

end # module
