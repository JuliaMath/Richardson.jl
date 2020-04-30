# Richardson package for Julia

[![Build Status](https://travis-ci.org/JuliaMath/Richardson.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Richardson.jl)

The `Richardson` package provides a function `extrapolate` that
extrapolates a given function `f(x)` to `f(x0)`, evaluating
`f` only  at a geometric sequence of points `> x0`
(or optionally `< x0`).

The key algorithm is [Richardson extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation) using a Neville–Aitken
tableau, which adaptively increases the degree of an extrapolation
polynomial until convergence is achieved to a desired tolerance
(or convergence stalls due to e.g. floating-point errors).  This
allows one to obtain `f(x0)` to high-order accuracy, assuming
that `f(x0+h)` is has a Taylor series or some other power
series in `h`.   (See e.g. [these course notes by Prof. Flaherty at RPI](http://www.cs.rpi.edu/~flaherje/pdf/ode4.pdf).)

## Usage

```jl
extrapolate(f, h; contract=0.125, x0=zero(h),
                  atol=0, rtol=atol>0 ? 0 : sqrt(ε), maxeval=typemax(Int))
```

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

## Examples

For example, let's extrapolate `sin(x)/x` to `x=0` (where the correct answer is `1`) starting at `x=1`, printing out the `x` value at each step so that we can see what the algorithm is doing.

(Since `f` is passed as the first argument to `extrapolate`, we
can use Julia's `do` syntax to conveniently define a multi-line
anonymous function to pass.)
```jl
extrapolate(1.0, rtol=1e-10) do x
    @show x
    sin(x)/x
end
```
giving the output:
```
x = 1.0
x = 0.125
x = 0.015625
x = 0.001953125
x = 0.000244140625
x = 3.0517578125e-5
(1.0000000000000002, 2.0838886172214188e-13)
```
That is, it evaluates our function `sin(x)/x` for 6 different values of `x` and returns `1.0000000000000002`, which is accurate to machine precision (the error is `≈ 2.2e-16`).  The returned error estimate of `2e-13` is conservative, which is typical for extrapolating well-behaved functions.

Since `sin(x)/x` is an *even* (symmetric) function around `x=0`,
its Taylor series contains only even powers of `x`.  We can
exploit this fact to *accelerate convergence for even functions* by
passing `power=2` to `extrapolate`:
```jl
extrapolate(1.0, rtol=1e-10, power=2) do x
    @show x
    sin(x)/x
end
```
gives
```
x = 1.0
x = 0.125
x = 0.015625
x = 0.001953125
x = 0.000244140625
(1.0, 0.0)
```
which converged to machine precision (in fact, the exactly rounded result) in only 5 function evaluations (1 fewer than above).

A classic use of Richardson extrapolation is accurately evaluating derivatives via [finite-difference approximations](https://en.wikipedia.org/wiki/Finite_difference) (although analytical derivatives, e.g. by automatic differentiation, are of course vastly more efficient when they are available); a classic description of this combination is found in [Ridders (1982)](https://www.sciencedirect.com/science/article/abs/pii/S0141119582800570).   In this example, we use Richardson extrapolation on the forward-difference approximation `f'(x) ≈ (f(x+h)-f(x))/h`, for which the error decreases as `O(h)` but a naive application to a very small `h` will yield a huge [cancellation error](https://en.wikipedia.org/wiki/Loss_of_significance) from floating-point roundoff effects.   We differentiate `f(x)=sin(x)` at `x=1`, for which the correct answer is `cos(1) ≈ 0.5403023058681397174009366...`, starting with `h=0.1`
```jl
extrapolate(0.1, rtol=0) do h
    @show h
    (sin(1+h) - sin(1)) / h
end
```
Although we gave an `rtol` of `0`, the `extrapolate` function will terminate after a finite number of steps when it detects that improvements are limited by floating-point error:
```
h = 0.1
h = 0.0125
h = 0.0015625
h = 0.0001953125
h = 2.44140625e-5
h = 3.0517578125e-6
(0.5403023058683176, 1.7075230118734908e-12)
```
The output `0.5403023058683176` differs from `cos(1)` by `≈ 1.779e-13`, so in this case the returned error estimate is only a little conservative.   Unlike the previous example, `extrapolate` is not able
to attain machine precision (the floating-point cancellation error in this function is quite severe for small `h`!), but it is able to get surprisingly close.

Using the `x0` keyword argument, we can compute the limit of `f(x)`
as `x ⟶ x0`.  In fact, you can pass `x0 = Inf` to compute a limit as
`x ⟶ ∞` (which is accomplished internally by a change of variables `x = 1/u` and performing Richardson extrapolation to `u=0`). For example:
```jl
extrapolate(1.0, x0=Inf) do x
    @show x
    (x^2 + 3x - 2) / (x^2 + 5)
end
```
gives
```
x = 1.0
x = 8.0
x = 64.0
x = 512.0
x = 4096.0
x = 32768.0
x = 262144.0
(1.0000000000000002, 1.2938539128981574e-12)
```
which is the correct result (`1.0`) to machine precision.