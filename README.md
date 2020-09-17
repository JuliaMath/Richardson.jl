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
that `f(x0+h)` has a Taylor series or some other power
series in `h`.   (See e.g. [these course notes by Prof. Flaherty at RPI](http://www.cs.rpi.edu/~flaherje/pdf/ode4.pdf).)

## Usage

```jl
extrapolate(f, h; contract=0.125, x0=zero(h),
                  atol=0, rtol=atol>0 ? 0 : sqrt(ε), maxeval=typemax(Int), breaktol=2)
```

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
oscillate much between `x0` and `x0+h`.  i.e. `h` should be a typical
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
which converged to machine precision (in fact, the exact result) in only 5 function evaluations (1 fewer than above).

### Infinite limits

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

#### Extrapolating series

One nice application of infinite limits is extrapolating infinite series.   If we start with an integer `x`, then the default `contract=0.125` will increase `x` by a factor of `8.0` on each iteration, so `x` will always be an exact integer and we can use it as the number of terms in a series.

For example, suppose we are computing the infinite series `1/1² + 1/2² + 1/3² + ⋯`.  This is the famous [Basel problem](https://en.wikipedia.org/wiki/Basel_problem), and it converges to `π²/6 ≈ 1.644934066848226…`.   If we compute it by brute force, however, we need quite a few terms to get high accuracy:
```jl
julia> sum(n -> 1/n^2, 1:100) - π^2/6
-0.009950166663333482

julia> sum(n -> 1/n^2, 1:10^4) - π^2/6
-9.99950001654426e-5

julia> sum(n -> 1/n^2, 1:10^9) - π^2/6
-9.999985284281365e-10
```
Even with 10⁹ terms we get only about 9 digits.   Instead, we can use `extrapolate` (starting at 1 term):
```jl
julia> val, err = extrapolate(1, x0=Inf, rtol=0) do N
           @show N
           sum(n -> 1/n^2, 1:Int(N))
       end
N = 1.0
N = 8.0
N = 64.0
N = 512.0
N = 4096.0
N = 32768.0
N = 262144.0
N = 2.097152e6
(1.644934066848228, 0.0)

julia> val - π^2/6
1.5543122344752192e-15
```
By `32768` terms, the extrapolated value is accurate to about 15 digits.

### Numerical derivatives

A classic application of Richardson extrapolation is the accurate evaluation of derivatives via [finite-difference approximations](https://en.wikipedia.org/wiki/Finite_difference) (although analytical derivatives, e.g. by automatic differentiation, are of course vastly more efficient when they are available).  In this example, we use Richardson extrapolation on the forward-difference approximation `f'(x) ≈ (f(x+h)-f(x))/h`, for which the error decreases as `O(h)` but a naive application to a very small `h` will yield a huge [cancellation error](https://en.wikipedia.org/wiki/Loss_of_significance) from floating-point roundoff effects.   We differentiate `f(x)=sin(x)` at `x=1`, for which the correct answer is `cos(1) ≈ 0.5403023058681397174009366...`, starting with `h=0.1`
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
The output `0.5403023058683176` differs from `cos(1)` by `≈ 1.779e-13`, so in this case the returned error estimate is only a little conservative.   Unlike the `sin(x)/x` example, `extrapolate` is not able
to attain machine precision (the floating-point cancellation error in this function is quite severe for small `h`!), but it is able to get surprisingly close.

Another possibility for a finite-difference/Richardson combination was suggested by [Ridders (1982)](https://www.sciencedirect.com/science/article/abs/pii/S0141119582800570), who computed both `f'(x)` and `f''(x)` (the first and second derivatives) simultaneously using a center-difference approximation, which requires two new `f(x)` evaluations for each `h`.  In particular, the center-difference approximations are `f'(x) ≈ (f(x+h)-f(x-h))/2h` and `f''(x) ≈ (f(x+h)-2f(x)+f(x-h))/h²`, both of which have errors that go as `O(h²)`.   We can plug both of these functions *simultaneously* into `extrapolate` (so that they share `f(x±h)` evaluations) by using a vector-valued function returning `[f', f'']`.   Moreover, since these center-difference approximations are even functions of `h` (identical for `±h`), we can pass `power=2` to `extrapolate` in order to exploit the even-power Taylor expansion.  Here is a function implementing both of these ideas:

```jl
# returns (f'(x), f''(x))
function riddersderiv2(f, x, h; atol=0, rtol=atol>0 ? 0 : sqrt(eps(typeof(float(real(x+h))))), contract=0.5)
    f₀ = f(x)
    val, err = extrapolate(h, atol=atol, rtol=rtol, contract=contract, power=2) do h
        f₊, f₋ = f(x+h), f(x-h)
        [(f₊-f₋)/2h, (f₊-2f₀+f₋)/h^2]
    end
    return val[1], val[2]
end
```
(This code could be made even more efficient by using [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) for the `[f', f'']` vector.)   The original paper by Ridders accomplishes something similar in `< 20` lines of [TI-59 calculator](https://en.wikipedia.org/wiki/TI-59_/_TI-58) code, by the way; so much for high-level languages!

For example,
```jl
julia> riddersderiv2(1, 0.1, rtol=0) do x
           @show x
           sin(x)
       end
x = 1
x = 1.1
x = 0.9
x = 1.05
x = 0.95
x = 1.025
x = 0.975
x = 1.0125
x = 0.9875
x = 1.00625
x = 0.99375
(0.5403023058681394, -0.841470984807975)
```
evaluates the first and second derivatives of `sin(x)` at `x=1` and obtains the correct answer `(cos(1), -sin(1))` to about 15 and 13 decimal digits, respectively, using 11 function evaluations.

### Handling problematic convergence

It is useful to consider a finite-difference approximation for the derivative of
the function `1/x` at some `x ≠ 0`: i.e. computing the limit of `f(h) = (1/(x+h) - 1/x) / h`
as `h` goes to zero similar to above.

This function `f(h)` has a [pole](https://en.wikipedia.org/wiki/Zeros_and_poles) at `h=-x`, i.e. `f(-x)` blows up.   This means
that the Taylor series of `f(h)` [only converges](https://en.wikipedia.org/wiki/Radius_of_convergence) for `h` values small enough to avoid this pole, and in
fact for `|h| < |x|`.   Since Richardson extrapolation is essentially
approximating the Taylor series, this means that the extrapolation process doesn't
converge if the starting `h` is too large, and `extrapolation` will give up and
halt with the wrong answer.

This lack of convergence is easily observed: set `x=0.01` (where the correct derivative of `1/x` is `-10000`) and consider what happens for a starting `h`
that is too large compared to `|x|`:
```jl
julia> extrapolate(1.0) do h
           @show h
           x = 0.01
           (1/(x+h) - 1/x) / h
       end
h = 1.0
h = 0.125
h = 0.015625
(-832.4165749908325, 733.4066740007335)
```
Before reached an `|h| < 0.01` where the power series could begin to converge, `extrapolate` gave up and returned a wrong answer (with a large error estimate to let you know that the result is garbage)!   In contrast, if we start with a small enough `h` then it converges just fine and returns the correct answer (`-10000`) to nearly machine precision:
```jl
julia> extrapolate(0.01) do h
           @show h
           x = 0.01
           (1/(x+h) - 1/x) / h
       end
h = 0.01
h = 0.00125
h = 0.00015625
h = 1.953125e-5
h = 2.44140625e-6
h = 3.0517578125e-7
(-10000.000000000211, 4.066770998178981e-6)
```

Of course, if you know that your function blows up like this, it is easy to choose good initial `h`, but how can we persuade `extrapolate` to do a better job automatically?

The trick is to use the `breaktol` keyword argument.  `breaktol` defaults to `2`,
which means that `extrapolate` gives up if the best error estimate increases by
more than a factor of 2 from one iteration to the next.  Ordinarily, this
kind of breakdown in convergence arises because you hit the limits of floating-point precision, and halting extrapolation is the right thing to do.  Here, however, it would converge if we just continued shrinking `h`.   So, we simply set `breaktol=Inf` to force extrapolation to continue, which works even for a large initial `h=1000.0`:
```jl
julia> extrapolate(1000.0, breaktol=Inf) do h
           @show h
           x = 0.01
           (1/(x+h) - 1/x) / h
       end
h = 1000.0
h = 125.0
h = 15.625
h = 1.953125
h = 0.244140625
h = 0.030517578125
h = 0.003814697265625
h = 0.000476837158203125
h = 5.9604644775390625e-5
h = 7.450580596923828e-6
h = 9.313225746154785e-7
h = 1.1641532182693481e-7
(-10000.000000029328, 5.8616933529265225e-8)
```
Not that it continues extrapolating until it reaches small `h` values where the power-series converges, and in the end it again returns the correct answer to nearly machine precision (and would reach machine precision if we set a smaller `rtol`).  (The `extrapolate` function *automatically* discards the initial points where the polynomial extrapolation fails.)