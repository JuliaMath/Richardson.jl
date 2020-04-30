using Richardson, Test, LinearAlgebra

@testset "basic tests" begin
    X = Float64[]
    val, err = extrapolate(1.0, rtol=1e-10, contract=0.1) do x
        push!(X, x)
        sin(x)/x
    end
    @test val ≈ 1 rtol=1e-10
    @test err < 1e-10
    @test all(x -> x > 0, X)
    @test length(X) == 6

    X2 = Float64[]
    val, err = extrapolate(1.0, rtol=1e-10, contract=0.1, maxeval=3) do x
        push!(X2, x)
        sin(x)/x
    end
    @test X2 == X[1:3]
    @test val ≈ 1 rtol=1e-2
    @test err < 1e-2

    empty!(X)
    val, err = extrapolate(-1.0, rtol=1e-10, contract=0.1) do x
        push!(X, x)
        sin(x)/x
    end
    @test val ≈ 1 rtol=1e-10
    @test err < 1e-10
    @test all(x -> x < 0, X)
    @test length(X) == 6

    empty!(X)
    val, err = extrapolate(1.0, x0=Inf) do x
        push!(X, x)
        (x^2 + 3x - 2) / (x^2 + 5)
    end
    @test val ≈ 1 rtol=1e-12
    @test err < 1e-10
    @test length(X) == 7

    # make sure this terminates rather than looping on NaN values
    val, err = extrapolate(x -> log(x)/sqrt(x-1), 0.2, x0=1.0)
    @test abs(val) < 1e-6 && abs(err) < 1e-6

    # vector-valued function support
    val, err = extrapolate(x -> [sin(x)/x, cos(x)], 0.1, rtol=1e-10)
    @test val ≈ [1,1] rtol=1e-10
    @test err < norm([1,1])*1e-10

    # Float32 specified via starting h
    val, err = extrapolate(x -> sin(x)/x, 1.0f0, rtol=1e-5)
    @test val ≈ 1.0 rtol=1e-5
    @test err < 1e-7
    @test val isa Float32
    @test err isa Float32

    # extrapolation in a vector space
    val, err = extrapolate([0.1,-0.2],x0=[0.0,1.0]) do x
        sin(x[1])/x[1]*log(x[2])/(x[2]-1)
    end
    @test val ≈ 1 rtol=1e-12
    @test err < 1e-8
end