using Richardson, Test

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

    empty!(X)
    val, err = extrapolate(-1.0, rtol=1e-10, contract=0.1) do x
        push!(X, x)
        sin(x)/x
    end
    @test val ≈ 1 rtol=1e-10
    @test err < 1e-10
    @test all(x -> x < 0, X)
    @test length(X) == 6
end
