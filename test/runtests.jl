using SafeTestsets

@time begin
    @time @safetestset "IO Tests" begin include("io_test.jl") end
    @time @safetestset "Full Runs Tests" begin include("full_run_test.jl") end
    @time @safetestset "Post-process Tests" begin include("postprocess_test.jl") end
    @time @safetestset "Quadrature Tests" begin include("quadrature_test.jl") end
    @time @safetestset "Struct Tests" begin include("struct_test.jl") end
    @time @safetestset "Tracer Tests" begin include("tracers_test.jl") end
end
