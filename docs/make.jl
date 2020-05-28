using Documenter, eFEM

push!(LOAD_PATH,"../src/")

makedocs(;
    modules=[eFEM],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pseastham/eFEM.jl/blob/{commit}{path}#L{line}",
    sitename="eFEM.jl",
    authors="Patrick Eastham",
    assets=String[],
)

deploydocs(;
    repo="github.com/pseastham/eFEM.jl",
)
