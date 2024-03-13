using ClimFluids
using Documenter

DocMeta.setdocmeta!(ClimFluids, :DocTestSetup, :(using ClimFluids); recursive=true)

makedocs(;
    modules=[ClimFluids],
    authors="The ClimFlows contributors",
    sitename="ClimFluids.jl",
    format=Documenter.HTML(;
        canonical="https://ClimFlows.github.io/ClimFluids.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ClimFlows/ClimFluids.jl",
    devbranch="main",
)
