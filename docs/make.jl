
# julia --project=docs docs/make.jl 
using PiecewiseLinearOpt
using Documenter

DocMeta.setdocmeta!(
    PiecewiseLinearOpt,
    :DocTestSetup,
    :(using PiecewiseLinearOpt, JuMP);
    recursive = true,
)

makedocs(;
    modules = [PiecewiseLinearOpt],
    sitename = "PiecewiseLinearOpt.jl",
    authors = "Joey Huchette and contributors",
    format = Documenter.HTML(;
        canonical = "https://jump-dev.github.io/PiecewiseLinearOpt.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Formulation Methods" => "methods.md",
        "API Reference" => "api.md",
    ],
    warnonly = true,
)

deploydocs(;
    repo = "github.com/jump-dev/PiecewiseLinearOpt.jl",
    devbranch = "master",
)
