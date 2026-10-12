# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Documenter
using DocumenterCodeBlocks
using PiecewiseLinearOpt

include("figures.jl")
univariate_figure(joinpath(@__DIR__, "src", "assets", "univariate.svg"))
triangulation_figure(joinpath(@__DIR__, "src", "assets", "triangulations.svg"))

makedocs(;
    sitename = "PiecewiseLinearOpt.jl",
    authors = "Joey Huchette and contributors",
    modules = [PiecewiseLinearOpt],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://jump.dev/PiecewiseLinearOpt.jl/stable",
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/functions.md", "manual/models.md"],
        "methods.md",
        "extending.md",
        "migration.md",
        "references.md",
        "api.md",
    ],
    plugins = [CodeBlocks()],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/jump-dev/PiecewiseLinearOpt.jl.git",
    push_preview = true,
)
