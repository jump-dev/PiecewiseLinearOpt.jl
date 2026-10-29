# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Documenter

Documenter.makedocs(;
    sitename = "PiecewiseLinearOpt.jl",
    authors = "Joey Huchette and contributors",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

Documenter.deploydocs(;
    repo = "github.com/jump-dev/PiecewiseLinearOpt.jl.git",
    push_preview = true,
)
