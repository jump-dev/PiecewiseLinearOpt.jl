# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Documenter
using DocumenterCodeBlocks
using Literate
using PiecewiseLinearOpt

# The plots of the tutorials are drawn by GR without a display
ENV["GKSwstype"] = "100"

include("figures.jl")
univariate_figure(joinpath(@__DIR__, "src", "assets", "univariate.svg"))
triangulation_figure(joinpath(@__DIR__, "src", "assets", "triangulations.svg"))

# The tutorials are Literate.jl scripts in src/tutorials, which are converted
# to Markdown pages next to them. The generated pages are not committed.
const TUTORIALS = ["transportation.jl", "hydropower.jl"]
const TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorials")
tutorial_pages = map(TUTORIALS) do file
    Literate.markdown(
        joinpath(TUTORIAL_DIR, file),
        TUTORIAL_DIR;
        documenter = true,
        credit = false,
    )
    return "tutorials/" * replace(file, ".jl" => ".md")
end

# A full build runs every example and tutorial, which takes several minutes.
# For edits to the text, set DOCS_DRAFT=true to build in draft mode, which
# renders the pages without running the code:
#
#     DOCS_DRAFT=true julia --project=docs docs/make.jl
#
# A tutorial under work can be run on its own, as a plain Julia script, with
# `julia --project=docs docs/src/tutorials/<name>.jl`.
makedocs(;
    sitename = "PiecewiseLinearOpt.jl",
    authors = "Joey Huchette and contributors",
    modules = [PiecewiseLinearOpt],
    draft = get(ENV, "DOCS_DRAFT", "false") == "true",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://jump.dev/PiecewiseLinearOpt.jl/stable",
        # The tutorial pages contain several plots
        size_threshold_ignore = tutorial_pages,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" =>
            ["manual/functions.md", "manual/models.md", "manual/methods.md"],
        "Tutorials" => tutorial_pages,
        "Background" => ["background/references.md"],
        "Developers" =>
            ["developers/extending.md", "developers/biclique_covers.md"],
        "API reference" => "api.md",
        "Release notes" => ["changelog.md", "migration.md"],
    ],
    plugins = [CodeBlocks()],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/jump-dev/PiecewiseLinearOpt.jl.git",
    push_preview = true,
)
