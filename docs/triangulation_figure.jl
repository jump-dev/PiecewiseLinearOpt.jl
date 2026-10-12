# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Draws the grid triangulations of a bivariate function as an SVG figure, with
# each triangle shaded by the mean function value at its corners.

using PiecewiseLinearOpt
import PiecewiseLinearOpt: segments, input_values, output_values

# The viridis color map, sampled at five points
const _COLORS = [
    (68, 1, 84),
    (59, 82, 139),
    (33, 144, 141),
    (93, 200, 99),
    (253, 231, 37),
]

function _color(t)
    s = clamp(t, 0, 1) * (length(_COLORS) - 1)
    k = min(floor(Int, s), length(_COLORS) - 2)
    a, b = _COLORS[k+1], _COLORS[k+2]
    c = map((u, v) -> round(Int, u + (s - k) * (v - u)), a, b)
    return "rgb($(c[1]),$(c[2]),$(c[3]))"
end

function triangulation_figure(
    path;
    f = (a, b) -> sin(3a) * cos(3b),
    breakpoints = 0:0.25:1,
    triangulations = [
        "K1Triangulation()" => K1Triangulation(),
        "UnionJackTriangulation()" => UnionJackTriangulation(),
        "UpperTriangulation()" => UpperTriangulation(),
        "LowerTriangulation()" => LowerTriangulation(),
        "BestFitTriangulation()" => BestFitTriangulation(),
        "RandomTriangulation()" => RandomTriangulation(),
    ],
    columns = 3,
)
    size, gap, label = 180, 40, 30
    rows = cld(length(triangulations), columns)
    width = columns * size + (columns + 1) * gap
    height = rows * (size + label) + (rows + 1) * gap / 2
    values = [f(a, b) for a in breakpoints, b in breakpoints]
    lo, hi = extrema(values)
    x_lo, x_hi = extrema(breakpoints)
    io = IOBuffer()
    println(
        io,
        """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $width $height" width="$width" font-family="sans-serif" font-size="13">""",
    )
    println(io, """<rect width="$width" height="$height" fill="white"/>""")
    for (k, (name, triangulation)) in enumerate(triangulations)
        row, column = divrem(k - 1, columns)
        x0 = gap + column * (size + gap)
        y0 = gap / 2 + row * (size + label + gap / 2)
        # The y axis points up, as in the documentation of the triangulations
        px(a) = x0 + size * (a - x_lo) / (x_hi - x_lo)
        py(b) = y0 + size * (1 - (b - x_lo) / (x_hi - x_lo))
        pwl = BivariatePWLFunction(
            breakpoints,
            breakpoints,
            f;
            triangulation = triangulation,
        )
        for segment in segments(pwl)
            points = join(
                ("$(px(a)),$(py(b))" for (a, b) in input_values(segment)),
                " ",
            )
            mean = sum(only, output_values(segment)) / 3
            println(
                io,
                """<polygon points="$points" fill="$(_color((mean - lo) / (hi - lo)))" stroke="white" stroke-width="1.5"/>""",
            )
        end
        println(
            io,
            """<text x="$(x0 + size / 2)" y="$(y0 + size + 20)" text-anchor="middle" fill="#333">$name</text>""",
        )
    end
    println(io, "</svg>")
    mkpath(dirname(path))
    write(path, take!(io))
    return path
end
