# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# SVG figures for the manual, drawn from the functions of the package so that
# they always match it.

using PiecewiseLinearOpt
import PiecewiseLinearOpt: segments, input_values, output_values

# The viridis color map, sampled at five points
const _COLORS =
    [(68, 1, 84), (59, 82, 139), (33, 144, 141), (93, 200, 99), (253, 231, 37)]

function _color(t)
    s = clamp(t, 0, 1) * (length(_COLORS) - 1)
    k = min(floor(Int, s), length(_COLORS) - 2)
    a, b = _COLORS[k+1], _COLORS[k+2]
    c = map((u, v) -> round(Int, u + (s - k) * (v - u)), a, b)
    return "rgb($(c[1]),$(c[2]),$(c[3]))"
end

# The grid triangulations of a bivariate function, with each triangle shaded by
# the mean function value at its corners
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

# Univariate functions, with their breakpoints and, if given, the function they
# interpolate as a dashed curve
function univariate_figure(
    path;
    panels = [
        "UnivariatePWLFunction([0, 1, 2, 4], [0, 2, 1, 3])" =>
            (UnivariatePWLFunction([0, 1, 2, 4], [0, 2, 1, 3]), nothing),
        "UnivariatePWLFunction(0:1.0:2π, sin)" =>
            (UnivariatePWLFunction(0:1.0:2π, sin), sin),
    ],
)
    width, height, gap, label = 340, 200, 50, 30
    total_width = length(panels) * (width + gap) + gap
    total_height = height + label + gap
    io = IOBuffer()
    println(
        io,
        """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $total_width $total_height" width="$total_width" font-family="sans-serif" font-size="13">""",
    )
    println(
        io,
        """<rect width="$total_width" height="$total_height" fill="white"/>""",
    )
    for (k, (name, (pwl, f))) in enumerate(panels)
        x0, y0 = gap + (k - 1) * (width + gap), gap / 2
        points = sort!(
            unique(
                (only(p), only(v)) for s in segments(pwl) for
                (p, v) in zip(input_values(s), output_values(s))
            ),
        )
        x_lo, x_hi = extrema(first, points)
        samples = range(x_lo, x_hi; length = 200)
        values = [last.(points); f === nothing ? Float64[] : f.(samples)]
        y_lo, y_hi = extrema(values)
        pad = 0.1 * (y_hi - y_lo)
        y_lo, y_hi = y_lo - pad, y_hi + pad
        px(a) = x0 + width * (a - x_lo) / (x_hi - x_lo)
        py(b) = y0 + height * (1 - (b - y_lo) / (y_hi - y_lo))
        # Frame, grid lines and tick labels
        println(
            io,
            """<rect x="$x0" y="$y0" width="$width" height="$height" fill="none" stroke="#bbb"/>""",
        )
        for a in ceil(Int, x_lo):floor(Int, x_hi)
            println(
                io,
                """<line x1="$(px(a))" y1="$y0" x2="$(px(a))" y2="$(y0 + height)" stroke="#eee"/>""",
                """<text x="$(px(a))" y="$(y0 + height + 16)" text-anchor="middle" fill="#555">$a</text>""",
            )
        end
        step = y_hi - y_lo > 3 ? 1.0 : 0.5
        for b in (ceil(y_lo/step)*step):step:y_hi
            println(
                io,
                """<line x1="$x0" y1="$(py(b))" x2="$(x0 + width)" y2="$(py(b))" stroke="#eee"/>""",
                """<text x="$(x0 - 6)" y="$(py(b) + 4)" text-anchor="end" fill="#555">$(isinteger(b) ? Int(b) : b)</text>""",
            )
        end
        if f !== nothing
            curve = join(("$(px(a)),$(py(f(a)))" for a in samples), " ")
            println(
                io,
                """<polyline points="$curve" fill="none" stroke="#999" stroke-width="1.5" stroke-dasharray="5,4"/>""",
            )
        end
        # The function, drawn segment by segment, and its breakpoints
        for s in segments(pwl)
            (a,), (b,) = input_values(s)
            (fa,), (fb,) = output_values(s)
            println(
                io,
                """<line x1="$(px(a))" y1="$(py(fa))" x2="$(px(b))" y2="$(py(fb))" stroke="rgb(59,82,139)" stroke-width="2.5"/>""",
            )
        end
        for (a, b) in points
            println(
                io,
                """<circle cx="$(px(a))" cy="$(py(b))" r="4.5" fill="rgb(33,144,141)" stroke="white" stroke-width="1.5"/>""",
            )
        end
        println(
            io,
            """<text x="$(x0 + width / 2)" y="$(y0 + height + label + 8)" text-anchor="middle" fill="#333">$name</text>""",
        )
    end
    println(io, "</svg>")
    mkpath(dirname(path))
    write(path, take!(io))
    return path
end
