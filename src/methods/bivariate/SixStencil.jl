# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    SixStencil(axis_method = Logarithmic())

Formulation of bivariate continuous piecewise linear functions on any grid
triangulation (Huchette and Vielma, 2023).

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with 6 binary variables to select the triangle in the chosen cell. Fewer are
needed on small grids, or when the diagonals of all cells have the same
direction, as levels that separate nothing are left out.
"""
struct SixStencil{M<:PWLMethod} <: TriangleSelectionMethod
    axis_method::M
end
SixStencil() = SixStencil(Logarithmic())

function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::SixStencil,
)
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)
    cover = _BicliqueCover()

    # Diagonal lines of cells running from SW to NE, grouped with an offset of
    # 3 into three levels. Along each line, the ends of the diagonals that are
    # not edges alternate between the two sides of the level.
    for o in 0:2
        A_o = Set{Tuple{Int,Int}}()
        B_o = Set{Tuple{Int,Int}}()
        for off_1 in o:3:(n_1-2)
            sw_in_A = true
            for i in (1+off_1):(n_1-1)
                j = i - off_1
                if !(1 ≤ j ≤ n_2 - 1)
                    continue
                end
                if triangle_direction[i, j]
                    if sw_in_A
                        push!(A_o, (i, j))
                        push!(B_o, (i + 1, j + 1))
                    else
                        push!(A_o, (i + 1, j + 1))
                        push!(B_o, (i, j))
                    end
                    sw_in_A = !sw_in_A
                end
            end
        end
        for off_2 in (3-o):3:(n_2-1)
            sw_in_A = true
            for j in (off_2+1):(n_2-1)
                i = j - off_2
                if !(1 ≤ i ≤ n_1 - 1)
                    continue
                end
                if triangle_direction[i, j]
                    if sw_in_A
                        push!(A_o, (i, j))
                        push!(B_o, (i + 1, j + 1))
                    else
                        push!(A_o, (i + 1, j + 1))
                        push!(B_o, (i, j))
                    end
                    sw_in_A = !sw_in_A
                end
            end
        end
        push!(cover, (sort!(collect(A_o)), sort!(collect(B_o))))
    end

    # The same for the lines running from SE to NW
    for o in 0:2
        A_o = Set{Tuple{Int,Int}}()
        B_o = Set{Tuple{Int,Int}}()
        for off_1 in o:3:(n_1-2)
            se_in_A = true
            for j in 1:(n_2-1)
                i = n_1 - j - off_1
                if !(1 ≤ i ≤ n_1 - 1)
                    continue
                end
                if !triangle_direction[i, j]
                    if se_in_A
                        push!(A_o, (i + 1, j))
                        push!(B_o, (i, j + 1))
                    else
                        push!(A_o, (i, j + 1))
                        push!(B_o, (i + 1, j))
                    end
                    se_in_A = !se_in_A
                end
            end
        end
        for off_2 in (3-o):3:(n_2-1)
            se_in_A = true
            for j in (off_2+1):(n_2-1)
                i = n_1 - j + off_2
                if !(1 ≤ i ≤ n_1 - 1)
                    continue
                end
                if !triangle_direction[i, j]
                    if se_in_A
                        push!(A_o, (i + 1, j))
                        push!(B_o, (i, j + 1))
                    else
                        push!(A_o, (i, j + 1))
                        push!(B_o, (i + 1, j))
                    end
                    se_in_A = !se_in_A
                end
            end
        end
        push!(cover, (sort!(collect(A_o)), sort!(collect(B_o))))
    end

    add_biclique_cover!(formulation, λ, cover)
    return nothing
end
