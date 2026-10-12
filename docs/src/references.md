# References

- E. Balas. Disjunctive programming: Properties of the convex hull of feasible
  points. _Discrete Applied Mathematics_ 89(1–3):3–44, 1998.
  [doi:10.1016/S0166-218X(98)00136-X](https://doi.org/10.1016/S0166-218X(98)00136-X).
  The formulation of a union of polyhedra used by [`ConvexPieces`](@ref).
- C. D'Ambrosio, A. Lodi and S. Martello. Piecewise linear approximation of
  functions of two variables in MILP models. _Operations Research Letters_
  38(1):39–46, 2010.
  [doi:10.1016/j.orl.2009.09.005](https://doi.org/10.1016/j.orl.2009.09.005).
  Compares bivariate formulations with the
  [separable reformulation](@ref "Separable reformulation of products").
- J. Huchette and J. P. Vielma. Nonconvex piecewise linear functions: Advanced
  formulations and simple modeling tools. _Operations Research_, 2023.
  [doi:10.1287/opre.2019.1973](https://doi.org/10.1287/opre.2019.1973),
  [arXiv:1708.00050](https://arxiv.org/abs/1708.00050).
  The zig-zag and grid formulations ([`ZigZagBinary`](@ref),
  [`ZigZagInteger`](@ref), [`SixStencil`](@ref), [`NineStencil`](@ref),
  [`OptimalTriangleSelection`](@ref), [`K1`](@ref)) and a computational
  comparison.
- J. Huchette and J. P. Vielma. A combinatorial approach for small and strong
  formulations of disjunctive constraints. _Mathematics of Operations
  Research_ 44(3):793–820, 2019.
  [arXiv:1607.04803](https://arxiv.org/abs/1607.04803).
  Independent branching formulations, including
  [`OptimalIndependentBranching`](@ref).
- B. Lyu, I. V. Hicks and J. Huchette. Building formulations for piecewise
  linear relaxations of nonlinear functions. _Operations Research_, 2025.
  [doi:10.1287/opre.2023.0187](https://doi.org/10.1287/opre.2023.0187),
  [arXiv:2304.14542](https://arxiv.org/abs/2304.14542).
  The Gray code and biclique cover formulations of
  [`IndependentBranching`](@ref).
- J. P. Vielma. Embedding formulations and complexity for unions of polyhedra.
  _Management Science_ 64(10):4721–4734, 2018.
  [arXiv:1506.01417](https://arxiv.org/abs/1506.01417).
  The logarithmic embedding formulation [`LogarithmicEmbedding`](@ref).
- J. P. Vielma, S. Ahmed and G. Nemhauser. Mixed-integer models for
  nonseparable piecewise-linear optimization: Unifying framework and
  extensions. _Operations Research_ 58(2):303–315, 2010.
  The convex combination, multiple choice, disaggregated logarithmic and
  incremental formulations ([`ConvexCombination`](@ref),
  [`MultipleChoice`](@ref), [`DisaggregatedLogarithmic`](@ref),
  [`Incremental`](@ref)).
- J. P. Vielma and G. Nemhauser. Modeling disjunctive constraints with a
  logarithmic number of binary variables and constraints. _Mathematical
  Programming_ 128(1–2):49–72, 2011.
  The logarithmic independent branching formulation
  [`LogarithmicIndependentBranching`](@ref) and the [`UnionJack`](@ref)
  formulation.
