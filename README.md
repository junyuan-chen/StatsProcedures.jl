# StatsProcedures.jl

*An interface framework for sharing intermediate steps across statistical methods*

[![CI-stable][CI-stable-img]][CI-stable-url]
[![codecov][codecov-img]][codecov-url]
[![PkgEval][pkgeval-img]][pkgeval-url]

[CI-stable-img]: https://github.com/junyuan-chen/StatsProcedures.jl/workflows/CI-stable/badge.svg
[CI-stable-url]: https://github.com/junyuan-chen/StatsProcedures.jl/actions?query=workflow%3ACI-stable

[codecov-img]: https://codecov.io/gh/junyuan-chen/StatsProcedures.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/junyuan-chen/StatsProcedures.jl

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/StatsProcedures.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/S/StatsProcedures.html

[StatsProcedures.jl](https://github.com/junyuan-chen/StatsProcedures.jl)
is a Julia package that provides interfaces for generic statistical methods
that allow sharing results from intermediate steps
across multiple model specifications that may require different methods.
It is originally a component of
[DiffinDiffsBase.jl](https://github.com/JuliaDiffinDiffs/DiffinDiffsBase.jl)
for reducing unnecessary repetitions of identical intermediate steps
when estimating multiple specifications.
It is not intended to be used as a standalone package.
