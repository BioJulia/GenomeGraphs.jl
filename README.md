# BioSequenceGraphs

| **Release**                                                     | **Documentation**                                                               | **Maintainers**                             |
|:---------------------------------------------------------------:|:-------------------------------------------------------------------------------:|:-------------------------------------------:|
| [![][release-img]][release-url] [![][license-img]][license-url] | [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | ![][maintainer-img] |

## Description

BioSequenceGraphs.jl provides a representation of sequence graphs for representing
genome assemblies and population graphs of genotypes/haplotypes and variation.

## Installation

BioSequenceGraphs is currently in pre-release development.
But you can clone BioSequenceGraphs from the Julia REPL:

```julia
julia> Pkg.clone("https://github.com/BioJulia/BioSequenceGraphs.jl.git")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

## Testing

BioSequenceGraphs.jl is tested against Julia `0.6` and current `0.7-dev` on Linux, OS X, and Windows.

| **PackageEvaluator**                                            | **Latest Build Status**                                                                                |
|:---------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------:|
| [![][pkg-0.6-img]][pkg-0.6-url] [![][pkg-0.7-img]][pkg-0.7-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url]        |


## Contributing and Questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.
Please go to the [Contributing Guidelines](https://biojulia.net/Contributing)
for more information.

If you have a question about
contributing or using this package, you are encouraged to use the
[Bio category of the Julia discourse
site](https://discourse.julialang.org/c/domain/bio).

[release-img]: https://img.shields.io/github/release/BioJulia/BioSequenceGraphs.jl.svg
[release-url]: https://github.com/BioJulia/BioSequenceGraphs.jl/releases/latest

[license-img]: https://img.shields.io/badge/license-MIT-green.svg
[license-url]: https://github.com/BioJulia/BioSequenceGraphs.jl/blob/master/LICENSE

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://biojulia.github.io/BioSequenceGraphs.jl/latest
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://biojulia.github.io/BioSequenceGraphs.jl/stable

[maintainer-img]: https://img.shields.io/badge/BioJulia%20Maintainer-Ward9250-orange.svg

[pkg-0.6-img]: https://pkg.julialang.org/badges/BioSequenceGraphs_0.6.svg
[pkg-0.6-url]: https://pkg.julialang.org/detail/BioSequenceGraphs
[pkg-0.7-img]: https://pkg.julialang.org/badges/BioSequenceGraphs_0.7.svg
[pkg-0.7-url]: https://pkg.julialang.org/detail/BioSequenceGraphs

[travis-img]: https://travis-ci.org/BioJulia/BioSequenceGraphs.jl.svg?branch=master
[travis-url]: https://travis-ci.org/BioJulia/BioSequenceGraphs.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/m4o4t50l3ih35jrt/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/Ward9250/biosequencegraphs-jl/branch/master

[codecov-img]: https://codecov.io/gh/BioJulia/BioSequenceGraphs.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/BioJulia/BioSequenceGraphs.jl
