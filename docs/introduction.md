# Introduction

## What is RH 1.5D

RH 1.5D is a modified version of the RH radiative transfer code that
runs through a 3D/2D/1D atmosphere column-by-column, in the same fashion
of rhf1d. It was developed as a way to efficiently run RH
column-by-column (1.5D) over large atmospheres, simplifying the output
and being able to run in supercomputers. It is MPI-parallel and scales
well to at least 10,000 processes.

While initially developed as another geometry on the RH tree, the
requirements of the parallel version required changes to the RH core
tree. Thus, it is more than a [wrapper]{.title-ref} over RH and is
distributed with a modified version of RH (see below for differences
from RH distributed by Han Uitenbroek).

## Acknowledging RH 1.5D

If you use RH 1.5D for your work, we would appreciate if you would
acknowledge it appropriately in a publication, presentation, poster, or
talk. For a publication, this is best done by citing the RH 1.5D
([Pereira & Uitenbroek
2015](http://adsabs.harvard.edu/abs/2015A%26A...574A...3P)) and RH
([Uitenbroek 2001](http://adsabs.harvard.edu/abs/2001ApJ...557..389U))
papers. In addition, if the journal allows it please include a link to
its [Github repository](https://github.com/ITA-Solar/rh).

## What this manual covers

Because there is much in common with RH, this manual should be seen as
an incremental documentation of the 1.5D parallel side. This manual
focuses on what is different from RH. Users should refer to the RH
documentation by Han Uitenbroek for more detailed information on RH.

## Comparison with RH

RH 1.5D inherits most of the code base from RH, but some features are
new. The code is organised in the same was as the RH source tree, with
the routines specific to the 1.5D version residing on a subdirectory
`rh15d` of the `rh` source tree. In this way, it works similarly to the
subdirectories in `rh` for different geometries. The compilation and
linking proceeds as for the other geometries: first the general
`librh.a` library should be compiled, and then the code in `rh15d` will
be compiled and linked to it. The run directory is very similar to that
of a given geometry in RH: most of the `*.input` files are used in the
same way.

The lists below show a comparison between RH 1.5D and RH for a 1D
plane-parallel geometry:

### Commonalities between RH 1.5D and RH

-   Core RH library
-   Structure and location of `*.input` files (some new options
    available)
-   Wavetable, Atom, Molecule, line list, and any other input files
    except the atmospheres
-   Directory-level compilation and general structure of run directory
    (new subdirectories needed)

### What is new in RH 1.5D

-   MPI-parallelism with dynamic load balancing and efficient I/O
-   Different format of atmosphere files
-   Different formats of output files
-   Select which quantities should be in output
-   New options for `keyword.input` and `atoms.input`
-   Hybrid angle-dependent PRD mode
-   PRD-switching
-   Option for using escape probability approximation as initial
    solution
-   Exclude from the calculations the higher parts of the atmosphere,
    above a user-defined temperature threshold
-   Depth scale optimisation
-   Option for cubic Hermite interpolation of source function in formal
    solver
-   Option for cubic Bézier and DELO Bézier interpolation of source
    function in formal solvers, both for polarised and unpolarised light
-   Support for more types of collisional excitations
-   Easy re-run of non-converged columns
-   Option for keeping background opacities in memory and not in disk
-   New analysis suite in Python

### What is not supported in RH 1.5D

-   Currently only a fraction of the RH output is written to disk (to
    save space), but more output can be added
-   The old IDL analysis suite does not currently support the new output
    format
-   `solveray` is no longer used
-   `backgrcontr` no longer works with the new output
-   Any other geometry aside from 1.5D
-   Continuing old run by reading populations and output files
-   Full Stokes NLTE iterations and background polarisation (might work
    with little effort, but has not been tested)
-   Thread parallelisation
