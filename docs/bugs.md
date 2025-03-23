# Known bugs and limitations

RH 1.5D is always evolving, and there are likely to be bugs and
limitations. If you find bugs, please file them as [issues on github](https://github.com/ITA-Solar/rh/issues). They will be dealt with as time permits.

## Current issues

-   Check the github [RH issues
    page](https://github.com/ITA-Solar/rh/issues) for an updated list.
-   If the `scratch` or `output` directories are not present, the code
    will crash. The error message is not very clear.
-   In `keyword.input`, if one sets a `SNAPSHOT` value to be more than
    what is in the atmosphere file, the code will stop with an error
    message: `Index exceeds dimension bound`. This error should be made
    more clear.
-   The atom files must not end with a blank line, otherwise `gencol`
    will fail and the program stops.
-   Line buffered or full buffered log options still require the user to
    change the source code.
-   Depth refinement fails in some cases due to problems caused by cubic
    interpolation artefacts.
-   Using more than 4000 cores and writing full output may cause I/O
    slowdowns and Lustre contention in some systems.

## Planned features

-   Support for multiple snapshots in the output files.
-   `pool` mode be more flexible, with the possibility of
    several `overlord` nodes, useful for running with more
    than 4000 processes.
-   More flexible control of what output is written.
