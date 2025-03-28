# Running the code

## Binaries and execution {#binaries-label}

Compilation should produce three executables: `rh15d_ray_pool`,
`rh15d_ray`, and `rh15d_lteray`. The latter is a special case for
running only in LTE. The other two are the main programs. They represent
two different modes of job distribution: *normal* and
*pool*.

In the *pool* mode there is a process that works as
*overlord*: its function is to distribute the work to other processes.
The other processes (*drones*) ask the
*overlord* for a work unit. When they finish their unit,
they go back and ask for more, until all tasks are completed. Because of
race conditions and because different columns will run at different
speeds, it is not possible to know which columns a given process will
run beforehand. Due to the *overlord*, `rh15d_ray_pool`
needs to run with two or more processes. The advantage of the
*pool* mode is that the dynamic load allocation ensures the
most efficient use of the resources. With the normal mode it may happen
that some processors will work on columns that take longer to converge
(especially as they are adjacent), and in the end the execution will
have to wait for the process that takes longer. In some cases
(especially with PRD) the *pool* mode can be 2-3 times
faster than the *normal* mode. When one runs with a large
number of processes (> 2000) and each column takes little time to
calculate, the *pool* mode can suffer from communication
bottlenecks and may be slower because a single *overlord*
cannot distribute the tasks fast enough. The only disadvantage of the
pool mode (so far) is that not all output is currently supported with
this mode.

!!! warning
    
    The *normal* mode is deprecated and will be removed in a
    later revision. Use only for single processor runs or if you know what
    you're doing!


In the *normal* mode the jobs (all the atmosphere columns
for which one wants to calculate) are divided by the number of processes
at the start of the execution. There is no communication between
processes, and each process knows from the start all the columns it is
going to run. These columns are adjacent. If the number of columns is
not a multiple of the number of processes, there will be some processes
with larger workloads. There is no minimum number of processes to run,
and `rh15d_ray` can also be run in a single process. Regions of an
atmosphere can take a lot longer to run than others, and the processes
that work on those will take longer to finish. In the
*normal* mode this means that the slowest process will set
the overall running time, and therefore in practice it can take more
than 10x longer than the *pool* mode (and is therefore not
recommended).

As an MPI program, the binaries should be launched with the appropriate
command. Some examples:

    mpirun -np N ./rh15d_ray_pool
    srun ./rh15d_ray_pool       # Use in Betzy or other system with SLURM
    mpiexec ./rh15d_ray_pool    # use in Pleiades
    aprun -B ./rh15d_ray        # use in Hexagon or other Cray

## The run directory

!!! warning

    Before running, make sure you have the sub-directories `scratch` and
`output` in the run directory.

The run directory contains the configuration files, the binaries, and
the `scratch` and `output` directories. As the names imply, temporary
files will be placed under `scratch` and the final output files in
`output`. No files under `scratch` will be used after the run is
finished (they are not read for re-runs).

The `scratch` directory contains different types of files. Most of them
are binary files write by RH 1.5D to save memory. Example files are the
`background_p*.dat` with the background opacities, files with PRD
weights, and the `rh_p*.log` log files. Each process creates one of
those files, and they will have the suffix `_pN.*`, where `N` is the
process number. The log files have the same format as in RH. The writing
of each process's log file is buffered by line. Because these are
updated often, when running with many processes this can be a drag on
some systems. Therefore, it is possible to run full buffering (meaning
log files are only written when the main program finishes). This option
is not exposed in the configuration files, so one needs to change the
file `parallel.c` in the following part:

``` c
/* _IOFBF for full buffering, _IOLBF for line buffering */
setvbuf(mpi.logfile, NULL, _IOLBF, BUFSIZ_MPILOG);
```

One should replace `_IOLBF` by `_IOFBF` to change from line buffering to
full buffering.

The `output` directory will contain the three output files:
`output_aux.hdf5`, `output_indata.hdf5`, and `output_ray.hdf5`. See
[output](output.md) for more details on the
structure of these files. If doing a re-run, these files must already
exist; they will be updated with the new results. Otherwise, if these
files are already in `output` before the execution, they will be
overwritten. At the start of the execution, the output files are written
with a special a fill value. This means that the disk space for the full
output must be available at the start of the run, and no CPU time will
be wasted if at the end of the run there is not enough disk space. The
files are usually written every time a process finishes work on a given
column. The fill value arrays are overwritten with the data. One
advantage of this method is that even if the system crashes or the
program stops, it is possible to recover the results already written
(and a re-run can be performed for just the missing columns).

All the processes write asynchronously to all the output files. In some
cases this can cause contention in the filesystem, with many processes
trying to access the same data at the same time. In the worst case
scenario, the contention can create bottlenecks which practically stop
the execution. Therefore, it is highly recommended that the users tune
their filesystem for the typical loads of RH. Many supercomputers make
use of Lustre, a parallel filesystem. With Lustre, resources such as
files can be divided in different stripes that can be placed in several
different machines (OSTs). For running RH with more than 500 processes,
one should use as many OSTs as available in the system, and select the
lustre stripe size to the typical amount of data written to a file per
simulation column. The stripe can set with the `lfs
setstripe` command:

    lfs setstripe -s stripe_size -c stripe_count -o stripe_offset directory|filename

It can be run per file (e.g. `output_ray.hdf5`), or for the whole
`output` directory. Using a stripe count of `-1` will ensure that the
maximum number of OSTs is used. For the typical files RH 1.5D produces,
it is usually ok to apply the same Lustre settings to the whole `output`
directory, and the following settings seem to reasonable (tested on
Betzy and Pleiades):

    lfs setstripe -S 8M -c -1 output/

Similarly, the directory where you will store the input atmospheres, and
the `scratch` directory can also benefit from Lustre striping. For
`scratch`, because most files there are small, it is recommended to use
a stripe count of 1.

!!! info

    You need a parallel filesystem if you are running RH 1.5D across more
    than one host. Most supercomputers have parallel filesystems, but if you
    are running in a smaller cluster this may not be the case. RH 1.5D will
    always run, but the HDF5 writes will not work and the results will be
    unreadable. NFS is **not** a parallel file system.


## Logs and messages

In addition to the logs per process saved to `scratch`, a much smaller
log will be printed in `stdout`. This log is a smaller summary of what
each process is doing. Here is an example of typical messages:

    Process   1: --- START task   1, (xi,yi) = (  0,156)
    Process 232: --- START task   1, (xi,yi) = (  0,159)
    Process  36: --- START task   1, (xi,yi) = (  0,162)
    Process  12: --- START task   1, (xi,yi) = (  0,171)
    (...)
    Process  12: *** END   task   1 iter, iterations = 121, CONVERGED
    Process   3: *** END   task   1 iter, iterations = 200, NO Convergence
    Process   4: *** SKIP  task   1 (crashed after 81 iterations)
    Process   3: --- START task   2, (xi,yi) = ( 23, 64)
    Process  12: --- START task   2, (xi,yi) = ( 23, 65)
    Process   4: --- START task   2, (xi,yi) = ( 23, 65)
    (...)
    *** Job ending. Total 262144 1-D columns: 262142 converged, 1 did not converge, 1 crashed.
    *** RH finished gracefully.

In this example one can see the three possible outputs for a
single-column calculation: convergence, non-convergence (meaning the
target `ITER_LIMIT` was not met in `N_MAX_ITER` iterations), or a crash
(many reasons). If there are singular matrices or other causes for a
column to crash, RH 1.5D will skip that column and proceed to the next
work unit. Such cases can be re-run with different parameters. In some
cases (e.g. inexistent files) it is not possible to prevent a crash, and
RH 1.5D will finish non-gracefully.

## Reruns and lack of convergence

Dynamic atmospheres often have large gradients in temperature, density,
velocity, etc. that cause problems when solving the non-LTE problem.
This may lead to some (or all) columns not converging, and is dependent
on the input atmosphere, model atoms, and run options. In RH terms,
non-converged or "crashed", represent the same problem. In some cases,
the iterations diverge strongly (crash), while in others they fail to
reach the target limit for convergence in the allocated maximum number
of iterations (non-convergence).

The output for atmosphere columns that did not converged or crashed is
not saved to disk (a fill value is used instead). The recommended
procedure in these cases is to rerun RH with different input options
(e.g. less agressive acceleration, a larger number of maximum
iterations). A special rerun mode is available to save time calculating
again the columns that have already converged. When `15D_RERUN = TRUE`
in `keyword.input`, RH will read the output and run again only for the
columns that did not converge.

The rerun mode requires all the previous output to be under `output/`.
The user has the possibility of changing options in `keyword.input` for
the rerun. Not all options can be changed, because this could lead to
non-sensical results (e.g. changing the input atmosphere or atom files).
Only the following options can be changed in `keyword.input`:

|  |
:-----------------------|
|  `15D_RERUN`
|  `15D_DEPTH_CUT`
|  `15D_TMAX_CUT`
|  `15D_DEPTH_REFINE`
|  `N_PESC_ITER`
|  `COLLRAD_SWITCH`
|  `COLLRAD_SWITCH_INIT`
|  `PRD_SWITCH`
|  `NRAYS`
|  `N_MAX_SCATTER`
|  `N_MAX_ITER`
|  `ITER_LIMIT`
|  `NG_DELAY`
|  `NG_ORDER`
|  `NG_PERIOD`
|  `S_INTERPOLATION`
|  `PRD_N_MAX_ITER`
|  `PRD_ITER_LIMIT`
|  `B_STRENGTH_CHAR`


All the other options are locked to the values used for the firs run.
Likewise, it is not possible to change `atoms.input`, the atom files, or
the line list files. When RH 1.5D is first run, nearly all input options
(except the input atmosphere and molecule files) are saved into the
output, and in a rerun these are read from the output and not from the
original files. This also means that a rerun can be performed even if
the original files are no longer available.

## Helper script

There is a Python script called `runtools.py` designed to make it easier
to run RH 1.5D for large projects. It resides in
`rh/python/runtools.py`. It requires [Python](http://www.python.org/)
with the [numpy](http://www.numpy.org/) and [h5py](http://www.h5py.org/)
(or [netCDF4](http://code.google.com/p/netcdf4-python/)) modules. It was
made to run a given RH 1.5D setup over many simulation snapshots,
spanning several atmosphere files. It supports a progressive rerun of a
given problem, and allows the of use different `keyword.input`
parameters for different columns, tackling columns harder to converge.

The first part of `runtools.py` should be modified for a users's need.
It typically contains:

``` python
atmos_dir = '/mydata_dir'
seq_file = 'RH_SEQUENCE'
outsuff = 'output/output_ray_mysim_CaII_PRD_s%03i.hdf5'
mpicmd = 'mpiexec'
bin = './rh15d_ray_pool'
defkey = 'keyword.save'
log = 'rh_running.log'
tm = 40
rerun = True
rerun_opt = [ {'NG_DELAY': 60, 'NG_PERIOD': 40, '15D_DEPTH_REFINE': 'FALSE',
               '15D_ZCUT': 'TRUE', 'N_MAX_ITER': 250, 'PRD_SWITCH': 0.002 },
              {'NG_DELAY': 120, 'NG_PERIOD': 100, '15D_DEPTH_REFINE': 'TRUE',
               'PRD_SWITCH': 0.001 } ]
```

The different options are:



  Name   |       Type  |   Description|
|:-----------|:-------|:---------------------------------------------------
  `atmos_dir` |  string  | Directory where the atmosphere files are kept.
  `seq_file`  |  string  | Location of sequence file. This file contains the names of the atmosphere files to be used (one file per line). The script will then run RH 1.5D for every snapshot in every file listed.
  `outsuff`   |  string  | Template to write the `output_ray.ncdf` files. The `%03i` format will be replaced with the snapshot number.
  `mpicmd`    |  string  | System-dependent command to launch MPI. The script knows that for aprun the -B option should be used. This option also activates system specific routines (e.g. how to kill the run in pleiades).
  `bin`       |  string |  RH 1.5D binary to use.
  `defkey`    |  string |  Default template for `keyword.input`. Because the of the rerun options, `keyword.input` is overwritten for every rerun. This file is used as a template it (i.e., most of its options will be unchanged, unless specified in `rerun_opt`).
  `log`       |  string |  File where to save the main log. Will be overwritten for each new snapshot.
  `tm`        |  int    |  Timeout (in minutes) to kill execution of code, if there is no message written to main log. Used to prevent code from hanging if there are system issues. After killed, program is relaunched. If `tm = 0`, program will never be killed.
  `rerun`     |  bool   |  If `True`, will re-run the program (with different settings) to achieve convergence if any columns failed. Number of reruns is given by size of `rerun_opt`.
  `rerun_opt`  | list   |  Options for re-run. This is a list made of dictionaries. Each dictionary contains the keywords to update `keyword.input`. Only the keywords that differ from the `defkey` file are necessary.


!!! info

    Only the first line of the sequence time is read at a time. The script
    reads the first line, deletes it from the file, and closes the file. It
    then reads the first line again and continues running, until there are
    no more lines in the file. This behaviour enables the file to be worked
    by multiple scripts at the same time, and allows one to dynamically
    change the task list at any time of the run.


!!! info

    The script also includes a tenaciously persistent wait and relaunch
    feature designed to avoid corruption if there are system crashes or
    problems. Besides the `tm` timeout, if there is any problem with the
    execution, the code will wait for some periods and try and relaunch the
    code. For example, if one of the atmosphere files does not exist,
    `runtools.py` will try three times and then proceed to the next file.

