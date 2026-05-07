# RH 1.5D — Parallel I/O knobs

Quick reference for the runtime-toggleable I/O optimisations added on top of
baseline commit `580f9118`. Use these to A/B test the impact of the recent
`[OPT]` series (pool-mode collective writes, node-shared atmosphere cache,
Lustre MPI-IO tuning) on performance and filesystem load.

All keywords below go into `keyword.input`. All env vars are read once in
`initParallel` and logged by rank 0 at startup.

## Keywords (keyword.input)

| Keyword | Default | Effect when `TRUE` | Effect when `FALSE` | Source |
|---|---|---|---|---|
| `15D_POOL_COLLECTIVE_WRITE` | `TRUE` | Phase-2 pool flush uses `H5FD_MPIO_COLLECTIVE` | Uses `H5FD_MPIO_INDEPENDENT` (each rank writes its own hyperslabs) | c6f3782 |
| `15D_USE_NODE_ATMOS_CACHE` | `FALSE` | `readAtmos_hdf5` serves from the node-shared cache window (one file read per node) | Every rank reads its columns from the atmosphere file directly | 3347388 … ede246c |
| `15D_ATMOS_CACHE_CYCLIC` | `FALSE` | Cyclic row decomposition across nodes (planned, not yet honored in v1) | Contiguous row decomposition | — |
| `15D_MPIIO_LUSTRE_HINTS` | `TRUE` | `create_hdf5_fapl` passes `mpi.info` (romio + stripe + cb_nodes) to `H5Pset_fapl_mpio` **and** enables collective HDF5 metadata ops on output | Passes `MPI_INFO_NULL`, no collective metadata ops (ROMIO defaults) | 42fc036, e1767d1 |
| `15D_FLUSH_INTERVAL` | `16` | Columns buffered per rank before phase-2 flush | same | 77fd7bd |

### Baseline reproduction

To approximate pre-`580f9118` I/O behavior for the runtime-toggleable OPTs,
set all four toggles off:

```
15D_POOL_COLLECTIVE_WRITE = FALSE
15D_USE_NODE_ATMOS_CACHE  = FALSE
15D_MPIIO_LUSTRE_HINTS    = FALSE
15D_FLUSH_INTERVAL        = 16
```

**Caveats — not runtime-toggleable.** The following OPTs are baked into the
code and would need an actual code rollback to isolate:

- `readAtmos_hdf5` per-column `H5Dread_multi` consolidation (1dc2f25, cef4ae0)
- pool-flush rate-dataset consolidation and removed barrier (77fd7bd)
- `MPI_Ssend` → `MPI_Send` in pool overlord dispatch (05f2c34)

## Environment variables (MPI-IO Lustre hints)

Read by `initParallel` in `rh15d/parallel.c`. Only meaningful when
`15D_MPIIO_LUSTRE_HINTS = TRUE`. Defaults match the reference production
job script (`lfs setstripe -c 12 -S 1M output/`).

| Env var | Default | MPI-IO hint | Notes |
|---|---|---|---|
| `RH_LUSTRE_STRIPE_COUNT` | `12` | `striping_factor` | Should match `lfs getstripe -c <outdir>` |
| `RH_LUSTRE_STRIPE_SIZE` | `1048576` (1 MB) | `striping_unit` | Bytes; match `lfs getstripe -S <outdir>` |
| `RH_MPIIO_CB_BUFFER_SIZE` | `16777216` (16 MB) | `cb_buffer_size` | Collective-buffer size per aggregator |
| `RH_MPIIO_CB_NODES` | `mpi.n_nodes` | `cb_nodes` | Number of aggregators; default = 1 per node |

Auto-track the real Lustre layout from the job script:

```bash
export RH_LUSTRE_STRIPE_COUNT=$(lfs getstripe -c output/)
export RH_LUSTRE_STRIPE_SIZE=$(lfs getstripe  -S output/)
srun ... rh15d_ray_pool
```

Unconditional MPI-IO hints (not env-overridable here): `romio_cb_write=enable`,
`romio_ds_write=disable`, `romio_cb_read=enable`, `romio_ds_read=disable`,
`cb_config_list=*:1`. If you need to experiment with these, set
`ROMIO_HINTS=/path/to/hints.file` — ROMIO only picks up hints from that file
for keys the application has not already set, so you'll need to remove the
corresponding `MPI_Info_set` line in `parallel.c` first.

Unconditional HDF5 FAPL tuning (always applied, not gated): 1 MB alignment,
8 MB metadata block size, file locking disabled. These are pure
HDF5-level tuning / correctness on Lustre, not performance knobs worth
toggling.

## Related files

- `rh15d/parallel.c` — `initParallel`, `create_hdf5_fapl`, `create_hdf5_fapl_indep`
- `rh15d/pool_buffer.c` — `writeCollective_pool` (phase-2 flush)
- `rh15d/hdf5atmos.c` — `readAtmos_hdf5` (node-cache read path)
- `readinput.c` — keyword registration
- `inputs.h` — `InputData` struct fields
