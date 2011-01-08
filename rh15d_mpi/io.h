/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 3-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Wed Dec 22 16:44:00 2010 --

       --------------------------                      ----------RH-- */


#ifndef __IO_H__
#define __IO_H__

#include <netcdf.h>



/* Definitions for brs_p */
#define  BRS_DOT_OUT  "brs.out"
#define  BRS_FILE     "scratch/brs_out_p"
#define  BRS_EXT      ".ncdf"
#define  HASLINE_VAR  "hasline"
#define  ISPOL_VAR    "ispolarised"
#define  BGREC_VAR    "backgrrecno"

/* Definitions for readj_p */
#define  J_FILE_TEMPLATE "%s.ncdf"

/* Definitions for readatmos_ncdf */
#define  TEMP_NAME "temperature"
#define  VZ_NAME   "velocity_z"
#define  NE_NAME   "electron_density"
#define  NH_NAME   "hydrogen_populations"
#define  Z_NAME    "z"

/* Definitions for background_p */
#define  FILE_EXT ".dat"



/* For keeping the netCDF file and variable IDs */
typedef struct {
  /* for the BRS file */
  int  brs_ncid,          brs_hl_var,        brs_ip_var,         brs_nrec_var;
  /* for the J file */
  int  j_ncid,            j_jlambda_var,     j_j20_var; 
  /* for the spectrum file*/ 
  int  spec_ncid,         spec_int_var,      spec_flux_var,     spec_wave_var, 
       spec_stokes_u_var, spec_stokes_q_var, spec_stokes_v_var; 
  /* for the input data file. Note: this netCDF file has several groups */
  int  in_ncid,           in_input_ncid,     in_atmos_ncid,     in_mpi_ncid;
  int  in_atmos_T,        in_atmos_ne,       in_atmos_vz,       in_atmos_vt,
       in_atmos_B,        in_atmos_gB,       in_atmos_chiB,     in_atmos_nh,
       in_atmos_ew,       in_atmos_ab,       in_atmos_eid,      in_atmos_mu,
       in_atmos_wmu,      in_atmos_z,        in_atmos_x,        in_atmos_y;
  int  in_mpi_xnum,       in_mpi_ynum,       in_mpi_tm,         in_mpi_tn,
       in_mpi_it,         in_mpi_conv,       in_mpi_dm,         in_mpi_ntsk,
       in_mpi_host,       in_mpi_ft;
  /* for the aux file */
  int  aux_ncid,         *aux_atom_ncid,     aux_op_ncid,      *aux_atom_pop,
      *aux_atom_poplte,  *aux_atom_RijL,    *aux_atom_RjiL,    *aux_atom_RijC,
      *aux_atom_RjiC,    *aux_atom_coll,    *aux_atom_damp,     aux_op_chi_ai,
       aux_op_chi_ad,     aux_op_eta_ai,     aux_op_eta_ad;
  /* for atom file positions */
  long *atom_file_pos;
} IO_data;



#endif /* !__IO_H__ */
              
/* ---------------------------------------- io.h -------------------- */
