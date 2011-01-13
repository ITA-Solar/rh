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

/* Definitions for the Aux file */
#define ARR_STRLEN  30
#define AUX_FILE    "output_aux.ncdf"
#define POP_NAME    "populations"
#define POPLTE_NAME "populations_LTE"
#define RIJ_L_NAME  "Rij_line"
#define RJI_L_NAME  "Rji_line"
#define RIJ_C_NAME  "Rij_continuum"
#define RJI_C_NAME  "Rji_continuum"
#define COLL_NAME   "collision_rates"
#define DAMP_NAME   "damping"
#define VBROAD_NAME "broadening_velocity"
#define NW_AD_NAME  "nwave_angle_dep"
#define NW_AI_NAME  "nwave_angle_ind"
#define CHI_AI_NAME "chi_angle_ind"
#define CHI_AD_NAME "chi_angle_dep"
#define ETA_AI_NAME "eta_angle_ind"
#define ETA_AD_NAME "eta_angle_dep"
#define WAVET_NAME  "wavelength_nm"
#define WAVE_AD_IDX_NAME "wave_angle_dep_indices"
#define WAVE_AI_IDX_NAME "wave_angle_ind_indices"


#define STOPFILE_TEMPLATE "scratch/STOP_FILE_p%d"


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
      *aux_atom_RjiC,    *aux_atom_coll,    *aux_atom_damp,    *aux_atom_vbroad,
       aux_op_chi_ai,     aux_op_chi_ad,     aux_op_eta_ai,     aux_op_eta_ad;
  /* for atom file positions */
  long *atom_file_pos;
} IO_data;



#endif /* !__IO_H__ */
              
/* ---------------------------------------- io.h -------------------- */
