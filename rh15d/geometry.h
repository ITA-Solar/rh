/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Tue Nov 23 22:09:58 2010 --
                       Graham S. Kerr: Fri Sept 3rd 2021.
                       Added IRRADIATED_INTP keyword. This will allow
                       the user to provide a grid of emissivities as 
                       functions of temperature, electron density, and
                       wavelength that is interpolated in each grid 
                       cell, summed to provide intensity, and injected
                       at the top boundary. 

       --------------------------                      ----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <hdf5.h>
#include <hdf5_hl.h>

/* --- Define geometric quantities for 1-D plane-parallel version --  */
enum boundcond  {ZERO, THERMALIZED, IRRADIATED, REFLECTIVE, IRRADIATED_INTP};
enum mass_scale {GEOMETRIC, COLUMN_MASS, TAU500};
enum atmos_type {HDF5, MULTI};
enum vertical   {TOP=0, BOTTOM};

typedef struct {
  enum     mass_scale  scale;
  enum     boundcond vboundary[2];
  enum     atmos_type atmos_format;
  int      Ndep, Nrays, save_Nrays;
  double   save_muz, save_mux, save_muy, save_wmu;
  double  *height, *cmass, *tau_ref, *mux, *muy, *muz, *wmu, *vel,
          *xscale, *yscale, **Itop, **Ibottom;
} Geometry;


/* For the input atmos file */
typedef struct {
  char     *file_name;
  hid_t     ncid,     nx_id,    ny_id,    nz_id,    nhyd_id,     z_varid;
  hid_t     T_varid,  ne_varid, vz_varid, nh_varid, vturb_varid;
  hid_t     H_popsin_varid, He_popsin_varid, Ca_popsin_varid;
  hid_t     Bx_varid, By_varid, Bz_varid;
  size_t    nx,       ny,       nz,       NHydr;
  double   *x, *y;
} Input_Atmos_file;

 
/* For the background stuff */
typedef struct {
  bool_t    do_fudge;
  int       Nfudge;
  // for now here, in the future perhaps write things in io.h:
  int       j_ncid, jlambda_var, j20_var;
  double   *lambda_fudge, **fudge;
  double **chi_c,**eta_c,*sca_c,**chip_c;
} BackgroundData;


/* --- Associated function prototypes --               -------------- */
void convertScales(Atmosphere *atmos, Geometry *geometry);
void getAngleQuad(Geometry *geometry);
void getBoundary(Geometry *geometry);
void MULTIatmos(Atmosphere *atmos, Geometry *geometry);
void writeGeometry(Geometry *geometry);

void init_atmos(Atmosphere *atmos, Geometry *geometry,
                Input_Atmos_file *infile);
void readAtmos(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
		       Input_Atmos_file *infile);
void close_atmos(Atmosphere *atmos, Geometry *geometry,
                 Input_Atmos_file *infile);
void init_hdf5_atmos(Atmosphere *atmos, Geometry *geometry,
                     Input_Atmos_file *infile);
void readAtmos_hdf5(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
		            Input_Atmos_file *infile);
void close_hdf5_atmos(Atmosphere *atmos, Geometry *geometry,
                      Input_Atmos_file *infile);
void readAtmos_multi(Atmosphere *atmos, Geometry *geometry,
                      Input_Atmos_file *infile);
void convertScales(Atmosphere *atmos, Geometry *geometry);
void setTcut(Atmosphere *atmos, Geometry *geometry, double Tmax);
void realloc_ndep(Atmosphere *atmos, Geometry *geometry);
void depth_refine(Atmosphere *atmos, Geometry *geometry, double tmax);
void readPopsin_hdf5(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
                Input_Atmos_file *infile, Atom *atom);
void readPopsin(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
                Input_Atmos_file *infile, Atom *atom);


/* --- Formal solution related --                      -------------- */
double Feautrier(int nspect, int mu, double *chi, double *S,
		 enum FeautrierOrder order, double *P, double *Psi);
void Piecewise_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi);
void PiecewiseStokes(int nspect, int mu, bool_t to_obs,
		     double *chi_I, double **S, double **I, double *Psi);
void Piecewise_Hermite_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi);
void PieceBezier_1D(int nspect, int mu, bool_t to_obs,
          double *chi, double *S, double *I, double *Psi);

#endif /* !__GEOMETRY_H__ */
/* ---------------------------------------- geometry.h -------------- */
