/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Authors:       Han Uitenbroek   (huitenbroek@nso.edu)
                      Jorrit Leenaarts (j.leenaarts@astro.uu.nl)

       Last modified: Thu Mar  4 08:22:25 2004 --

       --------------------------                     -----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__


/* --- Defines structures to contain the angle quadrature and
       interpolation mesh for 3-D Cartesian version. --  ------------ */

#define NCC 4

enum boundary  {FIXED, PERIODIC};
enum boundval  {IRRADIATED, ZERO, PLANCK};
enum vertical  {TOP=0, BOTTOM};
enum direction {UPWIND=0, DOWNWIND};
enum intersect {XY, XZ, YZ};

typedef struct Stencil Stencil;
typedef struct Longchar Longchar;

struct Stencil {
  enum intersect plane;
  int    quadrant, xbase[NCC], ybase[NCC], zbase[2];
  double ds, xkernel[NCC], ykernel[NCC], zkernel[2];
  Longchar *longchar;
};

struct Longchar {
  int      Nst;
  Stencil *stencil;
};

typedef struct {
  enum      boundary x_boundary, y_boundary;
  enum      boundval z_boundary_value[2]; 
  int       Nx, Ny, Nz, Nplane, Nrays, Nlongchar;
  double    dx, dy, *z, *mux, *muy, *muz, *wmu, *vx, *vy, *vz,
          **Itop, **Ibottom;
  Stencil **stencil_uw, **stencil_dw;
} Geometry;

typedef struct {
  enum intersect plane;
  double s, x, y, z;
} Intersect;

/* --- Associated function prototypes --               -------------- */

void   readAtmos(Atmosphere *atmos, Geometry *geometry);
void   ShortChar(Geometry *geometry, int nspect, int mu,
		 bool_t to_observer, double *chi, double *S, double *I, 
                 double *Psi);
void   PiecewiseStokes(Geometry *geometry, int nspect, int mu,
		       bool_t to_observer, double *chi,
		       double **S, double **I, double *Psi);
void   getAngleQuadr(Geometry *geometry);
void   fillMesh(Geometry *geometry);
void   getBoundary(Atmosphere *atmos, Geometry *geometry);
bool_t writeFlux(char *fileName);
void   writeGeometry(Geometry *geometry);
double Interpolate_3D(double *f, Geometry *geometry, Stencil *st,
		      int l, int m);

void   StokesK_3D(int nspect, Geometry *geometry, Stencil *st,
		  int l, int m, double chi_I, double K[4][4]);
double SolveLong(Geometry *geometry, Longchar *lc, int k, int l, int m,
		 double *chi, double *S, double *I);
void   SolveLongStokes(Geometry *geometry, Longchar *lc,
		       int nspect, int k, int l, int m,
		       double *chi, double **S, double **I,
		       double *I_uw);

#endif /* !__GEOMETRY_H__ */
 
/* ---------------------------------------- geometry.h -------------- */

