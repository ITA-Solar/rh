/* ------- file: -------------------------- trilinear_interp.c ------------------
 *
 *        Version:        rh_ita
 *        Author:         Graham S. Kerr (graham.s.kerr@nasa.gov; gskerr89@gmail.com)
 *        First modified: Mon May 20th 2019  --
 *
 * -----------------------------------------                     ----------RH-- */

/* --- Perform a trilinear interpolation
 *
 * If we define the vertices Cxyz of a rectangle as 
 *    P000, P100, P010, P110, P001, P101, P110, P111, 
 * such that the point of interest P lies within this cube, so within x0->x1
 * y0->y1, and z0->z1, then the value of this point can be written as 
 * 
 *   f(x,y,x) ~ C0 + C1dx + C2dy + C3dz + C4dxdy + C5dydz + C6*dx*dz + C7dxdydz
 *
 * The coeffcients are written below in the main routine 
 *
 *
 * ---------------------------------------------------------------------------- */

#include "rh.h"

/* ------- begin ------------------------- TrilinearInterp.c ------------------ */

void TrilinearInterp(int nX, int nY, int nZ, 
                     double *grid, double *xarr, double *yarr, double *zarr,
                     int x0, int y0, int z0, 
                     double xp, double yp, double zp, 
                     double *p)
{ 
  double p000, p100, p010, p110, p001, p101, p011, p111; 
  double c0, c1, c2, c3, c4, c5, c6, c7;
  double dx, dy, dz; 
/* Define the vertices locations
 * x0, y0, z0 are the the indices of the array grid[x,y,z], such that 
 * the point of interest lies between [x0,x0+1],[y0,y0+1],[z0,z0+1].
 * On the grid to search through these vertices are found by:
 *     grid[ x0 + y0*nX + z0*nX*nY]
 */
 p000 = grid[(x0)   + (y0)*nX   + (z0)*nX*nY ];
 p100 = grid[(x0+1) + (y0)*nX   + (z0)*nX*nY ];
 p010 = grid[(x0)   + (y0+1)*nX + (z0)*nX*nY ];
 p110 = grid[(x0+1) + (y0+1)*nX + (z0)*nX*nY ];
 p001 = grid[(x0)   + (y0)*nX   + (z0+1)*nX*nY];
 p101 = grid[(x0+1) + (y0)*nX   + (z0+1)*nX*nY];
 p011 = grid[(x0)   + (y0+1)*nX + (z0+1)*nX*nY];
 p111 = grid[(x0+1) + (y0+1)*nX + (z0+1)*nX*nY];

 c0 = p000;
 c1 = p100 - p000;
 c2 = p010 - p000; 
 c3 = p001 - p000; 
 c4 = p110 - p010 - p100 + p000; 
 c5 = p011 - p001 - p010 + p000;
 c6 = p101 - p001 - p100 + p000;
 c7 = p111 - p011 - p101 + p100 + p010 - p000;

 dx = (xp-xarr[x0]) / (xarr[x0+1]-xarr[0]);
 dy = (yp-yarr[y0]) / (yarr[y0+1]-yarr[0]);
 dz = (zp-zarr[z0]) / (zarr[z0+1]-zarr[0]);

 *p = c0 + c1*dx + c2*dy + c3*dz +c4*dx*dy + c5*dy*dz + c6*dz*dx + c7*dx*dy*dz; 

}
/* ------- end ----------------------------TrilinearInterp.c ----------------- */

