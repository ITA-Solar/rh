#include <string.h>
#include <stdlib.h>

#include "rh.h"
#include "atmos.h"
#include "geometry.h"
#include "error.h"
#include "spectrum.h"
#include "inputs.h"
#include "io.h"

/* --- Global variables --                             -------------- */
extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern IO_data io;
extern char messageStr[];


void readRayInput(void) {
    /* Reads ray.input file */
    const char routineName[] = "init_atmos";
    bool_t exit_on_EOF;
    int Nspect, Nread, Nrequired, checkPoint;
    double muz;
    FILE  *fp_ray;
    char  inputLine[MAX_LINE_SIZE];

    /* --- Read direction cosine for ray --              -------------- */
    if ((fp_ray = fopen(RAY_INPUT_FILE, "r")) == NULL) {
      sprintf(messageStr, "Unable to open inputfile %s", RAY_INPUT_FILE);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    getLine(fp_ray, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    Nread = sscanf(inputLine, "%lf", &muz);
    checkNread(Nread, Nrequired=1, routineName, checkPoint=1);
    if (muz <= 0.0  ||  muz > 1.0) {
      sprintf(messageStr,
  	    "Value of muz = %f does not lie in interval <0.0, 1.0]\n", muz);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    } else io.ray_muz = muz;
    /* --- read how many points to write detailed S, chi, eta, etc ---- */
    Nread = fscanf(fp_ray, "%d", &Nspect);
    checkNread(Nread, 1, routineName, checkPoint=2);
    io.ray_nwave_sel = Nspect;
    io.ray_wave_idx = NULL;
     /* --- Read wavelength indices for which chi and S are to be
         written out for the specified direction --    -------------- */
    if (Nspect > 0) {
      io.ray_wave_idx = (int *) malloc(Nspect * sizeof(int));
      Nread = 0;
      while (fscanf(fp_ray, "%d", &io.ray_wave_idx[Nread]) != EOF) Nread++;
      checkNread(Nread, Nspect, routineName, checkPoint=3);
      fclose(fp_ray);
    }
    /* --- Save geometry values to change back after --    ------------ */
    geometry.save_Nrays = atmos.Nrays;
    geometry.save_wmu = geometry.wmu[0];
    geometry.save_muz = geometry.muz[0];
    geometry.save_mux = geometry.mux[0];
    geometry.save_muy = geometry.muy[0];
}


void checkValuesRayInput(void) {
    /* Checks if wavelength indices make sense */
    const char routineName[] = "checkValuesRayInput";
    int i;

    for (i = 0;  i < io.ray_nwave_sel;  i++) {
        if (io.ray_wave_idx[i] < 0  ||  io.ray_wave_idx[i] >= spectrum.Nspect) {
            sprintf(messageStr, "Illegal value of wave_index[n]: %4d\n"
                    "Value has to be between 0 and %4d\n",
                    io.ray_wave_idx[i], spectrum.Nspect);
            Error(ERROR_LEVEL_2, routineName, messageStr);
        }
    }
}
