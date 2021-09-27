----------------------------------
----------------------------------
RH15D
Using fixed populations from input file
----------------------------------
----------------------------------

Graham Kerr
Sept 2021

graham.s.kerr@nasa.gov; grahamkerr.astro@gmail.com


--------
SUMMARY
--------

Describes the changes made to (RH 1.5D version) to use user provided populations for a selected species, rather than solving the population equations. That is, only the radiation transport part is solved, with the assumption that the populations provided are correct. These can either be from a previously run simulation of the same atmosphere but with additional/different active atoms, or from another simulation such as RADYN (the original intended purpose). 


The downside is that the population changes caused by, e.g. overlapping transitions or PRD etc., are not taken into account. The advantage is that non-equilibrium populations can be used.


I have not implemented this in a very elegant way. It was not clear if a second HDF5 file could be opened at read at the same time as the atmosphere file (I suspect not), so I was forced to add 'popsin' variables for each species individually. I have so only done this for H, He, Ca, which are the most likely candidates. This can be changed easily in the future. 

-------
TO DO
-------

* include the rhf1d version (should be rather straightforward)
* improve doc of routines in

-------
rh.h
-------

-- Added FIXED_POP_FROM_FILE case to "enum 'solution'" on line 26:
		      enum solution       {UNKNOWN=-1, LTE_POPULATIONS, ZERO_RADIATION,
		                           OLD_POPULATIONS, NEW_J, OLD_J, ESCAPE_PROBABILITY,
		                           FIXED_POP_FROM_FILE};

--------
readatom.c
--------

-- Added flag for the FIXED_POP_FROM_FILE initial solution:
      
		      else if (strstr(popsKey, "FIXED_POPS_FROM_FILE")) {
		         atom->initial_solution = FIXED_POPS_FROM_FILE;
		      }

-----------
statequil.c
-----------

-- Added the following so that the equations are not solved in the case that FIXED_POPS_FROM_FILE
keyword is set:

			int dosolvelineq = 1;
			    if (atom->initial_solution == FIXED_POPS_FROM_FILE) {
			       dosolvelineq = 0;
			    } else {
			       /* printf("\n Running standard version version..."); */
			        n_k[i_eliminate] = atom -> ntotal[k];
			        for (j = 0;  j < Nlevel;  j++) Gamma_k[i_eliminate][j] = 1.0;
			    }


			    /* --- Solve for new population numbers at location k -- -------- */
			    if (dosolvelineq){

			       SolveLinearEq(Nlevel, Gamma_k, n_k, TRUE);

			    }


--------
io.h      (rh15d)
--------

-- Added defn for H_POPSIN_NAME, HE_POPSIN_NAME, CA_POPSIN_NAME

				#define H_POPSIN_NAME "H_popsin"
				#define He_POPSIN_NAME "He_popsin"
				#define Ca_POPSIN_NAME "Ca_popsin"

--------
geometry.h (rh15d)
--------

-- Added to the hdf5 input atmos structure: 

  				hid_t     H_popsin_varid, He_popsin_varid, Ca_popsin_varid;

-- Added fn defns for:

				void readPopsin_hdf5(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
				                Input_Atmos_file *infile, Atom *atom);
				void readPopsin(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
				                Input_Atmos_file *infile, Atom *atom);


---------
readpopsin.c (rh15d)
---------

-- Added this script which calls the function to read the pops from the atmosphere file 
(probably redundant and can skip!... started life when I thought I would be using a seperate file.)
                 
                 readPopsin_hdf5(xi, yi, atmos, geometry, infile, atom);

--------
hdf5popsin.c (rh15d)
--------

-- Added this script to read in the pops from the atmos file. It is currently set up to read the atom ID
of the input atom (via initial_p.c), and search for the appropriate populations to assign to atom->n

------------
hdf5atmos.c (rh15d)
------------

-- Added the following to the init_hd5f_atmos routine, to make sure that if they are ever read, the code knows to close them:

					  infile->H_popsin_varid = -1;
					  infile->He_popsin_varid = -1;
					  infile->Ca_popsin_varid = -1;

-- Added the following to close_hdf5_atmos routine, to close the variables if they are read:

					  if (infile->H_popsin_varid != -1) ierror = H5Dclose(infile->H_popsin_varid);
					  if (infile->He_popsin_varid != -1) ierror = H5Dclose(infile->He_popsin_varid);
					  if (infile->Ca_popsin_varid != -1) ierror = H5Dclose(infile->Ca_popsin_varid);

--------
initial_p.c (rh15d)
--------

-- Added FIXED_POPS_FROM_FILE case. which calls readPopsin:

				readPopsin(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry,
                  			&infile, atom);

    Had to add these as globals:

			    Geometry geometry;
				Input_Atmos_file infile;

******** DO I NEED TO DO ANYTHING IN PARALLEL.C (updateatmosdep) ... or do I need to keep to keep track of zcut etc., 
For now get it working without refinement or cuts


