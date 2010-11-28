PRO write1datmos, atmos, FILENAME=filename

  CM_TO_M = 1.0E-2
  KM_TO_M = 1.0E+3

  openw, unit, filename, /GET_LUN

  printf, unit, atmos.atmosID + " (extended)"

  printf, unit, $
   FORMAT='(2X,"Height scale")'
  printf, unit, $
   FORMAT='("*",/"* lg g",/F6.2,/"*",/"* Ndep",/I4,/"*")', $
   alog10(atmos.gravitation), atmos.Ndep
  printf, unit, $
   FORMAT='("*height [km]   ", 5X,"Temperature",8X,"Ne",9X,"V",14X,"Vturb")'

  FOR k=0, atmos.Ndep-1 DO $
   printf, unit, FORMAT='(E17.8, 4E15.6)', $
   atmos.height[k], atmos.T[k], atmos.n_elec[k] * CM_TO_M^3 , $
   atmos.v[k]/KM_TO_M, atmos.vturb[k]/KM_TO_M

  printf, unit, FORMAT='("*",/"* Hydrogen populations (LTE)")'
 
  printf, unit, FORMAT=$
   '("*",4X,"nh(1)",7X,"nh(2)",7X,"nh(3)",7X,"nh(4)",7X,"nh(5)",7X,"np")'

  FOR k=0, atmos.Ndep-1 DO $
   printf, unit, FORMAT='(6E12.4)', atmos.nH[k, *]* CM_TO_M^3

  free_lun, unit
END
