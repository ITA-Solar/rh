;PRO solve1d

  ;; IDL driver for one-dimensional two-level atom solver with given
  ;; Planck function Bp and collisional destruction parameter epsil.
  ;;
  ;; Han Uitenbroek

  ;; Quantities:
  ;;            NmaxIter  - Maximim number OF iterations
  ;;            NgOrder   - Order of Ng convergence acceleration
  ;;            NgPeriod  - Period of Ng convergence acceleration
  ;;            Ndep      - Number of depth points
  ;;            Nrays     - Number of angles
  ;;            Nlambda   - Number of wavelength points
              
  NmaxIter  = 200L  &  iterLimit = 1.0D-4
  Ngdelay   =  3L   &  Ngorder   = 2L      &  Ngperiod = 4L     
       Ndep = 100L  &      Nrays = 3L      &  Nlambda  = 21L
  
  heightmin  = 0.0   & heightmax = 600.0     &  scaleHeight = 60.0
        
  ;;            chi       - Opacity

  height = heightmin + (heightmax - heightmin)*findgen(Ndep)/float(Ndep - 1)
  height = double(reverse(height))
  chi0  = 1.0D+06
  chi   = double(chi0 * 10.0^(-height/scaleHeight))
  
  ;;            Boundary  - upper: Iminus = Boundary[0]*Iplus  + Boundary[1]
  ;;                        lower: Iplus  = Boundary[2]*Iminus + Boundary[3]

  ;;  BoundCond = double([0.0, 0.0, 1.0, 0.0])  ;; ie reflective at bottom

  ;;            Adamp     - Damping parameter (set to 0.0 for Doppler profile,
  ;;                        but limit to 7 Doppler widths)
  ;;            lambda    - Wavelength in Doppler units

  Adamp  = double(0.01)
  lambda = 15.0 * dindgen(Nlambda)/double(Nlambda-1)
  
  ;;            Bp      - Planck function
  ;;            epsil   - Collisional destruction parameter

  Bp = dblarr(Ndep) + 1.0  
  epsil = dblarr(Ndep) + 1.0D-4

  tau = dblarr(Ndep)
  tau(0) = 0.0
  FOR k=0, Ndep-2 DO BEGIN
    dtau     = 0.5*(chi(k) + chi(k+1))*(height(k) - height(k+1))
    tau(k+1) = tau(k) + dtau
  ENDFOR
  Sth = (1.0 - (1.0 - sqrt(epsil)) * exp(-sqrt(epsil)*tau)) * Bp
  
  ;; Spawn C program and write input data to pipe

  spawn, './solve1d', /NOSHELL, UNIT=pipe
;;  openw, pipe, 'solve1d.dat', /GET_LUN

  writeu, pipe, NmaxIter, double(iterLimit), Ngdelay, Ngorder, Ngperiod
  writeu, pipe, Nrays, Nlambda, Ndep
  writeu, pipe, height
  writeu, pipe, Adamp, lambda
  writeu, pipe, chi, Bp, epsil

  ;; allocate space and read back results from pipe

  phi = dblarr(Nlambda)  &  wphi = 0.0D+0
  S   = dblarr(Ndep)     &  Ie   = dblarr(Nlambda, Nrays)

  readu, pipe, phi, wphi, S, Ie
  free_lun, pipe
  
  loadct, 1
  plot, tau(1:Ndep-1), S(1:Ndep-1), YTITLE='Source function', $
   xtitle='Optical depth', /XLOG, /YLOG, /NODATA

  oplot, tau(1:Ndep-1), S(1:Ndep-1), COLOR=225
  oplot, tau(1:Ndep-1), Sth(1:Ndep-1), COLOR=175
END
