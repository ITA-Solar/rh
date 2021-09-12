;+
;Graham Kerr
;NASA/GSFC
;graham.s.kerr@nasa.gov ; kerrg@cua.edu
;
;NAME:        makeemistabdat_RH.pro
;
;PURPOSE:     To create a look up table of emissivity  
;             using CHIANTI data, tabulated on a grid of temperature
;             and density. Several ions are excluded.
;
;INPUTS:      A file called chianti_gofnt.sav must be present 
;             in the working directory. This contains all of the 
;             necessary data from CHIANTI.
;
;OUTPUTS:     A binary file
;
;KEYWORDS:    filen_out -- The name to call the output file. Default is 
;                          coronal_rad.dat
;             savefile  -- Set to save the emissivities grid 
;
;NOTES:       This first pass will be hardcoded to the 'standard' set
;             that I inlude with RH. So, those transitions will be 
;             removed from the CHIANTI set. You can easily modify that
;             list depending on what active set you use with RH. 
;             You can also omit full ions. 
;
;             Based on mketab3.pro from Joel Allred
;
;             This will produce a grid of emissivities as fns of
;             [wavelegnth, temperature, density] that can then be
;             interpolated on given an atmosphere. 
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;     PROCEDURE CALL AND INITIAL SET UP    ;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO  makeemistabdat_RH, file_goft = file_goft,$
                        filen_in = filen_in,$
                        filen_out = filen_out,$
                        savefile = savefile, $
                        loadfile = loadfile
 
   if n_elements(filen_out) eq 0 then filen_out = 'emiss_tab' 

   ;Restore the chianti contribution fns, set up the wavelength arrays,
   ;and define some sizes
   if n_elements(file_goft) eq 0 then restore,'chianti_gofnt.sav' else restore, file_goft
   logT = chianti_gofnt[0].IONEQ_LOGT
   n_dens = n_elements(logD)
   n_temp = n_elements(logT)
   bwl = chianti_gofnt[0].wvl_limits[0]
   ewl = chianti_gofnt[0].wvl_limits[1]
   ;ewl = 40500d0
   ;bwl = 500d0
   ;ewl = 5000d0

   print, '>>>>> bwl = ',bwl
   print, '>>>>> ewl = ',ewl
   
   specres = 1.0d0
   nlambda = long((ewl - bwl) /specres +1 )
   lambda = dindgen(nlambda)/(nlambda-1) *(ewl-bwl) + bwl 
   

;; If you have already processed the CHIANTI data and produced a file containing 
;; sp, lambda etc., then you can load it here and skip the generation step. 
;; Otherwise the code will loop through denisty and temperature to produce the 
;; emissivity in each cell. 
if KEYWORD_SET(loadfile) then begin
    restore, filen_in, /v
endif else begin

   ;The specific transitions to ignore. 
   ;This Default set includes
   ; Hydrogen (some H I transitions not included)
   ; Calcium II
   ; Mg II (some Mg II transitions not included)
   ; Silicon I & II (there are no Si I)
   ; Carbon I & II (there are no C I b-b transitions in RH)
   ; Oxygen I 
   trans = [ {atom: 1, ion: 1, wl: 1215.6760},$
             {atom: 1, ion: 1, wl: 1025.7250},$
             {atom: 1, ion: 1, wl: 6564.7339},$
             {atom: 1, ion: 1, wl: 972.53900},$
             {atom: 1, ion: 1, wl: 4862.7432},$
             {atom: 1, ion: 1, wl: 18757.408},$
             {atom: 1, ion: 1, wl: 949.74500},$
             {atom: 1, ion: 1, wl: 4341.7290},$
             {atom: 1, ion: 1, wl: 12822.375},$
             {atom: 1, ion: 1, wl: 40524.550}];,$
         ;    {atom: 20, ion: 2, wl: 3969.5911},$
         ;    {atom: 20, ion: 2, wl: 3934.7771},$
         ;    {atom: 20, ion: 2, wl: 8664.5215},$
         ;    {atom: 20, ion: 2, wl: 8500.3584},$
         ;    {atom: 20, ion: 2, wl: 8544.4385},$
         ;    {atom: 12, ion: 2, wl: 2803.5300},$
         ;    {atom: 12, ion: 2, wl: 2796.3521},$
         ;    {atom: 12, ion: 2, wl: 1071.6840},$
         ;    {atom: 12, ion: 2, wl: 2929.4900},$
         ;    {atom: 12, ion: 2, wl: 2937.3689},$
         ;    {atom: 12, ion: 2, wl: 2798.7539},$
         ;    {atom: 12, ion: 2, wl: 2791.6001},$
         ;    {atom: 12, ion: 2, wl: 2798.7539},$
         ;    {atom: 12, ion: 2, wl: 10929.343},$
         ;    {atom: 12, ion: 2, wl: 4482.3828},$
         ;    {atom: 12, ion: 2, wl: 4482.5840},$
         ;    {atom: 12, ion: 2, wl: 9220.7734},$
         ;    {atom: 12, ion: 2, wl: 8236.9092},$
         ;    {atom: 12, ion: 2, wl: 7898.2168},$ 
         ;    {atom: 14, ion: 2, wl: 1533.4310},$
         ;    {atom: 14, ion: 2, wl: 1309.2760},$
         ;    {atom: 14, ion: 2, wl: 6348.8628},$
         ;    {atom: 14, ion: 2, wl: 1194.5000},$
         ;    {atom: 14, ion: 2, wl: 1816.9280},$
         ;    {atom: 14, ion: 2, wl: 3857.1111},$
         ;    {atom: 14, ion: 2, wl: 1264.7380},$
         ;    {atom: 6, ion: 2, wl: 2322.6860},$
         ;    {atom: 6, ion: 2, wl: 1334.5770},$
         ;    {atom: 6, ion: 2, wl: 2326.1130},$
         ;    {atom: 6, ion: 2, wl: 1335.7080},$
         ;    {atom: 6, ion: 2, wl: 1335.6630},$
         ;    {atom: 8, ion: 1, wl: 1355.5980},$
         ;    {atom: 8, ion: 1, wl: 1358.5120},$
         ;    {atom: 8, ion: 1, wl: 1302.1689},$
         ;    {atom: 8, ion: 1, wl: 1304.8580},$
         ;    {atom: 8, ion: 1, wl: 1306.0291},$
         ;    {atom: 8, ion: 1, wl: 1641.3051}]

    sp = dblarr(nlambda,n_temp,n_dens)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;     COMPUTE THE EMISSIVITIES    ;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;; Start with the first density value to get the various sizes 
   ;; required to collate the files later.
   ;; Here you can also define which species, and specific ions, to 
   ;; ignore in addition to those transitions defined above.
   dind = 0
   print, ' '
   print, ' '
   print, ' '
   print, '>>> Density ind = ',dind
   ind=replicate(1b,n_elements(CHIANTI_GOFNT[dind].lines))
   for i= 0, n_elements(trans)-1 do $
      ind=ind $
      	      and ((CHIANTI_GOFNT[dind].lines.iz ne trans[i].atom) or (CHIANTI_GOFNT[dind].lines.ion ne trans[i].ion) or (abs(trans[i].wl-CHIANTI_GOFNT[dind].lines.wvl) gt 3.0)) $
      	      and ((CHIANTI_GOFNT[dind].lines.iz ne 1)) $
              and ((CHIANTI_GOFNT[dind].lines.iz ne 2)) ;$
      	      ;and ((CHIANTI_GOFNT[dind].lines.iz ne 14) or (CHIANTI_GOFNT[dind].lines.ion ne 2)) $
      	      ;and ((CHIANTI_GOFNT[dind].lines.iz ne 14) or (CHIANTI_GOFNT[dind].lines.ion ne 3)) $
      	      ;and ((CHIANTI_GOFNT[dind].lines.iz ne 14) or (CHIANTI_GOFNT[dind].lines.ion ne 4)) $ 
              ;and ((CHIANTI_GOFNT[dind].lines.iz ne 6) or (CHIANTI_GOFNT[dind].lines.ion ne 1)) $ 
              ;and ((CHIANTI_GOFNT[dind].lines.iz ne 6) or (CHIANTI_GOFNT[dind].lines.ion ne 2)) $ 
              ;and ((CHIANTI_GOFNT[dind].lines.iz ne 6) or (CHIANTI_GOFNT[dind].lines.ion ne 3)) 

   ind=where(ind)
   ;; Create a structure containing only the transitions of interest. It is identical to gofnt but removes
   ;; the transitions defined above
   makespecstr,CHIANTI_GOFNT[dind], str_tmp, ind=ind
   str = str_tmp
   ;; Loop through temperature and grab the emissivity in each cell.
   ;; The outputs here are lambda and spec, the emissivity in cgs units
   for tind = 0, n_temp - 1 do begin
          print, '....... Temperature ind = ', tind
	  uvlines_etab, bwl, ewl, 10d0^(logT[tind]), 10d0^(logD[dind]), lambda, spec, $
                        nthreads = 5,$
                        ;nx = 101L,$
                        binsize = specres, contribfunc=1, chianti_gofnt = str_tmp, abundfile = '/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2009_asplund.abund'
      sp[*,tind,dind] = spec
   endfor
   ;sun_photospheric_2009_asplund.abund
   ;sun_coronal_2012_schmelz_ext.abund

   ;; Repeat the above for the remaining density values    
   for dind = 1, n_dens-1 do begin
       print, ' '
       print, ' '
       print, ' '
       print, '>>> Density ind = ',dind
       ind=replicate(1b,n_elements(CHIANTI_GOFNT[dind].lines))
   for i= 0, n_elements(trans)-1 do $
       ind=ind $
      	      and ((CHIANTI_GOFNT[dind].lines.iz ne trans[i].atom) or (CHIANTI_GOFNT[dind].lines.ion ne trans[i].ion) or (abs(trans[i].wl-CHIANTI_GOFNT[dind].lines.wvl) gt 3.0)) $
      	      and ((CHIANTI_GOFNT[dind].lines.iz ne 1)) $
              and ((CHIANTI_GOFNT[dind].lines.iz ne 2)) ;$
      	      ;and ((CHIANTI_GOFNT[dind].lines.iz ne 14) or (CHIANTI_GOFNT[dind].lines.ion ne 2)) $
      	      ;and ((CHIANTI_GOFNT[dind].lines.iz ne 14) or (CHIANTI_GOFNT[dind].lines.ion ne 3)) $
      	      ;and ((CHIANTI_GOFNT[dind].lines.iz ne 14) or (CHIANTI_GOFNT[dind].lines.ion ne 4)) $
              ;and ((CHIANTI_GOFNT[dind].lines.iz ne 6) or (CHIANTI_GOFNT[dind].lines.ion ne 1)) $ 
              ;and ((CHIANTI_GOFNT[dind].lines.iz ne 6) or (CHIANTI_GOFNT[dind].lines.ion ne 2)) $ 
              ;and ((CHIANTI_GOFNT[dind].lines.iz ne 6) or (CHIANTI_GOFNT[dind].lines.ion ne 3))  
 
   ind=where(ind)
   makespecstr,CHIANTI_GOFNT[dind], str_tmp, ind=ind
   str = [str, str_tmp]
   for tind = 0, n_temp - 1 do begin
          print, '....... Temperature ind = ', tind
	  uvlines_etab, bwl, ewl, 10d0^(logT[tind]), 10d0^(logD[dind]), lambda, spec, $
                             nthreads = 5,$
                             binsize = specres, contribfunc=1, chianti_gofnt = str_tmp, abundfile = '/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2009_asplund.abund'
      sp[*,tind,dind] = spec
    endfor
  endfor

  ;; Save the data either for later use, or to make it easier to inspect 
  if KEYWORD_SET(savefile) then begin
    filen_emiss = filen_out+'_grid.sav'   
    save, sp, lambda, logT, logD,str, file = filen_emiss
  endif

  delvar, str
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;     WRITE THE GRID    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Integrate over solid angle
sp*=4*!dpi 
;Divide by 2 since only half the radiation is directed downward.
sp/= 2 

;sp = sp*2


;; Cell-center the values
bins = lambda
n_bins=n_elements(bins)
ebin=dblarr(n_temp,n_dens,n_bins-1)
cwl=dblarr(n_temp,n_dens,n_bins-1)
wlbin=dblarr(n_bins-1)
dlam = [lambda[1:*] - lambda,0]
for i=0, n_bins-2 do begin
  wlbin[i]=(bins[i+1]+bins[i])/2.
  ind=where(lambda ge bins[i] and lambda lt bins[i+1])
  ebin[*,*,i] = total(sp[ind,*,*] * (replicate(dlam[ind],n_temp,n_dens)),1)
  cwl[*,*,i] = total(sp[ind,*,*] * (replicate((dlam[ind] * lambda[ind]),n_temp,n_dens)),1) / ebin[*,*,i]
  ebin[*,*,i] /= (bins[i+1]-bins[i])
endfor

ind=where(ebin lt 1d-99)
if (ind[0] ne -1) then ebin[ind]=1d-99

mwl=dblarr(n_elements(bins)-1)
for i=0,n_elements(bins)-2 do begin
  mwl[i]=total(ebin[*,*,i]*cwl[*,*,i],/nan)/total(ebin[*,*,i],/nan)
  if (mwl[i] eq 0) then mwl[i]=(bins[i+1]+bins[i])/2.
endfor
; Write:
;  nlambda
;  n_temp
;  n_dens
;  lambda
;  cwl
;  dlam
;  10^logT
;  10^logD
;  ebin
;  
get_lun, lu
openw, lu, filen_out+'.dat'
writeu, lu, nlambda-1, n_temp, n_dens, mwl, 10d0^(logT), 10d0^(logD), ebin
free_lun, lu

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
END



