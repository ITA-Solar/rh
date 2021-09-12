PRO uvlines_etab, bwl, ewl, tempr, edens, lambda, spec, $
                  nx=nx, istep=istep, binsize = binsize1, contribfunc = contribfunc1, los = los1, $
                  sngl_ion = sngl_ion, chianti_gofnt = chianti_gofnt, xlim = xlim1, abundfile = abundfile, trind = trind,$
                  exclude_ion = exclude_ion, linesout = linesout, minintensity = minintensity1, vbroad = vbroad1, nthreads = nthreads1

timet = [0d0]
z1t = [1d0]
tg1t = [tempr]
ne1t = [edens]
vz1t = [0d0]
dzt = [0d0]
d1t = 0.83*ne1t
grph = 1d0
ndep = long(1)

if (n_elements(binsize1) eq 0) then binsize = 1d0 else binsize = double(binsize1)
if (n_elements(istep) eq 0) then istep = lindgen(n_elements(timet))
if (n_elements(contribfunc1) eq 0) then contribfunc = 0L else contribfunc = long(contribfunc1)
if (n_elements(xlim1) eq 0) then xlim = 5d0 else xlim = double(xlim1)
if (n_elements(los1) lt ndep) then los = replicate(1d0,ndep) else los = double(los1)
if (n_elements(minintensity1) eq 0) then minintensity = -1d0 else minintensity = double(minintensity1)
if (n_elements(vbroad1) eq 0) then vbroad = dblarr(ndep,n_elements(istep)) else vbroad = double(vbroad1)
if (n_elements(nthreads1) eq 0) then nthreads = long(!cpu.tpool_nthreads) else nthreads = long(nthreads1)
if (n_elements(nx) eq 0) then nx = 101L

if (n_elements(chianti_gofnt) eq 0) then begin
  restore,'~/IDL_Programs/CHIANTI_Analysis/si4goft_chianit_ioneq.sav'
  print, '>>>>>>>> LOADING PREEXISTING GOFT'
endif else begin
  logd = alog10(chianti_gofnt.model_ne)
endelse

if (n_elements(abundfile) eq 0) then begin
  ;abundfile = '~/ssw/packages/chianti/dbase/abundance/sun_photospheric_2009_asplund.abund'
 ; abundfile = '/usr/local/ssw/packages/chianti/dbase/abundance/sun_coronal_1999_fludra_ext.abund'
  abundfile = '/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2009_asplund.abund'
 ;abundfile = '/usr/local/ssw/packages/chianti/dbase/abundance/sun_coronal_1992_feldman_ext.abund'
endif
readcol,abundfile,n, abund
abund = 10^(abund - 12d0)

dir = routine_filepath('uvlines_etab')
dir = strmid(dir,0,strpos(dir,'uvlines_etab.pro')-1)
incdir = !dir + '/external/include'

;cc='gcc -mcmodel=large -fopenmp -O3 -march=native -lgomp -fPIC -I"'+incdir+'" %c -c -o %o'
;ld='gcc -mcmodel=large -fopenmp -shared -o %L %O %X'
cc='gcc -mcmodel=medium -mlarge-data-threshold=1000 -fopenmp -O3 -march=native -lgomp -fPIC -I"'+incdir+'" %c -c -o %o'
ld='gcc -mcmodel=medium -mlarge-data-threshold=1000 -fopenmp -shared -o %L %O %X'
if (nthreads gt 1) then $
  make_dll,'uvlines.par', 'uvlines.par', ['uvlines_natural', 'uvlinesc'], input_directory = dir, output_directory = dir, dll_path = file, /reuse,/verb,cc=cc,ld=ld $
else $
  make_dll,'uvlines', 'uvlines', ['uvlines_natural', 'uvlinesc'], input_directory = dir, output_directory = dir, dll_path = file, /reuse,/verb,cc=cc,ld=ld

if n_elements(trind) ne 0 then begin
    lind = trind
    nlines = 1
    linesgoft = chianti_gofnt.lines[lind].goft
    lineswvl = chianti_gofnt[0].lines[lind].wvl
    izn = long(chianti_gofnt[0].lines[lind].iz-1)
    goto, nextstep
endif  

if (n_elements(sngl_ion) ne 0) then begin
  lind = bytarr(n_elements(chianti_gofnt[0].lines))
  for i = 0,n_elements(sngl_ion)-1 do begin
    convert_ion_to_num,sngl_ion[i],iz,izion
    lind = lind or (chianti_gofnt[0].lines.iz eq iz and chianti_gofnt[0].lines.ion eq izion)
  endfor
endif else begin
  lind = replicate(1b,chianti_gofnt[0].nlines)
endelse
if (n_elements(exclude_ion) ne 0) then begin
  for i = 0, n_elements(exclude_ion)-1 do begin
    convert_ion_to_num,exclude_ion[i],iz,izion
    lind = lind and (chianti_gofnt[0].lines.iz ne iz or chianti_gofnt[0].lines.ion ne izion)
  endfor
endif
lind = lind and ((chianti_gofnt[0].lines.wvl gt bwl-3d0) and (chianti_gofnt[0].lines.wvl lt ewl + 3d0)) 

lind = where(lind)
if (lind[0] eq -1) then begin
  print,'No lines in specified region.'
  return
endif


nlines = n_elements(lind)
linesgoft = chianti_gofnt.lines[lind].goft
lineswvl = chianti_gofnt[0].lines[lind].wvl
izn = long(chianti_gofnt[0].lines[lind].iz-1)

;nextstep:

logt = chianti_gofnt[0].ioneq_logt
nt = n_elements(istep)
ngdens = n_elements(logd)
ngtemp = n_elements(logt)
nabund = n_elements(abund)
nz= long(ndep)

nl = long((ewl - bwl)/binsize+1)
if ((nl-1)*binsize + bwl lt ewl) then nl+=1
lambda = dindgen(nl)*binsize + bwl

if (contribfunc) then begin
  spec = dblarr(nl, ndep, nt)
  if nt le 1 then linesint = dblarr(ndep,nlines) else linesint = dblarr(ndep,nt,nlines)
endif else begin
  spec = dblarr(nl, nt)
  linesint = dblarr(nt,nlines)
endelse

;linesint = reform(linesint)
;spec = reform(spec)

v = vz1t * los ; DOT V with LOS unit vector

res = call_external(file, 'uvlinesc', tg1t, ne1t, d1t, v, dzt, vbroad, linesgoft,  [logd], logt,$ 
                    lineswvl, lambda, spec, linesint, abund, izn, minintensity, grph, xlim, contribfunc, nz, nt, nlines, ngdens, $
                    ngtemp, nl, nabund,  nthreads, $
                    /i_value, value = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1])
;print, '>>>>>>>>>> MADE IT OUT OF UVLINESC'

if (contribfunc and nt gt 1) then itemp = dblarr(ndep,nt) else if (contribfunc and nt le 1) then itemp = dblarr(ndep) else if (contribfunc and nt le 1 and ndep le 1) then itemp = 0d0 else if (~contribfunc and nt gt 1) then itemp = dblarr(nt) else if (~contribfunc and nt le 1) then itemp = 0d0
linesout = replicate({int:itemp, snote:'', wvl:0d0, tmax:0.0, index:0L},nlines)

linesout.int = (linesint)
linesout.wvl = lineswvl
linesout.snote = chianti_gofnt[0].lines[lind].snote
linesout.tmax = chianti_gofnt[0].lines[lind].tmax
linesout.index = lind

if (n_elements(istep) eq 1) then spec = reform(spec)

END
