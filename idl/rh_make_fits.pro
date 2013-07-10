pro rh_make_fits,quantity,rayfile,input_dir,bifrost_dir,bifrost_snapname,runname,runid=runid,$
                 obstitle=obstitle,debug=debug,cp_to_rep=cp_to_rep
;
;+
; NAME:
;	RH_MAKE_FITS
;
; PURPOSE:
;       make FITS file version of rh15d output
;
; CATEGORY:
;       Bifrost, rh
;	
; CALLING SEQUENCE:
;	RH_MAKE_FITS,quantity,rayfile,runid,bifrost_dir,runname,rh_dir=rh_dir,$
;                  obstitle=obstitle,debug=debug,cp_to_rep=cp_to_rep
;
; INPUTS:
;       quantity    - either 'intensity' or 'zt1'
;       rayfile     - rh 15d ray output ncdf file
;	bifrost_dir - directory where bifrost atmosphere file is
;                     located
;  bifrost_snapname - name of the bifrost snapshots, e.g. 'cb24bih'
;	runname     - character description of the region
;                     modelled with Bifrost. Must be supplied to link
;                     this file to the corresponding bifrost atmosphere
;                     FITS files. e.g, 'en024048_hion'
;     
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;     obstitle - a short optional description of the run. e.g, 'Mg II 1.5D
;                PRD'
;     debug    - stop execution before returning.
;     runid    - directory where input files will be stored, if not
;                set, will be created automatically. The runid can be set
;                by hand but this is not recommended
;  cp_to_rep   - string, set to the subdirectory of
;                /Volums/SAM5/oslo_sim/[runname]
;                where the fits file will be
;                copied to. If set, the fits file will be copied there
;                and the content of input_dir will be copied to the
;                runid directory. If not set, the fits file will be
;                created in the calling directory, and the input data
;                will not be copied anywhere.
;
; OPTIONAL KEYWORD PARAMETERS:
;
; OUTPUTS:
;       writes data to fits file, and copies input data to the
;       oslo_sim repository.
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;       Procedures:  MKHDR, SXADDPAR, WRITEFITS,rh_idl2fits,rh_idl2fits_varrep,jncdf_varexist,jncdf_varread
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   v1.0 Jul  5 2013 jorritl - started
;   v1.1 Jul 10 2013 jorritl - added more options, seems to work nicely
;
;   $Id:$ 
;-
;

  if(n_params() lt 6) then begin 
     print,'pro rh_make_fits,quantity,rayfile,input_dir,bifrost_dir,bifrost_snapname,runname,runid=runid,'
     print,'                 obstitle=obstitle,debug=debug,cp_to_rep=cp_to_rep'
     return
  endif

; constants and conversions
  m_to_mm=1e-6

; check that quantity is correct
  quantities=['intensity','zt1']
  iw=where(quantity eq quantities)
  if (iw[0] lt 0) then begin
     message,'quantity not element of quantity list:',/cont
     message,'['+quantities[0]+','+quantities[1]+']',/cont
     message,'returning',/cont
     return
  endif

  if not keyword_set(obstitle) then obstitle=' '

;code revision
  rev=strsplit(jncdf_glattread(rayfile,'svn_id'),' ',/extract)
  instrumr=rev[1]
  
; select and read  variable
  case quantity of
     'intensity': begin
        
        if jncdf_varexist(rayfile,'intensity') then begin
           var=transpose(jncdf_varread(rayfile,'intensity'),[2,1,0])
           wave=jncdf_varread(rayfile,'wavelength')
           var=reverse(var,2)
           btype='intensity'
           bunit='W m^2 Hz^-1 ster^-1'
        endif else begin
           message,'variable "intensity" does not exist in rayfile',/info
           return
        endelse
     end
     'zt1': begin

        if jncdf_varexist(rayfile,'tau_one_depth') then begin
           var=transpose(jncdf_varread(rayfile,'tau_one_depth'),[2,1,0])
           var=reverse(var,2)*m_to_mm
           iw=jncdf_varread(rayfile,'wavelength_indices')
           wave=(jncdf_varread(rayfile,'wavelength'))[iw]
           btype='z(tau=1)'
           bunit='Mm'
        endif else begin
           message,'variable "tau_one_depth" does not exist in rayfile',/info
           return
        endelse

     end
  endcase

; replace netcdf undefined values ( in practive for columns that did
; not converge) with NaNs. I'm not sure how netcdf defines this
; value, but in IDL on my MacPro its 9.96921e+36. I'm afraid
; checking exact equivalence might lead to roundoff trouble, so I
; check data for being larger than 1e35. In practive I'm sure
; this really signifies unconverged columns or other problems.
  iw=where(var gt 1e35)
  var[iw]=!values.F_Nan


; get bifrost time and isnap and x and y axes
  isnap=jncdf_glattread(rayfile,'snapshot_number')
  cd,bifrost_dir,current=calling_dir
  d=br_obj(bifrost_snapname,isnap=isnap)
  u=d->getunits()
  x=d->getx()
  y=d->gety()
  t=(d->gettsnap())*u.ut
  obj_destroy,d
  cd,calling_dir

; x and y spacing
  get_xy,input_dir+'/keyword.input',x_start,x_end,x_step,y_start,y_end,y_step,vacuum_to_air
  x0=x[x_start]
  y0=y[y_start]
  dx=abs(x[x_start+x_step]-x[x_start])
  dy=abs(y[y_start+y_step]-y[y_start])

; air to vacuum conversion
  if strcmp(vacuum_to_air,'true') then wave=airtovacuum(wave)
  dl=wave[1]-wave[0]


; hard coded for now, might be read form file later
  muz=1.0

; construct filename
  get_atom_names,input_dir+'/atoms.input',atom_names 
  Natoms=n_elements(atom_names)
  filename='RH_'+quantity+'_BIFROST_'+runname+'_'+string3(isnap)
  for i=0,Natoms-1 do filename+=('_'+atom_names[i])
  filename+='.fits'

; get run id
  if not keyword_set(runid) then get_runid,'rh',runid
  print,'rh_make_fits: setting runid to "'+runid+'"'

  mkhdr,hdr,var,/extend
  sxaddpar,hdr,'INSTRUME','RH'        ,' Data generated by RH',before='DATE'
  sxaddpar,hdr,'INSTRUMR',instrumr    ,' RH code revision number',before='DATE'
  sxaddpar,hdr,'BTYPE   ',btype       ,'Data variable',before='DATE'
  sxaddpar,hdr,'BUNIT   ',bunit       ,'Data unit',before='DATE'
  sxaddpar,hdr,'CDELT1  ',dx          ,' [Mm] x-coordinate increment',before='DATE'
  sxaddpar,hdr,'CDELT2  ',dy          ,' [Mm] y-coordinate increment',before='DATE'
  sxaddpar,hdr,'CDELT3  ',dl          ,' [nm] (non-uniform) wavelength increment',before='DATE'
  sxaddpar,hdr,'CRPIX1  ',1           ,' Reference pixel x-coordinate',before='DATE'
  sxaddpar,hdr,'CRPIX2  ',1           ,' Reference pixel y-coordinate',before='DATE'
  sxaddpar,hdr,'CRPIX3  ',1           ,' Reference pixel wavelength-coordinate',before='DATE'
  sxaddpar,hdr,'CRVAL1  ',x0          ,' [Mm] Position pixel 1 x-coordinate',before='DATE'
  sxaddpar,hdr,'CRVAL2  ',y0          ,' [Mm] Position pixel 1 y-coordinate',before='DATE'
  sxaddpar,hdr,'CRVAL3  ',wave[0]     ,' [nm] Position pixel 1 lambda-coordinate',before='DATE'
  sxaddpar,hdr,'CTYPE1  ','x'         ,' [Mm] Label for x-coordinate',before='DATE'
  sxaddpar,hdr,'CTYPE2  ','y'         ,' [Mm] Label for y-coordinate',before='DATE'
  sxaddpar,hdr,'CTYPE3  ','WAVE'      ,' [nm] Label for wavelength-coordinate',before='DATE'
  sxaddpar,hdr,'CUNIT1  ','Mm'        ,' Unit for x-coordinate',before='DATE'
  sxaddpar,hdr,'CUNIT2  ','Mm'        ,' Unit for y-coordinate',before='DATE'
  sxaddpar,hdr,'CUNIT3  ','nm'        ,' Unit for wavelength-coordinate',before='DATE'
  sxaddpar,hdr,'ELAPSED ',t           ,' [s] Time of snapshot',before='DATE'
  sxaddpar,hdr,'ISNAP   ',isnap       ,' isnap from Bifrost file',before='DATE'
  sxaddpar,hdr,'DATA_LEV',2           ,' Data level',before='DATE'
  sxaddpar,hdr,'MUX     ',0           ,' cosine of ray angle with x axis',before='DATE'
  sxaddpar,hdr,'MUY     ',0           ,' cosine of ray angle with y axis',before='DATE'
  sxaddpar,hdr,'MUZ     ',muz         ,' cosine of ray angle with z axis',before='DATE'
  sxaddpar,hdr,'NATOM   ',Natoms      ,' number of atoms in non-LTE',before='DATE'
  for i=0,Natoms-1 do begin
     atm='ATOM'+string(i)+'   '
     sxaddpar,hdr,atm   ,atom_names[i],' name of first non-LTE atom',before='DATE'
  endfor
  sxaddpar,hdr,'OBSTITLE',obstitle    ,' short run description',before='DATE' 
  sxaddpar,hdr,'RUNID   ',runid       ,' unique run ID',before='DATE' 
  sxaddpar,hdr,'ORIGIN  ','ITA/Oslo'  ,' Origin of data',before='DATE'
                                ;
  rh_idl2fits,input_dir+'/keyword.input',hdr
  sxaddpar,hdr,'COMMENT',$
           "Non-uniform wavelength axis",before='DATE'

  writefits,filename,var,hdr

  mkhdr,hdr,wave,/image
  sxaddpar,hdr,'EXTNAME ','wavelength axis',' Extension name'
  sxaddpar,hdr,'BTYPE   ','WAVE',' Data variable'
  sxaddpar,hdr,'BUNIT   ','nm',' Unit for wavelength'
  writefits,filename,wave,hdr,/append

  if keyword_set(cp_to_rep) then begin

;     construct dir to copy to
     repdir='/Volumes/SAM5/oslo_sim/'+runname+'/'+cp_to_rep
     runiddir='/Volumes/SAM5/oslo_sim/input/'+runid
; make sure repdir exists and runid dir does not exist
     if not file_exist(repdir) then begin
        message,'repdir does not exist, not moving data',/info
        return
     endif
     if file_exist(runiddir) then begin
        message,'runiddir exists already, not copying data',/info
        return
     endif

; copy fits file
     file_move,filename,repdir+'/'+filename,/overwrite

; copy input directory
     file_copy,input_dir,runiddir,/recursive,/overwrite
     
  endif

  if keyword_set(debug) then stop

end

pro get_xy,input_file,x_start,x_end,x_step,y_start,y_end,y_step,vacuum_to_air

  qarr=['x_start','x_end','x_step','y_start','y_end','y_step']

  openr,lur,input_file,/get_lun
  while (not eof(lur)) do begin
     text=' '
     readf,lur,text
     text=strtrim(text,2)
     if(strlen(text) gt 0 and strmid(text,0,1) ne '#') then begin
        ip=strpos(text,'=')
        if(ip lt 0) then message,'text does not have equals sign'
        var=strlowcase(strtrim(strmid(text,0,ip),2))
        value=strtrim(strmid(text,ip+1,strlen(text)-ip),2)
        iw=where(var eq qarr,count)
        if count gt 0 then succes=execute(text)
        if var eq 'vacuum_to_air' then begin
           vacuum_to_air=strlowcase(value)
        endif
     endif
  endwhile 
  free_lun,lur

end

pro get_atom_names,input_file,atom_names

  line=0
  Nmetal=0
  Nactive=0
  openr,lur,input_file,/get_lun
  while (not eof(lur)) do begin
     text=' '
     readf,lur,text
     text=strtrim(text,2)
     if(strlen(text) gt 0 and strmid(text,0,1) ne '#') then begin
        line++
        if line eq 1 then begin
           reads,text,Nmetal
           atom_names=strarr(Nmetal)
        endif else begin
           words=strsplit(text,' ',/extract)
           if words[1] eq 'ACTIVE' then begin
                                ; split words 0, the substring after
                                ; the last /, if present, is the atom
                                ; name
              pn=strsplit(words[0],'/',/extract)
              atom_names[Nactive]=pn[n_elements(pn)-1]
; strip .atom extension
              atom_names[Nactive]=strmid(atom_names[Nactive],0,strpos(atom_names[Nactive],'.atom'))
              Nactive++
           endif
        endelse
        
     endif
  endwhile 
  free_lun,lur
  atom_names=atom_names[0:Nactive-1]
end

FUNCTION airtovacuum, wave, TO_VACUUM_LIMIT=to_vacuum_limit                   
;+
; NAME:
;	AIRTOVACUUM
; PURPOSE:
;	Convert air wavelengths to vacuum wavelengths, i.e. correct 
;	for the index of refraction of air under standard conditions.  
;	Wavelength values below 200.0 nm will not be altered.  Uses the IAU 
;	standard for conversion given in Morton (1991 Ap.J. Suppl. 77, 119)
;
; CALLING SEQUENCE:
;	W_VAC = AIRTOVACUUM(WAVE [, /TO_VACUUM_LIMIT)
;
; INPUT/OUTPUT:
;	WAVE - Wavelength in Angstroms, scalar or vector
;		WAVE should be input as air wavelength(s), it will be
;		returned as vacuum wavelength(s).
;
; EXAMPLE:
;	If the air wavelength is W = 605.6125 (a Krypton line), 
;	then AIRTOVACUUM yields a vacuum wavelength of W = 605.78019
;
; METHOD:
;	See Morton (Ap. J. Suppl. 77, 119) for the formula used
;
; REVISION HISTORY
;	Written W. Landsman                November 1991
;       Revised H. Uitenbroek, Jul 23 1996
;-
  On_error, 2

  IF N_params() EQ 0 THEN BEGIN
    print,'Syntax - w_vac = AIRTOVACUUM(WAVE [, /TO_VACUUM_LIMIT)'
    print,'WAVE (Input) is the air wavelength in nm'
    print,'On output WAVE contains the vacuum wavelength in nm'
    return, 0.0
  ENDIF

  IF (keyword_set(TO_VACUUM_LIMIT)) THEN $
   limit = to_vacuum_limit $
  ELSE $
   limit = 199.9352

  sigma2 = (1.0D+7 / wave)^2
  fact = 1.0000834213D+0 + $
   2.406030D+6/(1.30D+10 - sigma2) + 1.5997D+4/(3.89D+9 - sigma2)
  fact = fact * (wave GE limit) + 1.0 * (wave LT limit)
  
  return, wave * fact
END

pro get_runid,code,runid
;
; create unique directory to store input data
;
  if not (code eq 'rh' or code eq 'm3d') then begin
     message,'Code descriptor not valid',/info
     stop
  endif

  input_dir='/Volumes/SAM5/oslo_sim/input/'
  files=file_search(input_dir+code+'*',count=nfiles)

  if nfiles gt 0 then begin
     
     id=intarr(nfiles)
     for i=0,nfiles-1 do begin
        num=0
        words=strsplit(files[i],'_',/extract)
        nwords=n_elements(words)
        reads,words[nwords-1],num
        id[i]=num
     endfor
     
     nid=max(id)+1
     runid=code+'_'+strtrim(string(nid),2)
     
  endif else begin
     
     runid=code+'_1'

  endelse

end

function jncdf_varread,file,varstring,debug=debug,count=count,offset=offset,stride=stride

  if n_params() lt 1 then begin
     print,'function jncdf_varread,file,var,count=count,offset=offset,stride=stride,debug=debug'
     return,0
  endif
  
; parse variable string
  words=strsplit(varstring,'/',count=ct,/extract)
  varname=words[ct-1]

  ingroup=0
; set groupnames if present
  if ct gt 1 then begin
     ingroup=1
     ngrp=ct-1
     grp=words[0:ct-1]
  endif

  ncid = NCDF_OPEN(file)

  if ingroup eq 0 then begin
     ncdf_varget,ncid,varname,var,count=count,offset=offset,stride=stride
  endif else begin

     current_id=ncid

; descend group tree
     for i=0,ngrp-1 do begin

        group_id=ncdf_groupsinq(current_id)
        group_id=[group_id]
        ngr=n_elements(group_id)
        for k=0,ngr-1 do begin
           group_name=ncdf_groupname(group_id[k])
           if group_name eq grp[i] then begin
              current_id=group_id[k]
              break
           endif
        endfor

     endfor

; we now have the id of the group that contains var, so read
     ncdf_varget,current_id,varname,var,count=count,offset=offset,stride=stride

  endelse

  if keyword_set(debug) then stop

  ncdf_close,ncid

  return,var

end

function jncdf_varexist,file,varstring,debug=debug
;
; function returns 1 if variable exist in ncdf file, else returns 0.
;

  if n_params() lt 1 then begin
     print,'function jncdf_varexist,file,var,debug=debug'
     return,0
  endif
  
; parse variable string
  words=strsplit(varstring,'/',count=ct,/extract)
  varname=words[ct-1]

  ingroup=0
; set groupnames if present
  if ct gt 1 then begin
     ingroup=1
     ngrp=ct-1
     grp=words[0:ct-1]
  endif

  ncid = NCDF_OPEN(file)

  if ingroup eq 0 then begin
     ve=ncdf_varid(ncid,varname)
  endif else begin

     current_id=ncid

; descend group tree
     for i=0,ngrp-1 do begin

        group_id=ncdf_groupsinq(current_id)
        group_id=[group_id]
        ngr=n_elements(group_id)
        for k=0,ngr-1 do begin
           group_name=ncdf_groupname(group_id[k])
           if group_name eq grp[i] then begin
              current_id=group_id[k]
              break
           endif
        endfor

     endfor

; we now have the id of the group that contains var, so 
     ve=ncdf_varid(current_id,varname)

  endelse

  if keyword_set(debug) then stop

  ncdf_close,ncid

  if ve eq -1 then return,0 else return,1
 
end
