pro rh_delete_fits,f,force=force
;
; deletes an RH fits file from the oslo_sim repository and removes the
; corresponding runid directory
;
; Will ask for confirmation. Suppress confirmation by setting /force
;

if n_params() ne 1 then begin
   print,'pro rh_delete_fits,f,force=force'
   return
endif

if not file_exist(f) then begin
   message,'File "'+f+'" does not exist, returning.',/info
   return
endif

hdr=fheader(f)
runid=fxpar(hdr,'RUNID')

repdir='/Volumes/SAM5/oslo_sim/input/'+runid
if not file_exist(repdir) then begin
   message,'The RUNID directory "'+repdir+'" does not exist, returning.',/info
   return
endif

if not keyword_set(force) then begin
   message,'Do you want to delete the following file and directory:',/info
   message,'file:      "'+f+'"'
   message,'directory: "'+repdir+'"'
   message,'Press y to confirm, any other key to abort.',/info
   if get_kbrd() ne 'y' then begin
      message,'File deletion aborted',/info
      return
   endif
endif

file_delete,f
file_delete,runid,/recursive

end
