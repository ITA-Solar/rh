 
function starts_with_numeral,st
  fc=strmid(strtrim(st,2),0,1)
  iw=where(fc eq ['1','2','3','4','5','6','7','8','9','0'],count) 
  if count eq 1 then return, 1 else return,0
end


pro rh_idl2fits,inputfile,hdr ,debug=debug
;
;+
; NAME:
;	RH_IDL2FITS
;
; PURPOSE:
;       add RH input file variables to hdr
;
; CATEGORY:
;       RH15D
;	
; CALLING SEQUENCE:
;	RH_IDL2FITS,inputfile,hdr
;
; INPUTS:
;	inputfile - name of input file
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OPTIONAL KEYWORD PARAMETERS:
;	debug - add debug printout
;
; OUTPUTS:
;	adds input file variables to hdr
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
;
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
; $Id:$
;-
;
if(n_params() lt 2) then begin 
  message,'Syntax: rh_idl2fits,inputfile,hdr',/info
  return
endif

; argument checking

ff=file_search(inputfile,count=count)
if(count eq 0) then begin
  message,'file not found: '+inputfile,/info
  return
endif
if(n_elements(hdr) eq 0) then begin
  message,'hdr not given',/info
  return
endif

; set translation of variables with long names
rh_idl2fits_varrep,varrep

sxaddpar,hdr,'COMMENT',"Variables from multi3d.input file:",before='DATE'
openr,lur,inputfile,/get_lun
text=''
nid=0
id=strarr(1000)  ; array to store Id strings
while (not eof(lur)) do begin

  readf,lur,text
  text=strtrim(text,2)

  if(strlen(text) gt 0) then begin ; skip blank lines
   
     i=strpos(text,'; $Id:')
     if(i eq 0) then begin
        id[nid]=strmid(text,1,strlen(text)-1)
        nid=nid+1
     endif

     i=strpos(text,'#')
     if(i gt 0) or (i lt 0) then begin      ; remove comment line
        if(i gt 0) then text=strmid(text,0,i) ; remove comment
        ip=strpos(text,'=')
        if(ip lt 0) then message,'text does not have equals sign'
        var=strlowcase(strtrim(strmid(text,0,ip),2))
        value=strtrim(strmid(text,ip+1,strlen(text)-ip),2)
        if not starts_with_numeral(value) then value='"'+value+'"' ; make string of non-numerical values
        if(strpos(value,'"') eq 0) then value='"'+strtrim(strmid(value,1,strlen(value)-2),2)+'"' ; trim blanks from string value
        from=string(' ',format='(a30)') ; original variable name in .idl file
        if(strlen(var) gt 8) then begin
           iw=where(var eq strtrim(varrep[0,*],2))
           iw=iw[0]
           if(iw lt 0) then begin ; replacement not found
              message,'replacement not found for long variable name: '+var,/info
          if(keyword_set(debug)) then stop
        endif else begin
          var=strtrim(varrep[1,iw],2)
          strput,from,strtrim(varrep[0,iw],2),0
        endelse
      endif else begin
        strput,from,var,0
      endelse
      ip1=strpos(value,'"')
      ip2=strpos(value,"'")
      if(ip1 eq 0) or (ip2 eq 0) then begin

        if(keyword_set(debug)) then print,'string ',var,'=',value
        value=strmid(value,1,strlen(value)-2)   ; text value

      endif else begin
        ip=strpos(value,'.')
        if(ip gt 0) then begin
          if(keyword_set(debug)) then print,'float  ',var,'=',value
          value=float(value)
        endif else begin
          if(keyword_set(debug)) then print,'long   ',var,'=',value
          if(var ne 'niter_mg') then value=long(value)
        endelse
      endelse
     
        sxaddpar,hdr,var,value,' '+from,before='DATE'
   
     endif

  endif ; blank line if

endwhile ;eof loop

for i=0,nid-1 do begin
  sxaddpar,hdr,'COMMENT',strtrim(id[i],2),before='DATE'
endfor

free_lun,lur

end
   
