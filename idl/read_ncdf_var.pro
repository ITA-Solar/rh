FUNCTION read_ncdf_var, infile, varname, groupname=groupname
;
;+
; NAME:
;    read_ncdf_var
;
; PURPOSE:
;    Read a variable from a netCDF file.
;
; CALLING SEQUENCE:
;    data = read_ncdf_var(infile, varname, groupname)
;
; INPUTS:
;    infile - string with file name.
;    varname - string with variable name. Case sensitive.
;    
; KEYWORD PARAMETERS:
;    groupname - string with group name. Case sensitive. Will give an error
;                if group does not exist.
;
; OUTPUTS:
;    data - array with variable.
; -
;
  ; Ensure that a file is specified.
  IF N_PARAMS() LT 2 THEN $
    MESSAGE, 'Incorrect number of arguments.'
  
  ; Open the NetCDF file for reading.
  ncid = NCDF_OPEN(infile)
  IF (ncid EQ -1) THEN $
    MESSAGE, 'Error opening file: ' + infile
  if (n_elements(groupname) eq 0) then begin
    gid = ncid
  endif else begin
    groups = ncdf_groupsinq(ncid)
    if (groups[0] eq -1) then message, 'Group '+strtrim(groupname, 2)+' not found.'
    n_groups = n_elements(groups)
    gid = -1
    for i=0, n_groups -1 do begin
        if (ncdf_groupname(groups[i]) eq groupname) then gid = groups[i]
    endfor
    if (gid eq -1) then message, 'Group '+strtrim(groupname, 2)+' not found.'
  endelse
  
  varid = ncdf_varid(gid, varname)
  ncdf_varget, gid, varid, data
  
  ncdf_close, ncid
  
  return, data
  
END