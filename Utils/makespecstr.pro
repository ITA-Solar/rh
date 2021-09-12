PRO makespecstr, sourcestr, output, ind=ind


linestr = { $
         iz: 0     ,$
         ion: 0    ,$
         ident: '', $
         ident_latex: '', $
         snote:'', $
         lvl1: 0   ,$
         lvl2: 0   ,$
         tmax: 0.  ,$
         wvl:  0.d ,$
         flag: 0, $
         goft:  dblarr(n_elements(sourcestr.IONEQ_LOGT))  }
if (n_elements(ind) ne 0) then begin
  linestr=replicate(linestr,n_elements(sourcestr.lines[ind]))
  struct_assign,sourcestr.lines[ind],linestr
endif else begin
  linestr=replicate(linestr,n_elements(sourcestr.lines))
  struct_assign,sourcestr.lines,linestr
endelse

output = {lines:linestr, $
          ioneq_logt:sourcestr.ioneq_logt,$
          ioneq_name:sourcestr.ioneq_name,$ 
          ioneq_ref:sourcestr.ioneq_ref, $
          wvl_limits: sourcestr.wvl_limits, $
          model_file:sourcestr.model_file,$ 
          model_name:sourcestr.model_name,$
          model_ne:sourcestr.model_ne,$
          model_pe:sourcestr.model_pe,$
          model_te:sourcestr.model_te, $
          wvl_units: sourcestr.wvl_units, $
          int_units: sourcestr.int_units, $
          add_protons:sourcestr.add_protons, $
          date: sourcestr.date, $
          version:sourcestr.version, $
          photoexcitation:sourcestr.photoexcitation, $
          rphot: sourcestr.rphot, $
          radtemp: sourcestr.radtemp,$
          nlines: n_elements(linestr)}
END
