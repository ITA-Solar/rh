PRO chianti_gofnt_comp
  
bwl = 1d0
ewl = 10000d0
nd = 21
logd = dindgen(nd)*.5 + 6d0

;ioneq = 'radyn_qsslht_basicphot_t0_intpl.ioneq'
;lioneq = 'radyn_val3c_t0_intpl.ioneq'
;ioneq = '/Users/gskerr1/Applications/SSW/packages/chianti/dbase/ioneq/chianti.ioneq'
ioneq = '/data2/gskerr/SSW/packages/chianti/dbase/ioneq/chianti.ioneq'

goft = []
for i =0, nd-1 do begin

  print, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  print, ' '
  print, ' WORKING ON DENSITY #', i
  print, ' '
  print, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

  ch_synthetic, bwl, ewl, output=specstr, density= 10d0^logd[i], ioneq_name=ioneq,$ 
                          /goft, rphot = 1d0, /verbose
  save,file = 'chianti_gofnt_out.'+strng(i)+'.sav',specstr, logd  
  goft = [goft,specstr]

endfor

chianti_gofnt = goft
save, file = 'gofnt_chianti_out.sav',chianti_gofnt,logd

END
