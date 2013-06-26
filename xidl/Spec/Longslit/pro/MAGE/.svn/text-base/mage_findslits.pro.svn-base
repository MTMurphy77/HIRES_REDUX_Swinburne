function mage_findslits, traceflat, orderfile = orderfile

;  files = ['mage0114.fits','mage0115.fits','mage0116.fits','mage0117.fits','mage0118.fits']
  
;  for i=0, n_elements(files)-1 do begin
;     mage_proc, '../Raw/'+files[i], tmp, hdr=hdr
;     if i EQ 0 then begin
;        trcimg = tmp
;     endif else begin
;        trcimg += tmp
;     endelse;
;
;  endfor

  mage_proc, traceflat, trcimg, hdr=hdr
  tset_slits = mage_traceorders(trcimg, /chk)
  ordermask=mage_ordermask(tset_slits) 
   ;;slitmask = long_slits2mask(tset_slits)
   ;;ordermask = 0 * slitmask
   ;;ordermask[WHERE(slitmask GT 0)] = -slitmask[WHERE(slitmask GT 0)] + 21L

   if keyword_set(orderfile) then begin
      mwrfits, ordermask, orderfile, /create
      mwrfits, tset_slits, orderfile
   endif else begin
      file = "Orders.fits"
      mwrfits, ordermask, file, /create
      mwrfits, tset_slits, file
   endelse

   return, tset_slits

   print, "mage_findslits: all done!"

end
