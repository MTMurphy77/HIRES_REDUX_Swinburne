;+ 
; NAME:
; x_cumulative   
;    Version 1.1
;
; PURPOSE:
;  Calculates a histogram of sorts [defunct]
;
; CALLING SEQUENCE:
;   
; 
;
; INPUTS:
;   Array, binsize or nbin
;
; RETURNS:
;   cumulative array 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  XARR -- Value of xarray
;
; COMMENTS:
;
; EXAMPLES:
;   cumul = x_cumulative(array, 0.05)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Dec-2004 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_cumulative, arr, bin, NBIN=nbin, XMNX=xmnx, XARR=xarr

  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'cumul = x_cumulative(arr, [bin], NBIN=, XMNX=) [v1.1]'
      return, -1
  endif 

  if not keyword_set(XMNX) then begin
      xmnx = dblarr(2)
      xmnx[0] = min(arr, max=mx)
      xmnx[1] = mx
  endif

  ;; BIN or NBIN?
  if not keyword_set(bin) then begin
      if not keyword_set(NBIN) then begin
          print, 'x_cumlative: Must set bin or NBIN'
          return, -1
      endif
      ;; Set
      bin = (xmxn[1] - xmnx[0])/float(nbin)
  endif 
  if not keyword_set(NBIN) then nbin = (xmnx[1] - xmnx[0])/float(bin)
  xarr = dblarr(nbin)

  xarr = xmnx[0] + bin*dindgen(nbin)
  cumul = lonarr(nbin)
  for jj=0L,nbin-1 do begin
      ;; SDSS
      gd = where(arr GE xarr[jj], ngd)
      cumul[jj] = ngd 
  endfor
  return, cumul

end
