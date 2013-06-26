;+ 
; NAME: 
; igm_calc_lox
;    Version 1.1
;
; PURPOSE:
;    Calculate l(X) given an f(N,X) over an N_HI interval
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;  l(X)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   July-2011 Written by JXP
;-
;------------------------------------------------------------------------------
function igm_calc_lox, fn_strct, z, NHI_min, NHI_max, NEVAL=neval

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'lX = igm_calc_lox(fn_strct, z, NHI_min, [NHI_max]) [v1.0]'
    return, -1
  endif 

  if not keyword_set(NEVAL) then neval = 10000L
  if not keyword_set(NHI_MAX) then begin 
     NHI_MAX = 23.
     infinity=1
  endif

  ;; Brute force (should be good to ~0.5%)
  lgNHI = NHI_min + (NHI_MAX-NHI_MIN)*dindgen(neval)/(neval-1)
  dlgN = lgNHI[1]-lgNHI[0]
  
  ;; Evaluate f(N,X)
  lgfNX = eval_powerfn(fn_strct, lgNHI, z)

  ;; Sum
  lX = total(10.d^(lgfNX+lgNHI)) * dlgN * alog(10.)

  ;; Infinity?
  if keyword_set(INFINITY) then begin
     NEVAL2 = 1000L
     lgNHI = NHI_max + (99.-NHI_MAX)*dindgen(neval2)/(neval2-1)
     dlgN = lgNHI[1]-lgNHI[0]
     lgfNX = eval_powerfn(fn_strct, lgNHI, z)
     lX2 = total(10.d^(lgfNX+lgNHI)) * dlgN * alog(10.)
     ;; 
     lX = lX + lX2
  endif

  return, lX
end
