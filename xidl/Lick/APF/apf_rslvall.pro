;+ 
; NAME:
; apf_rslvall  
;     Version 1.1
;
; PURPOSE:
;   Resolves a number of key codes for HIRES redux
;
; CALLING SEQUENCE:
;   
;  apf_rslvall
;
; INPUTS:
;
; RETURNS:
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
;   apf_rslvall
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro apf_rslvall

  compile_opt strictarr
  ;; Subbias
  resolve_routine, 'apf_subbias', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'apf_proc', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Arcs
  resolve_routine, 'apf_procarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'apf_fitarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'apf_fit2darc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'apf_tracearc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'apf_fittrcarc', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'apf_mkaimg', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'apf_arcxyoff', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'apf_allarc', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'apf_arcalign', /no_recompile, /either, /COMPILE_FULL_FILE
  ;; Flats
  resolve_routine, 'x_trcflat', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_fntobj', /no_recompile, /either, /COMPILE_FULL_FILE
  resolve_routine, 'x_ordermask', /no_recompile, /either, /COMPILE_FULL_FILE
;  resolve_routine, 'slitflat_qa', /no_recompile, /either, /COMPILE_FULL_FILE

  ;; Other
  ;resolve_routine, 'apf_slitlen', /no_recompile, /either, /COMPILE_FULL_FILE

end
