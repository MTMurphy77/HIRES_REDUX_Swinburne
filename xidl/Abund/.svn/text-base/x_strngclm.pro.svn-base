;+ 
; NAME:
; x_strngclm
;   Version 1.1
;
; PURPOSE:
;    Convert column density structure to a string for Latex tables
;
; CALLING SEQUENCE:
;   clmc = x_strngclm(dla.ion[6].state[1], /LOG, /INDIV, MIN=)
;
; INPUTS:
;   Colm sturcture
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /LOG -- Input structure has linear wavelength values (convert to
;          LOG)
;  /INDIV -- Columns are ionic transitions, not combined values
;  MIN=  -- Minimum error for the column density
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   05-Oct-2004 Written by JXP
;   24-Aug-2011 Added extra flgclm option for CII* MR
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_strngclm, struct, LOG=log, INDIV=indiv, MIN=min

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'strng = x_strngclm(struct) [v1.0]'
    return, ''
  endif 

  ;; Log?
  if keyword_set(LOG) or struct.clm GT 1e2 then $
    x_logclm, struct.clm, struct.sigclm, logN, logS $
  else begin
      logN = struct.clm
      logS = struct.sigclm
  endelse

  ;; Minimum error
  if keyword_set(MIN) then logS = logS > MIN

  if keyword_set(INDIV) then begin
      case struct.flgclm of
          0: lin = ''
          1: lin = ''
          2: lin = '>'
          3: lin = '>'
          4: lin = '<'
          5: lin = '<'
          6: lin = '<'
          8: lin = ''
          9: lin = ''
          12: lin = ''
          13: lin = '<'
          else: stop
      endcase
  endif else begin
      case struct.flgclm of
          0: lin = ''
          1: lin = ''
          2: lin = '>'
          3: lin = '<'
          else: stop
      endcase
  endelse

  ;; Zero value
  if struct.clm EQ 0. then return, lin

  lin = lin + string(logN, FORMAT='(f6.2)' )
  if (struct.flgclm MOD 8) LE 1 then begin
      lin = lin + ' \pm '
      lin = lin + string( logS, format='(f4.2)')
  endif

  return, lin
end
