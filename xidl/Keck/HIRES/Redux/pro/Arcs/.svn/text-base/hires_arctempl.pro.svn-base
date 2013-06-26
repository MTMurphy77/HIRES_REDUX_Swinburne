;+ 
; NAME:
; hires_arctempl
;     Version 1.1
;
; PURPOSE:
;   Identify the best matching archived wavelength solution 
;
; CALLING SEQUENCE:
;   
;  templfil = hires_arctempl(hires, ordrs)
;
; INPUTS:
;   hires  -  HIRES structure
;   [ordrs]  -  Orders for fitting
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Feb-2005 Created by JXP 
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function hires_arctempl, hires, ordrs

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'templfil = hires_arctempl(hires, [ordrs]) [v1.1]'
      return,-1
  endif 
  
  ;;  Optional Keywords

  ;; Read the tempalte database
  readcol, getenv('XIDL_DIR')+ $
    '/Keck/HIRES/Redux/pro/Arcs/Templates/hires_templ.dat', $
    filnm, chip, xdisp, echa, xda, deck, rbin, cbin, ordi, ordf, $
    FORMAT='A,I,A,D,D,A,I,I,L,L', skipline=1

  ;; Cut on Cross
  gd = where(strcompress(hires.cross, /rem) EQ xdisp)
  chip = chip[gd]
  xdisp = xdisp[gd]
  echa = echa[gd]
  xda = xda[gd]
  deck = deck[gd]
  rbin = rbin[gd]
  cbin = cbin[gd]
  ordi = ordi[gd]
  ordf = ordf[gd]
  filnm = filnm[gd]

  ;; Max order
  maxo = max(ordf)

  ;; Closest XDANGL irrespective of binning
  if not keyword_set(ORDRS) then begin
      ;; Best = EDANGL LT 0.05 and XDANGL LT 0.2
      mtch = where(abs(hires.echangl-echa) LT 0.05 AND $
                   abs(hires.xdangl-xda) LT 0.2 AND $
                   hires.chip EQ chip, nclos)
      if nclos EQ 0 then begin
          mtch = where(abs(hires.echangl-echa) LT 0.05 AND $
                       abs(hires.xdangl-xda) LT 0.4 AND $
                       hires.chip EQ chip, nclos)
          if nclos EQ 0 then begin
              mtch = where(hires.chip EQ chip, nclos)
          endif
      endif 
  endif else begin
      ;; Specified orders
      gdo = where(min(ordrs) GE ordi AND $
                  max(ordrs) LE ordf, ngdo)
      if ngdo EQ 0 then begin
         mtch = where(hires.cross EQ xdisp, nmt)
         gdo = where(min(ordrs) GE ordi[mtch] AND $
                     (max(ordrs) < maxo) LE ordf[mtch], ngdo)
         if ngdo EQ 0 then begin
            print, 'hires_arctempl:  No archived wavelengths fitting your' + $
                   'setup.  Contact JXP ASAP (xavier@ucolick.org)!'
            stop
         endif
         mtch = mtch[gdo]
      endif else mtch = gdo 
   endelse

  ;; Improve the match from binning
;  gd = where(hires.colbin EQ cbin[mtch], ngd)
;  if ngd NE 0 then mtch = mtch[gd]
;  gd = where(hires.rowbin EQ rbin[mtch], ngd)
;  if ngd NE 0 then mtch = mtch[gd]
  
  ;; Closet ECHANGL
  mn = min(abs(hires.echangl-echa[mtch]),imn)
  allx = where(abs(echa[mtch]-echa[mtch[imn]]) LE 0.001, nax)

  if nax NE 1 then begin  ;; Close XDANGL
      mn2 = min(abs(hires.xdangl-xda[mtch[allx]]),imne)
      idx = mtch[allx[imne]]
  endif else idx = mtch[allx[0]]


  ;; Final file
  finfil = getenv('HIRES_CALIBS')+'/ARCS/'+filnm[idx]

  ;; Check
  if x_chkfil(finfil,/silent) EQ 0 then begin
      ;; 
      print, 'hires_arctempl: You need ', finfil
      print, 'hires_arctempl: Odds are you have not grabbed the Arc templates.'
      print, 'hires_arctempl: Please do so here'
      print, 'http://www.ucolick.org/~xavier/HIRedux/CALIBS/index.html'
      stop
  endif

  return, finfil
      
end
