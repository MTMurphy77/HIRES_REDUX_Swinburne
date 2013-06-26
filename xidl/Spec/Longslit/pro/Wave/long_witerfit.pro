;+
; NAME:
;   long_witerfit
;
; PURPOSE:
;   Iterate on wavelength solution
;
; CALLING SEQUENCE:
;
; INPUTS:
;   
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Mar-2005  Written by S. Burles (MIT), David Schlegel (LBL), and 
;                Joe Hennawi (UC Berkeley)
;-
;------------------------------------------------------------------------------
FUNCTION long_witerfit, arc1d, lines, fit0, wstruct, shft = shft1 $
                        , gdfit = gdfit, rejpt = rejpt $
                        , DEBUG = DEBUG, FIT_FLAG = FIT_FLAG

fit_flag = 0
pkwdth = wstruct.pkwdth
toler = wstruct.toler
npanic = wstruct.npanic
nord_panic = wstruct.nord_panic
psig = wstruct.PSIG
mxoff = wstruct.MXOFF
sigrej = wstruct.sigrej
nord = wstruct.nord
flg_qual = wstruct.FLG_QUAL
fweight=wstruct.FWEIGHT


IF NOT KEYWORD_SET(FORDR) THEN FORDR = 9L
goodind = WHERE(nord NE 0, niter)

fit_last = fit0
FOR kk = 0L, niter-1L DO BEGIN
    ;; Grab new lines
    lines.flg_plt = 0
    ;; Re-identify using this fit
    IF (kk EQ 0) AND KEYWORD_SET(shft1) THEN BEGIN
       IF abs(shft1) GT 0.01 THEN shft = shft1 
    ENDIF ELSE shft = 0
    x_templarc, arc1d, lines, fit_last,  FORDR = FORDR $
                , PKSIG = psig[kk], MXOFF = MXOFF[kk], FLG = FLG_TEMPL $
                , PKWDTH = pkwdth, toler = toler, shft = shft $
                , MAXQUAL = FLG_QUAL[kk], /THIN, FWEIGHT=FWEIGHT
    ;if flg_templ EQ -1 then begin
    ;    splog, 'Insufficient lines for AUTO!!'
    ;    stop
    ;endif
    ;; Check the number of good lines
    gdfit = where(lines.flg_plt EQ 1, ngd)
    ;; Perform a fit
    fit_now = fit_last
    fit_now.hsig = sigrej[kk]
    fit_now.lsig = sigrej[kk]
    fit_now.flg_rej = 1
;   fit_now.maxrej = ngd/2 ??
;   fit_now.minpt = ngd/2  ??
    fit_now.nord = nord[kk]
    IF (ngd GT 0 AND ngd LT npanic) AND (kk EQ niter-1L) THEN BEGIN
        splog, 'ngd=', strcompress(ngd, /rem), ' < npanic=' $
               , strcompress(npanic, /rem) 
        splog, 'Reverting to lower order fit nord_panic=' $
               , strcompress(nord_panic, /rem)
        stop
        fit_now.nord = nord_panic
    ENDIF
    IF ngd GE 3 THEN BEGIN
        wave_fit = x_fitrej(lines[gdfit].pix, lines[gdfit].wave $
                            , FITSTR = fit_now, REJPT = rejpt) 
        ;; don't crash on failure of Autoid JFH 05/08
        IF wave_fit[0] EQ -1 then begin
            splog, 'AUTO Failed!!'
            SPLOG, 'This could be a bad slit. Setting wavelenghts to zero'
;            stop
            *fit_now.ffit = 0.0*(*fit_now.ffit)
            fit_flag = 0
            rejpt = -1L
        ENDIF ELSE fit_flag = 1
    ENDIF ELSE BEGIN
        splog, 'AUTO Failed 2!!'
        SPLOG, 'This could be a bad slit. Setting wavelenghts to zero'
        stop
        *fit_now.ffit = 0.0*(*fit_now.ffit)
        fit_now.nord = n_elements(*fit_now.ffit)
        fit_flag = 0
        rejpt = -1L
        BREAK ;; break out of the loop and continue on to next slit
    END
    IF KEYWORD_SET(DEBUG) THEN BEGIN
        print, '*********FITTED**********'
        print, 'rms=', strcompress(string(fit_now.rms, '(F6.2)'), /rem)
        forprint, lines[gdfit].pix, lines[gdfit].WAVE, textout = 2
        IF rejpt[0] NE -1 THEN BEGIN
            print, '*********REJECTED**********'
            forprint, lines[gdfit[rejpt]].pix, lines[gdfit[rejpt]].WAVE $
                      , textout = 2 
        ENDIF
        IF fit_now.rms GT 1.0 THEN STOP
    ENDIF
    fit_last = fit_now
ENDFOR

ny = n_elements(arc1d)
nyby2 = ny/2L
;if fit_flag EQ 0 then stop
wave_vec = x_calcfit(dindgen(ny), fitstr = fit_now)
wave_cen = wave_vec[nyby2]
disp_cen = wave_vec[nyby2]-wave_vec[nyby2-1L]
;; Add central wavelength and central dispersion to fit struct
IF TAG_EXIST(fit_now, 'WAVE_CEN') EQ 0 THEN $
  fit_now = struct_addtags(fit_now, create_struct('WAVE_CEN', wave_cen)) $
ELSE fit_now.WAVE_CEN = wave_cen 
IF TAG_EXIST(fit_now, 'DISP_CEN') EQ 0 THEN $
  fit_now = struct_addtags(fit_now, create_struct('DISP_CEN', disp_cen)) $
ELSE fit_now.DISP_CEN = disp_cen

RETURN, fit_now
END
