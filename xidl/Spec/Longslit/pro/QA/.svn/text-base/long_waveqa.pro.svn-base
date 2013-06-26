;+
; NAME:
;   long_waveqa
;
; PURPOSE:
;   Generate QA for the wavelength solution.
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
PRO LONG_WAVEQA, wlines, fit, arc1d, rejpt, islit $
                 , fwhm, fwhmset, fwhmmask, qafile, panic = panic $
                 , BADFIT = BADFIT

angstrom = STRING(197B)
charsize = 1.1D
!p.multi = [0, 1, 2]
title_string = 'Slit#' + strcompress(string(islit+1), /rem)
IF KEYWORD_SET(PANIC) THEN title_string = title_string + '  PANIC'
IF KEYWORD_SET(BADFIT) THEN title_string = title_string + '  BADFIT'

ny = n_elements(arc1d)
nlines = n_elements(wlines)

wvfit  = x_calcfit(dindgen(ny), fitstr = fit)
linfit = x_calcfit(wlines.PIX, fitstr = fit)

ffit = 1.0*(*fit.ffit)
nrm = fit.nrm
func = fit.FUNC
nord = fit.NORD
; zero non-linear terms
IF func EQ 'POLY' AND nord GT 1 THEN ffit[0, 2:*] = 0.0 $
ELSE IF func EQ 'CHEBY' AND nord GT 2 THEN ffit[2:*] = 0.0

wvfit_lin = x_calcfit(dindgen(ny), func, ffit $
                      , nrm = nrm, nord = nord)
linfit_lin    = x_calcfit(wlines.PIX, func, ffit $
                          , nrm = nrm, nord = nord)
dwv = wvfit - shift(wvfit, 1)
dwv[0] = dwv[1]
dwv_med = djs_median(dwv)
rms = fit.rms/abs(dwv_med)

fwhmcolor = replicate('black', nlines)
IF (min(fwhmmask) EQ 0) THEN fwhmcolor[where(fwhmmask EQ 0)] = 'red'
traceset2xy, fwhmset, findgen(ny), fwhmfit


nl_wdev = wlines.WAVE-linfit_lin
all_nl_wdev = wvfit - wvfit_lin


yrange = [(0.05*djs_median(arc1d)) > 1.0, 50.0*max(arc1d)]
wrange = [min(wvfit), max(wvfit)]

djs_plot, wvfit, abs(arc1d), xrange = wrange $
          , yrange = yrange $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(' + angstrom + ')' $
          , ytitle = 'Arc Spectrum' $
          , title = title_string $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = charsize, /ylog
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'default' $
    ELSE this_color = 'red'
    w_string = strcompress(string(wlines[k].WAVE, FORMAT = '(F8.2)') $
                           , /rem)
    djs_xyouts, wlines[k].WAVE, 0.5*yrange[1] $
                , w_string, color = this_color $
                , charthick = 3.0, charsize = 0.7 $
                , orientation = 270.0
    djs_arrow, wlines[k].WAVE, 0.1*yrange[1], wlines[k].WAVE $
               , 0.03*yrange[1], /DATA, color = this_color, thick = 3.0 $
               , hsize = 0.0, /solid
ENDFOR

dwv_str = '\Delta\lambda = ' + strcompress(string(dwv_med, format = '(F7.4)') $
                                         , /rem) +  ' '  + angstrom
rms_str = 'RMS = ' +  $
  strcompress(string(rms, format = '(F7.4)'), /rem) + ' (pix)'
cen_str = '\lambda_{cen} = ' + string(wvfit[ny/2], FORMAT = '(F7.1)') $
  +  ' '  + angstrom
  
IF dwv_med LT 0 THEN xrange = [ny, 0] $
ELSE xrange = [0, ny]

djs_plot, findgen(ny), wvfit, xrange = xrange $
          , yrange = [0.98*min(wvfit), 1.02*max(wvfit)] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = 'pixel' $
          , ytitle =  '\lambda ' + '(' + angstrom + ')' $
          , title = 'Fit' $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = charsize

IF dwv_med LT 0 THEN xy_pos = 0.9 $
ELSE xy_pos = 0.1
djs_xyouts, xy_pos*ny, wrange[0] + (wrange[1]-wrange[0])*0.85 $
            ,  dwv_str, color = 'default' $
            , charthick = 3.0, charsize = charsize*1.2
djs_xyouts, xy_pos*ny,  wrange[0] + (wrange[1]-wrange[0])*0.75 $
            , rms_str, color = 'default' $
            , charthick = 3.0, charsize = charsize*1.2
djs_xyouts, xy_pos*ny,  wrange[0] + (wrange[1]-wrange[0])*0.65 $
            , cen_str, color = 'default' $
            , charthick = 3.0, charsize = charsize*1.2


FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'default' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].PIX], [wlines[k].WAVE] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = charsize, color = this_color, psym = 6
ENDFOR

djs_plot, wvfit, 0.0*findgen(ny), xrange = wrange $
          , yrange = [-6.0*rms, 6.0*rms] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(' + angstrom + ')' $
          , ytitle =  'Deviation [pix]' $ 
          , title = 'Residuals' $
          , color = 'default' $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = charsize
dwv_line = interpol(dwv, wvfit, wlines.WAVE)
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'default' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].WAVE], [(linfit[k]-wlines[k].WAVE)/dwv_line[k]] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = charsize, color = this_color, psym = 6
ENDFOR

;     Make plot of nonlinear part of wavelength solution
djs_plot, wvfit, all_nl_wdev/dwv $
          , xrange = wrange $
          , yrange = [min(all_nl_wdev/dwv)-2, max(all_nl_wdev/dwv) + 2] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = '\lambda ' + '(' + angstrom + ')' $
          , ytitle = 'x - x_{lin}' $
          , color = 'default', title = title_string $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 1.1
FOR k = 0L, nlines-1L DO BEGIN
    tind = WHERE(k EQ rejpt, nt)
    IF nt EQ 0 THEN this_color = 'default' $
    ELSE this_color = 'red'
    djs_oplot, [wlines[k].WAVE], [nl_wdev[k]/dwv_line[k]] $
               , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
               , charsize = charsize, color = this_color, psym = 6
ENDFOR

;     Plot fits for FWHM of arc lines
djs_plot, [wlines.PIX], [fwhm] $
          , xrange = xrange, yrange = [-0.1, 12.0] $
          , xstyle = 1, ystyle = 1 $
          , xtitle = 'Pixel', ytitle = 'FWHM' $
          , title = title_string, psym = 6 $
          , thick = 3.0, charthick = 3.0, xthick = 3.0, ythick = 3.0 $
          , charsize = 1.0, color = fwhmcolor
djs_oplot, findgen(ny), fwhmfit, thick = 3.0, color = 'blue'


!p.multi = 0
RETURN
END
