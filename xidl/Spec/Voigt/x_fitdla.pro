;+ 
; NAME:
; x_fitdla   
;   Version 1.11
;
; PURPOSE:
;    GUI used to fit DLA profiles interactively.  The user determines
;    the continuum at the same time.
;
; CALLING SEQUENCE:
;  x_fitdla, flux_fil, [err_fil], XSIZE=, YSIZE=, TITLE=, 
;               WAVE=, /BLOCK, FWHM=, INIFIT=, INFLG=
;
; INPUTS:
;   flux_fil   - FITS file (or array) containing flux
;   [err_fil]  - FITS error array
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   xsize      - Draw window xsize (pixels)
;   ysize      - Draw window ysize (pixels)
;   wave       - wavelength array 
;   FWHM=      - Resoltuion (FWHM in pixels)  [default: 2]
;   INFLG=     - Usual flag for reading the FITS file (see x_readspec)
;   INIFIT=    - IDL file containing a saved version of the fit
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_fitdla, 'spec.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  XGETX_PLT
;  XGETY_PLT
;  XGETXPIX_PLT
;  XGETYPIX_PLT
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;   02-Aug-2009 Mendez altered the matrix
;   07-Sep-2011 added keyword old_init to deal with old data (MN).
;   07-Sep-2011 added keyword smooth=# which smooths the data by # pixels (MN)
;-
;------------------------------------------------------------------------------

;;;;
; Events
;;;;

pro x_fitdla_event, ev


  WIDGET_CONTROL, ev.top, get_uvalue = state, /no_copy
  WIDGET_CONTROL, ev.id, get_uvalue = uval
  flg_plt = 0

  case uval of
      'LYB': begin
          x_fitdla, state.fx, state.sig, wave=state.wave, FWHM=state.fwhm, $
                    ZIN=state.lines[0].zabs, /BETA, /BLOCK, INFLG=4
      end
      'CRUDEHI': begin
          state.crude_val[1] = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CRUDELO': begin
          state.crude_val[0] = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CRUDE': begin
          state.crude = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'EXACT': begin
          state.exact = ev.value
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'ZABS' : begin
          widget_control, state.zabs_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].zabs = tmp
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'CONTI' : begin
          widget_control, state.conti_id, get_value=tmp
          state.conti = tmp
      end
      'BVAL' : begin
          widget_control, state.bval_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].b = tmp
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'NCOLM' : begin
          widget_control, state.Ncolm_id, get_value=tmp
          gdset = where(state.lines[0:state.nlin-1].set EQ state.curlin,ngd)
          state.lines[gdset].N = tmp
          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
      end
      'ERRORB' : widget_control, state.error_msg_id, set_value=''
      'DRAW_BASE' : begin
          if ev.enter EQ 0 then begin ; Turn off keyboard
              widget_control, state.text_id, sensitive = 0
          endif
          WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
          return
      end
      'DRAW' : begin
          widget_control, state.text_id, /sensitive, /input_focus ; Turn on keyb
          case ev.type of
              0 : begin ; Button press
                  case ev.press of
                      1 : begin     ; Left button = Zoom
                          x_speczoomreg, state
                          if state.flg_zoom NE 2 then begin
                              WIDGET_CONTROL, state.base_id, $
                                set_uvalue = state,  /no_copy
                              return
                          endif else state.flg_zoom = 0
                      end 
                      4 : begin
                          state.lines[state.curlin].zabs = $
                            (state.xpos / state.wrest) - 1.
                          x_fitline_updfit, state, FLG_PLT=flg_plt, EXACT=state.exact
                          widget_control, state.zabs_id, $
                            set_value=state.lines[state.curlin].zabs
                      end
                      else: 
                  endcase
              end
              1 : begin ; Button Release
                  WIDGET_CONTROL, state.base_id, set_uvalue = state,  /no_copy
                  return
              end
              2 : begin ; Motion event
                  state.xcurs = ev.x
                  state.ycurs = ev.y
                  state.xpos = xgetx_plt(state, /strct)
                  state.ypos = xgety_plt(state, /strct)
                  widget_control, state.xpos_id, set_value=state.xpos
                  widget_control, state.ypos_id, set_value=state.ypos
                  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
                  return
              end
          endcase
      end
      'TEXT' : begin
          eventch = string(ev.ch)
          if eventch EQ  'q' then begin
              widget_control, ev.top, /destroy
              return
          end
          ;; FIT LINE
          x_fitline, state, eventch, FLG_PLT=flg_plt
      end
      ;; BUTTONS
      'IDLOUT' : x_fitline_idlout, state, /ERR
      'SPECPLT' : begin
          if state.nlin NE 0 then zin = state.lines[state.curlin].zabs $
            else zin = 0
          x_specplot, state.fx, state.sig, wave=state.wave, inflg=4, $
            zin=zin, /lls, /block
      end
      'DONE' : begin
         widget_control, ev.top, /destroy
         return
      end
  endcase

  ;; PXMNX
  if state.xymnx[0] NE state.old_xymnx[0] OR $
    state.xymnx[2] NE state.old_xymnx[2] then begin
      state.old_xymnx = state.xymnx
      state.xpmnx = x_getxpmnx(state)
  endif

; Update Plot
  if (flg_plt MOD 2) EQ 1 then x_fitdla_UpdatePlot, state
  ; Return state to Base
  WIDGET_CONTROL, state.base_id, set_uvalue = state, /no_copy
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;  Plot
;;;;;;;;;;;;;;;;;;;;


pro x_fitdla_UpdatePlot, state
  
; Plot Data

  if state.psfile NE 1 then begin
      widget_control, state.draw_id, get_value=wind
      wset, wind
  endif

  clr = getcolor(/load)

;  if state.flg_dum EQ 1 then stop
;  state.flg_dum = 0
;  stop
  plot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.fx[state.xpmnx[0]:state.xpmnx[1]], psym=state.psym, $
    position=state.pos, xrange=[state.xymnx[0],state.xymnx[2]], $
    yrange=[state.xymnx[1],state.xymnx[3]], xstyle=1, ystyle=1, $
    xtitle='!17Wavelength', ytitle='Flux', $
    title=state.title, $
    background=clr.white, $
    color=clr.black, $
    xcharsize=1.7, $
    ycharsize=1.7

  ; Plot Error array
  if state.flg_sig EQ 1 then begin
    oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.sig[state.xpmnx[0]:state.xpmnx[1]], psym=state.psym, color=clr.red
	oplot, [-10000., 10000.], [0.,0.], color=clr.black, linestyle=5
  endif

  ;; FIT
  if state.nlin NE 0 then begin
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
        state.fit[state.xpmnx[0]:state.xpmnx[1]], color=clr.green
      ;; Mark all lines
      for i=0L,state.nlin-1 do begin
          oplot, replicate( (state.lines[i].zabs+1.)*$
                            state.lines[i].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.blue, linestyle=1
      endfor
      ;; Mark current set
      gdlin = where(state.lines[0:state.nlin-1].set EQ state.curlin, ngd)
      for i=0L,ngd-1 do begin
          oplot, replicate( (state.lines[gdlin[i]].zabs+1.)*$
                            state.lines[gdlin[i]].wrest, 2), $
            [state.xymnx[1], state.xymnx[3]], color=clr.red, linestyle=1
      endfor
  endif

; CONTINUUM
  ;; Line
  oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
    state.conti[state.xpmnx[0]:state.xpmnx[1]], color=clr.purple, $
    linestyle=1
  ;; Points
  if state.cstr.npts NE 0 then begin
      gdc = where(state.cstr.msk EQ 1)
      oplot, [state.cstr.xval[gdc]], [state.cstr.yval[gdc]], psym=1, $
        color=clr.orange, symsize=5
  endif

; Crude ERROR
  if state.crude NE 0 and state.nlin NE 0 then begin
      ;; HIGH
      tmp = state.lines[0:state.nlin-1]
      mx = max(tmp.N, imx)
      tmp[imx].N = tmp[imx].N + state.crude_val[1]
      crude_hi = x_allvoigt(state.wave, tmp, SIGMA=state.FWHM/2.35)
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
        crude_hi[state.xpmnx[0]:state.xpmnx[1]], color=clr.red, linestyle=2
      ;; LOW
      tmp = state.lines[0:state.nlin-1]
      tmp[imx].N = tmp[imx].N + state.crude_val[0]
      crude_lo = x_allvoigt(state.wave, tmp, SIGMA=state.FWHM/2.35)
      oplot, state.wave[state.xpmnx[0]:state.xpmnx[1]], $
        state.conti[state.xpmnx[0]:state.xpmnx[1]]* $
        crude_lo[state.xpmnx[0]:state.xpmnx[1]], color=clr.red, linestyle=2
  endif

end


;;;;;;;;;;;;;;;;;;;;
;  Reset
;;;;;;;;;;;;;;;;;;;;

pro x_fitdla_Reset, state


; Plotting
  state.xymnx = state.svxymnx
  state.old_xymnx = state.svxymnx

end

;;;;;;;;;;;;;;;;;;;;
;  INIT FIT
;;;;;;;;;;;;;;;;;;;;

pro x_fitdla_inifit, state, inifit


  ;; File?
  a = findfile(inifit, count=na)
  if na EQ 0 then begin
      print, 'x_fitdla: FILE ', inifit, ' does not exist!'
      stop
  endif

  ;; Restore
  restore, inifit, /relaxed_structure_assignment

  ;; Error
  if keyword_set( SIG ) then begin
      if sig[0] NE 0. then begin
          state.crude_val = sig
          state.crude = 1
          widget_control, state.crude_id, set_value=1
      endif
  endif

  ;; Continuum
  if keyword_set(conti) then state.conti[*] = conti else state.conti[*] = 1.

  ;; conti_str
  if keyword_set( CSTR ) then state.cstr = cstr

  ;; Lines
  a = where(lines.zabs GT 0., nlin)
  state.nlin = nlin
  tmp = state.lines[0:state.nlin-1]
  copy_struct, lines[0:nlin-1], tmp
  state.lines[0:state.nlin-1] = tmp
  state.curlin = 0L

  ;; Zoom in 
  cen = (1.+state.lines[state.curlin].zabs)*state.wrest
  gd = where(abs(state.wave-cen) LT 200.)
  mx = max(state.conti[gd])
  state.xymnx = [cen-200., -0.1*mx, cen+200., mx*2.]

  ;; Zabs, N
  widget_control, state.Ncolm_id, set_value=state.lines[state.curlin].N
  widget_control, state.bval_id, set_value=state.lines[state.curlin].b
  widget_control, state.zabs_id, set_value=state.lines[state.curlin].zabs

  ;; nset
  state.nset = state.lines[state.nlin-1].set

  ;; 
  x_fitline_updfit, state, EXACT=state.exact

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fitdla, yin, ysin, XSIZE=xsize, YSIZE=ysize, TITLE=title, $
              WAVE=wave, BLOCK=block, FWHM=fwhm, INIFIT=inifit, $
              INFLG=inflg, ZIN=zin, BETA=beta, OLD_INIT=OLD_INIT, $
              SMOOTH=smooth, LIMIT=limit, OUTFIT=outfit, OUTPS=outps

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_fitdla, fx, [ys_in] XSIZE=,YSIZE=, TITLE=, WAVE=, '
    print, '            INIFIT=, INFLG=, FWHM=, /BLOCK, ZIN=, /BETA) [v1.1]'
    return
  endif 

;  Optional Keywords

  device, get_screen_size=ssz
  if not keyword_set( XSIZE ) then begin
      if ssz[0] gt 2*ssz[1] then begin    ;in case of dual monitors
          ssz[0]=ssz[0]/2      
          ; force aspect ratio in case of different screen resolution,
          ; assumes widest resolution used is a 1.6 aspect ratio.
          if ssz[0]/ssz[1] lt 1.6 then ssz[1]=ssz[0]/1.6 
      endif
      xsize = ssz[0]-200
  endif
  if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
  if not keyword_set( FWHM ) then    fwhm = 2. ;; Pixels

; Read in the Data
  if not keyword_set( INFLG ) then inflg = 0
  ydat = x_readspec(yin, INFLG=inflg, head=head, $
                    NPIX=npix, WAV=xdat, FIL_SIG=ysin, SIG=ysig)

  if not keyword_set(YSIG) then ysig = replicate(1., npix)
  if keyword_set(smooth) then begin
     ydat=smooth(ydat, smooth, /NAN)
     ysig=ysig/sqrt(smooth)
  endif

  srt = sort(xdat)
  xdat = xdat[srt]
  ydat = ydat[srt]
  ysig = ysig[srt]

if keyword_set(old_init) then $
  tmp1 = { abslinstrct } $
else $
  tmp1 = { newabslinstrct }

  tmp2 = { contistrct }

;    STATE
  state = { fx: ydat, $
            wave: xdat, $
            sig: ysig, $
            npix: npix, $
            flg_sig: 0, $
            flg_zoom: 0, $
            pos: [0.1,0.1,0.95,0.95], $ ; Plotting
            FWHM: fwhm, $   ; FWHM of instrument (pix)
            wrest: 1215.6701d, $
            beta: 0, $
            nlin: 0, $ ; DLA flag
            EXACT: 1,$
            lines: replicate(tmp1,20), $ ; DLA flag
            fit: fltarr(n_elements(ydat)) + 1., $
            conti: replicate(1.,npix), $
            cstr: tmp2, $
            curlin: 0L, $
            nset: -1, $
            crude: 0, $  ; Crude Error
            crude_val: [-0.1, 0.1], $
            xpos: 0.d, $
            ypos: 0.d, $
            flg_dum: 1, $
            psfile: 0, $ ; Postscript
            psfilename: 'idl.ps', $  ; mendez
			outfilename:'fort.idl',$ ; mendez
            svxymnx: [0., $     ; xmin
                      min(ydat)-0.01*abs(max(ydat)-min(ydat)), $ ; ymin
                      float(n_elements(ydat)-1), $ ; xmax
                      max(ydat)+0.01*abs(max(ydat)-min(ydat))], $ ; ymax
            xymnx: fltarr(4), $
            old_xymnx:fltarr(4), $
            tmpxy: fltarr(4), $
            xpmnx: lonarr(2), $
            psym: 10, $
            title: '', $
            help: strarr(30), $
            size: lonarr(2), $
            xcurs: 0.0, $
            ycurs: 0.0, $
            base_id: 0L, $      ; Widgets
            lblordr_id: 0L, $
            draw_base_id: 0L, $
            draw_id: 0L, $
            text_id: 0L, $
            xpos_id: 0L, $
            ypos_id: 0L, $
            zabs_id: 0L, $
            bval_id: 0L, $
            Ncolm_id: 0L, $
            conti_id: 0L, $
            lines_id: 0L, $
            crude_id: 0L, $
            exact_id: 0L, $
            crudehi_id: 0L, $
            crudelo_id: 0L, $
            funclist_id: 0L, $
            error_msg_id: 0L $
          }

; WAVE
  if keyword_set(BETA) then begin
      state.wrest = 1025.7223d
      state.beta = 1
  endif
  if keyword_set(WAVE) then state.wave = wave
  if keyword_set( YSIG ) then state.flg_sig = 1

  resolve_routine, 'x_specplot', /NO_RECOMPILE
  resolve_routine, 'x_fitline', /NO_RECOMPILE

; Set svxymnx[0,2]
  state.svxymnx[0] = min(state.wave)
  state.svxymnx[2] = max(state.wave)

;    Title
  if size(yin, /type) EQ 7 then state.title = yin
  if keyword_set( TITLE ) then state.title = title

;  HEY ALEX MAKE SURE TO CORRECT THIS BEFORE CVS UPDATE
;  CHANGE FOR CVS UPDATE 
; this is same as a character, but I am leaving it in since I like to call it from the cmd line.
;  if keyword_set( LIMIT ) then begin   ; it should be this
  if (keyword_set( LIMIT )) OR (not keyword_set( LIMIT )) then begin
	state.svxymnx[1] = -0.2  ; reset the limits on the graph 1/0 error
	state.svxymnx[3] = 1.5   ; upper limit
  endif

;    WIDGET
  base = WIDGET_BASE( title = 'x_fitdla', /column)
  state.base_id = base
  
;      Toolbar
          
  toolbar = WIDGET_BASE( state.base_id, /row, /frame, /base_align_center,$
                           /align_center)
;        Version 
;  strlbl = strarr(10)
;  strlbl = ['x_fitdla', ' ', 'Ver 1.0']
;  verslabel = WIDGET_TEXT(toolbar, value=strlbl, xsize=15, ysize=3)

;        Help
  state.help[0] = '  :::Help Menu::: '
  state.help[1] = 'LMB/LMB -- Set region'
  state.help[2] = 's/s -- Set region'
  state.help[3] = 'l -- Set Left '
  state.help[4] = 'r -- Set Right '
  state.help[5] = 'b -- Set Bottom '
  state.help[6] = 't -- Set Top '
  state.help[7] = 'z -- Set ymin to 0.'
  state.help[8] = 'T -- Set ymax to 1.1'
  state.help[9] = 'w -- Reset the screen'
  state.help[11] = 'Z -- Set redshift by hand'
  state.help[12] = 'L -- Set redshift with a line'
  state.help[13] = 'i -- Zoom in'
  state.help[14] = 'o -- Zoom out'
  state.help[15] = '[ -- Pan left'
  state.help[16] = '] -- Pan right'
  state.help[17] = 'H -- Show this screen'
  state.help[18] = 'A -- Al III doublet'
  state.help[19] = 'C -- C IV doublet'
  state.help[20] = 'E -- EW measurement'
  state.help[21] = 'N -- AODM'
  state.help[22] = 'q -- Quit '


;;;;;;;;;
;  Toolbar

; zabs

  linbar = WIDGET_BASE( toolbar, /row, /frame, /base_align_center,$
                           /align_center)
  state.zabs_id = cw_field(linbar, title='zabs', value=0., /floating, $
                           /column, /return_events, xsize=10, uvalue='ZABS')
  state.Ncolm_id = cw_field(linbar, title='Ncolm', value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='NCOLM')
  state.bval_id = cw_field(linbar, title='bval', value=zabs, /floating, $
                           /column, xsize=10, /return_events, uvalue='BVAL')



; continuum
  state.conti_id = cw_field(toolbar, title='conti', value=state.conti, /floating, $
                           /column, xsize=12, /return_events, uvalue='CONTI')



; xy position

  xy_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.xpos_id = cw_field(xy_id, title='x:', value=0., xsize=13)
  state.ypos_id = cw_field(xy_id, title='y:', value=0., xsize=13)

; Crude error
  crude_err = widget_base(toolbar, /row, /align_center, frame=2)
  state.crude_id = CW_BGROUP(crude_err, ['No', 'Yes'], $
                              row=2, /exclusive, /no_release, $
                              set_value=0,  uvalue='CRUDE')
  crude_val = widget_base(crude_err, /column, /align_center, frame=2)
  state.crudehi_id = CW_FIELD(crude_val, value=state.crude_val[1], xsize=5, $
                              title='Crude High', /floating, $
                              /return_events, uvalue='CRUDEHI')
  state.crudelo_id = CW_FIELD(crude_val, value=state.crude_val[0], xsize=5, $
                              title='Crude Low', /floating, $
                              /return_events, uvalue='CRUDELO')
; Exact
  exbase_id = widget_base(toolbar, /column, /align_center, frame=2)
  state.exact_id = CW_BGROUP(exbase_id, ['No', 'Yes'], label_top='Exact?', $
                              row=2, /exclusive, /no_release, $
                              set_value=state.exact,  uvalue='EXACT')

  
  
;;;;;;;;;;;;
;      Drawing
  state.size[0] = xsize
  state.size[1] = ysize
  state.draw_base_id = widget_base(base, /column, /base_align_left, $
                                   uvalue='DRAW_BASE', frame=2, $
                                   /tracking_events)

  state.draw_id = widget_draw(state.draw_base_id, xsize=state.size[0], $
                              ysize=state.size[1], /frame, retain=2, $
                              /button_events, /motion_events, uvalue='DRAW')
;      Text
  state.text_id = widget_text(base, $
                              /all_events, $
                              scr_xsize = 1, $
                              scr_ysize = 1, $
                              units = 0, $
                              uvalue = 'TEXT', $
                              value = '')

;      Done
  buttbase = WIDGET_BASE( toolbar,  /column, /base_align_center,$
                           /align_center)
  specplt = WIDGET_BUTTON(buttbase, value='Specplt',uvalue='SPECPLT', $
                          /align_right)
  lyb = WIDGET_BUTTON(buttbase, value='Lyb fit',uvalue='LYB', /align_right)
  idlout = WIDGET_BUTTON(buttbase, value='IDL Out',uvalue='IDLOUT', /align_right)
  buttbase2 = WIDGET_BASE( toolbar,  /column, /base_align_center,$
                           /align_center)
  done = WIDGET_BUTTON(buttbase2, value='Done',uvalue='DONE', /align_right)
  
; Realize
  WIDGET_CONTROL, base, /realize
  
; Color Table 
;  if not keyword_set( NOCTB ) then loadct, 2, /silent

; HELP
  i = 0
  state.help[0] = '  :::Help Menu for x_fitdla::: '
  state.help[i++] = 'RMB -- Recenter line'
  state.help[i++] = 's/s -- Set region'
  state.help[i++] = 'lrbt -- Set Left Right Bottom Top'
  state.help[i++] = 'T -- Set ymax to 1.1'
  state.help[i++] = 'zz -- Zoom corners'
  state.help[i++] = 'io -- Zoom in/out'
  state.help[i++] = '{}[] -- Pan'
  state.help[i++] = 'w -- Reset the screen'
  state.help[i++] = 'H -- This widget'
  state.help[i++] = 'c -- Add new Lya line'
  state.help[i++] = 'L -- Add new LLS'
  state.help[i++] = 'B -- Add new Beta'
  state.help[i++] = 'd -- Delete current line'
  state.help[i++] = 'C -- Set continnuum' 
  state.help[i++] = 'nN -- Adjust colm'
  state.help[i++] = 'vV -- Adjust bvalue'
  state.help[i++] = '=- -- Loop through lines'
  state.help[i++] = 'P -- Print to postrscipt'
  state.help[i++] = 'I -- IDL output'
  state.help[i++] = '3,4 -- Add/move continuum point'

; Update
  x_fitdla_Reset, state

  ;; Init Fit and ZIN
  if keyword_set(INIFIT) then x_fitdla_inifit, state, inifit
  if keyword_set(INIFIT) and keyword_set(ZIN) then begin
     print, 'Do not set both INIFIT and ZIN at once!'
     stop
  endif
  if keyword_set( ZIN ) then begin
      state.xpos = (1.+zin)*state.wrest
      x_fitline_newlin, state, BETA=state.beta
  endif

  if keyword_set(OUTFIT) then state.outfilename = outfit
  if keyword_set(OUTPS) then state.psfilename = outps

  ;; Set pmnx
  state.xpmnx = x_getxpmnx(state)
  
  ;; Plot
  x_fitdla_UpdatePlot, state
  WIDGET_CONTROL, base, set_uvalue = state, /no_copy

  print, 'Press H for help'

; Send to the xmanager
  if not keyword_set(BLOCK) then xmanager, 'x_fitdla', base, /no_block $
  else xmanager, 'x_fitdla', base

return
end
