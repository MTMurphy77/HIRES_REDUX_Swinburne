;+ 
; NAME:
; hamspec_subbias   
;     Version 1.2
;
; PURPOSE:
;    Remove bias from all exposures listed in the "hamspec" structure.
;    The main driver is hamspec_suboscan which strips the image of the
;    overscan region (see that program for a full
;    description).  The hamspec_subbias routine will also remove an
;    archived bias image if requested.
;
; CALLING SEQUENCE:
;   
;  hamspec_subbias, hamspec, indx, /usebias, /nobiasrow, /clobber, /ARC,
;  /debug, BADROW=, RBIN=, CBIN=
;
; INPUTS:
;   hamspec  -  MIKE structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CBIN=  -- Column binning
;  RBIN=  -- Row binning
;  /NOBIASROW= if set, bias row is not used (normally, you should not use it)
;  /USEBIAS = use the bias image in addition to subtracting a fit to the
;             overscan columns (at right).  The code first generates a
;             smoothed version of the bias image.
;             If you do this, you want to be sure the same options
;             were used in hamspec_mkbias (i.e. /nobiasrow )
;             Presently, this step is not recommended.
;  /CLOBBER = Overwrite existing OV files 
;  /DEBUG -- Turn debug mode on
;  BADROW -- Rows identified as anomolous in the overscan region.
;            Generally the result of an electronics 'hiccup'.
;
; 
; OUTPUTS TO SNGL:
;   /NOFITS -- Do not write a fits file
;   OVIMG=    -- Bias subtracted image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;
; EXAMPLES:
;   hamspec_subbias, hamspec, indx, /usebias, /nobiasrow
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;   hamspec_subbias_sngl -- Subroutine under hamspec_subbias.  Uses most of
;                        the keywords described above.  Example:  
;                        rslt = hamspec_subbias_sngl('Raw/mb0020.fits')
;   hamspec_suboscan
;
; REVISION HISTORY:
;   01-Feb-2005 Written by JXP  (taken from mike_subbias)
;                                
;------------------------------------------------------------------------------

function hamspec_subbias_sngl, rawfil, $
                             USEBIAS=usebias, NOBIASROW = nobiasrow, $
                             CLOBBER=clobber, DEBUG=debug, IMTYP=imtyp, $
                             BADROW=badrow, OVIMG=ovimg, NOFITS=nofits, $
                             SILENT=silent, CBIN=cbin, RBIN=rbin, FRAME=frame

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 1  then begin 
      print,'Syntax:  ' + $
        'rslt = hamspec_subbias_sngl(rawfil, /noBIASROW, ' + $
        '/USEBIAS, OVROOT=, OVIMG=, /NOFITS, ' + $
        '/CLOBBER, /DEBUG, /SILENT [v1.2])'
      return, -1
  endif 
  
; Create OV directory (if necessary)
  a = findfile('OV/..', count=count)
  if count EQ 0 then file_mkdir, 'OV'
  
; Optional Keywords
  if not keyword_set( IMTYP ) then begin
      if keyword_set( ARC ) then imtyp = 'ARC' else imtyp = 'UNK'
  endif

  if not keyword_set( FRAME ) then begin
      pos1 = strpos(rawfil, 'HI.')
      if pos1 NE -1 then stop
      pos = strpos(rawfil, '.fits')
      frame = long(strmid(rawfil,pos-4,4))
  endif

  ;; Check for output
  outfil = hamspec_getfil('ov_fil', FRAME=frame, /name, CHKFIL=chkf)
  if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
      print,$
        'hamspec_subbias_sngl: File ', outfil, $
        ' found. Overscan already subtracted.'
      return, 0
  endif
          
  ;; Open Raw image
  raw = xmrdfits(strtrim(rawfil,2),0, head, /silent, /fscale)
  head0 = head
  if not keyword_set( SILENT ) then $
    print, 'hamspec_subbias: Subtracting the overscan for image ', rawfil
  sz = size(raw, /dimensions)
      
  ;; Set cbin, rbin
  if not keyword_set( CBIN ) then cbin = round(2048. / sz[0] )
  if not keyword_set( RBIN ) then rbin = round(4096. / sz[1] )
      
  ;; JXP :: Am saving the entire science image now
  cover = sxpar(head, 'COVER')
  n1 = sxpar(head,'NAXIS1')
  n2 = sxpar(head,'NAXIS2')
  dats = [0,n1-cover-1,0,n2-1]
  fincol = dats[1]
 
  ;; TRIM Vignetting
;  dats[3] = dats[3] < round(3990./rbin)

  ;; TRIM Red chip
;  if chip EQ 3L then dats[0] = dats[0] > round(110./cbin)
          
  if keyword_set (USEBIAS) then begin 
      bias = hamspec_getfil('bias_fil', SZ=[2048L/cbin,4096L/rbin],  FIL_NM=bias_fil)
      if not keyword_set( SILENT ) then $
        print, 'hamspec_subbias: Using BIAS file: ', bias_fil
      bias_size=size(bias,/dimensions)
      if not keyword_set( SILENT ) then $
        print, '   Smoothing bias image.'
      ;; DEBUG
      if keyword_set (debug) then begin
          for ii=10,15 do begin
              plot, findgen(100), bias[ii,100:199]
              smth =  smooth(bias[ii,*],59,/edge_truncate)
              oplot, findgen(100), smth[100:199]
              bias[ii,*] = smth
              oplot, findgen(100), bias[ii,100:199], $
                color=colors.red
          endfor
      endif

      for ii = 0,fincol-1 do begin 
          smth =  smooth(bias[ii,*],59,/edge_truncate)
          bias[ii,*] = smth
      endfor
      if keyword_set (debug) then begin
          xatv, bias, /block
          stop
      endif
              
  endif
          
  ;; Work
      
  if keyword_set (USEBIAS) then begin 
      ovimg = temporary(ovimg) - bias
      sxaddpar, head, 'BIAS', 'T', bias_fil
  endif

  ;; Ovesrcan
  x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=imtyp, $
    DEBUG= debug, SVBAD=badrow, CBIN=cbin, RBIN=rbin
      
  ;; Trim
;  dats = dats-1
  ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
  sxaddpar, head0, 'TRIM', 'T'
  sxaddpar, head0, 'ODSEC', strjoin(string(dats))
      
  ;; Rotate
  ovimg = transpose(ovimg)

  ;; Header
  mkhdr, main_head, ovimg
  sxdelpar, main_head, 'END'
  sxdelpar, head0, 'NAXIS'
  sxdelpar, head0, 'BITPIX'
  sxdelpar, head0, 'BZERO'
  sxdelpar, head0, 'SIMPLE'
  sxdelpar, head0, 'DATSUM'
  sxdelpar, head0, 'DATE'

  nhd = n_elements(head0)
  for qq=0L,nhd-1 do begin
      if strlen(strtrim(head0[qq],2)) GT 0 then $
        main_head = [main_head, head0[qq]]
  endfor
      
  ;; Write out bias-subtracted image
  if not keyword_set(NOFITS) then $
    mwrfits, ovimg, outfil, main_head, /create, /silent
      
  if not keyword_set( SILENT ) then $
    print, 'hamspec_subbias_sngl: All done!'

  return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_subbias, hamspec, indx, USEBIAS=usebias, NOBIASROW = nobiasrow, $
                  CLOBBER=clobber, DEBUG=debug, BADROW=badrow
              ; , VIEW=view

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 2  then begin 
      print,'Syntax:  ' + $
        'hamspec_subbias, hamspec, indx, /noBIASROW, /USEBIAS, OVROOT=, ' + $
        '/CLOBBER, /DEBUG, BADROW= [v1.1]'
;      print,'Recommended:  ' + $
;        'hamspec_subbias, hamspec, /noBIASROW, /USEBIAS, /CLOBBER '
      return
  endif 
  
  
  cbin = 0L
  rbin = 0L
  
  nindx = n_elements(indx)
  
  ;; Loop
  for q=0,nindx-1 do begin
      ;; Check for output
      outfil = hamspec_getfil('ov_fil', $
                   FRAME=hamspec[indx[q]].frame, /name, CHKFIL=chkf)
      if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
          print,$
            'hamspec_subbias: File ', outfil, $
            ' found. Overscan already subtracted.'
          hamspec[indx[q]].img_ov = outfil
          hamspec[indx[q]].flg_ov = 1
          continue
      endif
          
      ;; Open Raw image
      raw = xmrdfits(strtrim(hamspec[indx[q]].rootpth+hamspec[indx[q]].img_root,2), 0, head, $
                     /silent, /fscale)
      head0 = xheadfits(hamspec[indx[q]].rootpth+hamspec[indx[q]].img_root)
      print, 'hamspec_subbias: Subtracting the overscan for image ', $
        hamspec[indx[q]].rootpth+hamspec[indx[q]].img_root,  ' ', hamspec[indx[q]].type

      ;; FINCOL
      cover = sxpar(head, 'COVER')
      n1 = sxpar(head,'NAXIS1')
      n2 = sxpar(head,'NAXIS2')
      dats = [0,n1-cover-1,0,n2-1]
      fincol = dats[1]

      ;; TRIM Vignetting
;      dats[3] = dats[3] < round(3990./hamspec[indx[q]].rowbin)

      ;; TRIM Red chip
;      if chip EQ 3L then dats[0] = dats[0] > round(110./hamspec[indx[q]].colbin)

      if keyword_set (USEBIAS) then begin 
         stop ;; Not doing this so far
          bias = hamspec_getfil('bias_fil', SZ=[2048L/cbin,4096L/rbin], FIL_NM=bias_fil)
          print, 'hamspec_subbias: Using BIAS file: ', bias_fil
          bias_size=size(bias,/dimensions)
          print, '   Smoothing bias image.'
          ;; DEBUG
          if keyword_set (debug) then begin
              for ii=10,15 do begin
                  plot, findgen(100), bias[ii,100:199]
                  smth =  smooth(bias[ii,*],59,/edge_truncate)
                  oplot, findgen(100), smth[100:199]
                  bias[ii,*] = smth
                  oplot, findgen(100), bias[ii,100:199], $
                    color=colors.red
              endfor
          endif
          
          
          if keyword_set(SMOOTH) then begin
              for ii = 0,fincol do begin 
                  smth =  smooth(bias[ii,*],59,/edge_truncate)
                  bias[ii,*] = smth
              endfor
          endif
          
          if keyword_set (debug) then begin
              xatv, bias, /block
              stop
          endif
          
      endif
      

      
      if keyword_set (USEBIAS) then begin 
         stop ;; Not setup for this
          ovimg = temporary(ovimg) - bias
          hamspec[indx[q]].flg_anly = 3
          sxaddpar, head0, 'BIAS', 'T', bias_fil
      endif

      ;; Subtract overscan
      x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=hamspec[indx[q]].type, $ 
                  DEBUG= debug, SVBAD=badrow, CBIN=hamspec[indx[q]].colbin, $
                  RBIN=hamspec[indx[q]].rowbin
      
      ;; Trim
;      dats = dats-1
      ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
      head0 = head
      sxaddpar, head0, 'TRIM', 'T'
      sxaddpar, head0, 'ODSEC', strjoin(string(dats))

      ;; Rotate
      ovimg = transpose(ovimg)
      
      ;; Write out bias-subtracted image
      hamspec[indx[q]].img_ov = outfil
      hamspec[indx[q]].flg_ov = 1

      ;; Header
      mkhdr, main_head, ovimg
      sxdelpar, main_head, 'END'
      sxdelpar, head0, 'NAXIS'
      sxdelpar, head0, 'NAXIS1'
      sxdelpar, head0, 'NAXIS2'
      sxdelpar, head0, 'BITPIX'
      sxdelpar, head0, 'BZERO'
      sxdelpar, head0, 'SIMPLE'
      sxdelpar, head0, 'DATSUM'
      sxdelpar, head0, 'DATE'

      nhd = n_elements(head0)
      for qq=0L,nhd-1 do begin
          if strlen(strtrim(head0[qq],2)) GT 0 then $
            main_head = [main_head, head0[qq]]
       endfor
      
      mwrfits, ovimg, outfil, main_head, /create, /silent
      
  endfor

  print, 'hamspec_subbias: All done!'

  return
end
