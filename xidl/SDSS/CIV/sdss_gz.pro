;+ 
; NAME:
; sdss_gz.pro
;    Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;    2-Nov-2011  created by KLC (removed from sdss_functions)
;-
;------------------------------------------------------------------------------

function sdss_gz_mkzarr, zlim, zbinsize, header=header, nzbin=nzbin
  ;; consistent way of (re)constructing the global redshift array (binned)
  if n_params() ne 2 and not keyword_set(header) then begin
     print,'Syntax - sdss_gz_mkzarr( zlim, zbinsize, [header=, nzbin=])'
     return, -1
  endif 

  if keyword_set(header) then begin
     ;; Recreate necessary values
     zmin = sxpar(header,'ZMIN')
     zmax = sxpar(header,'ZMAX')
     zlim = [zmin,zmax]
     
     zbinsize = sxpar(header,'ZBINSZ')

     nzbin = sxpar(header,'NAXIS1')
  endif else nzbin = ceil((zlim[1] - zlim[0])/zbinsize) + 1L
  zglobal = zlim[0] + findgen(nzbin)*zbinsize
  
  return, zglobal

end                             ; sdss_gz_mkzarr()



function sdss_gz_mkewarr, ewlim, ewbinsize, header=header, newbin=newbin
  ;; consistent way of (re)constructing the global EWlim array (binned)
  if n_params() ne 2 and not keyword_set(header) then begin
     print,'Syntax - sdss_gz_mkzarr( ewlim, ewbinsize, [header=, newbin=])'
     return, -1
  endif 

  if keyword_set(header) then begin
     ;; Recreate necessary values
     ewmin = sxpar(header,'EWMIN')
     ewmax = sxpar(header,'EWMAX')
     ewlim = [ewmin,ewmax]
     
     ewbinsize = sxpar(header,'EWBINSZ')

     newbin = sxpar(header,'NAXIS2')
  endif else newbin = ceil((ewlim[1] - ewlim[0])/ewbinsize) + 1L
  ewglobal = ewlim[0] + findgen(newbin)*ewbinsize
  
  return, ewglobal

end                             ; sdss_gz_mkewarr()


pro sdss_gz_readgzfil, gzcum_fil, gzglobal, zarr, ewarr, $
                       zbinsize=zbinsize, nzbin=nzbin, $
                       ewbinsize=ewbinsize, newbin=newbin, $
                       header=header, zlim=zlim, ewlim=ewlim
  if n_params() ne 4 then begin
     print,'Syntax - sdss_gz_readgzfil, gzcum_fil, gzglobal, zarr, ewarr, '
     print,'                            [zbinsize=, nzbin=, ewbinsize=, newbin=,'
     print,'                            header=, zlim=, ewlim=]'
     return
  endif 

  gzglobal = xmrdfits(gzcum_fil,0,header,/silent) 
  zarr = sdss_gz_mkzarr(zlim,zbinsize,header=header,nzbin=nzbin)
  ewarr = sdss_gz_mkewarr(ewlim,ewbinsize,header=header,$
                          newbin=newbin)
end                             ; sdss_gz_readgzfil


function sdss_gz_rebin, gzcum_fil, zbinsize, zarr=zarr, $
                        nzbin=nzbin, zlim=zlim, $
                        header=header
  ;; Collapse gzglobal 2D grid to larger bins
  if n_params() ne 2 then begin
     print,'Syntax - '
     return, -1
  endif
  
  sdss_gz_readgzfil, gzcum_fil, gzglobal0, zarr0, ewarr0, $
                     zbinsize=zbinsize0, zlim=zlim0, nzbin=nzbin0, $
                     newbin=newbin0, header=header
  
  ;; Find nearest (lower) redshift bin
  if keyword_set(zlim) then begin
     izmin = floor((zlim[0] - zlim0[0]) / zbinsize0) 
  endif else begin
     izmin = 0
  endelse 
  sxaddpar,header,'ZMIN',zarr0[izmin]

  ;; Figure out number of z bins to sum over
  sumbin = floor(zbinsize / zbinsize0)
  zbinsize = sumbin * zbinsize0
  sxaddpar,header,'ZBINSZ',zbinsize

  ;; Figure out how many can fit this exactly
  nzbin = floor( (nzbin0 - izmin) / sumbin ) ; max out
;  if keyword_set(zlim) then $
;     nzbin = floor( (
  sxaddpar,header,'NAXIS1',nzbin
  sxaddpar,header,'ZMAX',(nzbin-1)*zbinsize

  ;; Instantiate (like in sdss_gz)
  gzglobal = fltarr(nzbin,newbin0,2) ; dz, dX
  for zz=0L,nzbin-1 do begin
     istart = izmin + sumbin*zz
     iend = istart + sumbin - 1
     gzglobal[zz,*,0] = total(gzglobal0[istart:iend,*,0],1)
     gzglobal[zz,*,1] = total(gzglobal0[istart:iend,*,1],1)

;     print,zz,istart,iend,gzglobal[zz,20,1],total(gzglobal0[istart:iend,20,1]),median(gzglobal0[istart:iend,20,1])
  endfor                        ; loop zz=nzbin
  zarr = sdss_gz_mkzarr(zlim,zbinsize,header=header)

;  stop
  return, gzglobal
end                             ; sdss_gz_rebin()


function sdss_gz_dx, gzcum_fil, zval, ewval, dz=dz, $
                     zarr=zarr, ewarr=ewarr
  ;; Interpolate 2D grid to return dX (or dz) for given zval and ewval
  ;; (which may either both be arrays or one array and one single value)
  if n_params() ne 3 then begin
     print,'Syntax - sdss_gz_dx(gzcum_fil, zval, ewval, [/dz, zarr=, ewarr=])'
     return, -1
  endif 

  if keyword_set(dz) then idelta = 0 else idelta = 1
  if size(gzcum_fil,/type) eq 7 then begin
     gzglobal = xmrdfits(gzcum_fil, 0, header, /silent)
     zarr = sdss_gz_mkzarr(zlim,zbinsize,header=header, nzbin=nzbin)
     ewarr = sdss_gz_mkewarr(ewlim,ewbinsize,header=header, newbin=newbin)
  endif else begin
     if not keyword_set(zarr) or not keyword_set(ewarr) then begin
        print,'sdss_gz_dx(): must set zarr= and ewarr='
        return, -1              ; EXIT
     endif 
     gzglobal = gzcum_fil       ; assume vector
     nzbin = (size(zarr,/dim))[0]
     zlim = [zarr[0],zarr[nzbin-1]] ; assume sorted
     zbinsize = zarr[1] - zarr[0]   ; assume uniform [THIS MAY BE WRONG LATER]

     newbin = (size(ewarr,/dim))[0]
     ewlim = [ewarr[0],ewarr[newbin-1]] 
     ewbinsize = ewarr[1] - ewarr[0]
  endelse

  nzval = (size(zval,/dim))[0] > 1 ; handle singularity
  newval = (size(ewval,/dim))[0] > 1 

  if nzval ne newval then begin
     if nzval eq 1 or newval eq 1 then grid = 1 $
     else begin
        print,'sdss_gz_dx(): cannot handle mismatched zval, ewval sizes ',$
              nzval, newval
        return, -1
     endelse
  endif 

  ;; Interpolate in 2D
  ;; First must identify the fractional "indices" to use with
  ;; IDL's interpolate()
  ;; Linear interpolation for mid-point (xm, ym) between points
  ;; (x1,y2) and (x2,y2) 
  ;; ym =  (y2-y1)/(x2-x1) * (xm - x1) + y1

  iz1 = floor((zval - zlim[0])/zbinsize)
  iz2 = iz1 + 1
  izm = (iz2-iz1)/(zarr[iz2]-zarr[iz1])  * (zval - zarr[iz1]) + iz1

  iw1 = floor((ewval - ewlim[0])/ewbinsize)
  iw2 = iw1 + 1
  iwm = (iw2-iw1)/(ewarr[iw2]-ewarr[iw1])  * (ewval - ewarr[iw1]) + iw1
  
  ;; use indices to interpolate in gzglobal slice and set out-of-bound
  ;; values to -1
  dx = interpolate(gzglobal[*,*,idelta], izm, iwm, missing=-1, grid=grid) 

;  stop

  return, dx

end                             ; sdss_gz_dX()


function sdss_gz_sanitycheck, dz=dz, dztot=dztot
  ;; Brute force way of making sure sdss_gz gets sensible results

  if not keyword_set(cosmology) then $
     cosmology = [71.9, 0.258, 0.742] ; WMAP5
  cosm_common, H0=cosmology[0], Omegavac=cosmology[2], $
               OmegaDM=cosmology[1],/silent 

  sdssdir = sdss_getsdssdir()
  snrstrct = xmrdfits(sdssdir+'inputs/dr7qso_noBAL_SNR.fit',1,/silent)
  gd = where(snrstrct.snr_civ[2] ge 4 and snrstrct.z_qso ge 1.7,ngd)

  dblt = dblt_retrieve('CIV')
  zmin = snrstrct[gd].wvobs_civ[0] / dblt.wvI - 1.
  zmax = snrstrct[gd].wvobs_civ[1] / dblt.wvI - 1.

  dxarr = cosm_xz(zmax,zmin=zmin,/silent,/exact,/noinit)
  dzarr = zmax - zmin

  dztot = total(dzarr)
  dxtot = total(dxarr)
  if keyword_set(dz) then begin
     tmp = dxtot
     dxtot = dztot
     dztot = tmp
  endif

  return,dxtot
     
end                             ; sdss_gz_sanitycheck()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_gz, list_fil, gzcum_fil, nsig=nsig, zbinsize=zbinsize, zlim=zlim,$
                   ewbinsize=ewbinsize, ewlim=ewlim,  cosmology=cosmology, $
                   clobber=clobber, debug=debug
  ;; Rudimentary completeness test results
  if n_params() ne 2 then begin
     print,'Syntax - sdss_concatgz, list_fil, gzcum_fil, [nsig=, zbinsize=, zlim=, '
     print,'                        ewbinsize=, ewlim=,  cosmology=, /clobber, /debug]'
     return
  endif

  if not keyword_set(cosmology) then $
     cosmology = [71.9, 0.258, 0.742] ; WMAP5
  cosm_common, H0=cosmology[0], Omegavac=cosmology[2], $
               OmegaDM=cosmology[1],/silent 
  
  sdssdir = sdss_getsdssdir()

  ;; Defaults to make something sensible
  if not keyword_set(nsig) then nsig = 3. ; ~99% c.l.
  if not keyword_set(zbinsize) then zbinsize = 0.005
  if not keyword_set(zlim) then zlim = [1., 6.]
  if not keyword_set(ewbinsize) then ewbinsize = 0.05 ; Ang
  if not keyword_set(ewlim) then ewlim = [0.05, 5.]   ; Ang
  dxsdsstot = sdss_gz_sanitycheck()

  ;; For reference, on the c.l. of given nsig
  ;; gsigma = 1. + findgen(10)*0.5
  ;; gprob = [68.26895123013139,86.63856338891091,$
  ;;          95.4499772097615,98.75806921717441,$
  ;;          99.7300213908383,99.9534744959942,$
  ;;          99.99366582300915,99.9993204775126,$
  ;;          99.99994267123438,99.99999620223669]

  ;; Bins are defined on the left-hand side by convention
  zglobal = sdss_gz_mkzarr(zlim, zbinsize, nzbin=nzbin)
  ewglobal = sdss_gz_mkewarr(ewlim, ewbinsize, newbin=newbin)

  ;; Create place to store EVERYTHING
  gzglobal = fltarr(nzbin,newbin,2) ; dz and dX

  ;; Read in data and get names
  readcol,list_fil,spec_fil,format='a',/silent,skip=1
  nspec = (size(spec_fil,/dim))[0]

  ;;  gz array will be [npix, [zabs, dz, EWlim, mask, nsig]]
  gz_fil = sdss_getname(spec_fil,/spec,/gz,dir=gzdir)

  ;; Loop over spectra and sum up dz and dX
  for ss=0L,nspec-1 do begin 
;  for ss=8519L,nspec-1 do begin 
     gz_los = xmrdfits(gzdir[ss]+gz_fil[ss],0,/silent)
     gd = where(gz_los[*,3] eq 1 and gz_los[*,2] gt 0.,npix)
     if npix ne 0 then begin 
        ;; Grep out pertinent arrays (to not be confused)
        zarr = gz_los[gd,0] 
        dz = gz_los[gd,1] 
        ewarr = gz_los[gd,2] * nsig / gz_los[gd,4]             ; scale
        dx = cosm_xz(zarr+dz,zmin=zarr,/silent,/exact,/noinit) ; calc

        ;; Minimize looping
        ;; focus on just redshift range needed to bin
        zstart = floor((zarr[0] - zlim[0])/zbinsize) > 0 
        zend = (floor((zarr[npix-1] - zlim[0])/zbinsize) + 1) < nzbin

        ;; First and last redshift bins may be partials
        for zz=zstart,zend do begin
;        for zz=214L,zend do begin
           sub = where(zarr ge zglobal[zz] and $
                       zarr lt zglobal[zz]+zbinsize,npix_sub)

           if npix_sub eq 0 then begin
              if keyword_set(debug) then $
                 print,'sdss_concatgz debug: zarr does not span zglobal[zz] = '$
                       + string(zglobal[zz],zglobal[zz]+zbinsize,$
                                format='(2(f7.5,2x))')
              continue          ; SKIP
           endif                ; npix_sub = 0

           ;; Minimize looping
           ;; take out sub arrays necessary for this zz redshift bin
           ewsrt = sub[sort(ewarr[sub])]
           zarr_sub = zarr[ewsrt]
           ewarr_sub = ewarr[ewsrt]
           dz_sub = dz[ewsrt]
           dx_sub = dx[ewsrt]

           ;; focus on just EW range needed to bin and loop
           ;; efficiently
           if ewarr_sub[0] gt ewlim[1]+ewbinsize then begin
              if keyword_set(debug) then $
                 print,'sdss_gz debug: all pixels out of EW bounds, ss = ',ss
              continue
           endif 

           ewstart = floor((ewarr_sub[0] - ewlim[0])/ewbinsize) > 0
           ewend = (floor((ewarr_sub[npix_sub-1] - ewlim[0])/ewbinsize) + 1) $
                   < newbin
           if keyword_set(debug) then wwarr = ewstart + lindgen(ewend-ewstart)

           iewmn_prev = 0
           ww = ewstart         ; map to gzglobal array
           while (ww lt ewend) and (iewmn_prev lt npix_sub) do begin
              ;; Looking for where ewarr is lower than current
              ;; ewglobal[ww] value because then all those pixels will be
              ;; sensitive enough
              ewmn = min(ewarr_sub[iewmn_prev:*]-(ewglobal[ww]+ewbinsize),$
                         iewmn,/abs)
              iewcurr = iewmn_prev + iewmn
              ;print,'iewcurr',iewcurr
              if ewarr_sub[iewcurr] gt ewglobal[ww]+ewbinsize then $
                 iewcurr--      ;  don't go into next bin
              ;print,'iewcurr',iewcurr

              ;; min() lists the first instance of minimum but could
              ;; be duplicates; so look ahead
              done = 0
              iewcurr++
              while not done do begin
                 if iewcurr eq npix_sub then begin
                    iewcurr--  ; revert
                    done = 1   ; out of bounds
                 endif else begin
                    if ewarr_sub[iewcurr] ge ewglobal[ww]+ewbinsize then begin
                       iewcurr-- ; revert
                       done = 1
                    endif else iewcurr++ ; increment
                 endelse 
                 ;print,'iewcurr',iewcurr
              endwhile
              ;print,'iewcurr',iewcurr


              if iewcurr lt 0 or iewcurr ge npix_sub then begin
                 ;; This probably should never happen
                 if keyword_set(debug) then $
                    stop,'sdss_concatgz debug: iewcurr out of bounds'
                 continue
              endif 
              
              ;; Since ewarr sparsely populated, figure out how much
              ;; can instantiate now
              if iewcurr lt npix_sub-1 then begin
                 wwmx = min(ewarr_sub[iewcurr+1]-ewglobal,iwwmx,/abs) 
                 ;print,'iwwmx',iwwmx
                 if iwwmx ne ww then begin
                    if (ewarr_sub[iewcurr+1] lt ewglobal[iwwmx]) then $
                       iwwmx--  ; don't go into next bin
                    ;print,'iwwmx',iwwmx
                    if iwwmx ne ewend then $
                       iwwmx--  ; this represents max range
                    ;print,'iwwmx',iwwmx
                 endif          ; iwwmx != ww
              endif else iwwmx = ww
                 ;print,'iwwmx',iwwmx

              if iwwmx lt ww then $
                 stop,'sdss_concatgz: indices messed up ',ww,iwwmx

              ;; Now can just add it up cumulatively 
              ;; What about fractional contributions?!!!
              gzglobal[zz,ww:*,0] = gzglobal[zz,ww:*,0] + $
                                        total(dz_sub[0:iewcurr])
              gzglobal[zz,ww:*,1] = gzglobal[zz,ww:*,1] + $
                                        total(dx_sub[0:iewcurr])
              
              if keyword_set(debug) then begin
                 ;; Watch it grow in this region
                 print,'sdss_concatgz debug: '+$
                       string(zz,ww,iwwmx,iewmn_prev,iewcurr,$
                              format='("zz=",i4,"; ww=",i4,":",i4,"; iew=",i4,":",i4)')
                 print,'ww','EWlim','dz','dX',format='(a4,2x,a6,2(2x,a10))'
                 printcol,wwarr,ewglobal[wwarr],gzglobal[zz,wwarr,0],$
                          gzglobal[zz,wwarr,1],format='(i4,2x,f6.4,2(2x,f10.4))'
                 print,''
              endif             ; /debug

              ;; Increment for next loop
              iewmn_prev = iewcurr + 1
              ww = iwwmx + 1

           endwhile             ; loop ww <= ewend 
        endfor                  ; loop zz=zend

     endif                      ; npix ne 0

     maxdx = max(gzglobal[*,*,1])
     if maxdx gt 0.5*dxsdsstot then $
        stop,'sdss_gz: exceeded 50% of total SDSS dX ',maxdx, dxsdsstot

     if (ss mod 500L) eq 0L then begin
        maxdz = max(gzglobal[*,*,0])
        print,'sdss_concatgz: ss = '+strtrim(ss,2),max(gzglobal)
        if keyword_set(debug) then stop
     endif 

  endfor                        ; loop ss=nspec


  ;; Ready to write out
  fxhmake, header, gzglobal

  ;; Add keywords
  sxaddpar, header, 'NSPEC', nspec, 'Number of spectra'
  sxaddpar, header, 'NSIG', nsig, 'Assumed nsigma in EWlim est'
  sxaddpar, header, 'ZBINSZ', zbinsize, 'Redshift bin size'
  sxaddpar, header, 'ZMIN', zlim[0], 'Minimum redshift'
  sxaddpar, header, 'ZMAX', zlim[1], 'Maximum redshift'
  sxaddpar, header, 'EWBINSZ', ewbinsize, 'EWlim bin size (Ang)'
  sxaddpar, header, 'EWMIN', ewlim[0], 'Minimum EWlim (Ang)'
  sxaddpar, header, 'EWMAX', ewlim[1], 'Maximum EWlim (Ang)'
  sxaddpar, header, 'H0', cosmology[0], 'Adopted Hubble constant'
  sxaddpar, header, 'OmDM', cosmology[1], 'Cosmology: Matter density'
  sxaddpar, header, 'Lambda', cosmology[2], 'Cosmology: Vacuum energy density'

  test = file_search(sdssdir+gzcum_fil+'*',count=ntest)
  if ntest eq 0 or keyword_set(clobber) then begin
     mwrfits,gzglobal,sdssdir+gzcum_fil,header,/create,/silent
     spawn,'gzip -f '+sdssdir+gzcum_fil
     print,'sdss_concatgz: created ',gzcum_fil
  endif 

  stop
end                             ; sdss_concatgz


