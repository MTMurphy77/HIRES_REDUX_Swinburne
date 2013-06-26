;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sdss_genprof.pro             
; Author: Kathy Cooksey                      Date: 23 Aug 2011
; Project: 
; Description: Generate civcandstrct of randomly generated
;              values of z, [N, b, dv, ncomp]
; Input: 
;   nabs - number of absorption lines to make for each wrest
;   zlim - bounds of randomly generated redshifts
;   nlim - bounds of randomly generated column densities (logN)
;   blim - bounds of randomly generated Doppler parameters (km/s)
;   nlimobs - bounds of observed column densities that define f(N)
;   nabsobs - number of absorbers found in nlimobs (normalize f(N))
;   aobs - power-law exponent f(N) = b*N^(-aobs)
; Optional:
;   seed - IDL randomu random number generation seed
;   dvabs - bounds of randomly generated velocity offsets for 
;           components (km/s)
;   ndv - number of components, used for each input wrest, and
;         cannot exceed number civcandstrct holds
;   debug - stop at opportune times for debugging purposes
; Output: 
;   strct_fil - name of output civcandstrct, with nabs elements,
;               with n_elements(wrest) times ndv each
; Example:
;   sdss_genprof,'randciv.fits',50,[0,1],[12,15],[10,40],
;      [12.7,14.5],21,1.83,,dvabs=[-200,200],ndv=5
;   [Danfort & Shull 2007ph: 12.7 < ncolm < 14.5, alpha = 1.83+/-0.32,
;                            delta(z) = 2.22, nabsobs = 21 for CIV]
; History:
;   23 Mar 2008 - created by KLC
;   23 Aug 2011 - adapted from civ_genprof, KLC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@spec_ewlin
pro sdss_genprof_plot,vpstrct,fwhmkms,savfil=savfil,fwhm=fwhm,$
                      _extra=extra
;; Display profiles
  if (N_params() lt 2) then begin
     print,'Syntax - '+ $
           'sdss_genprof_plot,strct_fil,fwhmkms,[savfil=savfil,'+$
           '   fwhm=fwhm,_extra=extra]'
     return
  endif

  if not keyword_set(fwhm) then fwhm = 2 ; assume STIS
  done = 0

  nstrct = n_elements(vpstrct)

  ;; Make (logarithmic) wavelength array with buffer
  bmx = max(vpstrct.b)
  nmx = max(vpstrct.n)
  zlo = min(vpstrct.zabs,max=zhi)
;     dv = (zhi-zlo)*2.998e5
  dv = 2.998e5 * ((1+zhi)/(1+zlo) - 1)
  wvlo = min(vpstrct.wrest,max=wvhi)*(1.+zlo-10*bmx/2.998e5)
  wvhi = wvhi*(1+zhi+10*bmx/2.998e5)
  cdelt1 = alog10(1.+fwhmkms/2.998e5)
  npix = ceil(alog10(wvhi/wvlo)/cdelt1)+1
  wave = 10^(alog10(wvlo)+dindgen(npix)*cdelt1)

  ;; Generate Voigt profile and plot
  vpfx = x_voigt(wave,vpstrct,fwhm=fwhm,_extra=extra) ;e.g., nosmooth
  unq = uniq(vpstrct.wrest,sort(vpstrct.wrest))
  ncomp = n_elements(vpstrct)/n_elements(unq)
  ttl = 'logN = '+string(nmx,format='(f5.2)')+'; b = '+$
        strtrim(round(bmx),2)+' km/s; ncomp = '+strtrim(ncomp,2) + $
        '; dv = '+strtrim(round(dv),2)+' km/s'
  x_splot,wave,vpfx,title=ttl,ymnx=[-0.1,1.1],/block

end                             ;sdss_genprof_plot


function sdss_genprof_config,config_fil
  ;; Read in and return structure from configuration file
  cfg = {nlim:fltarr(2), $
         blim:fltarr(2), $
         dvabs:fltarr(2), $
         ndv:0L, $
         dndz:0., $
         frac_rm:0., $
         pixscale:0., $
         dvelo:0., $
         nabs:0L, $
         zlim:fltarr(2) $
        }
  if not keyword_set(config_fil) then return,cfg

  tags = tag_names(cfg)
  ntags = (size(tags,/dim))[0]

  readcol,config_fil,val,descript,format='a,a',delimiter='|',/silent
  descript = strtrim(descript,2)
  nval = n_elements(val)

  for ii=0,ntags-1 do begin
     mtch = where(tags[ii] eq strupcase(descript),nmtch)

     if nmtch eq 0 then begin
        print,'sdss_genprof_config(): value not in config file ',tags[ii]
        continue
     endif 
     
     if nmtch gt 1 then $
        stop,'sdss_genprof_config(): multiple values in config file ',tags[ii]
     
     mtch = mtch[0]

     if n_elements(cfg.(ii)) eq 2 then begin
        prs = strsplit(val[mtch],/extract)
        cfg.(ii) = [float(prs[0]),float(prs[1])]
     endif else cfg.(ii) = float(val[mtch]) ; ndv should be long

  endfor                        ; loop ii=ntgas

  return,cfg
end                             ; sdss_genprof_config()


function sdss_genprof_setline, wrest
  if n_params() ne 1 then begin
     print,'sdss_genprof_setline(wrest)'
     return,-1
  endif 
  
  strct = create_struct(x_setline(wrest),$ ; includes fval
                        'ID_SYS',0L,$      ; 0 to 2*config_strct.nabs
                        'NCOMP',0) 
  return, strct 

end                             ; sdss_genprof_setline()

function sdss_genprof_mc, dblt_name, config_strct, seed=seed, oseed=oseed
  ;; Generate all the structures necessary
  if n_params() ne 2 then begin
     print,'Syntax - sdss_genprof_mc(dblt_name, config_strct, [seed=, oseed=])'
     return,-1
  endif 

  if keyword_set(seed) then oseed = seed
  if size(dblt_name,/type) eq 8 then $
     dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)

  rng = 2*lindgen(config_strct.nabs)         ; to get both doublets

  vp_tmplt = sdss_genprof_setline(dblt.wvI)        ; 0 to config_strct.ndv
  vpstrct = replicate(vp_tmplt,config_strct.nabs*2) ; each line
  vpstrct[rng+1] = sdss_genprof_setline(dblt.wvI)

  ;; Set System ID numbers
  vpstrct[rng].id_sys = lindgen(config_strct.nabs) ; "system" number
  vpstrct[rng+1].id_sys = vpstrct[rng].id_sys 

  ;; Generate and store random z, N, b for each doublet line
  num = randomu(oseed,config_strct.nabs*3,/double) ; want larger swath
  num2 = num[0:config_strct.nabs-1]
  vpstrct[rng].zabs = num2*(config_strct.zlim[1]-config_strct.zlim[0]) $
                      + config_strct.zlim[0]
  vpstrct[rng+1].zabs = vpstrct[rng].zabs

  num2 = num[config_strct.nabs:(2*config_strct.nabs-1)]   
  vpstrct[rng].n = num2*(config_strct.nlim[1]-config_strct.nlim[0]) $
                   + config_strct.nlim[0]
  vpstrct[rng+1].n = vpstrct[rng].n
  
  num2 = num[(2*config_strct.nabs):(3*config_strct.nabs-1)] 
  vpstrct[rng].b = num2*(config_strct.blim[1]-config_strct.blim[0]) $
                   + config_strct.blim[0]
  vpstrct[rng+1].b = vpstrct[rng].b


  if keyword_set(config_strct.ndv) then begin
     ;; Add components
     if not keyword_set(config_strct.dvabs) then $
        stop,'sdss_genprof(): must set config_strct.ndv and config_strct.dvabs[]'

     vpstrct[rng].ncomp = randomu(oseed,config_strct.nabs,$
                                  poisson=config_strct.ndv) < $
                          config_strct.ndv ; max out
     vpstrct[rng+1].ncomp = vpstrct[rng].ncomp
     
     sub = where(vpstrct[rng].ncomp gt 0,nsub) ; where to add?
     if nsub gt 0 then begin
        ;; Create array
        sub = rng[sub]
        ncomptot = total(vpstrct[sub].ncomp)
        vpstrct_sub = replicate(vp_tmplt,ncomptot*2) ; wvI, wvII
        rng = 2*lindgen(ncomptot)
        vpstrct_sub[rng+1] = sdss_genprof_setline(dblt.wvII)
        num = randomu(oseed,ncomptot*3,/double) ; N, b, z
        inum = 0

        ;; Loop over systems to add components too
        for ii=0,nsub-1 do begin
           ss = sub[ii]
           rng = 2*lindgen(vpstrct[ss].ncomp)
           if ii gt 0 then $
              rng = rng + 2*total(vpstrct[sub[0:ii-1]].ncomp)
           vpstrct_sub[rng].id_sys = ss / 2 ; map back to vpstrct.id_sys
           vpstrct_sub[rng+1].id_sys = vpstrct_sub[rng].id_sys

           ;; Set dv offset
           ;; z = dv/c*(1+z0) + z0
           num2 = num[inum:(inum+vpstrct[ss].ncomp-1)]
           inum++
           vpstrct_sub[rng].zabs = $
              (num2*(config_strct.dvabs[1]-config_strct.dvabs[0]) $
               + config_strct.dvabs[0])/2.998e5*$
              (1 + vpstrct[ss].zabs) + vpstrct[ss].zabs
           vpstrct_sub[rng+1].zabs = vpstrct_sub[rng].zabs
           
           ;; Set Doppler
           num2 = num[inum:(inum+vpstrct[ss].ncomp-1)]
           inum++
           vpstrct_sub[rng].b = num2*(config_strct.blim[1]-config_strct.blim[0]) $
                                + config_strct.blim[0]
           vpstrct_sub[rng+1].b = vpstrct_sub[rng].b
           
           ;; Total column originally selected is preserved; divided
           ;; between other components
           num2 = num[inum:(inum+vpstrct[ss].ncomp-1)]
           inum++
           jj = 0
           while jj lt vpstrct[ss].ncomp do begin
              vpstrct_sub[rng[jj]].n = $
                 num2[jj]*(vpstrct[ss].n-config_strct.nlim[0]) $
                 + config_strct.nlim[0]
              ;; Subtract it off to preserve total column
              if vpstrct_sub[rng[jj]].n lt vpstrct[ss].n then begin
                 vpstrct[ss].n = alog10(10.^vpstrct[ss].n - $
                                          10.^vpstrct_sub[rng[jj]].n)
                 ;; But stay above max
                 if vpstrct[ss].n lt config_strct.nlim[0] then begin
                    vpstrct_sub[rng[jj]].n = 0. ; exclude later
                    vpstrct[ss].n = vpstrct[ss+1].n ; undo
                    bad = 1
                 endif else bad = 0
              endif else bad = 1

              if bad then begin 
                 if keyword_set(debug) then $
                    print,'sdss_genprof() debug: truncating no. of components: ',$
                          string(vpstrct[ss].zabs,jj,vpstrct[ss].ncomp,$
                                 format='(f7.5,1x,i3,1x,i3)')
                 vpstrct[ss].ncomp = jj ; end while loop
                 vpstrct[ss+1].ncomp = vpstrct[ss].ncomp 
              endif
              jj++
           endwhile              ; loop jj < ncomp[ss]
           vpstrct[sub+1].n = vpstrct[sub].n
           vpstrct_sub[rng+1].n = vpstrct_sub[rng].n
           vpstrct_sub[rng].ncomp = vpstrct[ss].ncomp
           vpstrct_sub[rng+1].ncomp = vpstrct_sub[rng].ncomp
        endfor                  ; loop ii=config_strct.nabs

        ;; Append ones that have components instantiated
        gd = where(vpstrct_sub.n gt 0.,complement=bd)
        if gd[0] ne -1 then $
           vpstrct = [vpstrct,vpstrct_sub[gd]] ; don't sort!

     endif                      ; vpstrct.ncomp > 0

  endif                         ; config_strct.ndv > 0

  oseed = oseed[0]   ; not the array

  return, vpstrct
end                             ; sdss_genprof_mc()


function sdss_genprof_voigt, wave, vplin, config_strct, debug=debug
  ;; Take the input structure and make the actual profile
  if n_params() ne 3 then begin 
     print,'Syntax - sdss_genprof_voigt(wave, vplin, config_strct, [/debug])'
     return,-1
  endif 

  npix = (size(wave,/dim))[0]
  vpfx0 = replicate(1.,npix)

  nlin = (size(vplin.wrest,/dim))[0]
  wobs = vplin.wrest*(1+vplin.zabs)
  srt = sort(wobs)
  wobs = wobs[srt]
  wlo = wobs*(1 - 5*vplin[srt].b*3.33556e-6) ; auto-bounds
  whi = wobs*(1 + 5*vplin[srt].b*3.33556e-6)
  
  ;; Going to call x_voigt() only for the snippets that require it and
  ;; make sure no overlap
  istrt = 0
  iend = istrt
  for ii=0,nlin do begin                        ; must include very last line
     if ii eq 0 then if nlin gt 1 then continue ; must process if nlin=1
     if ii gt 0 and ii lt nlin-1 then begin
        if (whi[ii-1] ge wlo[ii]) then begin
           iend = ii
           if ii lt nlin-1 then continue ; else proceed
        endif                            ; keep buffering
     endif 

     dum = min(wave-wlo[istrt],imn,/absolute)
     dum = min(wave-whi[iend],imx,/absolute)
     if keyword_set(debug) then $
        print,'sdss_genprof_voigt() debug: ii = '+strtrim(ii,2)+$
              '; wlo[istrt] = '+string(wlo[istrt],format='(f7.2)')+$
              '; whi[iend] = '+string(whi[iend],format='(f7.2)')

     if imn eq imx then begin
        ;; Make sure enough pixels to not crash
        if (imn gt 0) and (imx lt npix-1) then begin
           imn = imn-1 
           imx = imx+1
        endif else begin
           if imn eq 0 then imx = imx+2
           if imx eq npix-1 then imn = imn-2
        endelse 
     endif
     ;; _extra include tau=
     vpfx0[imn:imx] = x_voigt(wave[imn:imx],vplin[srt[istrt:iend]],fwhm=2,$
                             /nosmooth,_extra=extra)

     ;; Reset for next batch
     istrt = iend+1
     iend = istrt
  endfor                        ; loop ii=nlin+1

  ;; Check
  bd = where(finite(vpfx0) ne 1)
  if bd[0] ne -1 then begin
     print,'sdss_genprof_voigt() WARNING: non-finite Voigt profile'
     vpfx0[bd] = 1.              ; just reset
  endif 
  
  ;; Better to smooth with convol, assuming constant resolution (not
  ;; quite right) and knowing constant in velocity space; 162 to 136
  ;; km/s 
  ;; In spSpec*.fit 6th extension, there is a dispersion
  ;; structure, which median resolution in pixels is 0.877945 for most
  ;; spectra (so that fraction times the 69 km/s log binning)
  nkpix = 31                                   ; some odd number, 15*69 km/s = 1035 km/s
  velo = findgen(nkpix)*config_strct.dvelo   ; km/s
  kernel = exp(-0.5*((velo - config_strct.dvelo*0.5*nkpix)/(140./2.354))^2)
  kernel = kernel/total(kernel) ; normalize
  vpfx = convol(vpfx0,kernel,/edge_truncate) 

  return,vpfx
end                             ; sdss_genprof_voigt()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function sdss_genprof,wave,dblt_name,config_strct,cleanspec,$ 
                      seed=seed,oseed=oseed,calcew=calcew,$
                      conti=conti, vpstrct=vpstrct, debug=debug
  if n_params() ne 4 then begin
     print,'Syntax - '+ $
           'sdss_genprof(wave, dblt_name, config_strct, cleanspec,'+$
           '            [seed=,oseed=,/calcew,conti=,/debug]'
     return,-1
  endif                         ; param prompt

  ;; Default Values
  if keyword_set(seed) then oseed = seed
  if size(dblt_name,/type) eq 8 then $
     dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)

  ;; Generate profiles
  if not keyword_set(vpstrct) then $
     vpstrct = sdss_genprof_mc(dblt, config_strct, seed=oseed, oseed=oseed) $
  else if (size(vpstrct,/dim))[0] lt 2*config_strct.nabs then $
     stop,'sdss_genprof(): too few profiles in input vpstrct'

  ;; Generate Voigt profile and convolve with un-normalize spectrum
  vpfx = sdss_genprof_voigt(wave, vpstrct, config_strct, debug=debug)

  ;; Save
  spectrum = cleanspec
  spectrum[*,0] = cleanspec[*,0]*vpfx
  spectrum[*,2] = cleanspec[*,2]*sqrt(vpfx)


  if keyword_set(debug) then begin
     ;; Print
     print,'sdss_genprof() debug: synthetic line summary'
     print,'#','Wobs','zabs','N','b','Ncomp',$
           format='(a3,1x,a7,1x,a7,1x,a5,1x,a5,1x,a5)'
     rng = 2*lindgen(config_strct.nabs)
     wobs = vpstrct.wrest*(1.+vpstrct.zabs)
     printcol,rng,wobs[rng],vpstrct[rng].zabs,vpstrct[rng].n,$
              vpstrct[rng].b,vpstrct[rng].ncomp,$
              format='(i3,1x,f7.2,1x,f7.5,1x,f5.2,1x,f5.2,1x,i3)'
     print,''

     ;; Figure out good number to replicate as indication of centroids
     gd = where(vpfx lt 1.)
     nlin = n_elements(vpstrct)/2
     fxabs = replicate(median(spectrum[gd,0]),nlin)
     wobs = wobs[2*lindgen(nlin)]
     ion = dblt.ion + ' '+strtrim(round(dblt.wvI),2)

     ;; Plot
     clr = getcolor(/load)
     subttl = string((1.+config_strct.zlim)*dblt.wvI,$
                     format='(2x,2(f7.2,1x))')
     if keyword_set(conti) then $
        x_splot,wave,spectrum[*,0],psym1=10,$
                ytwo=spectrum[*,2],psym2=10,$
                ythr=vpfx*conti,psym3=-3,color3=clr.blue,$
                xfou=wobs,yfou=fxabs,color4=clr.limegreen,psym4=1,$ ;cross
                ymnx=[min(spectrum[*,2]),max(spectrum[*,0])],$
                title='Synthetic Spectrum and Voigt Profile'+subttl,/block,$
                lgnd=['Flux','Error','Voigt',ion] $
     else x_splot,wave,cleanspec[*,0],psym1=10,color1=clr.gray,$
                  ytwo=cleanspec[*,2],psym2=10,color2=clr.magenta,$
                  ythr=spectrum[*,0],psym3=10,color3=clr.black,$
                  yfou=spectrum[*,2],psym3=10,color4=clr.red,$
                  xfiv=wobs,yfiv=fxabs,color5=clr.limegreen,psym5=1,$ ;cross
                  ymnx=[min(spectrum[*,2]),max(spectrum[*,0])],$
                  title='Cleaned and Synthetic Spectrum'+subttl,/block,$
                  lgnd=['Flux0','Error0','Flux','Error',ion]
  endif

  oseed = oseed[0]   ; not the array
  return,spectrum
end
