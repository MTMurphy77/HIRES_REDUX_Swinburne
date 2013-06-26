;+ 
; NAME:
; sdss_completeness
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
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   21 Sep 2011  Created by kLC
;-
;------------------------------------------------------------------------------
@sdss_fndlin                    ; resolve sdss_fndlin_srch()
@sdss_dblfitconti               ; resolve sdss_dblfitconti_fithybrid()
@sdss_civsearch                 ; resolve sdss_civsearch_srch()
@sdss_genprof                   ; resolve loads of things

function sdss_completeness_compare, mcstrct, civstrct

end                             ; sdss_completeness_compare()


pro sdss_completeness_los, dblt_name, spec_fil, config_fil, $
                           eigbasis=eigbasis, pca_fil=pca_fil, $
                           pca_head=pca_head, debug=debug, save_all=save_all, $
                           seed=seed, oseed=oseed, $
                           _extra=extra

  sdssdir = sdss_getsdssdir()
  ;; Real workhorse on a given sightline and parameters

  ;; Params
  if keyword_set(seed) then oseed = seed ; will return to same variable
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)
  if size(config_fil,/type) eq 8 then config_strct = config_fil $
  else config_strct = sdss_genprof_config(config_fil)

  if not keyword_set(eigbasis) then $ ; better not to read in every time
     eigbasis = xmrdfits(getenv("SDSSPATH")+"/eigenspectra/eigSpec_qso_all.fit",$
                         0, /silent)
  if not keyword_set(pca_fil) then $
     pca_fil = getenv('XIDL_DIR')+'/SDSS/PCA/pca_base2000.fits'
  if size(pca_fil,/type) eq 7 then $
     pca = xmrdfits(pca_fil, 0, pca_head, /silent) $
  else begin  
     pca = pca_fil 
     if not keyword_set(pca_head) then $
        stop,'sdss_completeness_los: must set pca_head' 
  endelse 

  cflg_eig = sdss_getcflg(/eig)
  cindx_eig = sdss_getcflg(/eig,/index)
  cflg_hyb = sdss_getcflg(/hyb)
  cindx_hyb = sdss_getcflg(/hyb,/index)

  if keyword_set(save_all) then begin
     ;; Set all the file names necessary
     cleanspec_fil = sdss_getname(spec_fil,/spec,dir=cleandir,/clean)
     cleanspec_fil = cleandir + cleanspec_fil
  endif 

  cstrct0 = xmrdfits(sdssdir+cstrct_fil,1,/silent)
  if cstrct0.cflg ne cflg_hyb then $
     stop,'sdss_completeness_los: completeness test geared toawards hybrid conti only'

  tmp = sdss_measuresnr('blank', wvlim_obs=wvlim_obs, /no_snr, dblt_name=dblt, $
                        zqso=cstrct0.z_qso) 


  ;; Read in spectrum, original abslin structure; prepare to pass around 
  parse_sdss, sdssdir+spec_fil, flux, wave, head=hdr, npix=npix, sig=sigma
  wave_rest = wave/(1.+cstrct0.z_qso)
  mask0 = replicate(1,npix)     ; should never change
  bdpix = where(sigma eq 0.)
  if bdpix[0] ne -1 then mask0[bdpix] = 0


  ;; Figure out median resolution in region sensitive to doublet
  resstrct = xmrdfits(sdssdir+spec_fil[ss],6,/silent) 
  gd = where(wave ge wvlim_obs[0] and wave le wvlim_obs[1], ngd) 
  if size(resstrct,/type) eq 8 then $
     config_strct.dvelo = median(resstrct[gd].dispersion) * config_strct.pixscale $
  else begin
     config_strct.dvelo = config_strct.pixscale
     print,'sdss_completeness_los: resolution ext=6 DNE ',spec_fil[ss]
  endelse 

  config_struct.zlim = wvlim_obs/dblt.wvI - 1. 
  config_struct.nabs = round(config_struct.dndz*$
                             (config_struct.zlim[1]-config_struct.zlim[0])) > 1
  pixmin = gd[0] - 10 > 0       ; buffer search zone
  dum = min(wave-dblt.wvII*(1+config_struct.zlim[1]), pixmax, /abs)
  pixmax = pixmax + 10 < npix
  

  ;; Generate a whole lot of profiles
  vpstrct = sdss_genprof_vpstrct(dblt, config_strct, seed=oseed, oseed=oseed, /max)

  ;; Loop over redshift
  uber_count = 0L
  istart = 0L
  istop = config_strct.nabs - 1L
  while uber_count lt config_strct.nabstot do begin
     ;; Select centroids to scrub based on observed EW; only do once 
     ;; in awhile (every 100 loops below)
     ;; _extra includes ewobs_lim=, /clobber
     cleanspec = sdss_cleanspec(wave, flux, sigma, cstrct0, config_strct.frac_rm, $
                                header=hdr, seed=oseed, oseed=oseed, debug=debug, $
                                cleanspec_fil=cleanspec_fil, _extra=extra) 

     ;; Randomly generate many profiles to be input into the same
     ;; cleaned spectrum
     for vv=0, 99 do begin
        ;; Take substrcture based on the fact that the systems in
        ;; vpstrct should be ordered consecutively
        if istop gt config_strct.nabstot then begin
           ;; New pool of many profiles
           vpstrct = sdss_genprof_vpstrct(dblt, config_strct, $
                                          seed=oseed, oseed=oseed, /max)
           ;; Increment ID number 
           vpstrct.id_sys = vpstrct.id_sys + uber_count
           istart = 0
           istop = config_strct.nabs - 1
        endif 
        gd = where(vpstrct.id_sys ge vpstrct_sub[istart].id_sys and $
                   vpstrct.id_sys le vpstrct_sub[istop].id_sys, nsub)
        vpstrct_sub = vpstrct[gd]
        uber_count = uber_count + config_strct.nabs 

        ;; Once in a great while throw in some huge absorber near
        ;; systemic of QSO 

        ;; Insert into spectrum Save spectrum in SDSS format [npix, 5] so
        ;; that 0th = flux, 1st = pixel mask (scrubbed lines); 2nd =
        ;; error, 3rd = voigt profile
        newspec = sdss_genprof(wave, dblt, config_struct, cleanspec, $
                               seed=oseed, oseed=oseed, $
                               conti=cstrct0.conti[0:cstrct0.npix-1, cindx_hyb], $
                               /calcew, vpstrct=vpstrct_sub, mcstrct=mcstrct_sub, $
                               debug=debug)
        snr = median(newspec[pixmin:pixmax, 0]/newspec[pixmin:pixmax, 2])

        ;; Fit eigen-conti and find centroids in just the region
        ;; sensitive to doublet
        contihdr = hdr
        neweig = eigqsoconti(wave_rest[cstrct0.ipix0:*], $
                             flux[cstrct0.ipix0:*], $
                             sigma[cstrct0.ipix0:*], eigbasis, $
                             finalmask=finalmask, header=contihdr, $
                             /silent)
        eigconti0 = dblarr(cstrct0.npix, 3, /nozero) ; MUST set values
        eigconti0[cstrct.ipix0:*, 0] = neweig[*, 0]
        eigconti0[cstrct.ipix0:*, 1] = finalmask
        eigconti0[cstrct.ipix0:*, 2] = neweig[*, 1]
        if cstrct0.ipix0 gt 0 then begin
           eigconti0[0:cstrct0.ipix0-1] = 0
        endif 
        cstrcte = cstrct0
        cstrcte.cflg = cflg_eig
        cstrcte.conti[0:cstrct0.npix-1, cindx_eig] = eigconti0[*, 0]
        cstrcte.sigconti[0:cstrct0.npix-1, cindx_eig] = eigconti[*, 2]
        ;; Update structure and find lines
        cstrcte = sdss_fndlin_srch(cstrcte, cstrcte.cflg, $
                                   wave=wave[pixmin:pixmax], $
                                   flux=newspec[pixmin:pixmax, 0], $
                                   sigma=newspec[pixmin:pixmax, 2], $
                                   mask=mask0[pixmin:pixmax], $
                                   debug=debug) 

        ;; Fit hybrid-conti and find centroids in just the region
        ;; sensitive to doublet; mimic the call in sdss_dblfitconti
        hybconti = $
           sdss_dblfitconti_fithybrid(wave, newspec[*, 0], newspec[*, 2], $
                                      snr, dblt, eigconti0, cstrcte, $
                                      ;; Options
                                      silent=silent, $
                                      ;; Things passed out
                                      eigconti=eigconti, $ ; new one
                                      ;; Things passed in (and deep)
                                      eigbasis=eigbasis, pca_fil=pca, $
                                      pca_head=pca_head, $
                                      contihdr=contihdr) ; to be modified
        ;; Update structure and find lines
        cstrcth = cstrcte
        cstrcth.cflg = cflg_hyb
        cstrcth.conti[0:cstrct0.npix-1, cindx_hyb] = hybconti[*, 0]
        cstrcth.sigconti[0:cstrct0.npix-1, cindx_hyb] = hybconti[*, 2]
        cstrcth = sdss_fndlin_srch(cstrcth, cstrcth.cflg, $
                                   wave=wave[pixmin:pixmax], $
                                   flux=newspec[pixmin:pixmax, 0], $
                                   sigma=newspec[pixmin:pixmax, 2], $
                                   mask=mask[pixmin:pixmax], $
                                   debug=debug) 

        ;; Find candidate doublets; mimic call in sdss_civsearch
        ;; Make sure cflg set appropriately before
        civstr = sdss_civsearch_srch(wave, newspec[pixmin:pixmax, 0], $
                                     newspec[pixmin:pixmax, 2], dblt, cstrcth, $
                                     debug=debug, silent=silent, $
                                     qstrct_tmplt=qstrct_tmplt, $
                                     civstr_tmplt=civstr_def, $
                                     count=count, _extra=extra)

        ;; Compare to input doublets
        ;; mcstrct_sub
        if count gt 0 then begin
           mcstrct_sub = sdss_completeness_compare(mcstrct_sub,civstr)
        endif else begin
        endelse 

        ;; Save
        if keyword_set(mcstrct) then mcstrct = [mcstrct, mcstrct_sub] $
        else mcstrct = mcstrct_sub

     endfor

     if count_rand ge nrandmin and $
        (count_rand mod nrandmin) eq 0 then begin
        ;; Test if any column density bin converged
        ;; If so, change column density range

     endif                      ; count_rant >= nrandmin

  endwhile                      ; ubercount < config_strct.nabstot

  
end                             ; sdss_completeness_los

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_completeness, list_fil, config_fil=config_fil, debug=debug, save_all=save_all, processor=processor, _extra=extra


  eigbasis = xmrdfits(getenv("SDSSPATH")+"/eigenspectra/eigSpec_qso_all.fit", $
                      0, /silent)
 
end                             ; sdss_completeness
