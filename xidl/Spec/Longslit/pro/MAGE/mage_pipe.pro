pro mage_pipe, mage, obj_id=obj_id, chk=chk ,illumflatfile=illumflatfile1, $
               clobber=clobber, sensfunc=sensfunc, singles=singles, NOCR=NOCR
      
  ;; setup directories
  IF NOT FILE_TEST('Object', /DIR) THEN FILE_MKDIR, 'Object'
  IF NOT FILE_TEST('Final', /DIR)  THEN FILE_MKDIR, 'Final'
  IF NOT FILE_TEST('FSpec', /DIR)  THEN FILE_MKDIR, 'FSpec'

  ;; ???? Bug. We need to be reading these in a slit specific manner
  tset_slits = mrdfits("Orders.fits",1)
  ordr_str=mrdfits("OStr_mage.fits", 1)
  pixflatfile = 'Flat/Pixflat.fits'


  ;; Pick out science targets
  isci = where(strtrim(mage.exptype) EQ 'SCIENCE' OR $
               strtrim(mage.exptype) EQ 'STD'     OR $
               strtrim(mage.exptype) EQ 'BRIGHT' AND $
               mage.obj_id GE 0, nscience)
  sciframes = mage[isci]
  obj_ids = sciframes[uniq(sciframes.obj_id)].obj_id
  ;;if we only are reducing one object
  if (keyword_set(OBJ_ID)) then begin
     obj_ids = [obj_id]
  endif

  ;; Loop over objects
  for iobj=0, n_elements(obj_ids)-1 do begin
     objframes = sciframes[where(sciframes.obj_id EQ obj_ids[iobj], nframes)]
     print, "MAGE_PIPE: Processing object ", objframes[0].object
     ;; add some header cards
     master_hdr = objframes[0].hdr
     sxdelpar, master_hdr, "NAXIS1"
     sxdelpar, master_hdr, "CDELT1"
     sxdelpar, master_hdr, "CD1_1"
     sxdelpar, master_hdr, "CRVAL1"
     sxaddpar, master_hdr, "COMMENT", "mage_pipe: v0.3, Build June 2009"
     sxaddpar, master_hdr, "COMMENT", "object reduced as: " + $
               objframes[0].exptype
     pipe_time = "mage_pipe: Reduced "+systime()
     sxaddpar, master_hdr, "COMMENT", pipe_time
     
     for iexp=0, nframes-1 do begin

        IF NOT keyword_set(CLOBBER) THEN BEGIN
           tmp = strsplit(objframes[iexp].fitsfile, 'mage', /extract)
           outfil = 'Object/ObjStr'+strtrim(tmp[0])
           
           if (file_test(outfil) EQ 1) then begin

              tt = xmrdfits(outfil, 1)
              ;reflux the reduced spectra here?
              
              IF NOT KEYWORD_SET(sensfunc) THEN BEGIN
                 print, 'Refluxing with ' + fileandpath(objframes[iexp].sensfunc)
                 mage_flux, objframes[iexp].sensfunc, tt, rej=0.05
                 sxaddpar, master_hdr, "COMMENT", "mage_pipe: Re-fluxed with "+  fileandpath(objframes[iexp].sensfunc)
              ENDIF
              
              print, "Already Reduced, adding in ", outfil
              
              if (iexp EQ 0) then begin
                 allframes = [ tt ]
              endif else begin
                 allframes = [ allframes, tt ]
              endelse
              continue
           endif
        endif

        ;; Is this a standard or bright object??
        STD = strcompress(objframes[iexp].exptype,/rem) EQ 'STD'
      
        BRIGHT = strcompress(objframes[iexp].exptype,/rem) EQ 'BRIGHT'
        
        IF KEYWORD_SET(ILLUMFLATFILE1) THEN illumflatfile=illumflatfile1 $
        ELSE illumflatfile = 'Flat/Illumflat_' + $
                             strcompress(objframes[iexp].SLIT,/rem) + '.fits'
       
        ;; Bias subtract and flat field
        mage_proc, objframes[iexp].rawpath+objframes[iexp].fitsfile $
                   ,sciimg, sciivar, pixflatfile=pixflatfile $
                   ,illumflatfile=illumflatfile, hdr=scihdr
        arcfile = objframes[iexp].rawpath+objframes[iexp].arcfile
        
        ;; Trace the arc + science data
        mage_proc, arcfile, arcimg
        piximg=mage_makescipix(arcimg,sciimg,tset_slits,chk=chk,std=std, bright=bright)
        
        ;; Generate a 2D skymodel
        print, "Generating the 2D sky model"
        skyimage = mage_skymodel(sciimg, sciivar, piximg=piximg $
                                 ,tset_slits=tset_slits)

        ;; Reject cosmic rays
        IF NOT KEYWORD_SET(NOCR) THEN BEGIN
           print, 'Cosmic ray rejection'
           ;;IF KEYWORD_SET(FWHMSET) THEN sigma_psf = $
           ;;   djs_median(fwhmset.median)/2.35482D $
           ;;ELSE 
           sigma_psf = 3.0D/2.35482D
           ;;   Take the PSF width to be that of the spectral direction.  
           ;;   This prevents the routine from rejecting sky lines
           ;;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
           ;;    
           ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1] 
           ;;    PSFVALS[0]          1.   PSFVALS[0]
           ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
           psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
           crmask  = psf_reject_cr(sciimg-skyimage, sciivar, psfvals $
                                   , satmask = (sciimg-skyimage) GT 8d4)
           sciivar = sciivar*(crmask EQ 0)
           ;; Do we need to sky-subtract before CR rejection???
           ;; Do we need to mask these out in the sciimg??
        ENDIF        





        ;; Generate a wavelenth solution for nearest arc
        print, "Generating wavelength solution"
        print, "Arc file: ", arcfile
        mage_arc, arcfile, ordr_str=ordr_str, outarc=arcimg_fil, /clobber
        waveimg = xmrdfits(arcimg_fil)
        
        ;; Find and trace objects
        IF NOT KEYWORD_SET(STD) THEN BEGIN
           filstd = 'Object/ObjStr'+ $
                    strcompress(strsplit(objframes[iexp].stdfile $
                                         , 'mage',/extract),/rem)
           if (filstd EQ 'Object/Objstr') then begin
              print, "ERROR: Go back and enter your standard star file name into the GUI box (Sens func tab)"
              return
           endif
           filstd=filstd[0]
        ENDIF ELSE filstd=0
        obj_strct=mage_findobj(sciimg-skyimage,sciivar,waveimg,tset_slits $
                               ,filstd=filstd,chk=chk)
        
        ;; Optimal extraction
        tmp = strsplit(objframes[iexp].fitsfile, 'mage', /extract)
        outfil = 'Final/f_'+objframes[iexp].fitsfile
        mage_echextobj,sciimg,sciivar,scihdr,skyimage $
                       ,piximg,waveimg,tset_slits $
                       ,obj_strct,outfil=outfil $
                       ,box_rad=box_rad,STD=STD, BRIGHT=BRIGHT, CHK=CHK 
        
        ;; Set some items in the obj_strct
        obj_strct.img_fil = objframes[iexp].rawpath+objframes[iexp].fitsfile
        obj_strct.arc_fil = objframes[iexp].rawpath+objframes[iexp].arcfile
        obj_strct.exp=float(sxpar(scihdr,'EXPTIME'))
        
        ;; Flux calibrate spectrum
        IF NOT KEYWORD_SET(STD) THEN mage_flux, objframes[iexp].sensfunc, obj_strct, rej=0.05
        sxaddpar, master_hdr, "COMMENT", "mage_pipe: fluxed with "+  objframes[iexp].sensfunc


        objfil = 'Object/ObjStr'+strcompress(tmp[0],/rem)
        ;write out object file
        mwrfits, obj_strct, objfil, /create
        
        if (iexp EQ 0) then begin
           allframes = [ obj_strct ]
        endif else begin
           allframes = [ allframes, obj_strct ]
        endelse
     endfor

     ;; ???? At present there is no way to coadd multiple std exposures 
     
     IF NOT KEYWORD_SET(STD) AND NOT KEYWORD_SET(sensfunc) THEN BEGIN
        ;; Combine exposures
        files = objframes.fitsfile
        sxaddpar, master_hdr, "COMMENT",  "mage_pipe: co-added "+strjoin(files, ', ')
        res = 299792.458/4100.*objframes[0].slit
        
        IF KEYWORD_SET(singles) AND N_ELEMENTS(allframes) GT 15 THEN BEGIN
           print, 'Extracting each spectrum before co-add'
           nexp = N_ELEMENTS(allframes)
           FOR kk=0, nexp-1, 15 DO BEGIN
              mage_combspec, allframes[kk:kk+14], fspec
              outflux = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_s_'+strtrim(string(kk/15), 2)+'_F.fits'
              outerr  = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_s_'+strtrim(string(kk/15), 2)+'_E.fits'
              combname = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_s_'+strtrim(string(kk/15), 2)+'_comb.fits'
              mage_1dspec, fspec, outflux, outerr, combname $
                           , hdr=master_hdr, resvel=res 
           ENDFOR

        ENDIF
 

        ;Always do the whole combine
        mage_combspec, allframes, fspec
        
        outflux = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_F.fits'
        outerr  = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_E.fits'
        combname = 'FSpec/'+strcompress(objframes[0].object,/rem)+'_comb.fits'
        mage_1dspec, fspec, outflux, outerr, combname $
                     , hdr=master_hdr, resvel=res 
     ENDIF

  ENDFOR

END
