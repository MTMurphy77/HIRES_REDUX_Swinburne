;+ 
; NAME:
; esi_echfltsct   
;     Version 1.1
;
; PURPOSE:
;    Subtracts a scattered light model from the flat
;    Normalizes the Image and Outputs to same file with
;    updated headers
;
; CALLING SEQUENCE:
;   
;  esi_echfltsct, esi, /IFLAT, GAIN_RTO=, /CHK
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  Image with scattered light removed (e.g. Flats/FlatECH##_D[I].fits)
;
; OPTIONAL KEYWORDS:
;   /IFLAT    - Use internal flats 
;   GAIN_RTO= - Value for gain ratio for two amp mode (default:
;      calculate from image)
;   /CHK      - Display final image and a few other steps
;   MINWID=   - Minimum width between orders for measuring scattered light
;   /PINH     - Indicates pinhole slit (used in early 2000)
;   /NOSCAT   - No scattered light subtraction
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only set for 1x1 binning -- Many 'hard wired' numbers
;
; EXAMPLES:
;   esi_echfltsct, esi, 0.75, /CHK
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Aug-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;   17-Feb-2004 Added kludge for order #2
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;pro esi_echfltsct, esi, slit, IFLAT=iflat, GAIN_RTO=gain_rto, CHK=chk $
;                   , MINWID = minwid, PINH = pinh, CBIN = cbin, RBIN = rbin $
;                   , KLUDGEOR = kludgeor, IFU = ifu, BSPLINE = BSPLINE $
;                   , NRMMED = NRMMED, NOSCAT = NOSCAT

;
;  if  N_params() LT 2  then begin 
;      print,'Syntax - ' + $
;        'esi_echfltsct, esi, slit, /IFLAT, /CHK, GAIN_RTO=, MINWID= '
;      print, '       /PINH [v1.1]'
;      return
;  endif 
  
;  Optional Keywords


  ;;skyflat = 1
  ;;noscat = 1
  ;; Take out scattered light for all orders except the last one for 
  ;; night sky flats. 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mage_superflat, filenames, slitfile, pixfile $
                         , order_vec = order_vec, skyflat = SKYFLAT $
                         , NOSCAT = NOSCAT, ILLUM = ILLUMFINAL
                         
  
  ;;scattfile = strtrim(archive_path, 2)+'scatt_template.fits'
  chk = 0
  nflat = n_elements(filenames)
  IF nflat LT 1 THEN message, 'No flat field files given in mage_superflat'

  ;; Irrelevant for only 1 or 2 files
  if (nflat LE 2) then sigrej = 1.0 $ 
  else if (nflat EQ 3) then sigrej = 1.1 $
  else if (nflat EQ 4) then sigrej = 1.3 $
  else if (nflat EQ 5) then sigrej = 1.6 $
  else if (nflat EQ 6) then sigrej = 1.9 $
  else sigrej = 2.0
  ;; dist from edge for scattered light
  IF NOT KEYWORD_SET(scatt_shift) THEN scatt_shift = 5L 

  ;; create ordermask from tset_slits structure
  tset_slits = xmrdfits(slitfile, 1)
  slitmask = long_slits2mask(tset_slits)
  ordermask = 0*slitmask
  ordermask[WHERE(slitmask GT 0)] = -slitmask[WHERE(slitmask GT 0)] + 21L
  slitmask_left  = long_slits2mask(tset_slits, xshift = -abs(scatt_shift))
  slitmask_right = long_slits2mask(tset_slits, xshift = abs(scatt_shift))
  scattmask = (slitmask_left EQ 0 AND slitmask_right EQ 0)
  ximg = long_slits2x(tset_slits, edgmask = edgmask, TOL_EDG = 3L)
  piximg = xmrdfits(pixfile, 0)
  ;; these are the same for all images
  ;;scatt_model = xmrdfits(scattfile, 0)
    
  ;; which orders are we flattening
  IF NOT KEYWORD_SET(order_vec) THEN BEGIN 
     order_vec = reverse(dindgen(15) + 6)
     norders = 15L
  ENDIF ELSE norders = n_elements(order_vec)
  
  IF KEYWORD_SET(SKYFLAT) THEN specsamp = 1.5D $
  ELSE specsamp = 5L ;; spectral direction sampling in pixels
  spatsamp = 0.05D ;; spatial sampling for illum func in fractional order widths
  
  ;; allocate an array 
  FOR ifile = 0L, nflat-1L DO BEGIN
     print, 'Now processing file: '+ $
            strtrim(filenames[ifile], 2)+', #'+strtrim(ifile+1, 2)+' of ' $
            +strtrim(nflat, 2)
     mage_proc, filenames[ifile], image
     ;; allocate image array
     IF ifile EQ 0 THEN BEGIN
        dims = size(image, /dim)
        nx = dims[0]
        ny = dims[1]
        xarr = findgen(nx) # replicate(1.0, ny)
        imgarr   = make_array(dimension = [dims, nflat], /float)
        illumarr = make_array(dimension = [dims, nflat], /float)
        ;; this mask has 1=good, 0=bad (for avsigclip)
        maskarr = make_array(dimension = [dims, nflat], /float)
        ;; create a scattered light template using the first flat
        IF NOT KEYWORD_SET(NOSCAT) THEN $
           scatt_model = mage_scattlight(image, tset_slits)
     ENDIF
     illum_flat = fltarr(nx, ny)
     maskarr[*, *, ifile] = (ordermask GT 0.0) AND EDGMASK EQ 0
     scattpix = WHERE(scattmask AND image GT 1.0 AND image LT 5d4 AND $
                      (xarr GT 0 AND xarr LT (nx-1)))
     ;; fit for the amplitude of scattered light for this flat using
     ;; the spatial template
     IF NOT KEYWORD_SET(NOSCAT) THEN BEGIN
        scatt_amp = total(image[scattpix]*scatt_model[scattpix])/ $
                    total(scatt_model[scattpix]^2)
        image = image - scatt_amp*scatt_model
     ENDIF
     ;; now normalize out the spectral and spatial depdendence
     FOR iorder = 0L, norders-1L DO BEGIN
        print, 'Generating flat field for order '+strtrim(order_vec[iorder], 2)
        inorder = where(ordermask EQ order_vec[iorder], nord)
        fitpix = where(ordermask  EQ order_vec[iorder]  $
                       AND finite(image)  $
                       AND abs(image) LE 6.0d4 $
                       AND NOT EDGMASK, nfit)
        ;;ximg GT 0.05  $
        ;;AND ximg LT 0.95, nfit)
        psort = sort(piximg[fitpix])
        med_spec_width = round(double(nfit)*double(specsamp)/double(ny))
        normspec1 = djs_median(image[fitpix[psort]], Width = med_spec_width $
                               , bound = 'reflect')
        ;; now smooth with a gausian 
        sig_res_spec = med_spec_width/15L
        nhalf =  long(sig_res_spec)*4L
        xkern = dindgen(2*nhalf+1)-nhalf
        kernel = gauss1(xkern, [0.0, sig_res_spec, 1.0])
        normspec = convol(normspec1, kernel, /edge_truncate)
        sset_spec = bspline_iterfit(piximg[fitpix[psort]], normspec $
                                    , upper = 3, lower = 3, yfit = yfit $
                                    , bkspace = 0.7d, maxiter = 20 $
                                    , maxrej = 10, /silent)
        image[inorder] = image[inorder]/bspline_valu(piximg[inorder], sset_spec)
        ;; now fit out average illuimination function across the order
        fitpix2 = where(ordermask EQ order_vec[iorder]  $
                        AND finite(image)  $
                        AND abs(image) LE 6.0d4, nfit2)
        psort2 = sort(ximg[fitpix2])
        med_spat_width = round(double(nfit2)*double(spatsamp))
        normspat1 = djs_median(image[fitpix2[psort2]], width = med_spat_width $
                               , bound = 'reflect')
        ;; now smooth with a gausian 
        sig_res_spat = med_spat_width/15L
        nhalf =  long(sig_res_spat)*4L
        xkern = dindgen(2*nhalf+1)-nhalf
        kernel = gauss1(xkern, [0.0, sig_res_spat, 1.0])
        normspat = convol(normspat1, kernel, /edge_truncate)
        statinds = WHERE(ximg[fitpix2[psort2]] GT 0.1 AND $
                         ximg[fitpix2[psort2]] LT 0.9)
        moms = moment(normspat[statinds])
        mean = moms[0]
        normspat = normspat/mean
        sset_spat = bspline_iterfit(ximg[fitpix2[psort2]], normspat $
                                    , upper = 3, lower = 3, yfit = yfit2 $
                                    , bkspace = spatsamp/50.0d $
                                    , maxiter = 20, maxrej = 10 $
                                    , /silent)
        ;; don't allow super big corrections to protect edges
        normbsp = (bspline_valu(ximg[inorder], sset_spat)) < 2.0
        normbsp = normbsp > 0.5
        image[inorder] = image[inorder]/normbsp
        illum_flat[inorder] = normbsp
     ENDFOR
     imgarr[*, *, ifile] = image
     illumarr[*, *, ifile] = illum_flat
     IF KEYWORD_SET(CHK) THEN BEGIN 
        thismask = fltarr(nx, ny)
        FOR kk = 0L, norders-1L DO BEGIN
           ii = where(ordermask EQ order_vec[kk] AND maskarr[*, *, ifile])
           thismask[ii] = 1L
        ENDFOR
        xatv, thismask*image, min = 0.9, max = 1.1, /block
     ENDIF
  ENDFOR
  IF size(imgarr, /n_dimensions) EQ 3 THEN BEGIN
     imgfinal = djs_avsigclip(imgarr, 3, sigrej = sigrej $
                              , maxiter = maxiter, inmask = (maskarr EQ 0) $
                              , outmask = outmask)
     illumfinal = djs_avsigclip(illumarr, 3, sigrej = sigrej $
                                , maxiter = maxiter, outmask = outmask)
  ENDIF ELSE BEGIN
     imgfinal = imgarr
     illumfinal = illumarr
  ENDELSE

  return, imgfinal
  
END
  


