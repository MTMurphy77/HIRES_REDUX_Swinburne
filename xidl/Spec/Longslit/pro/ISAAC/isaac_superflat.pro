;PRO isaac_superflat, flatfiles, flatdarkfiles,  pixflatfile $
;;, darkfiles = darkfiles, badpixmaskfile = badpixmaskfile $

, illumflatfile $



;                   , slitfile = slitfile $
;                    , objfile = objfile $
;                       , verbose = verbose, indir = indir $
;                       , tempdir = tempdir $
;                       , use_illum = use_illum $
;                       , use_pixel = use_pixel $
;                       , npoly = npoly, CHK = CHK $
;                       , _EXTRA = extra 

tempdir = 'Temp/'
path_dom = '/Users/joe/ISAAC_raw/nite3/'
dom_pre = 'ISAAC_SWS_MR_FLATS_' 

path_twi = path_dom
;path_twi = '/Users/joe/ISAAC_raw/nite2/'
twi_pre = 'ISAAC_SWS_SCI_'
idom = [7, 9, 11] ;; centeral wavelenght is: 1.6520
idrk = [8, 10, 12]
itwi = [2, 3, 4] ;; 1.5880 ;; 1st is 10s, last one we use as dark
domefiles = path_dom + dom_pre + string(idom, FORMAT = '(I4.4)') + '.fits'
darkfiles = path_dom + dom_pre + string(idrk, FORMAT = '(I4.4)') + '.fits'
ndome = n_elements(domefiles)
tempdome = tempdir + 'TempDome-' + fileandpath(domefiles)

twifiles  = path_twi + twi_pre + string(itwi, FORMAT = '(I4.4)') + '.fits'
darktwifile = path_twi + twi_pre + string(5, FORMAT = '(I4.4)') + '.fits'
ntwi = n_elements(twifiles)
temptwi = tempdir + 'TempTwi-' + fileandpath(twifiles)

;; Process the DARK + FLAT + FLAT + DARK sequence for Domeflats
FOR ii = 0L, ndome-1L DO BEGIN
   isaac_proc, domefiles[ii], flat1, ivar_flat1, hdr = scihdr
   isaac_proc, darkfiles[ii], dark1, ivar_dark1, hdr = hdr
   ;; Average darks and flats
   sig2_flat = (ivar_flat1 GT 0)/(ivar_flat1 + (ivar_flat1 EQ 0))
   sig2_dark = (ivar_dark1 GT 0)/(ivar_dark1 + (ivar_dark1 EQ 0))
   sig2_diff = sig2_dark + sig2_flat
   ivar_diff = (ivar_dark1 GT 0)*(ivar_flat1 GT 0)*(sig2_diff GT 0.0)/(sig2_diff + (sig2_diff EQ 0))
   mwrfits, float(flat1-dark1), tempdome[ii], scihdr, /create
   mwrfits, float(ivar_diff), tempdome[ii]
ENDFOR

dims = size(flat1, /dim)
nx = dims[0]
ny = dims[1]
tset_slits = isaac_slitset(nx, ny)
slitmask = long_slits2mask(tset_slits)
slitfile = 'slits-ISAAC.fits'
mwrfits, slitmask, slitfile, /create
mwrfits, tset_slits, slitfile

islit = where(slitmask EQ 1)
;; Use the last twiflat as the 'dark'
isaac_proc, darktwifile, darktwi, ivar_dtwi, hdr = hdr
FOR ii = 0L, ntwi-1L DO BEGIN
   isaac_proc, twifiles[ii], img1, ivar1, hdr = scihdr
   ;; Average darks and flats
   sig2_img1 = (ivar1 GT 0)/(ivar1 + (ivar1 EQ 0))
   sig2_darktwi = (ivar_dtwi GT 0)/(ivar_dtwi + (ivar_dtwi EQ 0))
   sig2_diff = sig2_img1 + sig2_darktwi
   ivar_diff = (ivar1 GT 0)*(ivar_dtwi GT 0)*(sig2_diff GT 0.0)/(sig2_diff + (sig2_diff EQ 0))
   diff = img1 - darktwi
   ;djs_iterstat, diff[islit], mean = mean, median = median, sigma = sigma
   ;ilines = where(abs(diff[islit]) GE (abs(mean) + 3.0*sigma))
   ;djs_iterstat, diff[islit[ilines]] $
   ;              , mean = mean2, median = median2, sigma = sigma2
   ;IF mean2 GT 0.0 THEN sign = 1.0 ELSE sign = -1.0
   mwrfits, diff, temptwi[ii], scihdr, /create
   mwrfits, float(ivar_diff), temptwi[ii]
ENDFOR


;; Create a wavefile for input to long_superflat
isaac_proc, darktwifile, arcimg, arcivar, hdr = hdr
scifile = 'sci-ISAAC.fits'
waveimg = isaac_waveimg(arcimg, arcivar, tset_slits, hdr, scifile $
                        , pixset = pixset)
wavefile = 'wave-ISAAC.fits'
mwrfits, waveimg, wavefile, /create 
mwrfits, pixset, wavefile 

;; output files
pixflatfile = 'pixflat-' + fileandpath(domefiles[0])
illumflatfile = 'illumflat-' + fileandpath(twifiles[0])


;tempfiles = tempdome
;use_pixel = [lonarr(ndome) + 1L]
;use_illum = [lonarr(ndome) + 1L]
;long_superflat, tempfiles, pixflatfile, illumflatfile $
;                , slitfile = slitfile, wavefile = wavefile $
;                , use_illum = use_illum, use_pixel = use_pixel $
;                , tempdir = tempdir, slitsamp = 5.0, /SOFI ;;, /CHK
RUN = 1
IF KEYWORD_SET(RUN) THEN BEGIN
   tempfiles = [tempdome, temptwi]
   use_pixel = [lonarr(ndome) + 1L, lonarr(ntwi)]
   use_illum = [lonarr(ndome), lonarr(ntwi) + 1L]
   long_superflat, tempfiles, pixflatfile, illumflatfile $
                   , slitfile = slitfile, wavefile = wavefile $
                   , use_illum = use_illum, use_pixel = use_pixel $
                   , tempdir = tempdir, slitsamp = 20.0, /SOFI ;;, /CHK
ENDIF

isaac_proc, twifiles[0], twiimg_illum, illumflatfile = illumflatfile
isaac_proc, darktwifile, darkimg_illum, illumflatfile = illumflatfile

isaac_proc, twifiles[0], twiimg_pix, pixflatfile = pixflatfile 
isaac_proc, darktwifile, darkimg_pix, pixflatfile = pixflatfile


isaac_proc, twifiles[0], twiimg_both $
            , pixflatfile = pixflatfile, illumflatfile = illumflatfile
isaac_proc, darktwifile, darkimg_both $
            , pixflatfile = pixflatfile, illumflatfile = illumflatfile

isaac_proc, twifiles[0], twiimg
isaac_proc, darktwifile, darkimg

diff_illum = twiimg_illum-darkimg_illum
diff_pix = twiimg_pix-darkimg_pix
diff_both = twiimg_both - darkimg_both
diff = twiimg-darkimg





END
