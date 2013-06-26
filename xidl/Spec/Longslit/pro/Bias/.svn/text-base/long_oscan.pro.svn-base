; NAME:
;   long_oscan
;
; PURPOSE:
;   Overscan subtract a standard CCD.  There are separate subroutines
;   for LRIS, DEIMOS, BCS, GMOS and KAST.
;
; CALLING SEQUENCE:
;   long_oscan, filename, [ rawsub, rawivar, hdr= , $
;    gain= , rnoise= , image=, /verbose ]
;
; INPUTS:
;   filename   -  name of a image file
;
; OPTIONAL INPUTS:
;   verbose    - If set, then verbose output
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   rawsub     - Bias subtracted image
;   rawivar    - Inverse variance of image
;   hdr        - Fits header of raw image
;   gain       - vector of gains for each amp
;   rnoise     - vector or readnoise for each amps
;   image      - original data image
;   bin        - Binning of the CCD [x,y]
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   mrdfits
;   djs_avsigclip
; REVISION HISTORY:
;   10-Mar-2005  Written by J. Hennawi (UCB) and Scott Burles (MIT)
;-
;------------------------------------------------------------------------------

function gmos_badpixelmask, ncol, nrow, ccd_num

  mask = bytarr(ncol, nrow) + 1B
  
  xbin = long(2048/ncol)
  ybin = long(4608/nrow)
  
  if ccd_num EQ 1 then begin
     mask[88/xbin:89/xbin, *] = 0
     mask[772/xbin:775/xbin, 708*4/ybin:*] = 0
     mask[986/xbin:987/xbin, 768*4/ybin:*] = 0
     return, mask
  endif
  
  if ccd_num EQ 2 then begin
     mask[1252/xbin:1255/xbin, 820*4/ybin:*] = 0
     mask[1280/xbin:1281/xbin, 315*4/ybin:*] = 0
     return, mask
  endif
  
  ;; no bad columns for CCD 3?
  return, mask
end 
      

; add support for the Bok/B&C Spectrograph
pro bok_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
bin_col = long(strmid(ccdsum, 0, 1))
bin_row = long(strmid(ccdsum, 1, 1))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)

data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region too close
; to data 
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

; generate the bad-pixel mask - assumes no spectral binning!!!
mask = rawsub*0+1
;mask[405-1:408-1,*] = 0
;mask[507-1:507-1,*] = 0
;mask[538-1:538-1,*] = 0
;mask[650-1:650-1,*] = 0
;mask[681-1:681-1,*] = 0
;mask[706-1:706-1,*] = 0
;mask[767-1:771-1,*] = 0
;mask[897-1:898-1,*] = 0

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

rawsub = transpose(rawsub)
rawivar = transpose(rawivar*mask) ; zero out!

sxaddpar, hdr, 'NAXIS1', (size(rawsub,/dim))[0] ; note!
sxaddpar, hdr, 'NAXIS2', (size(rawsub,/dim))[1]
sxdelpar, hdr, 'biassec'
sxdelpar, hdr, 'trimsec'
sxdelpar, hdr, 'datasec'
sxdelpar, hdr, 'ccdsec'
sxdelpar, hdr, 'origsec'
RETURN
END

pro gmos_oscan, filename, flux, invvar, hdr = hdr, thdr = thdr $
                , gain = gain, rnoise = rnoise_vec, image = image $
                , VERBOSE = VERBOSE, bin = bin, CCDONLY = CCDONLY


ccd = [1, 2, 3]
IF NOT KEYWORD_SET(CCDONLY) THEN i = ccd[0] $
ELSE i = ccd[ccdonly-1]
;  no overscan for rows...
;  look for filename
filelist = lookforgzip(filename)
if filelist[0] EQ '' then begin
    print, 'Could not find file named ', filename
    return
endif
hdr0 = xheadfits(filelist[0], exten = 0)
rawdata = xmrdfits(filelist[0], i, hdrt, /fscale)

if size(hdr0, /tname) EQ 'INT' then begin
    print, 'having trouble with ', filename
    return
endif

detsize = sxpar(hdr0, 'DETSIZE')
hccds   = sxpar(hdr0, 'NCCDS')
hamps   = sxpar(hdr0, 'NAMPS')

if hamps NE 1 then begin
    print, 'Only ready for single amplifier reads'
    return
endif

if detsize NE '[1:6144,1:4644]' OR hccds NE 3 then begin
    print, 'Expecting 3 CCDs of total size: [1:6144,1:4644]'
    return
endif

rnoise = sxpar(hdrt, 'RDNOISE')
gain    = sxpar(hdrt, 'GAIN')

rnoise_vec = fltarr(3)
rnoise_vec[i-1L] = rnoise

IF KEYWORD_SET(VERBOSE) THEN print, "CCD Number: ", i, " Gain:", gain $
  , " Read Noise:", rnoise, format = '(a,i2,a,f6.3,a,f5.2)'
    
rawcol = (size(rawdata))[1]
rawrow = (size(rawdata))[2]

biassec = strcompress(sxpar(hdrt, 'BIASSEC'), /rem)
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
oscan_pix=bias_arr[1]-bias_arr[0] + 1L
colbin = long(2L*1024/(rawcol-oscan_pix))
rowbin = long(4L*1152/rawrow)
bin = [colbin, rowbin]
bin = reverse(bin)

ncol = 2048L/colbin
nrow = rawrow

colsplit = lindgen(64L/colbin) + long(strmid(sxpar(hdrt, 'BIASSEC'), 1)) -1
datacols = lindgen(ncol) + long(strmid(sxpar(hdrt, 'DATASEC'), 1)) - 1
oscan = djs_avsigclip(rawdata[colsplit, *], 1)
suboscan = rebin(transpose(oscan), ncol, nrow)

flux1 = (rawdata[datacols, *] - suboscan) * gain

invvar1 = 1.0/(abs(flux1 - sqrt(2.0)*rnoise) + rnoise^2)
mask = gmos_badpixelmask(ncol, nrow, i)
invvar1 = invvar1*mask

n0 = n_elements(hdr0)
nt = n_elements(hdrt)
hdr = strarr(max([n0, nt]), 4)
hdr[0:n0-1L, 0] = hdr0[0:n0-1L]
hdr[0:nt-1L, 1] = hdrt[0:nt-1L]
nccd = n_elements(ccd)
IF KEYWORD_SET(CCDONLY) THEN BEGIN
    flux = transpose(reverse(flux1, 1))
    invvar = transpose(reverse(invvar1, 1))
    RETURN
ENDIF

flux1 = rebin(flux1, ncol, nrow, nccd)
invvar1 = rebin(invvar1, ncol, nrow, nccd)
for iccd = 1, nccd-1 do begin
    i = ccd[iccd]
    rawdata = xmrdfits(filename, i, hdrt, /silent)
    hdrt    = xheadfits(filename, exten = i)
    rnoise = sxpar(hdrt, 'RDNOISE')
    rnoise_vec[i-1] = rnoise
    gain    = sxpar(hdrt, 'GAIN')
    IF KEYWORD_SET(VERBOSE) THEN $
      print, "CCD Number: ", i, " Gain:", gain, " Read Noise:", rnoise $
      , format = '(a,i2,a,f6.3,a,f5.2)'
    
    colsplit = lindgen(64L/colbin) +  $
      long(strmid(sxpar(hdrt, 'BIASSEC'), 1)) -1
    datacols = lindgen(ncol) + long(strmid(sxpar(hdrt, 'DATASEC'), 1)) -1
    oscan = djs_avsigclip(rawdata[colsplit, *], 1)
    suboscan = rebin(transpose(oscan), ncol, nrow)
    
    f = (rawdata[datacols, *] - suboscan) *gain
    
    invv = 1.0/(abs(f-sqrt(2.0)*rnoise) + rnoise^2)
    
;    badpixels from Adam's pipeline:
    
    mask = gmos_badpixelmask(ncol, nrow, i)
    invv = invv*mask
    
    flux1[*, *, i-1] = f
    invvar1[*, *, i-1] = invv
      
    hdr[0:nt-1L, iccd+1L] = hdrt
endfor


flux = gmos_mosaic(flux1)
invvar = gmos_mosaic(invvar1)

return
end

pro mmt_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
 endif


;; Parse the header
gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))

      nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

case strtrim(detector,2) of
   '*ccd34*': BEGIN
      rawsub = transpose(rawsub)
      rawivar = transpose(rawivar) 
;  different trim for red channel
      rawsub = rawsub[80:610, 160:1029] ; 125:574
      rawivar = rawivar[80:610, 160:1029]
      dims = size(image, /dim)
      bin= [round(2700./ncol),1]  ;; Might want to calculate bin[1] too (not sure this was done right)
   END
   'mmtredchan': BEGIN
      ;; No trimming
      ;; Binning?
      dims = size(image, /dim)
      bin= [1,round(1032./dims[1])]  ;; This might not be right (JXP: 06/21/2010)
   end
   else: stop
endcase

;;  different trim for upgraded red channel (using transposed coordinates)
;;  keep commented for now since as it will require changing
;;  wavelength solution etc  24-04-2009
;IF strmatch(detector, 'mmtredchan*') THEN BEGIN
;    rawsub  =  rawsub[*, 7:1006] 
;    rawivar = rawivar[*, 7:1006] 
;ENDIF


RETURN
END

pro kpno_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
bin_col = long(strmid(ccdsum, 0, 1))
bin_row = long(strmid(ccdsum, 1, 1))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

rawsub = transpose(rawsub)
rawivar = transpose(rawivar) 
RETURN
END

pro caha_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin_col = long(sxpar(hdr, 'CCDBINX'))
bin_row = long(sxpar(hdr, 'CCDBINY'))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'CCDSENS'))
rnoise = double(sxpar(hdr, 'CCDRON'))
detector = strcompress(sxpar(hdr, 'CCDNAME'), /rem)
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[2]-bias_arr[0]

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = bias_arr[0] - data_arr[0] + 1L
imagerow = data_arr[3]-data_arr[1] + 1L

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;; Pitch the first four pixels for both the red and blue side as these
;; are hot pixels. 

rawsub  = rawsub[4:*, *]
rawivar = rawivar[4:*, *]
bin = reverse(bin)
IF strcmp(detector, 'SITe#22b_14') THEN BEGIN ;; BLUE SIDE
   rawsub = transpose(rawsub)
   rawivar = transpose(rawivar) 
ENDIF ELSE IF strcmp(detector, 'SITe#20b_12') THEN BEGIN ;; RED SIDE
   rawsub  = reverse(transpose(rawsub), 2)
   rawivar = reverse(transpose(rawivar), 2)
ENDIF
RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro caha_22m_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan when windowed!
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
 endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin_col = long(sxpar(hdr, 'CCDBINX'))
bin_row = long(sxpar(hdr, 'CCDBINY'))
bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'CCDSENS'))
rnoise = double(sxpar(hdr, 'CCDRON'))
detector = strcompress(sxpar(hdr, 'CCDNAME'), /rem)
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
if data_arr[0] LT bias_arr[0] then stop ;; You can do overscan subtraction and need to trim

temp_image = image * gain  ;; Electrons
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2) ;; This is wrong!!

RETURN
END


pro imacs_oscan, filename, rawsub, rawivar, hdr = hdr $
                 , gain = gain, rnoise = rnoise, image = image $
                 , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin= [round(2700./ncol),1]  ;; Might want to calculate bin[1] too

gain = double(sxpar(hdr, 'EGAIN'))
rnoise = double(sxpar(hdr, 'ENOISE'))
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;; trim the spectra 
rawsub  = rawsub[*, 1000:*]
rawivar = rawivar[*, 1000:*]
;; bad pixel mask (red camera)
rawsub[1688, 2504:*] = 0
rawivar[1688, 2504:*] = 0
rawsub[248, 2476:*] = 0
rawivar[248, 2476:*] = 0
rawsub[977, *] = 0
rawivar[977, *] = 0

RETURN
END



pro dis_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
bin_col = long(sxpar(hdr, 'CCDBIN1'))
bin_row = long(sxpar(hdr, 'CCDBIN2'))
bin = [1, 1]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
detector = sxpar(hdr, 'DETECTOR')
datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]
imagerow = data_arr[3]

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biascols = lindgen(nbias-oscan_buffer) + imagecol + oscan_buffer
oscan = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
osub = replicate(1, imagecol) # oscan
temp_image = (image[0:imagecol-1L, 0:imagerow-1L] - osub)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

If strmatch(detector, '*blue*') THEN BEGIN
    rawsub = transpose(rawsub)
    rawivar = transpose(rawivar) 
ENDIF ELSE IF strmatch(detector, '*red*') THEN BEGIN
    rawsub  = reverse(transpose(rawsub), 2)
    rawivar = reverse(transpose(rawivar), 2)
    ;; mask bad column on red side
    rawsub[212:*, 773] = 0.0D
    rawsub[212:262, 770:777] = 0.0
    rawivar[212:*, 773] = 0.0D
    rawivar[212:262, 770:777] = 0.0
ENDIF ELSE message, 'Unrecognized detector'

RETURN

END

pro p200_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
endif

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]
ccdsum = strcompress(sxpar(hdr, 'CCDSUM'), /rem)
bin_col = long(strmid(ccdsum, 0, 1))
bin_row = long(strmid(ccdsum, 1, 1))

bin = [bin_col, bin_row]

gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RON'))
detector = sxpar(hdr, 'FPA')
ccd = sxpar(hdr, 'DETNAM')
datasec = strcompress(sxpar(hdr, 'TSEC1'), /rem)
biassec = strcompress(sxpar(hdr, 'BSEC1'), /rem)
data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[1]-bias_arr[0] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 5L
imagecol = data_arr[1]-data_arr[0] + 1L
icol = data_arr[0] -1L
imagerow = data_arr[3]-data_arr[2] + 1L
biascol = bias_arr[0]-1L
rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)
biascols = lindgen(nbias-oscan_buffer) + biascol + oscan_buffer
oscan1 = djs_avsigclip(image[biascols, 0:imagerow-1L], 1)
oscan = djs_median(oscan1, width = 200, boundary = 'reflect')
osub = replicate(1, imagecol) # oscan
temp_image = (image[icol:icol+imagecol-1L, *] - osub)*gain 
rawsub[0:imagecol-1L, *] = temp_image
rawivar[0:imagecol-1L, *] = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

If strmatch(detector, 'DBSP_BLUE') THEN BEGIN
   rawsub  = reverse(rawsub, 2)
   rawivar = reverse(rawivar, 2) 
   dims = size(rawsub, /dim)
   case ccd of 
      'CCD44-82': 
      else: begin
         IF dims[1] NE 1880 THEN BEGIN
            IF dims[1] GT 1880 THEN BEGIN 
               rawsub1 = fltarr(dims[0], 1880)
               rawivar1 = fltarr(dims[0], 1880)
               rawsub1 = rawsub[*, 550L:(1880+550-1L)]
               rawivar1 = rawivar[*,550L:(1880+550-1L)]
;          rawsub1 = rawsub[*, 788L:(788+1880-1L)]
;          rawivar1 = rawivar[*, 788:(788+1880-1L)]
               rawsub = rawsub1
               rawivar = rawivar1
               splog, 'WARNING: Truncating P200 blue side arrays to 1880'
            ENDIF ELSE BEGIN
               rawsub1 = fltarr(dims[0], 1880)
               rawivar1 = fltarr(dims[0], 1880)
               rawsub1[*, 0:dims[1]-1L] = rawsub
               rawivar1[*, 0:dims[1]-1L] = rawivar
               rawsub = rawsub1
               rawivar = rawivar1
               splog, 'WARNING: Padding P200 blue side arrays to 1880'
            ENDELSE
         endif
      end
   endcase
ENDIF ELSE IF strmatch(detector, 'DBSP_RED') THEN BEGIN
   rawsub  = transpose(rawsub)
   rawivar = transpose(rawivar)
ENDIF ELSE message, 'Unrecognized detector'

RETURN

END



pro mars_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

bin = [1,1]

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
endif

ncol = (size(image))[1]
nrow = (size(image))[2]

biassec = [long(strmid(sxpar(hdr, 'BIASSEC'), 1))-1L $
           , long(strmid(sxpar(hdr, 'BIASSEC'), 6))-1L]
postcol = biassec[1]-biassec[0] + 1L

rowtrim = [215, 544]
gain = double(sxpar(hdr, 'GAIN'))
rnoise = double(sxpar(hdr, 'RDNOISE'))
imagecol = ncol-postcol 
imagerow = rowtrim[1]-rowtrim[0] + 1L
post_buffer = 4L
rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biaspix = 5L
biascols = reverse(ncol-post_buffer-lindgen(biaspix))
goodrows = rowtrim[0]+ lindgen(imagerow)
biasrows = [rowtrim[0]+150L, rowtrim[0]+199L]
oscan = djs_median(image[biascols, biasrows[0]:biasrows[1]])
temp_image = (image[0:imagecol-1L, rowtrim[0]:rowtrim[1]] - oscan)*gain 
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

rawsub = transpose(rawsub)
rawivar = transpose(rawivar)
bin = reverse(bin)

RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro deimos_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin

bin = [1,1]

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
endif

ncol = (size(image))[1]
nrow = (size(image))[2]

datasec = xregtovec(sxpar(hdr,'DATASEC'))
postcol = ncol-datasec[1]

rowtrim = [datasec[2]-1,datasec[3]-1]
biasrows = rowtrim
;gain = double(sxpar(hdr, 'GAIN'))
;rnoise = double(sxpar(hdr, 'RDNOISE'))
gain = 1.25
rnoise = 2.4

imagecol = datasec[1]-datasec[0]+1L  ; ncol-postcol 
imagerow = rowtrim[1]-rowtrim[0] + 1L
post_buffer = 4L

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

biaspix = (ncol - datasec[1] - 10)
biascols = reverse(ncol-post_buffer-lindgen(biaspix))
goodrows = rowtrim[0]+ lindgen(imagerow)

;; Get the overscan
oscan = djs_median(image[biascols, biasrows[0]:biasrows[1]],1)
;; Subtract
temp_image = (image[0:imagecol-1L, rowtrim[0]:rowtrim[1]] - $
              (replicate(1.,imagerow)#oscan))*gain 
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;rawsub = transpose(rawsub)
;rawivar = transpose(rawivar)
;bin = reverse(bin)

RETURN
END


pro lris_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin = bin
;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    ;; Images are packed differently when Red side upgraded
    if strmid(sxpar(hdr,'DATE'),0,10) GT '2009-07-01' then begin
        image = x_readmhdufits(filelist[0], header=hdr, /notrim, $
                               /nobias)
        ;; Flip in y for blue side (back to the old way)
        if strcmp(strcompress(sxpar(hdr, 'INSTRUME'), /rem),'LRISBLUE') then $
          image = rotate(image,7) else begin
            ;; Trim top and bottom
            PRELINE = sxpar(hdr, 'PRELINE')
            POSTLINE = sxpar(hdr, 'POSTLINE')
            szi = size(image, /dimen)
            image = image[*,preline:szi[1]-preline-postline-1]
        endelse
    endif else $
      image = xmrdfits(filelist[0], 0, hdr, /fscale) 
endif

ncol = (size(image))[1]
nrow = (size(image))[2]
bin = long(strsplit(sxpar(hdr, 'BINNING'), ',', /extract))
xbin = bin[0]
ybin = bin[1]
postpix = long(sxpar(hdr, 'POSTPIX'))
if (keyword_set(verbose)) then $
  splog, filelist[0], ': Binning is ', xbin, ' x', ybin
; are we red or blue?
instrument = strcompress(sxpar(hdr, 'INSTRUME'), /rem)
IF strcmp(instrument, 'LRISBLUE') THEN BEGIN
;   prepix is in different headers for red versus blue (old red chip)
    prepix = long(sxpar(hdr, 'PRECOL'))
    ;; Sometimes only 2 amps for Longslit data
    namps = sxpar(hdr,'NUMAMPS')
    ;; JXP -- Likely bug in headers for TWOAMPTOP
    if namps EQ 2 then begin
        prepix = 100/xbin
        postpix = 160/xbin
    endif
    gain   = ([1.55,1.54471,1.63,1.65991])[4-namps:*]
    rnoise = ([3.9, 4.2, 3.6, 3.6])[4-namps:*]
ENDIF ELSE IF strcmp(instrument, 'LRIS') THEN BEGIN
    if strmid(sxpar(hdr,'DATE'),0,10) GT '2009-07-01' then begin
       if sxpar(hdr, 'NUMAMPS') NE 4 then stop
       namps = 4
       if strmid(sxpar(hdr,'DATE'),0,10) GT '2010-12-03' then begin
          ;; Second upgrade
          gain = [1.25, 1.18, 1.19, 1.17]
          rnoise = [4.65, 4.7, 4.5, 4.65]
          print, 'long_oscan: LRISr re-upgrade :: Using gain = ', gain
       endif else begin ;; LRIRs upgrade (the first)
          gain = [0.955, 1.022, 0.916, 0.877]
          rnoise = [3.97, 4.57, 3.38, 3.39]
       endelse
       prepix = long(sxpar(hdr, 'PRECOL'))
    endif else begin
        namps = 2
        gain = [1.98, 2.17]
        rnoise = [6.1, 6.3]
        prepix = long(sxpar(hdr, 'PREPIX'))
    endelse
ENDIF ELSE message, 'Cant figure out what instrument this is'

; These are small buffers to avoid using the overscan region to close to data
post_buffer1 = 4L
post_buffer2 = 8L
; The factor of namps is because of the two different amps

;;Begin - MF Jan 2011: Original line
; imagecol = ncol-namps*(prepix+postpix)/xbin
; Updated line: 
if (strcmp(instrument, 'LRIS')) and (strmid(sxpar(hdr,'DATE'),0,10) GT '2010-12-10') then $
  imagecol = ncol-namps*(prepix+postpix-1)/xbin $
else imagecol = ncol-namps*(prepix+postpix)/xbin
;Note that this might be a problem even before 2010-12-10 for binned data.
;;End - MF Jan 2011
imagerow = nrow

rawsub = fltarr(imagecol, imagerow)
rawivar = fltarr(imagecol, imagerow)

datacol =  namps* (prepix/xbin) + lindgen(namps)*1024/xbin
postcol =  datacol[namps-1L] + (1024+post_buffer1)/xbin
FOR iamp = 0, namps-1L DO BEGIN
    biascols = lindgen((postpix-post_buffer2)/xbin) $
      + (iamp*postpix)/xbin + postcol
    oscan = djs_avsigclip(image[biascols, *], 1)
    osub = replicate(1, imagecol/namps) # oscan
    imagecols = lindgen(1024/xbin)+ iamp*1024/xbin
    temp_image = (image[imagecols + namps*(prepix/xbin), *] - osub)*gain[iamp] 
    rawsub[imagecols, *] = temp_image
    rawivar[imagecols, *] = 1.0/(abs(temp_image - $
                                     sqrt(2.0)*rnoise[iamp]) +rnoise[iamp]^2)
ENDFOR

IF strcmp(instrument, 'LRIS') and $
  strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
    ;; bad pixel mask
    rawsub = transpose(rawsub)
    rawivar = transpose(rawivar)
    bin = reverse(bin)
    IF nrow EQ 2048 THEN BEGIN
        rawivar[0:722, 822] = 0.0
        rawivar[0:1482, 956] = 0.0
        rawivar[0:1485, 960] = 0.0
        rawivar[*, 2047] = 0.0
    ENDIF ELSE IF nrow EQ 1000 THEN BEGIN
        rawivar[0:224, 822] = 0.0
        rawivar[0:444, 956] = 0.0
    ENDIF ELSE message $
      , 'Problem in lris_oscan readout geometry or binning not supported'
ENDIF

RETURN
END

;;;;;;;;;;;;;;;;;;;;;;;
pro kast_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin


;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
endif

;mjd_newcontrol = 2454534L ;; New controller
mjd_newcontrol = 2454892L ;; New controlor (March, 2009)

;; Blue or Red?
if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' OR $
  strmid(sxpar(hdr,'VERSION'),0,5) EQ 'kastb' then side = 1 else side = 2
;if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' then side = 1 else side = 2

dims = size(image, /dim)
ncol = dims[0]
nrow = dims[1]

dateobs = sxpar(hdr, 'DATE-OBS')
mjd = x_setjdate(strmid(dateobs,0,10))

if mjd LT mjd_newcontrol then begin
    ;; Old controllers 
    bin= [round(2700./ncol),1]  ;; Might want to calculate bin[1] too
    
    if side EQ 1 then begin  ;; Blue side
        case strtrim(sxpar(hdr,'CLKSPEED'),2) of
            'Slow': begin
                dobs = sxpar(hdr,'DATE-OBS')
                jdate = x_setjdate(strmid(dobs,0,10))
                if jdate GE 2454526L then begin 
                   print, 'long_oscan:  Assuming CCD#9 for the Kast blue CCD'
                   if jdate GE 2454721L then begin
                      gain = 1.2 ;; Dewar#9 on Sep 11, 2008 by E. Gates
                      rnoise = 3.7
                   endif else begin
                      gain = 2.3 ;; Dewar#9 on Feb 29, 2008 as measured by E. Gates on Mar 1, 2008
                      rnoise = 8
                   endelse
                endif else begin
                   gain = 3.8 ;; Original Kast CCD
                   rnoise = 6.
                endelse
             end
            else: stop
        endcase
    endif else begin  ;; Red side
        case strtrim(sxpar(hdr,'CLKSPEED'),2) of
            'Slow': begin
                gain = 3.0 
                rnoise = 12.5
            end
            else: stop
        endcase
    endelse 
    rawsub = image * gain
    rawivar = 1. / ( (image>3.)*gain + rnoise^2)
    ;; Trim
    if side EQ 1 then begin
        rawsub = rawsub[*,10:*]
        rawivar = rawivar[*,10:*]
    endif else begin
        rawsub = rawsub[*,0:144]
        rawivar = rawivar[*,0:144]
    endelse
endif else begin
  ;; New controller (installed -- )
    bin = [sxpar(hdr,'CBIN'), sxpar(hdr,'RBIN')]
    case round(sxpar(hdr,'READ-SPD')) of 
        80: begin ;; Slow
            if side EQ 1 then begin
                gain = 1.2
                rnoise = 4.
            endif else begin
                gain = 3.
                rnoise = 12.5
            endelse
        end
        else: stop
    endcase

    if side EQ 2 then begin ;; Red side (easier)
        ;; Subtract overscan
        ;RdS 2011
        col_sci = sxpar(hdr, 'DNAXIS1') > (sxpar(hdr, 'NAXIS1')-sxpar(hdr, 'COVER'))
        ;col_sci = sxpar(hdr, 'DNAXIS1')
        col_origin = sxpar(hdr, 'CRVAL1')
        ncol = sxpar(hdr, 'NAXIS1')
        col_oscan = col_sci-col_origin  ;; This might not be right
        
        ;; Get the overscan
        oscan = djs_median(image[col_oscan:*, *], 1)
        temp_image = (image[0:col_oscan-1L, *] - $
                  (replicate(1.,ncol)#oscan))*gain 
        rawsub = temp_image
        rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)
    endif else begin  ;; BLUE SIDE
        ncol_over = sxpar(hdr, 'COVER')
        namps = sxpar(hdr, 'AMPSROW')
        ncol = sxpar(hdr, 'NAXIS1')
        if namps NE 2 then stop

        ;RdS 2011
        col_sci = sxpar(hdr, 'DNAXIS1') > (sxpar(hdr, 'NAXIS1')-2*sxpar(hdr, 'COVER'))
        ;col_sci = sxpar(hdr, 'DNAXIS1')
        sci1 = [0L, col_sci/2 - 1]
        sci2 = [sci1[1]+1,col_sci-1]
        bias1 = [col_sci+2,col_sci+ncol_over-2]
        bias2 = [col_sci+ncol_over+2,ncol-1]
        oscan1 = djs_median(image[bias1[0]:bias1[1], *], 1)
        oscan2 = djs_median(image[bias2[0]:bias2[1], *], 1)
        ;; Save
        temp_image = fltarr(col_sci,nrow)
        temp_image[sci1[0]:sci1[1],*] = (image[sci1[0]:sci1[1], *] - $
                  (replicate(1.,sci1[1]-sci1[0]+1)#oscan1))*gain 
        temp_image[sci2[0]:sci2[1],*] = (image[sci2[0]:sci2[1], *] - $
                  (replicate(1.,sci2[1]-sci2[0]+1)#oscan2))*gain 
        rawsub = temp_image
        rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)
    endelse
endelse
    
;; Spectra running up and down
rawsub = transpose(rawsub)
rawivar = transpose(rawivar)
    
RETURN

END

pro fire_oscan, filename, rawsub, rawivar, hdr = hdr $
               , gain = gain, rnoise = rnoise, image = image $
               , verbose = verbose, bin = bin

  tmp = mrdfits(filename, 0, hh)

  hdr = hh
  gain = 1.2
  rnoise = 17.0/sqrt(4)
  bin = 1

  rawsub = (tmp)
  rawivar = 1.0/(abs(tmp - sqrt(2.0)*rnoise) +rnoise^2)
  rawivar = (rawivar)


END

pro fors2_oscan, filename, rawsub, rawivar, hdr = hdr $
                 , gain = gain, rnoise = rnoise, image = image $
                 , verbose = verbose, bin = bin

;  no overscan for rows...
if NOT keyword_set(image) then begin
;  look for filename
    filelist = lookforgzip(filename)
    if filelist[0] EQ '' then begin
        print, 'Could not find file named ', filename
        return
    endif
    image = xmrdfits(filelist[0], 0, hdr, /fscale) 
 endif

dims = size(image, /dim)
nx = dims[0]
ny = dims[1]
ind_binx = WHERE(stregex(hdr, 'HIERARCH ESO DET WIN1 BINX', /bool))
binx = double(strmid(hdr[ind_binx], 30, 14))
ind_biny = WHERE(stregex(hdr, 'HIERARCH ESO DET WIN1 BINY', /bool))
biny = double(strmid(hdr[ind_binx], 30, 14))
bin = [binx, biny]

ind_gain = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 GAIN', /bool))
gain = (double(strmid(hdr[ind_gain], 30, 14)))[0]
ind_ron = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 RON', /bool))
rnoise = (double(strmid(hdr[ind_ron], 30, 14)))[0]
ind_det =  WHERE(stregex(hdr, 'HIERARCH ESO DET CHIP1 NAME', /bool))
detector = strcompress(strmid(hdr[ind_det], 30, 14))

ind_prscx = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 PRSCX', /bool))
prscx = (long(strmid(hdr[ind_prscx], 30, 14)))[0]
ind_ovscx = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 OVSCX', /bool))
ovscx = (long(strmid(hdr[ind_ovscx], 30, 14)))[0]

ind_prscy = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 PRSCY', /bool))
prscy = (long(strmid(hdr[ind_prscy], 30, 14)))[0]
ind_ovscy = WHERE(stregex(hdr, 'HIERARCH ESO DET OUT1 OVSCY', /bool))
ovscy = (long(strmid(hdr[ind_ovscy], 30, 14)))[0]

data_arr = [prscx+1L, prscy+1L, nx-ovscx, ny-ovscy]
bias_arr = [prscx+1L, ny-ovscy+1L, nx-ovscx, ny]
;datasec = strcompress(sxpar(hdr, 'DATASEC'), /rem)
;biassec = strcompress(sxpar(hdr, 'BIASSEC'), /rem)
;data_arr = long(strsplit(datasec, '[*:*,*:*]', /extract))
;bias_arr = long(strsplit(biassec, '[*:*,*:*]', /extract))
nbias = bias_arr[3]-bias_arr[1] + 1L

; These are small buffers to avoid using the overscan region to close to data
oscan_buffer = 1L
imagex = data_arr[2]-data_arr[0] + 1L
imagey = data_arr[3]-data_arr[1] + 1L

rawsub = fltarr(imagex, imagey)
rawivar = fltarr(imagex, imagey)

biascols = lindgen(nbias-oscan_buffer) + imagey + oscan_buffer + prscy[0]
oscan1 = djs_avsigclip(image[0:imagex-1L, biascols], 0)
sig_res = 7
nhalf =  long(sig_res)*4L
xkern = dindgen(2*nhalf+1)-nhalf
kernel = gauss1(xkern, [0.0, sig_res, 1.0])
oscan = convol(oscan1, kernel, /edge_truncate)
osub = oscan # replicate (1, imagey)
temp_image = (image[prscx:prscx + imagex-1L $
                    , prscy: prscy + imagey -1L] - osub)*gain 
rawsub = temp_image
rawivar = 1.0/(abs(temp_image - sqrt(2.0)*rnoise) +rnoise^2)

;; Now transpose so that spectral direction vertical with blue down
rawsub = transpose(rawsub)
rawivar = transpose(rawivar) 
bin = reverse(bin)

RETURN
END



;;;;;;;;;;;;;;;;;;;;;;;
pro long_oscan, filename, rawsub, rawivar, hdr = hdr $
                , gain = gain, rnoise = rnoise, image = image $
                , verbose = verbose, bin=bin, CCDONLY = CCDONLY $
                , TRANSFORM = TRANSFORM


hdr1 = xheadfits(filename)

if (size(hdr1, /tname) NE 'STRING') then begin
    splog, 'Invalid FITS header for file ', filename
    flux = 0
    invvar = 0
    return
endif

;  Is this Keck or Gemini?
telescope = strcompress(sxpar(hdr1, 'TELESCOP'), /rem)
instrument = strcompress(sxpar(hdr1, 'INSTRUME'), /rem)
detector = strcompress(sxpar(hdr1, 'DETECTOR'), /rem)
telid =  strcompress(sxpar(hdr1, 'TELID'), /rem)

IF strcmp(telescope, 'KeckI') THEN $
  lris_oscan, filename, rawsub, rawivar, hdr = hdr $
  , gain = gain, rnoise = rnoise, image = image $
  , verbose = verbose, bin = bin $
ELSE IF (strcmp(telescope, 'Gemini-North') OR $
         strcmp(instrument, 'GMOS-N')) THEN BEGIN
   IF KEYWORD_SET(TRANSFORM) THEN $
      gmos_trnimg3to1, filename, rawsub, rawivar, hdr = hdr $
                    , gain = gain, rnoise = rnoise $
                    , verbose = verbose, bin = bin $
   ELSE gmos_oscan, filename, rawsub, rawivar, hdr = hdr, thdr = thdr $
                    , gain = gain, rnoise = rnoise, VERBOSE = VERBOSE $
                    , bin = bin, CCDONLY = CCDONLY
ENDIF ELSE IF strcmp(telescope, 'mmt') THEN $
   mmt_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
              , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF  strcmp(detector, 'mars') THEN $
   mars_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF  strcmp(telescope, 'kp4m') THEN $
   kpno_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, 'DEIMOS*') THEN $
  deimos_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF (stregex(instrument, 'KAST*',/bool,/fold_case) EQ 1) OR $
         (stregex(sxpar(hdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then $
  kast_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, 'DIS*') THEN  $
  dis_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(telid, '200') THEN $
  p200_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin $ 
ELSE IF strmatch(instrument, '*IMACS*') THEN $ 
  imacs_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin  $
ELSE IF strmatch(instrument, 'FIRE') THEN $ 
  fire_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin  $
ELSE IF  strcmp(telescope, 'CA-3.5') THEN $
   caha_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
               , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF  strcmp(telescope, 'CA-2.2') THEN $
   caha_22m_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
                   , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE IF strmatch(instrument, '*bcspeclamps*') THEN $ ; jm11jun08ucsd
  bok_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
  , rnoise = rnoise, image = image, verbose = verbose, bin = bin  $
ELSE IF strmatch(instrument, '*FORS2*') THEN $
   fors2_oscan, filename, rawsub, rawivar, hdr = hdr, gain = gain $
                , rnoise = rnoise, image = image, verbose = verbose, bin = bin $
ELSE message, 'Not sure what instrument you want here'

RETURN
END
