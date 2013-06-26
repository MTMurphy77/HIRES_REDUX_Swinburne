path = '/Users/joe/DATA/SOFI_DATA/2011-09-22/'
init = 152
nseq = 2
even = 2*lindgen(nseq)
odd = 2*lindgen(nseq) + 1
anum = init + even
bnum = init + odd
afiles = path + 'SOFI_' + string(anum, FORMAT = '(I4.4)') + '.fits'
bfiles = path + 'SOFI_' + string(bnum, FORMAT = '(I4.4)') + '.fits'

FOR ii = 0L, nseq-1L DO BEGIN
   IF file_test(afiles[ii]) THEN sofi_proc, afiles[ii], aimg ELSE stop
   IF file_test(bfiles[ii]) THEN sofi_proc, bfiles[ii], bimg ELSE stop
   IF ii EQ 0 THEN BEGIN
      dim = size(aimg, /dim)
      a_stack = fltarr(dim[0], dim[1], nseq)
      b_stack = fltarr(dim[0], dim[1], nseq)
   ENDIF
   a_stack[*, *, ii] = aimg
   b_stack[*, *, ii] = bimg
   diff = aimg-bimg
   djs_iterstat, diff, median = median, mean = mean, sigma = sigma, sigrej = 3.0
   atv, aimg - bimg, min = -10.0*sigma, max = 10.0*sigma, /block
ENDFOR

IF nseq GT 1 THEN BEGIN
   a_mas =  djs_avsigclip(a_stack, 3, sigrej = 3.0, outmask = amask)
   b_mas =  djs_avsigclip(b_stack, 3, sigrej = 3.0, outmask = bmask)
   aden = total((amask EQ 0), 3) 
   bden = total((bmask EQ 0), 3) 
   a_bar = (aden GT 0.0)*total(a_stack*(amask EQ 0), 3)/(aden + (aden EQ 0.0))
   b_bar = (bden GT 0.0)*total(b_stack*(bmask EQ 0), 3)/(bden + (bden EQ 0.0))
ENDIF ELSE BEGIN
   a_bar = a_stack
   b_bar = b_stack
ENDELSE
diff = a_bar - b_bar
djs_iterstat, diff, median = median, mean = mean, sigma = sigma, sigrej = 3.0
xatv, diff, min = -10.0*sigma, max = 10.0*sigma





END
