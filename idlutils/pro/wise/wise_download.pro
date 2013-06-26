;+
; NAME:
;   wise_download
;
; PURPOSE:
;   Download commands for only WISE data (all bands) within an RA,Dec box
;
; CALLING SEQUENCE:
;   wise_download, rarange=, decrange=
;
; INPUTS:
;   rarange    - RA range for field centers; default to [209.65,220.08] deg
;   decrange   - DEC range for field centers; default to [50.63,54.74] deg
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; REVISION HISTORY:
;   28-Apr-2013  Written by D. Schlegel, LBL
;-
;------------------------------------------------------------------------------
pro wise_download, rarange=rarange1, decrange=decrange1

   if (keyword_set(rarange1)) then rarange = rarange1 $
    else rarange = [209.65,220.08]
   if (keyword_set(decrange1)) then decrange = decrange1 $
    else decrange = [50.63,54.74]

   splog, file='wise_download.sh', /noname

; Directories with file lists
setenv,'WISE_IMAGE_DIR=/Users/schlegel/wise1;/Users/schlegel/wise1ext'
wdir = str_sep(getenv('WISE_IMAGE_DIR'), ';')

; Directories at IPAC to get frames
remotedir = ['http://irsa.ipac.caltech.edu/ibe/data/wise/allsky/4band_p1bm_frm/', $
 'http://irsa.ipac.caltech.edu/ibe/data/wise/prelim_postcryo/p1bm_frm/']

   for idir=0L, n_elements(wdir)-1 do begin
      ixfile = filepath('WISE-index-L1b.fits', root_dir=wdir[idir])
      ixlist = mrdfits(ixfile, 1, /silent)
      inear = where(ixlist.ra GE rarange[0] AND ixlist.ra LE rarange[1] $
       AND ixlist.dec GE decrange[0] AND ixlist.dec LE decrange[1], nnear)
      for i=0L, nnear-1 do begin
         subdirs = [ixlist[inear[i]].scangrp, $
          ixlist[inear[i]].scan_id, $
          string(ixlist[inear[i]].frame_num,format='(i3.3)')]
         subdir3 = subdirs[0]+'/'+subdirs[1]+'/'+subdirs[2]+'/'

         thiscmd = 'wget -r -N -nH -np -nv --cut-dirs=3' $
          + ' "'+remotedir[idir]+subdir3+'"
         splog,thiscmd,/noname
      endfor
if (idir EQ 0) then $
 mwrfits, ixlist[inear], 'WISE1-WISE-index-L1b.fits', /create $
else $
 mwrfits, ixlist[inear], 'WISE1EXT-WISE-index-L1b.fits', /create
   endfor

   splog, /close

   return
end
;------------------------------------------------------------------------------
