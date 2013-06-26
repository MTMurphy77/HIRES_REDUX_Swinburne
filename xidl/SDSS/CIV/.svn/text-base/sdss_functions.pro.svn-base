;+ 
; NAME:
; sdss_functions.pro
;    Version 1.0
;
; PURPOSE:
;   Many functions that retrieve all little information for/from SDSS
;   CIV/SiIV/CaII stuff
;
; CALLING SEQUENCE:
;   
;   sdss_functions will print all the possible functions to screen
;   (and compile everything).
;   Otherwise, compile at beginning of codes with:
;   @sdss_functions or by calling sdss_functions,/compile
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
;   27-May-2011  created by KLC
;   30-Jul-2011  updated to handle new structure, KLC
;-
;------------------------------------------------------------------------------
;; Conventions:
;; *Every function/procedure must have n_params() check or /help
;; keyword. 
;; *Every function/procedure must have description in cade and in
;; sdss_functions print-out.
;; *Comments!  Including associating final end with
;; function/procedure. 
;; *Function/procedure names:  sdss_get[], sdss_prnt[], sdss_plt[],
;; sdss_cp[]2[]. 

function sdss_getsdssdir,help=help
  ;; Set the default data structure (with appropriate /'s)
  if keyword_set(help) then begin
     print,'Syntax - sdss_setsdssdir([/help])'
     return,-1
  endif 
     
  return,getenv('SDSSPATH')+'/'+getenv('SDSSDR')+'/'
end                             ; sdss_getsdssdir()


function sdss_getname,sdsstab,spec=spec,hdr=hdr,eig=eig,spl=spl,abslin=abslin, $
                      hyb=hyb, clean=clean, gz=gz, $
                      plate=plate, root=root, dir=dir
  ;; Return the spSpec-jjjjj-pppp-fff*.fits for given structure or
  ;; header 
  if n_params() ne 1 then begin
     print,'Syntax - sdss_getname(sdsstab,[/spec,/hdr,/eig,/spl,/abslin,plate='
     print,'                      /hyb,/cleanspec, root=, dir=])'
     return,-1
  endif 
  
  ;; Parse pieces of names
  nfil = (size(sdsstab,/dim))[0]
  if keyword_set(spec) then begin
     ;; Assumes all names are of equal length and nature; avoids looping
     prs = strsplit(sdsstab[0],'/',count=nprs)
     name = strmid(sdsstab,prs[nprs-1])
     prs = strsplit(name[0],'.',count=nprs)
     name = strmid(name,0,prs[1]-1) ; remove .fit*
     prs = strsplit(name[0],'-',count=nprs)
     mjd = strmid(name,prs[1],prs[2]-prs[1]-1)
     plate = strmid(name,prs[2],prs[3]-prs[2]-1)
     fiber = strmid(name,prs[3])
  endif else begin

     if keyword_set(hdr) then begin
        ;; sdsstab is actually a header; and only 1
        nm = sxpar(sdsstab,'NAME')
        fiber = strtrim(sxpar(sdsstab,'FIBERID'),2)
        plate = strtrim(sxpar(sdsstab,'PLATEID'),2)
        prs = strsplit(nm,'-',/extract,count=nprs)
        mjd = prs[1]
     endif else begin
        ;; Use structure information
        if size(sdsstab,/type) ne 8 then $
           stop,'sdss_getname(): input must be structure'
        mjd = strtrim(sdsstab.smjd,2)
        fiber = strtrim(sdsstab.fiber,2)
        plate = strtrim(sdsstab.plate,2)
     endelse 
     
     ;; Set appropriate lengths
     bd = where(strlen(fiber) lt 3,nbd)
     while nbd ne 0 do begin
        fiber[bd] = '0'+fiber[bd]
        bd = where(strlen(fiber) lt 3,nbd)
     endwhile
     bd = where(strlen(plate) lt 4,nbd)
     while nbd ne 0 do begin
        plate[bd] = '0'+plate[bd]
        bd = where(strlen(plate) lt 4,nbd)
     endwhile
  endelse                       ; header or structure
  
  ;; Concatenate
  root = mjd + '-' + plate + '-' + fiber
  if not keyword_set(name) then name = 'spSpec-' + root

  ;; dir is relative to $SDSSPATH/$SDSSDR
  dir = 'spectro/'              ; may be overwritten

  ;; Suffix
  if keyword_set(eig) and keyword_set(spl) and keyword_set(abslin) and keyword_set(hyb) then $
     stop,'sdss_getname(): cannot have all keywords (/eig, /spl, /abslin, /hyb) set'
  if keyword_set(eig) then begin
     dir = 'conti/'
     name = name+'-eigconti'
  endif 
  if keyword_set(spl) then begin
     dir = 'conti/'
     name = name+'-splconti'
  endif 
  if keyword_set(abslin) then begin
     dir = 'abslin/'
     name = name+'-abslin'
  endif 
  if keyword_set(hyb) then begin
     dir = 'conti/'
     name = name+'-hybconti'
  endif 
  if keyword_set(clean) then begin
     dir = 'cleanspec/'
     name = name+'-clean'
  endif 
  if keyword_set(gz) then begin
     dir = 'gz/'
     name = name+'-gz'
  endif 

  ;; All have same subdirectory
  dir = dir + '1d_26/'+plate+'/1d/' 

  ;; Default
  return,name+'.fit'
  
end                             ; sdss_getname()


function sdss_getqsostrct,dr=dr, name=name, bal=bal, nobal=nobal, help=help
  ;; Return the appropriate data release catalog structure or file
  ;; name to be used 
  if keyword_set(help) then begin 
     print,'Syntax - sdss_getqsostrct([dr=, /name, /bal, /nobal, /help])'
     return, -1 
  endif 

  ;; Set Parameters
  if not keyword_set(dr) then dr = 7 ;data release
  
  case dr of 
     7: if keyword_set(bal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_BAL.fit' $
     else if keyword_set(nobal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_noBAL.fit' $
     else fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.fit'
     else: begin
        print,'sdss_getqsostrct(): DR='+strtrim(dr,2)+' not supported'
        return, -1
     end
  endcase
        
  ;; Check existence
  test = file_search(fil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_getqsostrct(): file DNE ',strtrim(fil,2)+'*'

  ;; Return either name or read
  if keyword_set(name) then return,fil $
  else begin
     strct = xmrdfits(fil,1,/silent)
     return, strct
  endelse 
  
end                             ; sdss_getqsostrct()


function sdss_getqsolist,dr=dr, read=read, help=help, count=count, bal=bal, $
                         nobal=nobal, abslin_dir=abslin_dir 
  ;; Return the appropriate data release catalog spectra list to be used 
  if keyword_set(help) then begin 
     print,'Syntax - sdss_getqsolist([dr=, /read, /help, count=, /bal, '
     print,'                         /nobal, abslin_dir=])'
     return, -1 
  endif 

  ;; Set Parameters
  if not keyword_set(dr) then dr = 7 ;data release
  
  case dr of 
     7: if keyword_set(bal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_BAL.list' $
     else if keyword_set(nobal) then $
        fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_noBAL.list' $
     else  fil = getenv('XIDL_DIR')+'/SDSS/CIV/dr7qso_srt.list'
     else: begin
        print,'sdss_getqsostrct(): DR='+strtrim(dr,2)+' not supported'
        return, -1
     end
  endcase
        
  ;; Check existence
  test = file_search(fil,count=ntest)
  if ntest eq 0 then $
     stop,'sdss_getqsostrct(): file DNE ',strtrim(fil,2)+'*'

  ;; Return either name or read
  if keyword_set(read) then begin
     readcol,fil,nam,format='a',/silent
     count = (size(nam,/dim))[0]
     abslin_dir = nam[0]
     count = count-1
     return,nam[1:count]
  endif else begin
     return, fil
  endelse 
  
end                             ; sdss_getqsolist()


pro sdss_wrqsolist,sdsssum,list_fil,absdir=absdir,silent=silent
  ;; Write the specifically formatted files to work with
  ;; e.g. sdss_fndlin
  if n_params() ne 2 then begin
     print,'Syntax - sdss_wrqsolist, sdsssum, list_fil, [absdir=,/silent]'
     return
  endif
  if not keyword_set(absdir) then absdir = 'abslin/' 
  if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)

  spec_fil = sdss_getname(sdsstab,plate=plate,dir=dir)

  openw,1,list_fil
  printf,1,absdir
  writecol,list_fil,dir+spec_fil,filnum=1
  close,1
  if not keyword_set(silent) then $
     print,'sdss_wrqsolist: created ',list_fil

end                             ; sdss_wrqsolist()


pro sdss_getqsoinlist,list_fil,sdsssum,silent=silent,snrstrct=snrstrct
  ;; Read list and find the matching QSOs
  if n_params() ne 2 then begin
     print,'Syntax - sdss_getqsoinlist'
     return
  endif 

  readcol,list_fil,spec,format='a',skip=1,/silent
  nspec = (size(spec,/dim))[0]
  if keyword_set(snrstrct) then begin
     if size(snrstrct,/type) eq 7 then sdsstab = xmrdfits(snrstrct,1,/silent) $
     else sdsstab = snrstrct
     qso_name = sdsstab.qso_name
  endif else begin
     sdsstab = sdss_getqsostrct()
     ;; Make name
     plate = strtrim(sdsstab.plate,2)
     bd = where(strlen(plate) ne 4)
     while bd[0] ne -1 do begin
        plate[bd] = '0'+plate[bd]
        bd = where(strlen(plate) ne 4)
     endwhile
     fiber = strtrim(sdsstab.fiber,2)
     bd = where(strlen(fiber) ne 3)
     while bd[0] ne -1 do begin
        fiber[bd] = '0'+fiber[bd]
        bd = where(strlen(fiber) ne 3)
     endwhile
     qso_name = strtrim(sdsstab.smjd,2)+'-'+plate+'-'+fiber
  endelse 

  ;; Loop
  istrt = (strsplit(spec[0],'spSpec-',count=nprs))[1]
  istop = strlen(spec[0]) - istrt - 4
  subqso_name = strmid(spec,istrt,istop)
  for ii=0L,nspec-1 do begin
     mtch = where(subqso_name[ii] eq qso_name,nmtch)
     if nmtch ne 1 then $
        stop,'sdss_getqsoinlist: multiple matches!'
     if ii eq 0 then subsdsstab = sdsstab[mtch[0]] $
     else subsdsstab = [subsdsstab,sdsstab[mtch[0]]]
  endfor                        ; loop ii=nspec

  if size(sdsssum,/type) eq 7 then begin
     mwrfits,subsdsstab,sdsssum,/create,/silent
     spawn,'gzip -f '+sdsssum
     if not keyword_set(silent) then $
        print,'sdss_getqsoinlist: created ',sdsssum
  endif else sdsssum = subsdsstab

end                             ; sdss_getqsoinlist


function sdss_getcflg, eig=eig, spl=spl, hyb=hyb, custom=custom, help=help, $
                       index=index
  ;; Return convention for whether absorption system defined with
  ;; eigen continuum or spline continuum.
  ;; eigen will be default
  ;; To map to index, take fix(alog(cflg)/alog(2))
  if keyword_set(help) then begin
     print,'Syntax - sdss_getcflg([/eig,/spl,/custom,/help]); default /eig'
     return,-1
  endif 
  cflg = 1                             ; /eig and default; index = 0
  if keyword_set(custom) then cflg = 8 ; index = 3
  if keyword_set(hyb) then cflg = 4    ; index = 2
  if keyword_set(spl) then cflg = 2    ; index = 1

  if keyword_set(index) then return,fix(alog(cflg)/alog(2)) $
  else return, cflg
end                             ; sdss_getcflg()


function sdss_getrating, definite=definite, good=good, maybe=maybe, bad=bad, $
                         unrated=unrated, help=help
  ;; Return convention sdss_chkciv evaluation
  if keyword_set(help) then begin
     print,'Syntax - sdss_getrating([/definite,/good,/maybe,/bad,/unrated])'
     print,'                  Default: /bad'
     return,-1
  endif 
  flg = 0
  
  if keyword_set(definite) and (keyword_set(good) or keyword_set(maybe) or $
                                keyword_set(unrated)) then begin
     print,'sdss_getrating(): multiple keywords set; defaulting'
     return,-1
  endif

  if not (keyword_set(definite) or keyword_set(good) or keyword_set(maybe) or $
          keyword_set(unrated)) then return, 0
  
  if keyword_set(definite) then return, 3
  if keyword_set(good) then return, 2
  if keyword_set(maybe) then return, 1
  if keyword_set(unrated) then return, -1
  
end                             ; sdss_getrating()


function sdss_srtcivstrct, civstrct_fil, by_rating=by_rating, index=index
  ;; Sort entries by QSO name (jjjjj-pppp-fff) and zabs
  if n_params() ne 1 then begin
     print,'Syntax - sdss_srtcivstrct( civstrct_fil, [/by_rating, index=])'
     return,-1
  endif 
  
  if size(civstrct_fil,/type) eq 8 then civstrct = civstrct_fil $
  else civstrct = xmrdfits(civstrct_fil,1,/silent)

  if keyword_set(by_rating) then begin
     index = civstrct.rating[0] * 0 - 1 ; quick instantiation
     rtg = civstrct[uniq(civstrct.rating[0],sort(civstrct.rating[0]))].rating[0]
     nrtg = (size(rtg,/dim))[0] > 1 ; foils singularity
     for rr=0,nrtg-1 do begin
        sub = where(civstrct.rating[0] eq rtg[rr])
        tmp = sdss_srtcivstrct(civstrct[sub],index=isub)
        if rr eq 0 then newcivstr = tmp $ 
        else newcivstr = [newcivstr,tmp]
        index[sub] = sub[isub]
     endfor                     ; loop rr=nrtg
     civstrct = newcivstr
  endif else begin
     civstrct.qso_name = strtrim(civstrct.qso_name,2)
     index = sort(civstrct.qso_name)
     civstrct = civstrct[index]
     unq = uniq(civstrct.qso_name)
     nunq = (size(unq,/dim))[0] 
     
     rng = lindgen(unq[0]+1)
     srt = sort(civstrct[rng].zabs_orig[0])
     civstrct[rng] = civstrct[rng[srt]] ; sort by QSO and z
     index[rng] = index[rng[srt]]
     
     for qq=1,nunq-1 do begin
        rng = unq[qq-1] + 1 + lindgen(unq[qq]-unq[qq-1])
        srt = sort(civstrct[rng].zabs_orig[0])
        civstrct[rng] = civstrct[rng[srt]] ; sort by QSO and z
        index[rng] = index[rng[srt]]
     endfor                                ; loop qq=nunq
  endelse

  return, civstrct
end                             ; sdss_srtcivstrct()


function sdss_rmcivduplicates, civstrct_fil, count=count, index=index, $
                               complement=complement, icomplement=icomplement, $
                               ncomplement=ncomplement, check=check
  ;; Remove duplicate candidates based on qso_name, zabs_orig, and
  ;; cflg as in sdss_chkciv
  if n_params() ne 1 then begin
     print,'Syntax - sdss_rmcivduplicates( civstrct_fil, [count=, index=,'
     print,'                               complement=, ncomplement=, '
     print,'                               icomplement=, check=])'
  endif 

  if size(civstrct_fil,/type) eq 8 then civstrct = civstrct_fil $
  else civstrct = xmrdfits(civstrct_fil,1,/silent)

  civstrct = sdss_srtcivstrct(civstrct)
  nciv = (size(civstrct,/dim))[0] > 1 ; defeats singularity prob
  mask = intarr(nciv)                 ; 0 = keep; 1 = duplicate
  check = mask - 1L                    ; to hold indices of matches

  unq = uniq(civstrct.qso_name) ; defined by the last per block
  nqso = (size(unq,/dim))[0]
  nciv_per_qso = [unq[0] + 1,(shift(unq,-1)-unq)[0:nqso-2]]
  
  for qq=0L,nqso-1 do begin
     if qq gt 0 then rng = unq[qq-1] + lindgen(nciv_per_qso[qq]) $
     else rng = lindgen(nciv_per_qso[qq])
     ;; printcol,rng,'  '+civstrct[rng].qso_name,civstrct[rng].zabs_orig[0],civstrct[rng].cflg 

     ;; redshifts are already sorted
     zunq = uniq(civstrct[rng].zabs_orig[0])
     nzunq = (size(zunq,/dim))[0] > 1 ; defeats singularity problem
     
     if nzunq ne nciv_per_qso[qq] then begin
        ;; Then there's duplicate to trace down so want to
        ;; modify rng and instantiate complement
        nciv_per_z = [zunq[0] + 1,(shift(zunq,-1)-zunq)[0:nzunq-2]]
        zunq = rng[zunq]
        bd = where(nciv_per_z ne 1,nbd)
        if nbd ne 0 then begin
           bd2 = where(nciv_per_z[bd] eq 2,nbd2,$ ; easier to handle
                      complement=bd3,ncomplement=nbd3)
           if nbd2 ne 0 then begin
              ;; Since just 2 then can do quick global where() grep
              ;; without looping; uniq() returns the *last* incidence
              ;; of a sorted list
              ;; qso_name and zabs_orig[0] already matched
              bd2 = zunq[bd[bd2]]
              subbd2 = where(civstrct[bd2].cflg eq $
                            civstrct[bd2-1].cflg,nsubbd2,$
                           complement=subgd,ncomplement=nsubgd) 
              if nsubbd2 ne 0 then begin
                 mask[bd2[subbd2]]++ ; duplicate

                 ;; Make sure to take most favorable rating for one
                 ;; being kept
;                 stop,'sdss_rmcivduplicates(): check this subbd2 is right'
                 ;; printcol,bd2,'  '+civstrct[bd2].qso_name,civstrct[bd2].zabs_orig[0],civstrct[bd2].cflg        
                 ;; printcol,bd2-1,'  '+civstrct[bd2-1].qso_name,civstrct[bd2-1].zabs_orig[0],civstrct[bd2-1].cflg
                 ;; print,bd2[subbd2],'  '+civstrct[bd2[subbd2]].qso_name,civstrct[bd2[subbd2]].zabs_orig[0],civstrct[bd2[subbd2]].cflg
                 ;; print,bd2[subbd2]-1,'  '+civstrct[bd2[subbd2]-1].qso_name,civstrct[bd2[subbd2]-1].zabs_orig[0],civstrct[bd2[subbd2]-1].cflg


                 civstrct[bd2[subbd2]-1].rating[0] = $
                    civstrct[bd2[subbd2]-1].rating[0] > $
                    civstrct[bd2[subbd2]].rating[0]

                 check[bd2[subbd2]-1] = bd2[subbd2]
              endif             ; nsubbd2 ne 0
           endif                ; nbd2 ne 0

           if nbd3 ne 0 then begin
              ;; Have to loop (slow!)
              bd3 = zunq[bd[bd3]]
              stop,'sdss_rmcivduplicates(): do not know how to handle triplets+'
           endif                ; nbd3 ne 0
        endif                   ; nbd ne 0

     endif                      ; nzunq ne nciv_per_qso[qq]

  endfor                        ; loop qq=nqso

  index = where(mask eq 0,count,complement=icomplement,ncomplement=ncomplement)
  subcivstr = civstrct[index]
  if ncomplement eq 0 then complement = -1 $
  else complement = civstrct[icomplement]
  
  return, subcivstr
end                             ; sdss_rmcivduplicates()


function sdss_getrateddblt, strct_fil, count=count, _extra=extra
  ;; Return the subsample with the desired rating[0]
  if n_params() ne 1 then begin
     print,'Syntax - sdss_getrateddblt(strct_fil, [count=, _extra=])'
     return,-1
  endif 

  if size(strct_fil,/type) eq 8 then strct = strct_fil $
  else strct = xmrdfits(strct_fil,1,/silent)

  gd = where(strct.rating[0] eq sdss_getrating(_extra=extra),count)
  if count eq 0 then return, -1 $
  else return, strct[gd]

end                             ; sdss_getrateddblt()


function sdss_getlimflg, upper=upper, lower=lower, help=help
  ;; Report upper (4) or lower (2) limit flags
  if keyword_set(help) then begin
     print,'Syntax - sdss_getlimflg([/upper,/lower,/help])'
     return, -1
  endif 
  if keyword_set(upper) then return,4
  if keyword_set(lower) then return,2
  if not keyword_set(upper) and not keyword_set(lower) then return,1 ; good!
end                             ; sdss_getlimflg()


function sdss_setlimflg, inflg, upper=upper, lower=lower
  ;; Binary combine upper (4) and/or lower (2) limit flags to inflg
  ;; For example, a good (flg=1) line that is a lower limit (2) has
  ;; binary combination 3. 
  if n_params() lt 1 then begin
     print,'Syntax - sdss_setlimflg(inflg,[/upper,/lower])'
     return, -1
  endif 
  outflg = inflg
  if keyword_set(upper) then outflg = outflg or sdss_getlimflg(/upper)
  if keyword_set(lower) then outflg = outflg or sdss_getlimflg(/lower)

  return, outflg
end                             ; sdss_setlimflg()


function sdss_getewflg, boxcar=boxcar, custom=custom, $
                        orig=orig, help=help
  ;; Report EW flags
  ;; 1: orig; 8: tau; 16: boxcar; 32: gauss; 64: custom
  if keyword_set(help) then begin
     print,'Syntax - sdss_getewflg([/boxcar,/custom,/orig,/help])'
     return, -1
  endif 
  if keyword_set(orig) then return,1
  if keyword_set(boxcar) then return,16
  if keyword_set(custom) then return,64
end                             ; sdss_getewflg()


function sdss_setewflg, inflg, boxcar=boxcar, custom=custom, $
                        orig=orig, bad=bad
  ;; Binary combine EW flags to inflg
  ;; For example, a orig (flg=1) line that is a lower limit (2) has
  ;; binary combination 3. 
  if n_params() lt 1 then begin
     print,'Syntax - sdss_setewflg(inflg,[/boxcar,/custom,/orig])'
     return, -1
  endif 
  outflg = inflg
  if keyword_set(boxcar) then outflg = outflg or sdss_getewflg(/boxcar)
  if keyword_set(custom) then outflg = outflg or sdss_getewflg(/custom)
  if keyword_set(orig) then outflg = outflg or sdss_getewflg(/orig)
  if keyword_set(bad) then outflg = 0 ; don't do anything!

  return, outflg
end                             ; sdss_setewflg()


function sdss_getbalflg, civbal=civbal, mgiibal=mgiibal, eig=eig, $
                         help=help
  ;; Report Shen et al. (2011) catalog (1: CIV BAL, 2: MgII BAL, 3:
  ;; both) or eig-selected (16) BAL flags
  ;; Shen et al. (2011), ApJS, in press, arXiv:1006.5178
  ;; http://das.sdss.org/va/qso_properties_dr7/dr7.htm
  if keyword_set(help) then begin
     print,'Syntax - sdss_getbalflg(/civbal or /mgiibal or /eig,[/help])'
     return, -1
  endif 
  if keyword_set(civbal) then return,1
  if keyword_set(mgiibal) then return,2
  if keyword_set(eig) then return,16
  return,4                      ; unknown
end                             ; sdss_getbalflg()


function sdss_fndbalqso, sdsssum, count=count, balflg=balflg, full=full, $
                         flgonly=flgonly
  ;; Match BAL QSOs from Shen et al (2011) to input QSO structure by
  ;; returning the where() and optionally, the full list of flags
  if n_params() ne 1 then begin
     print,'Syntax - sdss_fndbalqso(sdsssum, [count=, balflg=, /full, /flgonly]'
     return, -1
  endif 
  sdssdir = sdss_getsdssdir()

  if keyword_set(full) then begin
     ;; Restore save file (ibal and balflg); faster!
     restore,getenv('XIDL_DIR')+'/SDSS/CIV/sdss_balqso_map.sav'
     count = (size(ibal,/dim))[0]
     if keyword_set(flgonly) then return, balflg $
     else return,ibal
  endif 

  ;; Otherwise must search
  if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)
  nsdss = (size(sdsstab,/dim))[0]
  balflg = replicate(0,nsdss)
  
  ;; Load BAL data
  baltab = sdss_getqsostrct(/bal)
  nbal = (size(baltab,/dim))[0]
  
  ;; Compare with minimal looping
  nfil = nbal < nsdss
  for qq=0L,nfil-1 do begin
     ;; Search two ways
     if nbal lt nsdss then begin
        ibal = qq
        mtch = where(baltab[qq].mjd eq sdsstab.smjd and $
                     baltab[qq].plate eq sdsstab.plate and $
                     baltab[qq].fiber eq sdsstab.fiber,nmtch)
        if nmtch gt 1 then $
           stop,'sdss_instantbalflg: multiple matches!'
        if nmtch ne 0 then mtch = mtch[0]
     endif else begin
        mtch = where(baltab.mjd eq sdsstab[qq].smjd and $
                     baltab.plate eq sdsstab[qq].plate and $
                     baltab.fiber eq sdsstab[qq].fiber,nmtch)
        if nmtch gt 1 then $
           stop,'sdss_instantbalflg: multiple matches!'
        if nmtch ne 0 then begin
           ibal = mtch[0]
           mtch = qq            ; want sdsstab reference
        endif 
     endelse 

     ;; If no match then move on
     if nmtch eq 0 then continue
     balflg[mtch] = baltab[ibal].bal_flag
  endfor                        ; loop qq=nfil

  if keyword_set(flgonly) then return, balflg $
  else return, where(balflg gt 0,count)
end                             ; sdss_fndbalqso()


pro sdss_instantbalflg, sdsssum, debug=debug
  ;; Instantiate sdsscontistrct.balflg from Shen et al. (2011) BAL QSO
  ;; catalog 
  if n_params() ne 1 then begin
     print,'Syntax - sdss_instantbalflg, sdsssum, [/debug]'
     return
  endif 
  sdssdir = sdss_getsdssdir()

  if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)
  nsdss = (size(sdsstab,/dim))[0]
  
  ibal = sdss_fndbalqso(sdsstab,count=nbal,balflg=balflg)
  if nbal eq 0 then return

  for bb=0L,nbal-1 do begin 
     cfil = sdss_getname(sdsstab[ibal[bb]],/abslin,dir=cdir)
     cfil = cfil[0]
     cstrct = xmrdfits(sdssdir+cdir+cfil,1,/silent)
     cstrct.balflg = balflg[ibal[bb]]
     mwrfits,cstrct,sdssdir+cdir+cfil,/create,/silent
     spawn,'gz -f '+sdssdir+cdir+cfil
  endfor                        ; loop bb=nbal

  print,'sdss_instantbalflg: Number of BALs flagged',nbal

end                             ; sdss_instantbalflg


function sdss_getspecwave
  ;; Retun SDSS wavelength range
  return, [3820.d,9200.d]
end                             ; sdss_getspecwave


function sdss_getspecpixscale, loglam=loglam
  ;; Retun SDSS pixel scale
  val = alog(10) * 0.0001 
  if keyword_set(loglam) then return, val $ ; dwv = val * wave
  else return, val * 299792.458 ; 69 km/s
end                             ; sdss_getspecpixscale


pro sdss_getrandqsosmpl, nrand, ostrct_fil, list_fil=list_fil, $
                         zlim=zlim, rmaglim=rmaglim, bal=bal, _extra=extra
  ;; Draw random sampe of QSOs with distribution of whole
  if n_params() ne 2 then begin
     print,'Syntax - sdss_getrandqsosmpl, nrand, ostrct_fil, [list_fil=,'
     print,'                  /bal, zlim=,rmaglim=,_extra=]'
     return
  endif 
  
  sdsstab = sdss_getqsostrct()
  if keyword_set(bal) then begin
     ibal = sdss_fndbalqso(sdsstab,count=nsdss,/full) ; full is short cut
     sdsstab = sdsstab[ibal]
  endif else nsdss = (size(sdsstab,/dim))[0]

  ;; Subsets
  if keyword_set(zlim) then begin
     gd = where(sdsstab.z ge zlim[0] and sdsstab.z le zlim[1],nsdss)
     if nsdss eq 0 then begin
        print,'sdss_getrandqsosmpl: no QSOs in z cut'
        return
     endif 
     sdsstab = sdsstab[gd]
  endif 
  if keyword_set(rmaglim) then begin
     gd = where(sdsstab.psf_r ge rmaglim[0] and sdsstab.psf_r le rmaglim[1],nsdss)
     if nsdss eq 0 then begin
        print,'sdss_getrandqsosmpl: no QSOs in R-magnitude cut'
        return
     endif 
     sdsstab = sdsstab[gd]
  endif 


  ;; Get a random sample
  ;; See http://www.idlcoyote.com/code_tips/randperm.html
  x = lindgen(nsdss)
  y = randomu(dseed, nsdss)     ; uniformly scrambled #s
  z = x[sort(y)]
  if nrand lt nsdss then $
     substrct = sdsstab[z[0:nrand-1]] $; take first bit 
  else begin
     stop,'sdss_getrandqsosmpl: With cuts, not enough QSOs ',nsdss
     substrct = sdsstab
  endelse 

  ;; Sort by plate name
  spec_fil = sdss_getname(substrct,plate=plate)
  srt = sort(plate)
  plate = plate[srt]
  spec_fil = spec_fil[srt]
  substrct = substrct[srt]

  if keyword_set(list_fil) then $
     sdss_wrqsolist,substrct,list_fil,_extra=extra ; absdir=, /silent

  if size(ostrct_fil,/type) eq 7 then begin
     mwrfits,substrct,ostrct_fil,/create,/silent
     print,'sdss_getrandqsosmpl: created ',ostrct_fil
  endif else ostrct_fil = substrct
  
end                             ; sdss_getrandqsosmpl



function sdss_calcnormerr, flux, error, conti_fil, cflg=cflg, unnorm=unnorm, $
                           baderrval=baderrval
  ;; Compute the total error for spectrum, including continuum and
  ;; flux errors.
  ;; Can parse conti_fil as file name, sdsscontistrct format, or as
  ;; normal [npix,3] array
  if n_params() ne 3 then begin
     print,'Syntax - sdss_calcnormerr(flux,error,conti_fil,[cflg=,/unnorm,baderrval=])'
     return,-1
  endif 
  sdssdir = sdss_getsdssdir()
  if keyword_set(unnorm) then power = 1. else power = 2.
  if not keyword_set(baderrval) then baderrval = 0. ; for badpix
  npix = (size(flux,/dim))[0]
  bdpix = where(error eq 0.)

  ;; Read in many options
  case size(conti_fil,/type) of 
     7: begin                   ; file name
        conti = xmrdfits(sdssdir+conti_fil,0,/silent)
        if size(conti,/n_dimen) ne 2 then begin
           cstrct = xmrdfits(sdssdir+conti_fil,1,/silent) ; next extension
           if size(cstrct,/type) ne 8 then $
              stop,'sdss_calcnormerr(): non-sensical input file ',conti_fil
           if not keyword_set(cflg) then $
              stop,'sdss_calcnormerr(): must set cflg'
           ;; Use structure
           conti = dblarr(cstrct.npix,3)
           cindx = fix(alog(cflg)/alog(2))
           conti[*,0] = cstrct.conti[0:cstrct.npix-1,cindx]
           conti[*,2] = cstrct.sigconti[0:cstrct.npix-1,cindx]
           ipix0 = cstrct.ipix0
        endif 
     end 
     8: begin                   ; structure
        cstrct = conti_fil
        if keyword_set(cflg) then cstrct.cflg = cflg 
        if not keyword_set(cstrct.cflg) then $
           stop,'sdss_calcnormerr(): must set abslin cflg or cflg keyword'
        ;; Use structure
        conti = dblarr(cstrct.npix,3)
        cindx = fix(alog(cstrct.cflg)/alog(2))
        conti[*,0] = cstrct.conti[0:cstrct.npix-1,cindx]
        conti[*,2] = cstrct.sigconti[0:cstrct.npix-1,cindx]
        ipix0 = cstrct.ipix0
     end
     else: conti = conti_fil    ; better be right already
  endcase
  if not keyword_set(ipix0) then begin
     ;; Figure out ipix0
     bd = where(conti[*,0] eq 0. or finite(conti[*,0],/nan))
     if bd[0] eq -1 then ipix0 = 0 $
     else begin
        istrt = where(bd ne shift(bd,1)+1,ngap)
        istop = where(bd ne shift(bd,-1)-1,ntest)
        ipix0 = istop[0]        ; just the beginning
     endelse 
  endif 

  if npix ne (size(conti,/dim))[0] then $
     stop,'sdss_calcnormerr(): input error and continuum not same size'
  normerr = error * 0.
  gd = where(conti[*,0] ne 0. and finite(conti[*,0]),complement=bd)
  normerr[gd] = sqrt(error[gd]^2*conti[gd,0]^2 + $
                     conti[gd,2]^2*flux[gd]^2)/conti[gd,0]^power
  if ipix0 ne 0 then begin
     normerr[0:ipix0] = error[0:ipix0] 
     if not keyword_set(unnorm) then begin
        gd = where(conti[0:ipix0,0] ne 0. and finite(conti[0:ipix0]),complement=bd)
        if gd[0] ne -1 then normerr[gd] = normerr[gd]/conti[gd,0]
        if bd[0] ne -1 then normerr[bd] = baderrval
     endif 
  endif 
  if bdpix[0] ne -1 then normerr[bdpix] = baderrval

  return, normerr
end                             ;  sdss_calcnormerr()


pro sdss_pltconti, spec_fil, all=all, norm=norm, zabs=zabs, $
                   _extra=extra
  ;; Read in and plot the spectrum with continuum
  ;; if /all set, still have to put /eig or /spl
  if n_params() ne 1 then begin
     print,'Syntax - sdss_pltconti, spec_fil, [/all, /norm, zabs=, _extra=]'
     return
  end 
  sdssdir = sdss_getsdssdir()
  if keyword_set(all) and keyword_set(norm) then begin
     print,'sdss_pltconti: cannot display all spline and eigen-normalized spectra.'
     return
  endif 

  ;; Find file
  sfil = sdss_getname(spec_fil,/spec,dir=dir)
  test = file_search(sdssdir+dir+sfil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_pltconti: spec_fil DNE as given or in full path'

  ;; Read in files (use test b/c has full path and any .gz)
  parse_sdss,test,fx,wv,sig=er0,zqso=zqso ; zqso may be bad!

  ;; Find continuum (_extra includes /spl or /eig or /abslin)
  cfil = sdss_getname(spec_fil,/spec,dir=dir,_extra=extra)

  ;; Find continuum file
  test = file_search(sdssdir+dir+cfil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_pltconti: continuum file DNE ',dir+cfil
  
  ;; Read in from explicitly stored continuum
  ;; Test kind of file
  cstrct = xmrdfits(test,0,/silent) ; array or header
  if size(cstrct,/n_dimen) eq 2 then begin
     conti = cstrct[*,0]
     ;; er = sqrt(er0^2*cstrct[*,0]^2 + cstrct[*,2]^2*fx^2)/cstrct[*,0]
     er = sdss_calcnormerr(fx,er0,cstrct,/unnorm)

     if keyword_set(all) then begin
          cfil2 = sdss_getname(spec_fil,/spec,/spl)
          if cfil2 eq cfil then $
             cfil2 = sdss_getname(spec_fil,/spec,/eig) ;  try other
          tmp = xmrdfits(sdssdir+dir+cfil2,0,/silent)
          altconti = tmp[*,0]
     endif 
     print,'sdss_pltconti: Using spectrum header redshift, which may be wrong'
  endif else begin
     ;; Figure out which version
     cstrct = xmrdfits(test,1,/silent) ; structure
     zqso = cstrct.z_qso               ; better than the header!
     cindx = fix(alog(cstrct.cflg)/alog(2))
     conti = cstrct.conti[0:cstrct.npix-1,cindx]
     ;; er = sqrt(er0^2*cstrct.conti[0:cstrct.npix-1,cindx]^2 + $
     ;;           cstrct.sigconti[0:cstrct.npix-1,cindx]^2*fx^2)/$
     ;;           cstrct.conti[0:cstrct.npix-1,cindx]
     er = sdss_calcnormerr(fx,er0,cstrct,/unnorm,cflg=cstrct.cflg)
     if cstrct.cflg eq sdss_getcflg(/eig) or $
        cstrct.cflg eq sdss_getcflg(/hyb) then begin
        altconti = cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/spl,/index)]
     endif else begin
        conti = cstrct.conti[0:cstrct.npix-1,cindx]
        altconti = cstrct.conti[0:cstrct.npix-1,sdss_getcflg(/eig,/index)]
     endelse 
  endelse
  if keyword_set(all) then begin
     psym3 = -3 
     print,'sdss_pltconti: blue is other conti'
  endif else begin 
     altconti = er0
     psym3 = 10
     print,'sdss_pltconti: blue is original error'
  endelse 

  if keyword_set(norm) then begin
     fx = fx/conti
     er = er/conti
     conti = 0
  endif 

  ;; Plot
  if keyword_set(zabs) then zin = zabs else zin = zqso
  x_specplot,fx,er,wav=wv,inflg=4,ytwo=conti,psym2=-3,$
             ythr=altconti,psym3=psym3,zin=zin,/lls,$
             title=cfil+' zqso = '+string(zqso,format='(f8.5)'), /block
end                             ; sdss_pltconti


pro sdss_pltsnr, spec_fil, lsnr=lsnr, cfil2=cfil2, debug=debug,$
                 _extra=extra
  ;; Plot convolved S/N used to find absorption lines
  ;; _extra= passed to x_specplot and includes zin=, /lls, etc.
  if n_params() ne 1 then begin
     print,'Syntax - sdss_pltsnr, spec_fil, [lsnr=, _extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()
  if not keyword_set(lsnr) then lsnr = 3.5 ; to match sdss_fndlin
  
  cfil = sdss_getname(spec_fil,/spec,/abslin,dir=cdir)    
  
  cstrct = xmrdfits(sdssdir+cdir[0]+cfil[0],1,/silent)        
  
  parse_sdss,sdssdir+cstrct.sdss_obs[0],fx,wv
  
  cindx = fix(alog(cstrct.cflg)/alog(2))
  gd = where(finite(cstrct.snr_conv[0:cstrct.npix-1,cindx]),ngd)      
  title = strtrim(spec_fil,2)+$
          ' zqso = '+string(cstrct.z_qso,format='(f8.5)')

  if keyword_set(cfil2) then begin
     if size(cfil2,/type) eq 8 then cstrct2 = cfil2 $
     else cstrct2 = xmrdfits(sdssdir+cfil2,1,/silent)
     cindx2 = fix(alog(cstrct2.cflg)/alog(2))
     gd2 = where(finite(cstrct2.snr_conv[0:cstrct2.npix-1,cindx2]),ngd)
     x_specplot,cstrct.snr_conv[gd,cindx],replicate(lsnr,ngd),wave=wv[gd],inflg=4,$
                /block,title=title,_extra=extra,$
                two_wave=wv[gd2], ytwo=cstrct2.snr_conv[gd2,cindx2]
  endif else $
     x_specplot,cstrct.snr_conv[gd,cindx],replicate(lsnr,ngd),wave=wv[gd],inflg=4,$
                /block,title=title,_extra=extra

  if keyword_set(debug) then $
     stop,'sdss_pltsnr debug: stopping before exiting'
end                             ; sdss_pltsnr


function sdss_contiqual, spec_fil, wvlim=wvlim, _extra=extra
  ;; Measure the continuum chi^2
  if n_params() ne 1 then begin
     print,'Syntax - sdss_contiqual(spec_fil, [wvlim=,_extra=])'
     return,[-1,-1]
  end 
  sdssdir = sdss_getsdssdir()

  ;; Find file
  sfil = sdss_getname(spec_fil,/spec,dir=dir)
  test = file_search(sdssdir+dir+sfil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_contiqual(): spec_fil DNE as given or in full path'

  ;; Read in files (use test b/c has full path and any .gz)
  parse_sdss,test,fx,wv,sig=er,npix=npix

  ;; Find continuum (_extra includes /spl or /eig or /abslin)
  cfil = sdss_getname(spec_fil,/spec,dir=dir,_extra=extra)

  ;; Find continuum file
  test = file_search(sdssdir+dir+cfil+'*',count=ntest)
  if ntest eq 0 then $
     stop,'sdss_contiqual(): continuum file DNE ',dir+cfil
  
  ;; Read in from explicitly stored continuum
  ;; Test kind of file
  cstrct = xmrdfits(test,0,/silent) ; array or header
  if size(cstrct,/n_dimen) eq 2 then begin
     conti = cstrct[*,0]
     ;; Test what error should be becasue eigconti and splconti stored
     ;; differently; original er = 0. indicates gap
     test = where(er lt cstrct[*,2] and er ne 0.,ntest)
     if ntest eq 0 then $
        ;; er = sqrt(er^2*cstrct[*,0]^2 + cstrct[*,2]^2*fx^2)/conti^2 $
                er = sdss_calcnormerr(fx,er,cstrct) $
     else er = cstrct[*,2]
  endif else begin
     ;; Figure out which version
     cstrct = xmrdfits(test,1,/silent) ; structure
     cindx = fix(alog(cstrct.cflg)/alog(2))
     conti = cstrct.conti[0:cstrct.npix-1,cindx]
     er = cstrct.sigconti[0:cstrct.npix-1,cindx]
  endelse
  ;; Normalize 
  fx = fx/conti
  
  ;; Trim
  if keyword_set(wvlim) then begin
     gd = where(wv ge wvlim[0] and wv le wvlim[1],npix)
     if npix eq 0 then begin
        print,'sdss_contiqual(): no spectrum in wavelength limits'
        return,[-1,-1]
     endif 
     fx = fx[gd]
     er = er[gd]
  endif 

  ;; RMS and ...
  rms = sqrt(total(fx^2)/npix)
;  sigrms = sqrt( total((fx*er)^2) / (npix * total(fx^2)) ) ; weighted?
  sigrms = sqrt(total(er^2)/npix)

  return,[rms,sigrms]
end                             ; sdss_contqual()


function sdss_measuresnr, spec_fil, list=list, wvlim_obs=wvlim_obs, $
                          dvgal=dvgal, dvqso=dvqso, dblt_name=dblt_name, $
                          zqso=zqso, sdsssum=sdsssum, no_snr=no_snr, $
                          gz_fil=gz_fil, nsig=nsig, clobber=clobber, $
                          silent=silent
  ;; Measure S/N for unnormalized input spectrum
  if n_params() ne 1 then begin
     print,'Syntax - sdss_measuresnr(spec_fil, [/list, wvlim_obs=, dvgal='
     print,'                         dvqso=, zqso=, sdsssum=, dblt_name=, /no_snr'
     print,'                         gz_fil=, nsig=, /clobber, /silent])'
     return,-1
  endif 

  if not keyword_set(dblt_name) then dblt_name = 'CIV'
  if not keyword_set(dvgal) then dvgal = 5000.  ; km/s
  if not keyword_set(dvqso) then dvqso = -3000. ; km/s

  if size(dblt_name,/type) eq 7 then $
     dblt = dblt_retrieve(dblt_name) $
  else dblt = dblt_name
  wvspecmnx = sdss_getspecwave()
  cinv = 1./2.998e5             ; km^-1 s
  wvr_lim = [1250.,1.e5]        ; matches sdss_fndciv

  if dblt.ion eq 'CIV' then wvr_lim[0] = 1310. ; exclude OI/SiII forest

  sdssdir = sdss_getsdssdir()
  if keyword_set(list) then begin
     readcol,spec_fil,file,format='(a)',skip=1,/silent
     sv_spec_fil = spec_fil     ; will revert later
     spec_fil = file
  endif  
  nfil = (size(spec_fil,/dim))[0] > 1 ; foils singularity

  if keyword_set(nsig) then begin
     pixscale = sdss_getspecpixscale(/loglam) 

     if keyword_set(gz_fil) then begin
        if ((size(gz_fil,/dim))[0] > 1) ne nfil then $
           stop,'sdss_measuresnr(): gz_fil= size does not match spec_fil size'
     endif 
  endif else if keyword_set(gz_fil) then $
     stop,'sdss_measuresnr(): must sent nsig if setting gz_fil='

  if not keyword_set(sdsssum) and not keyword_set(zqso) then begin
     print,'sdss_measuresnr(): sdsssum and zqso not set so measure S/N for full spectrum.'
     zqso = replicate(100.,nfil) ; something large
  endif else begin
     if keyword_set(sdsssum) then begin
        if size(sdsssum,/type) eq 8 then sdsstab = sdsssum $
        else sdsstab = xmrdfits(sdsssum,1,/silent)
        if (size(sdsstab,/dim))[0] ne nfil then $
           stop,'sdss_measuresnr(): sdsstab not same size as input file list'
        ;; Assume everything in right order
        zqso = sdsstab.z
     endif
  endelse 

  ;; Set wavelength bounds for comparison
  ;; and it's ALLOWED to have wvlim_obs[*,1] < wvlim_obs[*,0]
  ;; but S/N will equal -9.99
  wvlim_obs = fltarr(nfil,2,/nozero)
  wvlim_obs[*,0] = (wvspecmnx[0] > (wvr_lim[0]*(1.+zqso))) > $
                   (dblt.wvI*(1.+dvgal*cinv))
  wvlim_obs[*,1] = wvspecmnx[1] < (dblt.wvII*(1.+zqso)*(1.+dvqso*cinv))
  
  ;; Exit point; very fast
  if keyword_set(no_snr) then begin
     if nfil eq 1 then $
        ;; Return normal looking arrays
        wvlim_obs = wvlim_obs[0,*]
     if keyword_set(sv_spec_fil) then spec_fil = sv_spec_fil
     return,-1
  endif 


  ;; ;;;;;;;
  snr = fltarr(nfil,3,/nozero)  ; median flux, error, median(flux/error)
  for ss=0L,nfil-1 do begin 
     parse_sdss, sdssdir+spec_fil[ss], flux, wave, npix=npix, sig=error, $
                 head=hdr
     ;; gz array will be [npix, [zabs, dz, EWlim, mask, nsig]]
     gz_los = fltarr(npix,5,/nozero) ; instantiate efficiently
     gz_los[*,0] = wave/dblt.wvI - 1.
     gz_los[*,1] = pixscale * wave     ; dlambda (later dz)
     gz_los[*,2] = !values.f_infinity   
     gz_los[*,3] = 0                   ; mask (0 = bad; 1 = useful)
     gz_los[*,4] = nsig                ; useful

     sub = where(wave ge wvlim_obs[ss,0] and wave le wvlim_obs[ss,1] and $
                 error ne 0.,nsub)
     if nsub ne 0 then begin
        snr_arr = flux[sub]/error[sub]
        snr[ss,*] = [median(flux[sub],/even),$
                     median(error[sub],/even),$
                     median(snr_arr,/even)] 

        if keyword_set(nsig) then begin
           ;; Estimate g(z) using Danforth & Shull (2008) EW_lim
           ;; EW_lim = nsig * lambda / (R_lambda * (S/N)_lambda)
           dispstrct = xmrdfits(sdssdir+spec_fil[ss],6,/silent)
           if size(dispstrct,/type) eq 8 then $
              dispersion = dispstrct.dispersion $
           else begin
              print,'sdss_measuresnr(): dispersion structure DNE ',spec_fil[ss]
              dispersion = replicate(1,npix)
           endelse
           gz_los[sub,2] = nsig * gz_los[sub,1] / $ ; why gz_los[*,1] dlambda
                           (dispersion[sub] * snr_arr)
           gz_los[sub,3] = 1                  ; useful, qualified pixels

           bd = where(dispersion[sub] eq 0,complement=gd)
           if bd[0] ne -1 then begin
              ;; Don't leave these as infinite though don't know why
              ;; would ever get here
              ;; likely never will
              gz_los[sub[bd],2] = interpol(gz_los[sub[gd],2],$
                                           gz_los[sub[gd],0],$
                                           gz_los[sub[bd],0])
              print,'sdss_measuresnr(): interpolated some pixels ',spec_fil[ss]
           endif                ; interpolate dispersion

           bd = where(gz_los[sub,2] le 0.)
           if bd[0] ne -1 then $
              gz_los[sub[bd],3] = 0 ; back to bad

        endif                   ; keyword_set(nsig)
     endif else $
        snr[ss,*] = -9.99

     gz_los[*,1] = gz_los[*,1]/dblt.wvI ; now dz
     ;; I think I have to correct for the fact that the resolution
     ;; element is larger/smaller than the pixel scale and I'm
     ;; double counting in cases
     if keyword_set(dispersion) then begin
        gd = where(dispersion gt 0.)
        if gd[0] ne -1 then $
           gz_los[gd,1] = gz_los[gd,1]/dispersion[gd]
     endif


     if keyword_set(nsig) then begin
        if size(gz_fil,/type) eq 7 then begin
           ;; Write file
           test = file_search(sdssdir+gz_fil[ss]+'*',count=ntest)
           if ntest eq 0 or keyword_set(clobber) then begin
              if ntest eq 0 then begin
                 ;; Check if need to make directory
                 dir = strmid(gz_fil[ss],0,strpos(gz_fil[ss],'/',$
                                                  /reverse_search))
                 dtest = file_search(sdssdir+dir,count=ndtest)
                 if ndtest eq 0 then $
                    spawn,'mkdir -p '+dir
              endif 
              mwrfits,gz_los,sdssdir+gz_fil[ss],hdr,/create,/silent
              spawn,'gzip -f '+sdssdir+gz_fil[ss]
              if not keyword_set(silent) then $
                 print,'sdss_measuresnr(): created ',gz_fil[ss]
           endif                   ; created file
           
        endif else gz_fil = gz_los ; return array
     endif                         ; keyword_set(nsig)

  endfor                        ; loop ss=nfil

  if nfil eq 1 then begin
     ;; Return normal looking arrays
     snr = snr[0,*]
     wvlim_obs = wvlim_obs[0,*]
  endif 
  
  if keyword_set(sv_spec_fil) then spec_fil = sv_spec_fil
  
  return, snr
end                             ; sdss_measuresnr()



pro sdss_mksnrtab, list_fil, sdsssum, out_fil, _extra=extra
  ;; Create structure of all S/N per doublet
  if n_params() ne 3 then begin
     print,'Syntax - sdss_mksnrtab, list_fil, sdsssum, out_fil, [_extra=]'
     return
  endif 
  sdssdir = sdss_getsdssdir()
  dblt_arr = ['SiIV','CIV','MgII','CaII']
  ndblt = (size(dblt_arr,/dim))[0]

  readcol,list_fil,spec_fil,skip=1,format='a',/silent
  nfil = (size(spec_fil,/dim))[0]

  if size(sdsssum,/type) eq 7 then sdsstab = xmrdfits(sdsssum,1,/silent) $
  else sdsstab = sdsssum
  if nfil ne (size(sdsstab,/dim))[0] then $
     stop,'sdss_mksnrtab: list_fil and sdsssum unequal sizes.'

  strct_tmp = {qso_name:'', $
               z_qso:0., $
               SNR_SiIV:fltarr(3,/nozero), $
               wvobs_SiIV:fltarr(2,/nozero), $
               SNR_CIV:fltarr(3,/nozero), $
               wvobs_CIV:fltarr(2,/nozero), $
               SNR_MgII:fltarr(3,/nozero), $
               wvobs_MgII:fltarr(2,/nozero), $
               SNR_CaII:fltarr(3,/nozero), $
               wvobs_CaII:fltarr(2,/nozero) $
              }
  tags = tag_names(strct_tmp)

  snr_strct = replicate(strct_tmp,nfil)
  istrt = strpos(spec_fil[0],'-')
  istop = strpos(spec_fil[0],'.',/reverse_search) 
  snr_strct.qso_name = strmid(spec_fil,istrt+1,istop-istrt-1)
  snr_strct.z_qso = sdsstab.z

  for ii=0,ndblt-1 do begin
     ;; _extra includes dvgal=, dvqso= 
     snr_arr = sdss_measuresnr(spec_fil,wvlim_obs=wvlim_obs,$
                               dblt_name=dblt_arr[ii],zqso=snr_strct.z_qso, $
                              _extra=extra)
     snrtag = where(tags eq 'SNR_'+strupcase(dblt_arr[ii]))
     wvtag = where(tags eq 'WVOBS_'+strupcase(dblt_arr[ii]))
     snr_strct.(snrtag)[0] = snr_arr[*,0]
     snr_strct.(snrtag)[1] = snr_arr[*,1]
     snr_strct.(snrtag)[2] = snr_arr[*,2]
     snr_strct.(wvtag)[0] = wvlim_obs[*,0]
     snr_strct.(wvtag)[1] = wvlim_obs[*,1]
  endfor                        ; for ii=ndblt
  

  ;; Write out
  if size(out_fil,/type) eq 7 then begin
     mwrfits,snr_strct,out_fil,/create,/silent
     spawn,'gzip -f '+out_fil
     print,'sdss_mksnrtab: created ',out_fil
  endif else out_fil = snr_strct

end                             ; sdss_mksnrtab


pro sdss_prntsnrtab, snrstrct_fil, out_fil, snrsumm=snrsumm, zbin=zbin, $
                     dblt_name=dblt_name, zmin=zmin
  ;; Print out a variety of cuts based on S/N
  ;; if dblt_name not set, input all 4 usuals
  ;; if zmin not set, then use all (zmin = 0)
  ;; Default table is QSO ID, zqso, S/N for each dblt_name
  ;; if snrsumm = 1, then print out a by-redshift breakdown for all
  ;; dblt_names for S/N >= 3, 4, and 5
  ;; if snrsumm > 1, then print out list of spectra with S/N >= snrsumm
  if n_params() ne 2 then begin
     print,'Syntax - sdss_prntsnrtab, snrstrct_fil, out_fil, [snrsumm=, zbin=,'
     print,'                          dblt_name=, zmin=]'
     return
  endif 
  if not keyword_set(zbin) then zbin = 0.2
  if not keyword_set(zmin) then zmin = 0.0

  if size(snrstrct_fil,/type) eq 7 then $
     snr_strct = xmrdfits(snrstrct_fil,1,/silent) $
  else snr_strct = snrstrct_fil
  tags = tag_names(snr_strct[0])

  close,1
  openw,1,out_fil

  if keyword_set(snrsumm) then begin
     if keyword_set(dblt_name) then dblt_arr = dblt_name $
     else dblt_arr = ['SiIV','CIV','MgII','CaII'] 
     ndblt = n_elements(dblt_arr) ; do not use size(/dim)


     if snrsumm gt 1 then begin
        ;; Print list
        for ii=0,ndblt-1 do begin
           snrtag = where(tags eq 'SNR_'+strupcase(dblt_arr[ii]))
           gd = where(snr_strct.(snrtag)[2] ge snrsumm and $
                      snr_strct.z_qso ge zmin,ngd)
           prs = strsplit(snr_strct[0].qso_name,'-')
           spec_fil = 'spectro/1d_26/'+strmid(snr_strct[gd].qso_name,prs[1],prs[2]-prs[1]-1)+$
                      '/1d/spSpec-'+snr_strct[gd].qso_name+'.fit'

           printf,1,'abslin/'
           
           writecol,out_fil,spec_fil,snr_strct[gd].(snrtag)[2],fmt='(a,5x,f7.2)',$
                    filnum=1
        endfor                  ; loop ii=ndblt

        close,1
        print,'sdss_prntsnrtab: created spectra list for S/N >= '+$
                 string(snrsumm,format='(f4.1," ")')+out_fil
        
     endif else begin
        ;; Print S/N cut summary
        for ii=0,ndblt-1 do begin

           snrtag = where(tags eq 'SNR_'+strupcase(dblt_arr[ii]))
           
           sub = where(snr_strct.(snrtag)[2] gt 0. and $
                       snr_strct.z_qso ge zmin,nsub) ; has coverage
           zmn = min(snr_strct[sub].z_qso,max=zmx)
           

           gd3 = sub[where(snr_strct[sub].(snrtag)[2] ge 3.,ngd3)]
           hist3 = histogram(snr_strct[gd3].z_qso,loc=loc3,min=zmn,max=zmx,binsize=zbin)

           gd4 = sub[where(snr_strct[sub].(snrtag)[2] ge 4.,ngd4)]
           hist4 = histogram(snr_strct[gd4].z_qso,loc=loc4,min=zmn,max=zmx,binsize=zbin)

           gd5 = sub[where(snr_strct[sub].(snrtag)[2] ge 5.,ngd5)]
           hist5 = histogram(snr_strct[gd5].z_qso,loc=loc5,min=zmn,max=zmx,binsize=zbin)

           printf,1,'' 
           printf,1,dblt_arr[ii]," Total: ",string(nsub,format='(i5)') 
           printf,1,'zbin','S/N=3','4','4/3','5','5/3',$
                  format='(a6,2x,a5,2x,2(a5,1x,"(",a4,")",2x))'
           writecol,out_fil,loc3,hist3,hist4,hist4/float(hist3),hist5,hist5/float(hist3),$
                    fmt='(f6.2,2x,i5,2x,2(i5,1x,"(",f4.2,")",2x))' ,filnum=1
           printf,1,"Total",total(hist3),total(hist4),total(hist5), $
                  format='(a6,2x,i5,2x,2(i5,9x))' 
           
        endfor                  ; loop ii=ndblt
        close,1
        print,'sdss_prntsnrtab: created S/N-z summary tables ',out_fil

     endelse 

  endif else begin
     ;; Print all S/N information
     printf,1,"DR7 QSO S/N Summary (for regions where ion's detectable)"
     printf,1,'QSO ID','zqso','SiIV','CIV','MgII','CaII',$
            format='(a15,2x,a7,4(1x,a6))'

     writecol,out_fil,snr_strct.qso_name,snr_strct.z_qso,$
              snr_strct.SNR_SiIV[2],snr_strct.SNR_CIV[2],$
              snr_strct.SNR_MgII[2],snr_strct.SNR_CaII[2],$
              fmt='(a15,2x,f7.5,4(1x,f6.2))',filnum=1
     close,1
     print,'sdss_prntsnrtab: created all-S/N summary table ',out_fil
  endelse 

end                             ; sdss_prntsnrtab


pro sdss_reorganize, help=help
  ;; Move old EIGCONTI/, SPLCONTI/, and ABSLIN/ stuff to new directory structure
  ;; Be in the right directory (also tests that $SDSSDR set)
  if keyword_set(help) then begin
     print,'Syntax - sdss_reorganize [/help]'
     return
  endif 
  pwd, cur_dir
  cd,sdss_getsdssdir()

  ;; Do it for all
  sdsstab = sdss_getqsostrct() 
  nsdss = (size(sdsstab,/dim))[0] 
  
  ;; Generate appropriate names
  spl_fil = sdss_getname(sdsstab,/spl,plate=plate,dir=cdir) 
  eig_fil = sdss_getname(sdsstab,/eig) 
  fnd_fil = sdss_getname(sdsstab,/abslin,dir=absdir) 

  ;; Loop by plates (since have to make directories)
  unq = uniq(plate,sort(plate))
  nunq = (size(unq,/dim))[0]

  for pp=0L,nunq-1 do begin 
     mtch = where(plate eq plate[unq[pp]],nmtch) 
     
     if nmtch eq 0 then begin 
        print,'sdss_reorganizeo: no match!' 
        return 
     endif 

     spawn,'mkdir -p '+cdir[unq[pp]]
     odir = '.gz '+cdir[unq[pp]] + '/.' 

     for ff=0,nmtch-1 do spawn,'rsync -ravu EIGCONTI/'+spl_fil[mtch[ff]]+odir 
     for ff=0,nmtch-1 do spawn,'rsync -ravu SPLCONTI/'+eig_fil[mtch[ff]]+odir 

     spawn,'mkdir -p '+absdir[unq[pp]]
     odir = '.gz '+absdir[unq[pp]] + '/.' 

     for ff=0L,nmtch-1 do spawn,'rsync -ravu ../ABSLIN/'+fnd_fil[mtch[ff]]+odir 

  endfor                        ; loop pp=nunq
  cd, cur_dir
end                             ; sdss_reorganize


pro sdss_srtbalqso, strct_fil, bal_fil, qso_fil
  ;; Take sdsscivstrct and sort between candidates in BAL QSO
  ;; sightlines and those in regular QSOs
  if n_params() ne 3 then begin
     print,'Syntax - sdss_srtbalqso, strct_fil, bal_fil, qso_fil'
     return
  endif 

  if size(strct_fil,/type) eq 8 then begin
     ;; Input is structure; sort and either return as structure or
     ;; write out to files
     strct = strct_fil 

     if tag_exist(strct,'BALFLG') then $
        qso = where(strct.balflg eq 0,complement=bal) $
     else qso = where(strct.flg_bal eq 0,complement=bal) 

     if qso[0] ne -1 then begin
        if size(qso_fil,/type) eq 7 then begin
           mwrfits,strct[qso],qso_fil,/create,/silent
           print,'sdss_srtbalqso: sorted QSOs into ',qso_fil
        endif else qso_fil = strct[qso]
     endif 
     if bal[0] ne -1 then begin
        if size(bal_fil,/type) eq 7 then begin
           mwrfits,strct[bal],bal_fil,/create,/silent
           print,'sdss_srtbalqso: sorted BAL QSOs into ',bal_fil
        endif else bal_fil = strct[bal]
     endif 
  endif else begin
     ;; Input is file which may have two extensions 
     done = 0
     ext = 1
     while not done do begin
        strct = xmrdfits(strct_fil,ext,/silent)
        if size(strct,/type) ne 8 then begin
           done = 1             ; EOF
           continue
        endif 

        if tag_exist(strct,'BALFLG') then $
           ;; Either {sdsscivstrct}
           qso = where(strct.balflg eq 0,complement=bal) $
           ;; Else {qalcharstrct}
        else qso = where(strct.flg_bal eq 0,complement=bal)  
        
        if qso[0] ne -1 then begin
           mwrfits,strct[qso],qso_fil,create=(ext mod 2),/silent
           print,'sdss_srtbalqso: sorted QSOs into ',qso_fil,$
                 ' ext = '+strtrim(ext,2)
        endif 
        if bal[0] ne -1 then begin
           mwrfits,strct[bal],bal_fil,create=(ext mod 2),/silent
           print,'sdss_srtbalqso: sorted BAL QSOs into ',bal_fil,$
                 ' ext = '+strtrim(ext,2)
        endif 

        ext = ext + 1
     endwhile
  endelse 

end                             ; sdss_srtbalqso



pro sdss_prntdr7summ, help=help
  ;; Print out some facts
  if keyword_set(help) then begin
     print,'Syntax - sdss_prntdr7summ [/help]'
     return
  endif 

  sdsstab = sdss_getqsostrct()
  file = sdss_getname(sdsstab,/abslin) ; dummy name
  wvspecmnx = sdss_getspecwave()
  nsdss = (size(sdsstab,/dim))[0]
  print,strtrim(nsdss,2)+' (100%) -- Total number of DR7 QSOs '
  nsdss = float(nsdss); for fractions
  
  balflg = sdss_fndbalqso(sdsstab,/full,count=nbal,/flgonly)
  print,' '+strtrim(nbal,2)+' ('+string(nbal/nsdss*100.,format='(f4.1)')+'%)'+$
        ' -- Number of BAL QSOs '
  
  ;; various redshift cuts reflect sdss_civsearch/sdss_fndciv 
  ;; ;;;;;;;;;
  ;; CIV
  ;; ;;;;;;;;;
  civ = dblt_retrieve('CIV')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_civ,dblt_name=civ,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_civ[*,0] ge wvspecmnx[0] and $
             wvobs_civ[*,1] le wvspecmnx[1] and $
             wvobs_civ[*,0] lt wvobs_civ[*,1],ngd)
  zmin = min(wvobs_civ[gd,0]/civ.wvI-1.)
  zmax = max(wvobs_civ[gd,1]/civ.wvII-1.)
  print,strtrim(ngd,2)+$
        ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zciv  < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+$
        ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; CIV with Lya
  ;; ;;;;;;;;;
  lya = dblt_retrieve('Lya')  
  lya.wvII = civ.wvI
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_civ_lya,dblt_name=lya,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_civ_lya[*,0] ge wvspecmnx[0] and $
             wvobs_civ_lya[*,1] le wvspecmnx[1] and $
             wvobs_civ_lya[*,0] lt wvobs_civ_lya[*,1],ngd)
  zmin = min(wvobs_civ_lya[gd,0]/civ.wvI-1.)
  zmax = max(wvobs_civ_lya[gd,1]/civ.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and Lya coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zciv  < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and Lya coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; CIV & SiIV: this seems wrong...
  ;; ;;;;;;;;;
  siiv = dblt_retrieve('SiIV')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_siiv,dblt_name=siiv,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_siiv[*,0] ge wvspecmnx[0] and $
             wvobs_siiv[*,1] le wvspecmnx[1] and $
             wvobs_siiv[*,0] lt wvobs_siiv[*,1] and $
             wvobs_civ[*,0] ge wvspecmnx[0] and $
             wvobs_civ[*,1] le wvspecmnx[1] and $
             wvobs_civ[*,0] lt wvobs_civ[*,1],ngd)
  zmin = min(wvobs_civ[gd,0]/civ.wvI-1.)
  zmax = max(wvobs_civ[gd,1]/civ.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and SiIV coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zciv  < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CIV and SiIV coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; SiIV
  ;; ;;;;;;;;;
  gd = where(wvobs_siiv[*,0] ge wvspecmnx[0] and $
             wvobs_siiv[*,1] le wvspecmnx[1] and $
             wvobs_siiv[*,0] lt wvobs_siiv[*,1],ngd)
  zmin = min(wvobs_siiv[gd,0]/siiv.wvI-1.)
  zmax = max(wvobs_siiv[gd,1]/siiv.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with SiIV coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zsiiv < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with SiIV coverage (excl. BAL QSOs) '


  ;; ;;;;;;;;;
  ;; SiV with Lya
  ;; ;;;;;;;;;
  lya.wvII = siiv.wvI
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_siiv_lya,dblt_name=lya,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_siiv_lya[*,0] ge wvspecmnx[0] and $
             wvobs_siiv_lya[*,1] le wvspecmnx[1] and $
             wvobs_siiv_lya[*,0] lt wvobs_siiv_lya[*,1],ngd)
  zmin = min(wvobs_siiv_lya[gd,0]/siiv.wvI-1.)
  zmax = max(wvobs_siiv_lya[gd,1]/siiv.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with SIIV and Lya coverage '+$
        string(min(zmin[gd]),zmax,format='(f5.3," < zsiiv < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)',+$
        ' -- Number with SIIV and Lya coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; CaII
  ;; ;;;;;;;;;
  caii = dblt_retrieve('CaII')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_caii,dblt_name=caii,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_caii[*,0] ge wvspecmnx[0] and $
             wvobs_caii[*,1] le wvspecmnx[1] and $
             wvobs_caii[*,0] lt wvobs_caii[*,1],ngd)
  zmin = min(wvobs_caii[gd,0]/caii.wvI-1.)
  zmax = max(wvobs_caii[gd,1]/caii.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CaII coverage '+$
        string(zmin,zmax,format='(f5.3," < zcaii < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+ ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with CaII coverage (excl. BAL QSOs) '

  ;; ;;;;;;;;;
  ;; MgII
  ;; ;;;;;;;;;
  mgii = dblt_retrieve('MgII')  
  dum = sdss_measuresnr(file,wvlim_obs=wvobs_mgii,dblt_name=mgii,$
                        zqso=sdsstab.z,/no_snr)
  gd = where(wvobs_mgii[*,0] ge wvspecmnx[0] and $
             wvobs_mgii[*,1] le wvspecmnx[1] and $
             wvobs_mgii[*,0] lt wvobs_mgii[*,1],ngd)
  zmin = min(wvobs_mgii[gd,0]/mgii.wvI-1.)
  zmax = max(wvobs_mgii[gd,1]/mgii.wvII-1.)
  print,strtrim(ngd,2)+' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with MgII coverage '+$
        string(zmin,zmax,format='(f5.3," < zmgii < ",f5.3)')
  gd = gd[where(balflg[gd] eq 0,ngd)]
  print,strtrim(ngd,2)+ ' ('+string(ngd/nsdss*100.,format='(f5.1)')+'%)'+$
        ' -- Number with MgII coverage (excl. BAL QSOs) '


  print,'CIV cut excludes OI/SiII "forest."'
  print,'SiIV cut excludes Lya forest.'
  print,'CaII cut exlcudes HVCs (up to 5000 km/s).'
  print,'MgII for z > 0.'

end                             ; sdss_prntdr7summ


function sdss_calcparalleljob, sdsssum, processor, count=count
  ;; Paired with sdss_runparallelsrch.sh can divide big structure and
  ;; run sdss_fndlin and sdss_civsearch on multiple processors
  if n_params() ne 2 then begin
     print,'Syntax - sdss_calcparalleljob(sdsssum, processor, [count=]'
     return, -1
  endif 

  if size(sdsssum,/type) eq 8 or (size(sdsssum,/dim))[0] gt 0 then sdsstab = sdsssum $
  else sdsstab = xmrdfits(sdsssum,1,/silent)
  nsdss = (size(sdsstab,/dim))[0]   ; if full DR7, 105783
  nsub = nsdss/processor[1]      ; if processor[1]=4, nsub = 26445
  istrt = (processor[0]-1)*nsub 
  ;; Stuff the remainder into the last processor job
  if processor[0] eq processor[1] then iend = nsdss - 1 $
  else iend = processor[0]*nsub - 1
  count = iend - istrt + 1L
  
  return,[istrt,iend]                        ; make sure to exit for script
end                             ; sdss_calcparalleljob()


pro sdss_catparalleljob, cand_fil, processor
  ;; Paired with sdss_runparallelsrch.sh can combine sub-structures  
  ;; from sdss_civsearch into one
  if n_params() ne 2 then begin
     print,'Syntax - sdss_catparalleljob, cand_fil, processor'
     return
  endif
  
  ;; Load conventional name from sdss_runparallelsrch.sh
  sub_fil = cand_fil+'.'+strtrim((lindgen(processor[1])+1),2)
  for pp=1,processor[1] do begin
     strct = xmrdfits(sub_fil[pp-1],1,/silent)

     if pp eq 1 then begin
        if size(strct,/type) eq 8 then allciv = strct ; may be nothing
     endif else begin
        if size(strct,/type) eq 8 then begin
           if keyword_set(allciv) then allciv = [allciv,strct] $
           else allciv = strct
        endif 
     endelse 

     strct = xmrdfits(sub_fil[pp-1],2,/silent)

     if pp eq 1 then begin
        if size(strct,/type) eq 8 then allbroadciv = strct ; may be nothing
     endif else begin
        if size(strct,/type) eq 8 then begin
           if keyword_set(allbroadciv) then allbroadciv = [allbroadciv,strct] $
           else allbroadciv = strct
        endif 
     endelse 
  endfor                        ; loop pp=processor[1]
  
  if keyword_set(allciv) then begin
     allciv = sdss_srtcivstrct(allciv)
     mwrfits,allciv,cand_fil,/create,/silent
  endif else mwrfits,-1,cand_fil,/create,/silent
  if keyword_set(allbroadciv) then begin
     allbroadciv = sdss_srtcivstrct(allbroadciv)
     mwrfits,allbroadciv,cand_fil,/silent ; ext=2
  endif else mwrfits,-1,cand_fil,/silent ; ext=2
  spawn,'gzip -f '+cand_fil
  print,'sdss_catparalleljob: created ',cand_fil
  
  return                        ; make sure to exit for script
end                             ; sdss_catparalleljob


pro sdss_mtchcivstrct, strct1_fil, strct2_fil, mtch_fil, diff_fil, $
                       log_fil=log_fil,final=final, force=force, dvtol=dvtol
  ;; Take sdsscivstrct and match on qso_name and wavelength limits...
  if n_params() ne 4 then begin
     print,'Syntax - sdss_mtchcivstrct, strct1_fil, strct2_fil, mtch_fil, diff_fil'
     print,'                            [log_fil=,/final,/force,dvtol=]'
     return
  endif 
  ;; Defaults
  if not keyword_set(dvtol) then dvtol = 0.d ;  km/s fuzzines
  cinv = 1./3.e5                             ; km^-1 s
  if not keyword_set(log_fil) then $
     log_fil = strmid(mtch_fil,0,strpos(mtch_fil,'.',/reverse_search))+'.log'
  close,1
  openw,1,log_fil

  ;; Read in
  if size(strct1_fil,/type) eq 7 then $
     strct1 = xmrdfits(strct1_fil,1,/silent) $
  else strct1 = strct1_fil
  strct1.qso_name = strtrim(strct1.qso_name,2)

  if size(strct2_fil,/type) eq 7 then $
     strct2 = xmrdfits(strct2_fil,1,/silent) $
  else strct2 = strct2_fil
  strct2.qso_name = strtrim(strct2.qso_name,2)

  ;; Plates
  plate1 = strmid(strct1.qso_name,strpos(strct1[0].qso_name,'-')+1,4)
  strct1 = strct1[sort(plate1)]
  plate2 = strmid(strct2.qso_name,strpos(strct2[0].qso_name,'-')+1,4)
  strct2 = strct2[sort(plate2)]

  nstrct1 = (size(strct1,/dim))[0]
  nstrct2 = (size(strct2,/dim))[0]
  if nstrct1 gt nstrct2 then begin
     ;; Swap
     tmp = strct1
     strct1 = strct2
     strct2 = tmp
     print,'sdss_mtchcivstrct: swapping order of input files.'
     printf,1,'sdss_mtchcivstrct: swapping order of input files.'
     if size(srct2_fil,/type) eq 7 then $
        printf,1,'strct1 = ',strct2_fil
     if size(srct1_fil,/type) eq 7 then $
        printf,1,'strct2 = ',strct1_fil
  endif else begin
     if size(srct1_fil,/type) eq 7 then $
        printf,1,'strct1 = ',strct1_fil
     if size(srct2_fil,/type) eq 7 then $
        printf,1,'strct2 = ',strct2_fil
  endelse 
  nstrct1 = (size(strct1,/dim))[0]
  nstrct2 = (size(strct2,/dim))[0]

  tags = tag_names(strct1)
  if keyword_set(final) then begin
     ;; Use final zabs and wvlim
     ztag = where(tags eq 'ZABS_FINAL')
     wvlimtag = where(tags eq 'WVLIM_FINAL')
  endif else begin
     ;; Use original zabs and wvlim
     ztag = where(tags eq 'ZABS_ORIG')
     wvlimtag = where(tags eq 'WVLIM_ORIG')
  endelse 

  ;; Instantiate useful arrays
  mask1 = replicate(-1,nstrct1)
  mask2 = replicate(-1,nstrct2)
  wvobs1 = strct1.wrest[0]*(1+strct1.(ztag)[0])
  wvobs2 = strct2.wrest[0]*(1+strct2.(ztag)[0])

  ;; Match on wavelength bounds
  for ii=0L,nstrct1-1 do begin
     ;; Use wavelength bounds to maximize inclusion
     ;; Allow boundaries to be "fuzzy" if think bounds aren't
     ;; naturally big enough (default: dvtol = 0.)
     mtch = where(strct1[ii].qso_name eq strct2.qso_name and $
                  ((wvobs1[ii] ge strct2.(wvlimtag)[0,0]+(1-dvtol*cinv) and $
                    wvobs1[ii] le strct2.(wvlimtag)[0,1]+(1+dvtol*cinv)) or $
                   (wvobs2 ge strct1[ii].(wvlimtag)[0,0]+(1-dvtol*cinv) and $
                    wvobs2 le strct1[ii].(wvlimtag)[0,1]+(1+dvtol*cinv))) and $
                  mask2 eq -1,nmtch)
     if nmtch gt 1 then begin
        printf,1,strct1[ii].qso_name,wvobs1[ii],$
               format='("Multiple matches to ",a15,2x,"wvobs = ",f7.5)'
        ;; Take closest
        mn = min(wvobs1[ii]-wvobs2[mtch],imn,/abs)
        name = strct2[mtch].qso_name
        name[imn] = '*'+name[imn] ; indidcate
        writecol,log_fil,name,wvobs2[mtch],filnum=1,$
                 fmt='(5x,a15,2x,"wvobs = ",f7.5)'
        ;; Set
        mtch = mtch[imn]
     endif else mtch = mtch[0]  ; which may be -1

     if mtch ne -1 then begin
        mask1[ii] = mtch
        mask2[mtch] = ii
     endif 
  endfor                        ; loop ii=nstrct1

  mtch1 = where(mask1 ge 0,nmtch1,complement=unmtch1,ncomplement=nunmtch1)
  mtch2 = where(mask2 ge 0,nmtch2,complement=unmtch2,ncomplement=nunmtch2)
  if keyword_set(force) and (nunmtch1 ne 0 or nunmtch2 ne 0) then begin
     ;; Secondary forced match by closest in sightline
     ;; First status
     printf,1,''
     printf,1,'sdss_mtchcivstrct: initial matches ',nmtch1
     printf,1,'sdss_mtchcivstrct: initial non-matches ',nunmtch1,nunmtch2

     if nunmtch1 ne 0 then begin
        printf,1,''
        printf,1,'Forced pairing between unmatched strct1 and all strct2.'
        frc_mask1 = replicate(' ',nstrct1)
        for ii=0L,nunmtch1-1 do begin
           los = where(strct1[unmtch1[ii]].qso_name eq strct2.qso_name)
           if los[0] eq -1 then begin
              printf,1,'sdss_mtchcivstrct: LOS not in strct2: ',$
                    strct1[unmtch1[ii]].qso_name

              if not keyword_set(qso1notin2) then $
                 qso1notin2 = strct1[unmtch1[ii]].qso_name $
              else qso1notin2 = [qso1notin2,strct1[unmtch1[ii]].qso_name]
              continue
           endif 
           mn = min(wvobs2[los]-wvobs1[unmtch1[ii]],imn,/abs)
           dv = 3.e5*mn/strct1[unmtch1[ii]].(ztag)[0]
           if abs(dv) lt 1.e3 then begin
              frc_mask1[unmtch1[ii]] = frc_mask1[unmtch1[ii]] + ' ' + $
                                       string(los[imn])
              printf,1,strct2[los[imn]].qso_name,wvobs2[los[imn]],$
                 round(dv),format='(5x,a15,2x,"wvobs2 = ",f7.5,2x,"dv21 = ",i6)'
           endif 
        endfor                  ; loop ii=nstrct1
        gd = where(frc_mask1 ne ' ')
        if gd[0] ne -1 then $
           writecol,log_fil,gd,' '+strct1[gd].qso_name,frc_mask1[gd],filnum=1
     endif                      ; nunmtch != 0

     if nunmtch2 ne 0 then begin
        printf,1,''
        printf,1,'Forced pairing between unmatched strct2 and all strct1.'
        frc_mask2 = replicate(' ',nstrct2)
        for ii=0L,nunmtch2-1 do begin
           los = where(strct2[unmtch2[ii]].qso_name eq strct1.qso_name)
           if los[0] eq -1 then begin
              printf,1,'sdss_mtchcivstrct: LOS not in strct1: ',$
                    strct2[unmtch2[ii]].qso_name
              if not keyword_set(qso2notin1) then $
                 qso2notin1 = strct2[unmtch2[ii]].qso_name $
              else qso2notin1 = [qso2notin1,strct2[unmtch2[ii]].qso_name]
              continue
           endif 
           mn = min(wvobs1[los]-wvobs2[unmtch2[ii]],imn,/abs)
           dv = 3.e5*mn/strct2[unmtch2[ii]].(ztag)[0]
           if abs(dv) lt 1.e3 then begin
              frc_mask2[unmtch2[ii]] = frc_mask2[unmtch2[ii]] + ' ' + $
                                       string(los[imn])
              printf,1,strct1[los[imn]].qso_name,wvobs1[los[imn]],$
                    round(dv),format='(5x,a15,2x,"wvobs1 = ",f7.5,2x,"dv12 = ",i6)'
           endif 
        endfor                  ; loop ii=nstrct1
        gd = where(frc_mask2 ne ' ')
        if gd[0] ne -1 then $
           writecol,log_fil,gd,' '+strct2[gd].qso_name,frc_mask2[gd],filnum=1
     endif                      ; nunmtch != 0

     ;; New conditionals
     mtch1 = where(mask1 ge 0 or frc_mask1 ne ' ',nmtch1,$
                   complement=unmtch1,ncomplement=nunmtch1)
     mtch2 = where(mask2 ge 0 or frc_mask2 ne ' ',nmtch2,$
                   complement=unmtch2,ncomplement=nunmtch2)

  endif                         ; /force
  
  
  ;; Output
  printf,1,''
  if nmtch1 ne nmtch2 and not keyword_set(force) then $
     stop,'sdss_mtchcivstrct: matched indices not aligned.'

  if nmtch1 ne 0 then begin
     mwrfits,strct1[mtch1],mtch_fil,/create,/silent ; ext = 1
     ;; Keep it in the same order
     mwrfits,strct2[mask1[mtch1]],mtch_fil,/silent ; ext = 2
     spawn,'gzip -f '+mtch_fil
     print,'sdss_mtchcivstrct: created 2-extensions of ',mtch_fil,nmtch1
     printf,1,'sdss_mtchcivstrct: created 2-extensions of ',mtch_fil,nmtch1
  endif
  if nunmtch1 ne 0 then begin
     ;; Sort
     tmp = sdss_srtcivstrct(strct1[unmtch1])
     mwrfits,tmp,diff_fil,/create,/silent ; ext = 1
     diff = tmp
  endif else mwrfits,[0],diff_fil,/create,/silent ; ext = 1
  if nunmtch2 ne 0 then begin
     tmp = sdss_srtcivstrct(strct2[unmtch2])
     mwrfits,tmp,diff_fil,/silent            ; ext = 2
     diff2 = tmp
     
     ;; Other info
     dvabs2 = (diff2.zabs_orig[1]-diff2.zabs_orig[0])/(1+diff2.zabs_orig[0])*3.e5
     gd = where(abs(dvabs2) le 150.,ngd)
     print,'sdss_mtchcivstrct: diff2 with dvabs <= 150 km/s',ngd
     printf,1,'sdss_mtchcivstrct: diff2 with dvabs <= 150 km/s',ngd
  endif else mwrfits,[0],diff_fil,/silent ; ext = 2
  spawn,'gzip -f '+diff_fil
  print,'sdss_mtchcivstrct: created 2-extensions of ',$
        diff_fil,nunmtch1,nunmtch2
  printf,1,'sdss_mtchcivstrct: created 2-extensions of ',$
         diff_fil,nunmtch1,nunmtch2

  save,/all,filename=strmid(mtch_fil,0,strpos(log_fil,'.',/reverse_search))+'.sav'

  close,1
  print,'sdss_mtchcivstrct: created ',log_fil
  

end                             ; sdss_mtchcivstrct



pro sdss_prntcandsumm, cand_fil, final=final, outfil=outfil, logarr=logarr, $
                       by_rtg=by_rtg, silent=silent
  ;; Print # of candidates, # LOS, z and EW median and range (for *_orig or
  ;; *_final) and division by rating
  if n_params() ne 1 then begin
     print,'Syntax - sdss_prntcandsumm, cand_fil, [/final, outfil=, logarr=, /by_rtg, /silent'
     return
  endif 

  logarr = strarr(100)
  ii = 0

  if size(cand_fil,/type) eq 7 then cand = xmrdfits(cand_fil,1,/silent) $
  else cand = cand_fil
  ncand = (size(cand,/dim))[0]
  getfnam,cand[0].wrest[0],fval,ion
  prs = strsplit(ion,/extract)
  ion = prs[0]

  tags = tag_names(cand)
  if keyword_set(final) then begin
     ztag = where(tags eq 'ZABS_FINAL')
     ewtag = where(tags eq 'EW_FINAL')
  endif else begin
     ztag = where(tags eq 'ZABS_ORIG')
     ewtag = where(tags eq 'EW_ORIG')
  endelse 
  
  unq = uniq(cand.qso_name,sort(cand.qso_name))
  nqso = (size(unq,/dim))[0]

  logarr[ii] = '# '+ion+' candidates = '+string(ncand,format='(i6)')+$
               '; # LOS = '+string(nqso,format='(i5)')
  ii++

  logarr[ii] = 'zabs median, min, max: '+$
               string(median(cand.(ztag)[0],/even),$
                      min(cand.(ztag)[0],max=mx),mx,$
                      format='(3(f7.5,2x))')
  ii++

  logarr[ii] = 'EW median, min, max: '+$
               string(median(cand.(ewtag)[0],/even),$
                      min(cand.(ewtag)[0],max=mx),mx,$
                      format='(3(f8.3,2x))')
  ii++  

  if keyword_set(by_rtg) then begin
     ;; Will this program recursively
     unq = uniq(cand.rating[0],sort(cand.rating[0]))
     nrtg = (size(unq,/dim))[0]
     for rr=0,nrtg-1 do begin
        logarr[ii] = ion+' with rating = '+strtrim(cand[unq[rr]].rating[0],2)
        ii++
        sub = where(cand.rating[0] eq cand[unq[rr]].rating[0])
        sdss_prntcandsumm,cand[sub],final=final,logarr=subarr,/silent
        gd = where(subarr ne '',ngd)
        logarr[ii:(ii+ngd-1)] = subarr[gd]
        ii = ii+ngd
     endfor                     ; loop rr=nrtg
  endif                         ; /by_rtg


  if keyword_set(outfile) then begin
     openw,1,outfile,append=ss
     printf,1,logarr[0:ii-1],format='(a)'
     close,1
     print,'sdss_prntcandsumm: created ',outfile
  endif else if not keyword_set(silent) then $
     print,logarr[0:ii-1],format='(a)' 

end                             ;  sdss_prntcandsumm


function sdss_mtchcandlist, cand_fil, list_fil, silent=silent, outfil=outfil
  ;; Given a candidate structure and a list of sightlines, find the
  ;; candidates in the list
  if n_params() ne 2 then begin
     print,'Syntax - sdss_mtchcandlist(cand_fil, list_fil, [/silent,outfil=])'
     return,-1
  endif 

  if size(cand_fil,/type) eq 7 then cand = xmrdfits(cand_fil,1,/silent) $
  else cand = cand_fil
  cand.qso_name = strtrim(cand.qso_name,2)
  cand = cand[sort(cand.qso_name)]

  ncand = (size(cand,/dim))[0]
  eqsowcand = uniq(cand.qso_name) ; ending index
  nqsowcand = (size(eqsowcand,/dim))[0]
  bqsowcand = [0,eqsowcand[0:nqsowcand-2]+1]; beginning index


  readcol,list_fil,spec_fil,skip=1,/silent,format='a'
  nspec = (size(spec_fil,/dim))[0]
  dum = sdss_getname(spec_fil,/spec,root=qsoinlist)

  nloop = nqsowcand < nspec     ; minimal looping
  for qq=0L,nloop-1 do begin
     ;; Search two ways
     if nqsowcand lt nspec then begin
        mtch = where(cand[eqsowcand[qq]].qso_name eq qsoinlist,nmtch)
        ;; Save candidates
        if nmtch eq 0 then begin
           if not keyword_set(silent) then $
              print,'sdss_mtchcandlist(): QSO w/ candidate not in list: ',$
                    cand[eqsowcand[qq]].qso_name
        endif else begin
           if keyword_set(subcand) then $
              subcand = [subcand,cand[bqsowcand[qq]:eqsowcand[qq]]] $
           else subcand = cand[bqsowcand[qq]:eqsowcand[qq]]
        endelse 
     endif else begin
        mtch = where(cand[eqsowcand].qso_name eq qsoinlist[qq],nmtch)
        if nmtch eq 0 then begin
           if not keyword_set(silent) then $
              print,'sdss_mtchcandlist(): QSO in list_fil has no candidate: ',$
                    qsoinlist[qq]
        endif else begin
           if nmtch ne 1 then $
              stop,'sdss_mtchcandlist(): should not be here...'

           if keyword_set(subcand) then $
              subcand = [subcand,cand[bqsowcand[mtch[0]]:eqsowcand[mtch[0]]]] $
           else subcand = cand[bqsowcand[mtch[0]]:eqsowcand[mtch[0]]]
        endelse 
     endelse 
  endfor                        ; loop qq=nloop
  
  subcand = sdss_srtcivstrct(subcand)
  if keyword_set(outfil) then begin
     mwrfits,subcand,outfil,/create,/silent
     spawn,'gzip -f '+outfil
     if not keyword_set(silent) then $
        print,'sdss_mtchcandlist(): created ',outfil
  endif 

  return,subcand

end                             ; sdss_mtchcandlist()


pro sdss_wrnotefil, civstrct_fil, note_fil, notes=notes, silent=silent
  if n_params() ne 2 then begin
     print,'Syntax - sdss_wrnotefil, civstrct_fil, note_fil, [notes=]'
     return
  endif 

  if size(civstrct_fil,/type) eq 7 then $
     civstr = xmrdfits(civstrct_fil,1,/silent) $
  else civstr = civstrct_fil
  ncivstrct = (size(civstr,/dim))[0]

  if keyword_set(notes) then begin
     ;; Add it on
     num = (size(notes,/dim))[0]
     if num eq 0 then $
        ;; One message for all
        tmp = replicate(notes,ncivstrct) $
     else if num ne ncivstrct then $
        stop,'sdss_wrnotefil: notes must have same size as structure' $
     else tmp = notes

     blank = where(strtrim(civstr.notes,2) eq '',complement=filled)
     if blank[0] ne -1 then civstr[blank].notes = tmp[blank]
     if filled[0] ne -1 then $
        civstr[filled].notes = strtrim(civstr[filled].notes,2)+'; '+tmp[filled]
  endif  
  civstr.notes = strtrim(civstr.notes,2)

  gd = where(civstr.zabs_final[0] gt 0.)
  zabs = civstr.zabs_orig[0]
  if gd[0] ne -1 then zabs[gd] = civstr[gd].zabs_final[0]
  writecol,note_fil,lindgen(ncivstrct),civstr.qso_name,$
           zabs,civstr.rating[0],civstr.notes,$
           fmt='(i7,1x,a16,1x,f7.5,1x,i2,":",2x,a)'
  if not keyword_set(silent) then print,'sdss_wrnotefil: saved ',note_fil

end                             ; sdss_wrnotefil


pro sdss_rdnotefil, note_fil, index, qso_name, zabs, rating, notes, $
                    quick=quick
  if n_params() ne 6 and not keyword_set(quick) then begin
     print,'Syntx - sdss_rdnotefil, note_fil, index, qso_name, zabs, '
     print,'                        rating, notes'
     print,' or     sdss_rdnotefil, note_fil, notes, /quick'
     return 
  endif 

  readcol, note_fil, dum1, notes, format='a,a', delimiter=':', /silent
  nqso = (size(dum1,/dim))[0]

  if keyword_set(quick) then begin
     ;; Just read the notes
     index = notes
     return
  endif 

  tmpnote_fil = 'sdss_rdnotefil_tmp.notes'
  close,1
  openw,1,tmpnote_fil
  printf,1,dum1,format='(a)'
  close,1

  readcol, tmpnote_fil, index, qso_name, zabs, rating, $
           format='i,a,f,i', /silent
  spawn,'\rm '+tmpnote_fil
 
end                             ; sdss_rdnotefil


function sdss_cpstrct, oldstrct_fil, newstrct_tmplt, excl_tag=excl_tag, verbose=verbose
  ;; Copy all elements of one civstrct to another (primarily used to
  ;; add new NOTES tag)
  if n_params() ne 2 then begin
     print,'Syntax - sdss_cpstrct(oldstrct_fil, newstrct_tmplt, [excl_tag=, /verbose])'
     return,-1
  endif 
  
  if size(oldstrct_fil,/type) eq 8 then oldstrct = oldstrct_fil $
  else oldstrct = xmrdfits(oldstrct_fil,1,/silent)
  nstrct = (size(oldstrct,/dim))[0] ; won't work if oldstrct is 1
  oldtags = tag_names(oldstrct[0])
  noldtags = (size(oldtags,/dim))[0]

  if size(newstrct_tmplt,/type) ne 8 then $
        stop,'sdss_cpstct(): newstrct_tmplt must be structure'
  newstrct = replicate(newstrct_tmplt,nstrct)

  newtags = tag_names(newstrct[0])
  nnewtags = (size(newtags,/dim))[0]

  for tt=0L,nnewtags-1 do begin
     mtch = where(newtags[tt] eq oldtags)
     if mtch[0] eq -1 then begin
        if keyword_set(verbose) then $
           print,'sdss_cpstrct: tag '+newtags[tt]+' DNE in old structure'
        continue
     endif else begin
        test = where(newtags[tt] eq strupcase(excl_tag))
        if test[0] eq -1 then begin
           newstrct.(tt) = oldstrct.(mtch[0]) ; actual copy
        endif else if keyword_set(verbose) then $
           print,'sdss_cpstrct: excluding tag '+excl_tag[test[0]]
     endelse  
  endfor                        ; loop tt=ntags

  return, newstrct
end
                    


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro sdss_functions, compile=compile
  ;; Print everything in this package
  ;; if /compile set, then run silent
  if keyword_set(compile) then return

  print,'SDSS_GETSDSSDIR() -- return concatenated $SDSSPATH/$SDSSDIR/'
  print,'SDSS_GETNAME() -- return spSpec-jjjjj-pppp-fff* name by convention.'
  print,'SDSS_GETQSOSTRCT() -- return DR X catalog structure or file name.'
  print,'SDSS_GETQSOLIST() -- return DR X catalog spectra list or file name.'
  print,'SDSS_WRQSOLIST -- write formatted QSO list from structure.'
  print,'SDSS_GETQSOINLIST -- find and save subset of SDSS structure based on list.'
  print,'SDSS_GETCFLG() -- return continuum flag 0 (spline) or 1 (eigen; default).'
  print,'SDSS_GETRATING() -- return convenction for rating tag (-1, 0, 1, 2, 3).'
  print,'SDSS_SRTCIVSTRCT() -- sort elements in sdsscivstrct by QSO name and zabs.'
  print,'SDSS_RMCIVDUPLICATES() -- remove duplicate qso_name, zabs_orig[0], cflg.'
  print,'SDSS_GETRATEDDBLT() -- return sample of doublets with desired rating.'
  print,'SDSS_GETLIMFLG() -- report the analyze, upper, or lower limit flags (default: 1).'
  print,'SDSS_SETLIMFLG() -- combine upper and/or lower limit flags to input.'
  print,'SDSS_GETEWFLG() -- report the EW type flags (default: 1 or orig).'
  print,'SDSS_SETEWFLG() -- combine EW type flags to input.'
  print,'SDSS_GETBALFLG() -- return the BAL flag (default: 4).'
  print,'SDSS_FNDBALQSO() -- return indices of BAL QSOs in input structure.'
  print,'SDSS_INSTANTBALFLG -- instantiate BALFLG tag in sdsscontistrct.'
  print,'SDSS_GETSPECWAVE() -- return SDSS wavelength range [3820, 9200].'
  print,'SDSS_GETSPECPIXSCALE() -- return SDSS standard pixel scale [69 km/s].'
  print,'SDSS_GETRANDQSOSMPL -- select random QSOs from DR X catalog.'
  print,'SDSS_CALCNORMERR() -- calculate the total error including continuum error.'
  print,'SDSS_PLTCONTI -- plot spectrum and eigen or spline continuum.'
  print,'SDSS_PLTSNR -- plot convolved S/N from abslin structure.'
  print,'SDSS_CONTIQUAL() -- meeasure eigen or spline continuum fit quality.'
  print,'SDSS_MEASURESNR() -- measure S/N for input spectra and/or compute dblt bounds.'
  print,'SDSS_MKSNRTAB -- create structure of all S/N and observed wvlim per doublet.'
  print,'SDSS_PRNTSNRTAB -- print out a variety of S/N-based tables/lists.'
  print,'SDSS_REORGANIZE -- move old EIGCONTI/, SPLCONTI/, ABSLIN/ files to new directories.'
  print,'SDSS_SRTBALQSO -- divide input into BAL sightlines and else.'
  print,'SDSS_PRNTDR7SUMM -- print some doublet cuts and statistics.'
  print,'SDSS_CALCPARALLELJOB() -- divide input for given processor.'
  print,'SDSS_CATPARALLELJOB -- concatenate sdss_civsearch-in-parallel outputs.'
  print,'SDSS_MTCHCIVSTRCT -- match elements in two sdsscivstrct files.'
  print,'SDSS_PRNTCANDSUMM -- print basic information from sdsscivstrct.'
  print,'SDSS_RDNOTEFIL -- read notes file from sdss_chkciv.'
  print,'SDSS_WRNOTEFIL -- write notes file.'
  print,'SDSS_CPSTRCT() -- copy one structure (array) to another.'
  return
end 
