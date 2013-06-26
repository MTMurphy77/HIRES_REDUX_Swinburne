PRO mage_sensfunc, mage, stdfile, fluxtable,sensfile $
                   ,illumflatfile=illumflatfile
  
  istd=WHERE(mage.FITSFILE EQ fileandpath(stdfile) or $
             mage.FITSFILE+'.gz' EQ fileandpath(stdfile))
  mage_pipe,mage[istd],illumflatfile=illumflatfile, /sensfunc
  tmp = strsplit(fileandpath(stdfile), 'mage', /extract, /regex)
  objfile = 'Object/ObjStr'+strtrim(tmp[0])
  test = file_search(objfile,count=ntest)
  if ntest eq 0 then $
     objfile = strmid(objfile,0,strpos(objfile,'.gz'))
  obj_strct=xmrdfits(objfile,1)
  mage_fitstd, fluxtable, obj_strct, sensfile
  RETURN
END
