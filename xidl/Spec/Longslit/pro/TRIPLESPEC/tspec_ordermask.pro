FUNCTION TSPEC_ORDERMASK, tset_slits
  nx = tset_slits[0].dims[0]
  ny = tset_slits[0].dims[1]
  yarr = replicate(1.0, nx) # findgen(ny)
  ordermask = long_slits2mask(tset_slits)
  ordermask[WHERE(ordermask GT 0)] = -ordermask[WHERE(ordermask GT 0)] + 8L
  RETURN, ordermask
END
