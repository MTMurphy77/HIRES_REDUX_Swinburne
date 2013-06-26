pro airtovac,wave                   
;+

;;;; ALM: READ NOTE AFTER THE CODE

; NAME:
;       AIRTOVAC
; PURPOSE:
;       Convert air wavelengths to vacuum wavelengths 
; EXPLANATION:
;
;       OLD: Wavelengths are corrected for the index of refraction of
;       air under standard conditions.  Wavelength values below 2000 A
;       will not be altered.  Uses the IAU standard for conversion
;       given in Morton (1991 Ap.J. Suppl. 77, 119)
;
;       CORRECTION (Michael Murphy 21 Apr. 2010): ThAr atlases are
;       (for historical reasons) given as air wavelengths. The
;       original measurements used in most atlases, e.g. the Palmer &
;       Engleman (1983) atlas were in vacuum and have been converted
;       to air wavelengths using the Edlen (1966) formula. Therefore,
;       to return an air wavelength back to the correct vacuum
;       wavelength, the inverse of the Edlen (1996) should be
;       used. However, Morton (1991) uses the Edlen (1953) formula and
;       that's what was use previously in this routine. Also, the
;       formula was not inverted, probably because this requires an
;       iterative solution whereby one calculates the refractive index
;       using the air wavelength, determines a first approximation to
;       the vacuum wavelength, then uses that value to re-approximate
;       the refractive index etc. until convergence to some tolerance
;       is reached. This routine now provides vacuum wavelengths to
;       within an absolute wavelength tolerance of 10^-14 Ang. of that
;       expected from the 1966 Edlen formula.
;
; CALLING SEQUENCE:
;       AIRTOVAC, WAVE
;
; INPUT/OUTPUT:
;       WAVE - Wavelength in Angstroms, scalar or vector
;               WAVE should be input as air wavelength(s), it will be
;               returned as vacuum wavelength(s).  WAVE is always converted to
;               double precision upon return.
;
; EXAMPLE:

;       OLD: If the air wavelength is W = 6056.125 (a Krypton line),
;       then AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
;
;       AFTER CORRECTION: If the air wavelength is W = 6056.125 (a
;       Krypton line), then AIRTOVAC, W yields an vacuum wavelength of
;       W = 6057.80185191001328
;
;
; METHOD:
;       OLD: See Morton (Ap. J. Suppl. 77, 119) for the formula used
;       AFTER CORRECTION: See above.
;
; REVISION HISTORY
;       Written W. Landsman                November 1991
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Corrected by Michael Murphy        April 2010
;-
   On_error,2

  if N_params() EQ 0 then begin
      print,'Syntax - AIRTOVAC, WAVE'
      print,'WAVE (Input) is the air wavelength in Angstroms'
      print,'On output WAVE contains the vacuum wavelength in Angstroms'
      return
  endif

; OLD ROUTINE
;  sigma2 = (1d4/double(wave) )^2.              ;Convert to wavenumber squared
;
; Compute conversion factor
;
;  fact = 1.D + 6.4328D-5 + 2.94981D-2/(146.D0 - sigma2) + $
;                            2.5540D-4/( 41.D0 - sigma2)
;    
;  fact = fact*(wave GE 2000.) + 1.0*(wave LT 2000.0)
;
;  wave = wave*fact              ;Convert Wavelength

; CORRECTED ROUTINE
  print,'Using corrected airtovac by M. Murphy'
  tol = 1.D-14
; old_wave = 0.D
  new_wave = double(wave)
; while (ABS(old_wave-new_wave) GE tol) do begin
; while total(abs(old_wave-new_wave) LT tol) LT n_elements(wave) do begin
; print, n_elements(wave)
  for i=0L, (n_elements(wave)-1) do begin
  old_wave = 0.D
    j = 0
    while (ABS(old_wave-new_wave(i)) GE tol) do begin 
      old_wave = new_wave(i)
      sigma2 = 1.D8/(new_wave(i)*new_wave(i))
      refindm1 = 1.D-8*( 8342.13D + 2406030.D/(130.D - sigma2) + $
                        15997.D/(38.9D - sigma2) )
      new_wave(i) = wave(i)*(1.D + refindm1)
      j = j + 1
      ; upper limit of 10 iterations
      if j GT 9 then break
      ; print, new_wave, format='(F20.15)'
    endwhile
  endfor
  wave = new_wave

  return            
  end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ALM: Note that the routine below is the updated version used in the 
;;;; idlutils package. It is likely more accurate than what was there before, 
;;;; as shown in the edits above.
;;;; I'm using the modified version based on Michael's algorithm here, for
;;;; consistency with old xidl results.
;pro airtovac,wave_air, wave_vac                  
;;+
;; NAME:
;;       AIRTOVAC
;; PURPOSE:
;;       Convert air wavelengths to vacuum wavelengths 
;; EXPLANATION:
;;       Wavelengths are corrected for the index of refraction of air under 
;;       standard conditions.  Wavelength values below 2000 A will not be 
;;       altered.  Uses relation of Ciddor (1996).
;;
;; CALLING SEQUENCE:
;;       AIRTOVAC, WAVE_AIR, [ WAVE_VAC]
;;
;; INPUT/OUTPUT:
;;       WAVE_AIR - Wavelength in Angstroms, scalar or vector
;;               If this is the only parameter supplied, it will be updated on
;;               output to contain double precision vacuum wavelength(s). 
;; OPTIONAL OUTPUT:
;;        WAVE_VAC - Vacuum wavelength in Angstroms, same number of elements as
;;                 WAVE_AIR, double precision
;;
;; EXAMPLE:
;;       If the air wavelength is  W = 6056.125 (a Krypton line), then 
;;       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
;;
;; METHOD:
;;	Formula from Ciddor 1996, Applied Optics 62, 958
;;
;; NOTES: 
;;       Take care within 1 A of 2000 A.   Wavelengths below 2000 A *in air* are
;;       not altered.       
;; REVISION HISTORY
;;       Written W. Landsman                November 1991
;;       Use Ciddor (1996) formula for better accuracy in the infrared 
;;           Added optional output vector, W Landsman Mar 2011
;;       Iterate for better precision W.L./D. Schlegel  Mar 2011
;;-
;   On_error,2
;   compile_opt idl2
;
;  if N_params() EQ 0 then begin
;      print,'Syntax - AIRTOVAC, WAVE_AIR, [WAVE_VAC]'
;      print,'WAVE_AIR (Input) is the air wavelength in Angstroms'
;       return
;  endif
;
;    wave_vac = double(wave_air)
;    g = where(wave_vac GE 2000, Ng)     ;Only modify above 2000 A
;    
;    if Ng GT 0 then begin 
; 
;  for iter=0, 1 do begin
;  sigma2 = (1d4/double(wave_vac[g]) )^2.     ;Convert to wavenumber squared
;
;; Compute conversion factor
;  fact = 1.D +  5.792105D-2/(238.0185D0 - sigma2) + $
;                            1.67917D-3/( 57.362D0 - sigma2)
;    
;
;  wave_vac[g] = wave_air[g]*fact              ;Convert Wavelength
;  endfor
;  if N_params() EQ 1 then wave_air = wave_vac
;  endif
;  
;  return            
;  end
