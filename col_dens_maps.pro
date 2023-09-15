pro col_dens_maps

;-------------------------------------------Load SPICE Data-----------------------------------------------------

  ; Clean any lingering kernels out of memory here:
  cspice_ktotal, 'all', count
  Print, 'Deleting ', strtrim(string(count),2), ' old SPICE kernels from memory'
  i=0
  while i lt count do begin
    cspice_kdata, 0, 'all', file, type, source, handle, found
    cspice_unload, file
    i = i + 1
  endwhile

  ; Load New Kernels
  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\lsk\naif0010.tls')         ; leap seconds kernel
  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\pck\pck00010.tpc')         ; Planet rotational states
  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\Jupiter_System\jup310.bsp')
  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\planets\de421.bsp')    ; SPK (ephemeris kernel) for planets
  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\generic_kernels\spk\satellites\sat319.bsp'); SPK (ephemeris kernel) for satellites

  cspice_ktotal, 'all', count
  Print, 'Loaded ', strtrim(string(count),2), ' new Spice kernels'
  
;-------------------------------------------SPICE all loaded---------------------------------------------------
  
;; first, want to match HIRES guider frames to spectra. do this using the time stamps. howeverrrr, the HIRES time stamps for the
;; guider frames are WRONG, so I need to use 'DATE' and 'UTC' from those headers and stitch them together.
  
  dir              = 'Z:\DATA\Keck\Europa Na\HIRES_20220928'
  guiders          = FILE_SEARCH(dir+'\MAGIQ files', 'hiresslit*'+'*.fits')
  spectra          = FILE_SEARCH(dir+'\just RAW data', 'hires0*'+'*.fits')
  guider_times     = []
  spectra_times    = []
  exp_times        = []
  
; first, get the time stamps for the guider frames, convert to et
  FOR i = 0, N_elements(guiders)-1 DO BEGIN
    guider_header  = headfits(guiders[i])
    cspice_UTC2ET,   sxpar(guider_header, 'DATE-OBS') + ' '+ sxpar(guider_header, 'UTC'), guider_ET
    guider_times   = [guider_times, guider_ET]
    ;print, guiders[i] + ' ' + sxpar(guider_header, 'DATE-OBS') + ' '+ sxpar(guider_header, 'UTC')
  ENDFOR
  
; now i'm going to unpack the last one just so that i can get image dimensions and put them into an array later
  guider_img = mrdfits(guiders[0], 0, guider_header, /fscale)
  s          = size(guider_img)
  
; now, get time stamps for spectra, convert to et AND get exposure times
  FOR j = 0, N_elements(spectra)-1 DO BEGIN
    spectra_header = headfits(spectra[j])
    cspice_UTC2ET,   sxpar(spectra_header, 'DATE-OBS') + ' ' + sxpar(spectra_header, 'UTC'), spectra_ET
    exposure       = sxpar(spectra_header, 'EXPTIME')
    exp_times      = [exp_times, exposure]
    spectra_times  = [spectra_times, spectra_ET]
    ;print, spectra[j] + ' ' + sxpar(spectra_header, 'DATE-OBS') + ' '+ sxpar(spectra_header, 'UTC')
  ENDFOR
  
;; this next part gathers ALL magiq files for EACH spectra (aka taken within the exposure time). 

  START_times_et        = spectra_times - exp_times
  magiq_per_frame = fltarr(s[1], s[2], n_elements(guiders))
  
  magiq_per_frame = fltarr(s[1], s[2], n_elements(guiders))
  
  FOR h = 165, 195-1 DO BEGIN ;   N_elements(start_times_et)-1 DO BEGIN                                                                   ; files 116 - 164 are of europa with the sodium filter in
     cspice_ET2UTC, start_times_et[h], "ISOC", 2, start_times
     print, spectra[h] + ' ' + sxpar(headfits(spectra[h]), 'DATE-OBS') + ' START: '+  STRMID(start_times, 11) + '      END: '+ sxpar(headfits(spectra[h]), 'UTC')
     print, 'has MAGIQ files..............................'
     
     FOR k = 0, N_elements(guider_times)-1 DO BEGIN
      IF guider_times[k] LT spectra_times[h] AND guider_times[k] GT START_times_et[h] then begin
        print, guiders[k]+ ' ' + sxpar(headfits(guiders[k]), 'DATE-OBS') + ' '+ sxpar(headfits(guiders[k]), 'UTC')+'   '
        print, k
        
        magiq_per_frame[*,*,k] = mrdfits(guiders[k], 0, guider_header, /fscale)                                                           ; saves all the magiq files per spectral observation
        if strmid(guiders[k], 59, 1) eq '_' then magiq_per_frame[*,*,k] =    rotate(magiq_per_frame[*,*,k], 1)
        if strmid(guiders[k], 50, 12) eq 'hiresSlit000' then magiq_per_frame[*,*,k] = REVERSE(mrdfits(guiders[k], 0, guider_header, /fscale), 2)
        
        window, 0, xs=512, ys=512, title=STRMID(spectra[h], 52, 9)+' --> '+strmid(guiders[k], 50, 17)
        cgimage, magiq_per_frame[*,*,k], minv=0.75*median(magiq_per_frame[*,*,k]), maxv=1.5*median(magiq_per_frame[*,*,k])
      ENDIF
     ENDFOR
   goodframes = []
     layered = TOTAL(magiq_per_frame, 3)
     window, 1, xs=512, ys=512, title='TOTALED MAGIQ FRAMES FOR THIS SPECTRUM'
     cgimage, layered, minv=0.75*median(layered), maxv=1.5*median(layered)
     print, '    '
  ENDFOR
  
  
  stop







end