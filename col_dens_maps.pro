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
    
;    guider_frame   = mrdfits(guiders[i], 0, guider_header, /fscale)
;    window, 4, xs=512, ys=512, title=strmid(guiders[i], 50)
;    cgimage, guider_frame, minv=0.75*mean(guider_frame), maxv=1.5*mean(guider_frame)
;    stop
    
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
  NSslits         = fltarr(7, 162) + !Values.F_Nan
  EWslits         = fltarr(162, 7) + !Values.F_Nan 
  
  FOR h = 165, 195-1 DO BEGIN ;   N_elements(start_times_et)-1 DO BEGIN   *****165 is where gg475 filt begins!                                                                ; files 116 - 164 are of europa with the sodium filter in
     cspice_ET2UTC, start_times_et[h], "ISOC", 2, start_times
     print, spectra[h] + ' ' + sxpar(headfits(spectra[h]), 'DATE-OBS') + ' START: '+  STRMID(start_times, 11) + '      END: '+ sxpar(headfits(spectra[h]), 'UTC')
     print, 'has MAGIQ files..............................'
     
     FOR k = 200, N_elements(guider_times)-1 DO BEGIN
      IF guider_times[k] LT spectra_times[h] AND guider_times[k] GT START_times_et[h] then begin
        print, strmid(guiders[k],50)+ ' ' + sxpar(headfits(guiders[k]), 'DATE-OBS') + ' '+ sxpar(headfits(guiders[k]), 'UTC')+'   ', k
;        print, k
        
        magiq_per_frame[*,*,k] = mrdfits(guiders[k], 0, guider_header, /fscale)                                                           ; saves all the magiq files per spectral observation
;        if strmid(guiders[k], 59, 1) eq '_' then magiq_per_frame[*,*,k] =    rotate(magiq_per_frame[*,*,k], 1)
;        if strmid(guiders[k], 50, 12) eq 'hiresSlit000' then magiq_per_frame[*,*,k] = REVERSE(mrdfits(guiders[k], 0, guider_header, /fscale), 2)
        print, round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.)
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 344 then begin          ; NS Oriented
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 1)
        endif
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 244 then begin          ; EW Oriented
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 2)
        endif
        
        window, 0, xs=512, ys=512, title=STRMID(spectra[h], 52, 9)+' --> '+strmid(guiders[k], 50, 17)
        cgimage, magiq_per_frame[*,*,k], minv=0.75*median(magiq_per_frame[*,*,k]), maxv=1.5*median(magiq_per_frame[*,*,k])
        ;wait, 1
        
; let's do centroids now
        
        maxes = []
        ylocs = []

        for cols = 0, s[2]-1 do begin
          column = magiq_per_frame[*,cols,k]
          maxes = [maxes, max(column)]
        endfor
        
        centroid_loc  = max(maxes, yloc)
        xloc = WHERE(magiq_per_frame[*,yloc,k] eq max(magiq_per_frame[*,yloc,k]))
        
        CNTRD, magiq_per_frame[*,*,k], xloc, yloc, xcen, ycen, 100
        
        if xcen[-1] eq -1. or ycen[-1] eq -1. then begin
          magiq_per_frame[*,*,k] = fltarr(s[1], s[2])
          continue
        endif
        
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 244 then magiq_per_frame[158:320,267:274,k] = !values.F_Nan
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 344 then magiq_per_frame[267:274,158:320,k] = !values.F_Nan
        
        shift_in_x = s[1]/2. - xcen
        shift_in_y = s[2]/2. - ycen
        
        magiq_per_frame[*,*,k] = shift(magiq_per_frame[*,*,k], shift_in_x, shift_in_y)
        
        window, 3, xs=512, ys=512, title='centroid'
        cgimage, magiq_per_frame[*,*,k], minv=0.75*median(magiq_per_frame[*,*,k]), maxv=1.5*median(magiq_per_frame[*,*,k])
        
      ENDIF
     ENDFOR
   goodframes = []
     layered = TOTAL(magiq_per_frame, 3)
     window, 1, xs=512, ys=512, title='TOTALED MAGIQ FRAMES FOR THIS SPECTRUM'
     cgimage, layered, minv=0.75*median(layered), maxv=1.5*median(layered)
     print, '    '
  ENDFOR
  
  window, 2, xs=512, ys=512, title='TOTALED MAGIQ FRAMES FOR THIS SPECTRUM'
  cgimage, layered, minv=0.75*median(layered), maxv=1.5*median(layered)
  
  
  stop







end