pro col_dens_maps, filt

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

  START_times_et                    = spectra_times - exp_times
  magiq_per_frame                   = fltarr(s[1], s[2], n_elements(guiders))
  
  NSslit_locations                  = fltarr(s[1], s[2])
  EWslit_locations                  = fltarr(s[1], s[2])
  NSslit_locations[267:274,158:320] = !Values.F_Nan
  EWslit_locations[158:320,267:274] = !Values.F_Nan
  slit_location_Juno                = ROT(NSslit_locations, 316)
  slit_location_344                 = ROT(NSslit_locations, 350)
  
  
  
  
if filt eq 'gg475' then begin
  FOR h = 164, 196-1 DO BEGIN ;   N_elements(start_times_et)-1 DO BEGIN   *****165 is where gg475 filt begins!

;;;  this next part gathers ALL magiq files for EACH spectra (aka taken within the exposure time)

     cspice_ET2UTC, start_times_et[h], "ISOC", 2, start_times
     print, spectra[h] + ' ' + sxpar(headfits(spectra[h]), 'DATE-OBS') + ' START: '+  STRMID(start_times, 11) + '      END: '+ sxpar(headfits(spectra[h]), 'UTC')
     print, 'has MAGIQ files..............................'
     
     FOR k = 0, N_elements(guider_times)-1 DO BEGIN
      IF guider_times[k] LT spectra_times[h] AND guider_times[k] GT START_times_et[h] then begin
        print, strmid(guiders[k],50)+ ' ' + sxpar(headfits(guiders[k]), 'DATE-OBS') + ' '+ sxpar(headfits(guiders[k]), 'UTC')+'   ', k

        magiq_per_frame[*,*,k] = mrdfits(guiders[k], 0, guider_header, /fscale)                                                           ; saves all the magiq files per spectral observation
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
        
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 244 then magiq_per_frame[*,*,k] = EWslit_locations   + magiq_per_frame[*,*,k] 
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 344 then magiq_per_frame[*,*,k] = slit_location_344  + magiq_per_frame[*,*,k] 
        
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
  
  cgPS_Open, filename = dir+'\Figures\gg475_guider_map.eps', $
    /ENCAPSULATED, xsize = 10, ysize = 10
    !P.font=1
    loadct, 3
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
  
    title = 'HIRES 2022-09-29 Slit Orientations Around Europa'
  
    cgimage, layered, minv=0.75*median(layered), maxv=1.5*median(layered)
  
  cgPS_Close
  stop
endif     ; gg475 filter data


  
;;;; this is for Na filtered data  
if filt eq 'Na' then begin  
  FOR h = 126, 165-1 DO BEGIN                                                             ; these are the Na filter files
    cspice_ET2UTC, start_times_et[h], "ISOC", 2, start_times
    print, spectra[h] + ' ' + sxpar(headfits(spectra[h]), 'DATE-OBS') + ' START: '+  STRMID(start_times, 11) + '      END: '+ sxpar(headfits(spectra[h]), 'UTC')
    print, 'has MAGIQ files..............................'
    
    ;if (h GT 136) and (h LT 141) then continue                                           ; comment this line out if you WANT to include juno flyby data

    FOR k = 0, N_elements(guider_times)-1 DO BEGIN
      IF guider_times[k] LT spectra_times[h] AND guider_times[k] GT START_times_et[h] then begin
        print, strmid(guiders[k],50)+ ' ' + sxpar(headfits(guiders[k]), 'DATE-OBS') + ' '+ sxpar(headfits(guiders[k]), 'UTC')+'   ', k
        
        if k eq 158 then continue                                                         ; idk what this frame is but it's not 10 S like the log says
        
        magiq_per_frame[*,*,k] = mrdfits(guiders[k], 0, guider_header, /fscale)           ; saves all the magiq files per spectral observation
        print, round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.)
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 334 then begin          ; NS Oriented
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 1)
        endif
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 244 then begin          ; EW Oriented
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 2)
        endif
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 245 then begin          ; EW Oriented
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 2)
        endif
;        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq  44 then begin          ; juno Oriented
;          magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
;          magiq_per_frame[*,*,k] =    ROT(magiq_per_frame[*,*,k], 44, /INTERP)
;        endif

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

        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 245 then magiq_per_frame[*,*,k] = EWslit_locations   + magiq_per_frame[*,*,k]         ; for some EW frames, we used 244.5 and some were 244
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 244 then magiq_per_frame[*,*,k] = EWslit_locations   + magiq_per_frame[*,*,k] 
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq 334 then magiq_per_frame[*,*,k] = NSslit_locations   + magiq_per_frame[*,*,k]
        if round(sxpar(headfits(spectra[h]), 'ROTPOSN') - 90.) eq  44 then magiq_per_frame[*,*,k] = slit_location_Juno + magiq_per_frame[*,*,k]         ; juno flyby byebye

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
  
  cgPS_Open, filename = dir+'\Figures\Na_guider_map_Juno.eps', $
      /ENCAPSULATED, xsize = 10, ysize = 10
    !P.font=1
    loadct, 3
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
  
    title = 'HIRES 2022-09-29 Na filter Slit Orientations Around Europa'
  
    cgimage, layered, minv=0.75*median(layered), maxv=1.5*median(layered)

  cgPS_Close
  
  stop
endif    ; Na filter data



end