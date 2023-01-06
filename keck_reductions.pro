FUNCTION match_scattered_sunlight, p, x=x, y=y, err=err, fit=fit
  common sunlight_fit_common, fitindices

  fit = P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]

  return, abs(y - fit)/err
end

FUNCTION scale_fit_sunlight, p, x
  return, P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]
end

pro keck_reductions, part=part

; written on dec. 15, 2022
; reduce europa juno flyby keck hires data taken on sept. 28, 2022
  
  directory = 'C:\DATA\HIRES_20220928\Processed\'

  ; ========================================================== LOAD SPICE ===================================================================
  kernel_directory    = 'C:\SPICE\'
  CSPICE_KTOTAL, 'all', count
  PRINT, STRCOMPRESS('Cleaning ' + STRING(count) + ' lingering kernels out of memory . . .')
  i = 0
  WHILE i LT count DO BEGIN
    CSPICE_KDATA, 0, 'all', file, type, source, handle, found
    CSPICE_UNLOAD, file
    i = i + 1
  ENDWHILE
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'pck00010.tpc')) ; Body rotational states
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'naif0010.tls')) ; Leap seconds
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'de421.bsp'))    ; SPK (ephemeris kernel) for planets MOST CURRENT
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'sat319.bsp'))   ; SPK (ephemeris kernel) for satellites (most applicable data)
  CSPICE_FURNSH, STRCOMPRESS(FILE_SEARCH(kernel_directory, 'jup310.bsp'))   ; SPK (ephemeris kernel) for satellites (most applicable data)
  CSPICE_KTOTAL, 'all', count
  Print, 'Loaded ', strcompress(count), ' new SPICE kernels'
; ====================================================== SPICE KERNELS ALL LOADED =========================================================


; all basic reductions (bias, flats, thar) are done in keck_europa_flyby.pro part=0.

; ------------------------------------------------------- cosmic ray subtraction ----------------------------------------------------------
  if part eq 0 then begin
    
    ;la_cosmic, directory+'Flat_BS.fits', outsuff = "CR", sigclip = 4.5, statsec = imaging_statsec, readn=RDnoise, gain=Gain
    Na_flatBSCR = MRDFITS(directory+'Flat_BSCR.fits', 0, header, /Dscale, /silent )
    Na_flat_spectra = Na_flatBSCR
    s_spectra = size(Na_flat_spectra, /dim)

; Process the spectral frames
   ; window, 0, xs = s_spectra[0], ys = s_spectra[1]
    
    filenames = file_search(directory+'hires0*') ; these are all bias subtracted and flat corrected; need to do cosmic rays
    
    for i = 0, n_elements(filenames)-1 do begin
      filename = strcompress(filenames[i])
      
      Gain = 1.
      Print, filenames[I], (sxpar(header, 'EXPOSURE'))
      new_filename = STRMID(filename, strlen(directory))
      new_filename = STRMID(new_filename, 0, strpos(new_filename,'.fits'))
      
      la_cosmic, directory + new_filename + '.fits', outsuff = "CR", sigclip = 8.0, statsec = statsec, readn=RDnoise, gain=Gain
      
      BSCR = MRDFITS(directory+new_filename+'CR.fits', 0, header, /Dscale, /silent )
      Dummy = BSCR
      
      Xshift       = 0.
      Yshift       = 0.                                                                                   ; checked & looking good for europa
      flat         = Na_flat_spectra
      shifted_flat = INTERPOLATE(flat, findgen(s_spectra[0]) + xshift, findgen(s_spectra[1]) + yshift, Missing = !values.F_NaN, /grid)
      
      ; Flat field the spectral channels
      ;BSCR = BSCR / shifted_flat        ; this line is screwing things up >:(
      
      if filename eq directory+'hires0127.Cleaned.fits' then stop
      ; Inspect
      wset, 0
      tv, bytscl(BSCR)

      ; Calculate the Doppler shift between absorption and emission line centers, i.e. Europa's heliocentric velocity
      UTC_String = sxpar(header, 'DATE')
      cspice_UTC2ET, UTC_string, ET
      cspice_spkezr, 'Europa', ET, 'J2000', 'LT', 'Sun', Europa_Sun_State, ltime
      theta  = cspice_vsep(Europa_Sun_State[0:2], Europa_Sun_State[3:5])
      Dopplershift = cos(theta) * norm(Europa_Sun_State[3:5])                  ; km/s heliocentric velocity: scalar projection of the relative velocity along the line of sight
      sxaddpar, raw_header, 'R_dot', Dopplershift, 'Europa approx heliocentric velocity [km/s]'
      Print, UTC_string, ' Europa approximate heliocentric velocity =', Dopplershift, ' km/s'

      MWRFITS, BSCR, directory+'Cosmic Rays\' + new_filename + '.FF.fits', raw_header, /silent, /create
      
    ENDFOR
    stop
  ENDIF


; ============================================ STRAIGHTENING IMAGE IN Y AND RECTIFYING IN X ===================================================

  if part eq 1 then begin

; --------------------------------------------------------- Straightening in Y ----------------------------------------------------------------

    files     = FILE_SEARCH(directory+'\Cosmic Rays\', "*FF.fits", count=n_Files)
    newimg    = fltarr(s_coords[1] - s_coords[0] + 1, s_coords[3] - s_coords[2] + 1)
    sunimg    = fltarr(s_coords[1] - s_coords[0] + 1, s_coords[3] - s_coords[2] + 1)
    pixelsofsun = 4.
    COADDED = make_array(1024, 480)

    FOR file = 0, 12 DO BEGIN
      array = MRDFITS(files[file], 0, header, /Fscale, /silent)
      array = float(array)
      specchan = array[s_coords[0]:s_coords[1], s_coords[2]:s_coords[3]]  ; Eliminates imaging channel
      sunlightspec = []
      sunmax = []                                                         ; Array that holds sun row for each column
      edges = make_array(1024., 720)
      specchan = [[edges], [specchan], [edges]]
      straightened_spectrum = specchan                                    ; Will write the straightened spectrum into this array

      FOR i = 0, 1024 - 1 DO BEGIN
        column = specchan[i,*]
        suncol = WHERE(column EQ MAX(column),count)                       ; Finds the sunlight spectrum in each column and puts that row into a 1D array
        suncol = suncol[0]
        sunmax = [sunmax, suncol[0]]
      ENDFOR

      newpos = n_elements(specchan[0,*])/4                                ; Shifting the spectrum to match up halfway in ALL Europa files --> coadd
      polyn = poly_fit(findgen(1024), sunmax,2)                           ; Finds 2nd order polynomial of the sunlight spectrum
      polynomial = poly(findgen(1024), polyn)                             ; Puts polynomial into a form polyn[0] + polyn[1]*X + polyn[2]*X^2
      shifted = polynomial - newpos                                       ; finds the difference between actual position and desired position

      FOR i = 0, 1024 - 1 DO BEGIN
        column = specchan[i,*]
        shifted_array = INTERPOLATE(column, findgen(n_elements(column)) $ ; Shifts each column so that it is positioned at the last suncol
          + shifted[i], Missing = !values.F_NaN)
        straightened_spectrum[i,*] = shifted_array                        ; This is the straightened image in the y-direction
      ENDFOR
      stop

      ; ----------------------------------------------------- Rectification Using Moon File -----------------------------------------------------------

      moonfile = Directory + 'Moon_Surface_DriiftScan.fits'
      moon     = float(MRDFITS(moonfile, 0, Moon_header, /Fscale, /silent))
      moonspecchan = moon[s_coords[0]:s_coords[1], s_coords[2]:s_coords[3]]
      moonspecchan = [[edges], [moonspecchan], [edges]]

      ;     straighten the moonspectrum
      straightened_moonspectrum = moonspecchan
      FOR i = 0, 1024 - 1 DO BEGIN
        column = moonspecchan[i,*]
        shifted_array = INTERPOLATE(column, findgen(n_elements(column)) $ ; Shifts each column so that it is positioned at the last suncol
          + shifted[i], Missing = !values.F_NaN)
        straightened_moonspectrum[i,*] = shifted_array                    ; This is the straightened image in the y-direction
      ENDFOR

      D1cen = 390.
      D2cen = 631.
      windowwidth = 15.
      D1 = straightened_moonspectrum[D1cen - windowwidth:D1cen + windowwidth,*]
      D2 = straightened_moonspectrum[D2cen - windowwidth:D2cen + windowwidth,*]
      D1_center = fltarr(n_elements(straightened_moonspectrum[0,*]))
      D2_center = fltarr(n_elements(straightened_moonspectrum[0,*]))

      ;     D1 LINE FIRST
      for i = 0, n_elements(straightened_moonspectrum[0,*]) - 1 do begin
        result = mpfitpeak(findgen(windowwidth*2. + 1), straightened_moonspectrum[D1cen - windowwidth:D1cen + windowwidth,i], a, STATUS = STATUS)
        if STATUS GT 1 then D1_center[i] = A[1] else D1_center[i] = !values.F_nan
      endfor

      y = findgen(n_elements(straightened_moonspectrum[0,*]))
      real = where(finite(D1_center), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D1_center[real], 2)
      D1location = poly(findgen(n_elements(straightened_moonspectrum[0,*])), coeff) + D1cen-windowwidth

      ;     D2 LINE SECOND
      for i = 0, n_elements(straightened_moonspectrum[0,*]) - 1 do begin
        result = mpfitpeak(findgen(windowwidth*2. + 1), straightened_moonspectrum[D2cen - windowwidth:D2cen + windowwidth,i], a, STATUS = STATUS)
        if STATUS GT 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
      endfor

      y = findgen(n_elements(straightened_moonspectrum[0,*]))
      real = where(finite(D2_center), /NULL)
      COEFF = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
      D2location = poly(findgen(n_elements(straightened_moonspectrum[0,*])), coeff) + D2cen-windowwidth

      ;     Take the average deviation from the mean location, then averaged between both D1 & D2 lines
      mean_shift = MEAN([[D1location - mean(D1location)], $                                        ; Seems good enough, more inspection would be a good idea though!
        [D2location - mean(D2location)]], dim=2)

      ; ---------------------------------------------------------- Apply Rectification to Europa ----------------------------------------------------

      rectified_europa = specchan                                                                  ; dummy array, about to overwrite
      FOR i = 0, s_coords[3] - s_coords[2] DO BEGIN
        row = straightened_spectrum[*,i]
        rectified_array = INTERPOLATE(row, findgen(n_elements(row)) + mean_shift[i], Missing = !values.F_NaN)
        rectified_europa[*,i] = rectified_array
      ENDFOR
      TV, bytscl(rectified_europa)
      stop

      ; ----------------------------- SAVE all rectified Europa files into a new folder: 'rips_europa\Rectified' WITH big header! -------------------

      new_filename = STRMID(files[file], strlen(directory+'rips_europa\New_Processed\'))
      new_filename = new_filename.REMOVE(-8)
      MWRFITS, rectified_europa, Directory+'rips_europa\Rectified\' + new_filename + '.RECT.fits', header, /silent, /create

      COADDED = [coadded + rectified_europa]

    ENDFOR
    TV, bytscl(COADDED)
    stop
  endif

  ; ========================================================= SUNLIGHT SUBTRACTION ==========================================================

  if part eq 2 then begin

    rectfiles = Directory + 'rips_europa\Rectified'
    files     = FILE_SEARCH(rectfiles, "*RECT.fits", count=n_Files)
    rectified_europa = make_array(1024, 960)

    FOR file = 0, 12 DO BEGIN
      array = MRDFITS(files[file], 0, header, /Fscale, /silent)
      arrays = float(array)
      rectified_europa = [rectified_europa + arrays]                          ; Co-adds all rectified Europa files together
    ENDFOR

    ; --------------------------------------------- Find the Dispersion & Sunlight vs Exosphere Indicies --------------------------------------
    D2Cen               = 408
    D1Cen               = 650
    windowwidth         = 15.
    pixelsofsun         = 20.
    spec_1D             = total(rectified_europa, 2, /Nan)
    result              = mpfitpeak(findgen(windowwidth*2. + 1), spec_1D[D2cen - windowwidth:D2cen + windowwidth], a, STATUS = STATUS)
    D2_Solar            = D2cen - windowwidth + a[1]
    result              = mpfitpeak(findgen(windowwidth*2. + 1), spec_1D[D1cen - windowwidth:D1cen + windowwidth], a, STATUS = STATUS)
    D1_Solar            = D1cen - windowwidth + a[1]
    dispersion          = (5895.92424 - 5889.95095) / ( D1_Solar - D2_Solar ) ; A/pixel
    Dopplershift        = sxpar(header, 'R_dot')
    Europa_D2_pixel     = float(D2_Solar - Dopplershift / (dispersion*cspice_clight() / 5889.95095))  ; converted from km/s to pixel units
    Europa_D1_pixel     = float(D1_Solar - Dopplershift / (dispersion*cspice_clight() / 5895.92424))  ; converted from km/s to pixel units

    ;     Get the pixel indices where the spectrum consists of just scattered sunlight
    fitindices = where( (abs(findgen(1024) - Europa_D1_pixel) gt 5) and $
      (abs(findgen(1024) - Europa_D2_pixel) gt 5), /null)                         ; Excludes the sodium emission from Europa

    fitindices = fitindices[where((fitindices gt 300) and (fitindices lt 750))]                       ; omit the edges where the data are noisier

    ; -------------------------------------------------------- Sunlight Subtraction -----------------------------------------------------------
    sunmax    = []                                                         ; Array that holds sun row for each column
    newimg    = fltarr(s_coords[1] - s_coords[0] + 1, s_coords[3] - s_coords[2] + 1)
    sunimg    = fltarr(s_coords[1] - s_coords[0] + 1, s_coords[3] - s_coords[2] + 1)

    FOR i = 0, 1024 - 1 DO BEGIN
      column = rectified_europa[i,*]
      suncol = WHERE(column EQ MAX(column),count)                          ; Finds the sunlight spectrum in each column and puts that row into a 1D array
      suncol = suncol[0]
      sunmax = [sunmax, suncol[0]]
    ENDFOR

    newpos      = n_elements(rectified_europa[0,*])/4                      ; Shifting the spectrum to match up with the mean sunlight spectrum location
    sunlight    = rectified_europa[*, newpos - pixelsofsun : newpos + pixelsofsun]
    TV, bytscl(sunlight)
    sunlight_1d = TOTAL(sunlight, 2, /Nan)
    P_returned  = fltarr(4,s_coords[3] - s_coords[2] + 1)                  ; Three coefficients + MPFIT's "Status"
    P_guessed   = Fltarr(3,s_coords[3] - s_coords[2] + 1)                  ; Initial Guess that we throw at MPFIT

    FOR i = 0, s_coords[3] - s_coords[2] DO BEGIN

      ; generate an initial guess for multipliciative scaling
      row         = rectified_europa[*,i]
      guess_scale = median(row[fitindices] / sunlight_1D[fitindices])
      row_err     = sqrt(abs(row))

      ; Fit a y = A*Gauss_smooth(x,C) + B function to the spectrum, where x is the reference solar spectrum
      p0 = [guess_scale, 0.0, 0.5]                                           ; Guess at initial coefficients
      parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p0))
      parinfo.value         = p0
      ;parinfo[1].fixed      = 1
      ;parinfo[2].fixed      = 1
      parinfo[2].limited    = [1, 1]
      parinfo[2].limits     = [0.0, 20.]

      WEIGHTS = 1./(abs(findgen(1024) - Europa_D2_pixel))^.4 + 1./(abs(findgen(1024) - Europa_D1_pixel))^.4

      fa = {x:sunlight_1d[fitindices], y:row[fitindices], err:1./weights[fitindices]}
      p = mpfit('match_scattered_sunlight', p0, PERROR = err_P, functargs=fa, status=status, parinfo=parinfo)
      P_guessed[*,i]  = p0
      p_returned[*,i] = [p, status]
      scaled_sunlight = scale_fit_sunlight(P, sunlight_1d)             ; Puts it into y = A*shift(Gauss_smooth(x,D),C) + B form

      sunimg[*,i] = scaled_sunlight
      sub         = guess_scale * sunlight_1d                          ; If you JUST want the multiplicative correction (no scattered sunlight accounted for w mpfit)
      totsubtrd   = row  - scaled_sunlight                             ; Change back to row - sub to get just the multiplicative factor
      newimg[*,i] = totsubtrd

    ENDFOR

    ; Window 1 shows the rectified lines (should be straight, and second plot should be ~ 0)
    ; Window 2 shows 2d rectified rebinned image
    ; Window 3 shows the bias/dark/flat/CR corrected, shifted, rectified, sunlight-subtracted 1D light from Europa

    ;  window, 1
    ;    !P.Multi = [0, 1, 2]
    ;    cgplot, findgen(n_elements(newimg[0,*])), D1_center + D1cen-windowwidth, YRANGE=[350.,700.], XRANGE=[0,544], title='Europa D1 and D2 Line Curvature', $
    ;      xtitle='Row (Pixel)', ytitle='Centroid of Absorption Line'
    ;    cgplot, findgen(n_elements(newimg[0,*])), D2_center + D2cen-windowwidth, color='black', /overplot
    ;    cgplot, findgen(n_elements(newimg[0,*])), D2_center - D1_center, /ynozero

    window, 2, xs=1500, ys = 600
    loadct, 0
    p = cglayout([1,3], xgap = 0, ygap = 2)
    p[2,*] = .8
    cgimage, rebin(rectified_europa, 1024, 120), /keep_aspect, /axes, pos = p[*,0]
    cgColorbar, /vertical, minrange = 0, maxrange = max(rebin(rectified_europa, 1024, 60)), pos = [P[2,0]+.05,P[1,0],P[2,0]+.12,P[3,0]], title = 'Data'

    cgimage, rebin(sunimg, 1024, 120), /keep_aspect, /axes, pos = p[*,1], /noerase
    cgColorbar, /vertical, minrange = 0, maxrange = max(rebin(sunimg, 1024, 60)), pos = [P[2,1]+.05,P[1,1],P[2,1]+.12,P[3,1]], title = 'Fitted Reflectance'

    loadct, 3
    cgimage, rebin(newimg, 1024, 120), minvalue = -5., maxvalue = 10., /keep_aspect, /axes, /noerase, pos = p[*,2]
    cgColorbar, /vertical, minrange = -5, maxrange = 10., pos = [P[2,2]+.05,P[1,2],P[2,2]+.12,P[3,2]], title = 'Na Residual'

    window, 3, xs = 1500, ys = 600
    pos = cglayout([2,2], xgap = 0, ygap = 0, iymargin = 0)
    pos[Where(pos eq 0.53356647)] = .7
    cgplot, fitindices, TOTAL(rectified_europa[fitindices,*], 2, /NaN), $
      title = 'COADDED', ytitle="Total Counts", $
      /xstyle, pos = pos[*,0], xtickformat = '(A1)', color = 'yellow', psym=14
    cgplot, findgen(1024), TOTAL(rectified_europa, 2, /NaN), /overplot
    cgplot, findgen(1024), TOTAL(sunimg, 2, /NaN), /overplot, color ='red'

    cgplot, [D1_Solar, D1_Solar], [-1.e6, 1.e6], linestyle = 2, color = 'blue', /overplot
    cgplot, [D2_Solar, D2_Solar], [-1.e6, 1.e6], linestyle = 2, color = 'blue', /overplot
    cgplot, findgen(1024), TOTAL(newimg, 2, /NaN), xtitle = "Column (pixels)", ytitle="Total Counts", xr=minmax(fitindices), pos = pos[*,2], /noerase
    cgplot, [Europa_D1_pixel, Europa_D1_pixel], [-1.e5, 1.e5], linestyle = 2, color = 'red', /overplot
    cgplot, [Europa_D2_pixel, Europa_D2_pixel], [-1.e5, 1.e5], linestyle = 2, color = 'red', /overplot
    cgtext, .4, 0.4, 'Reflectance Subtracted Spectrum', align = 0.5, /normal, charsize = 1.5
    cgtext, .35, 0.5, 'Data', color = 'black', /normal
    cgtext, .35, 0.55, 'Fit Region', color = 'yellow', /normal
    cgtext, .35, 0.6, 'Fitted Reflectance', color = 'red', /normal

    cgplot, TOTAL(rectified_europa, 1, /NaN), findgen(N_elements(TOTAL(rectified_europa, 1, /NaN))), /xstyle, pos = pos[*,1],/noerase, ytickformat = '(A1)'
    cgplot, [0, 1.e7], [newpos, newpos], linestyle = 2, color = 'red', /overplot
    cgplot, TOTAL(newimg, 1, /NaN), findgen(N_elements(TOTAL(newimg, 1, /NaN))), /xstyle, pos = pos[*,3],/noerase, ytickformat = '(A1)'
    cgplot, [-1.e6, 1.e6], [newpos, newpos], linestyle = 2, color = 'red', /overplot

    ; -------------------------------------------------------- Getting Angular Size & Pixel # ----------------------------------------------------------

    body     = 'Europa'
    frame    = 'J2000'
    abcorr   = 'LT+S'
    observer = 'Earth'

    obsdate = strcompress(sxpar(header, 'DATE'),/remove_all)                       ; gets dates of all observations
    cspice_str2et, obsdate, et                                                     ; converts dates to ephemeris time

    cspice_bodvrd, body, 'RADII', 3, radii                                         ; Get Body shape constants
    cspice_spkpos, body, et, frame, abcorr, observer, ptarg, ltime
    R_M = 206264.806 * atan(radii[0] / norm(ptarg))                                ; Radius of the body in arcsec
    sxaddpar, header, 'RADII', 2.*R_M, 'Europa angular diameter [arcsec]'          ; Adds ang. diameter to header


    ;    ENDFOR
    stop
  endif






end