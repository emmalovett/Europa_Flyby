;+
; NAME:
;       Keck_Europa_Flyby
; PURPOSE:
;       Analyze Keck HIRES data od Europa's Na and K emissions during the Sept 28, 2022 Juno Flyby
; EXPLANATION:
;       The output X and Y coordinates are scaled to be between
;       -90 and +90 to go from equator to pole to equator. Output map points
;       can be centered on the north pole or south pole.
;
; CALLING SEQUENCE:
;       Keck_Europa_Flyby, part = part
;
; INPUTS:
;       Part = 0 Bias subtration, cosmic ray
;              1 Flatten, extract orders and get the wavelength solutions
;
; OUTPUTS:
;
; KEYWORDS:
;
; REVISION HISTORY:
;       Skeleton adaptation from prior HIRES code C. Schmidt & E. Lovett December 2022
;-

FUNCTION MX_plus_B, X, P
  ; multiply P[0] and offset P[1]
  return, P[0]*X+ P[1]
end

FUNCTION Gaussian_For_MPFIT, p, x=x, y=y, err=err, fit=fit
  z = (x - p[1])/p[2]
  fitted = p[0]*exp(-1.*z^2); + abs(p[3])
  return, (y - fitted)/err
;  return, fitted
  ;RETURN, P[0] + GAUSS1(X, P[1:2])
end

FUNCTION med_filter, X, P ; Sigma filter in 1-dimension, if an array differs from it's median by P[1] standard deviations replace it with the median
  ; p[0] = width to evaluate the median
  ; p[1] = standard deviation about which to reject
  x_out = x
  med = median(x, P[0], /even)
  replace_pixels = where(abs(med - x) gt p[1]*stddev(x), /null)
  x_out[replace_pixels] = med[replace_pixels]
  return, x_out
end

FUNCTION match_scattered_sunlight, p, x=x, y=y, err=err, fit=fit
  common sunlight_fit_common, fitindices
  fit = P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]
;  fit = P[0]*shift(x,P[3]) + P[1]
;  cgplot, y
;  cgplot, fit, /overplot, color='red'
  return, abs(y - fit)/err
end ; matched with "scale_fit_sunlight" below

FUNCTION scale_fit_sunlight, p, x
  return, P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]
end

PRO Keck_Europa_Flyby_copy, part = part, dir = dir, filt = filt

  case dir of
    'Z:\DATA\Keck\Europa Na\HIRES_20220928': begin
      ;Europa           = '1003517' ; aka C/2017 K2 (PANSTARRS)
      biases          = string(indgen(13)+4, format='(I4.4)')
      Flats           = string(indgen(10)+30, format='(I4.4)')
      Lamps           = string(indgen(5)+40, format='(I4.4)')
      Star_Frames     = string(127, format='(I4.4)')                                   ; use a europa frame bc we have no star frame
      Europa_frames   = string(indgen(37)+127, format='(I4.4)')                        ; these are the frames for the Juno flyby, Na filter only though...
      Europa_frames   = [Europa_frames, string(indgen(31)+165, format='(I4.4)')]      ; this is for K and Na data; for just Na, use string(indgen(38)+127, format='(I4.4)')
      Jupiter_frames  = string(113, format='(I4.4)')                                   ; Jupiter disk center post eclipse
    end
  endcase

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
  ;  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\Small_Bodies\2017 K2.bsp')                 ; SPK (ephemeris kernel) for 2017 K2 PANSTARRS
  ;  CSPICE_FURNSH, STRCOMPRESS('C:\SPICE\Small_Bodies\Didymos.bsp')                 ; SPK (ephemeris kernel) for Didymos

  cspice_ktotal, 'all', count
  Print, 'Loaded ', strtrim(string(count),2), ' new Spice kernels'

  ;-------------------------------------------Load Constants-----------------------------------------------------
  ; Define Rest wavelengths
  Na = [5889.95095, 5895.92424, 8183.256, 8194.824]
  K  = [7664.89913, 7698.96456]
  O  = [5577.330, 6300.304, 6363.776, 7771.944, 7774.166, 7775.388, 8446.25, 8446.36, 8446.76]
  S  = [9212.865, 9228.092, 9237.538]                         ; See Ajello et al. 2008
  SO = [9549.18, 9626.21]                                     ; 0-0 and 1-1 band heads. See Setzer et al. Journal of Molecular Spectroscopy 198, 163–174 (1999), converted to Air wavelength
  C  = [8335.15, 9405.73]                                     ; Worth a check but nothing here
  Cl = [8085.56, 8086.67, 8375.94, 8585.97, 9121.15, 9592.22] ; Worth a check but nothing here
  O2_plus = [5631.9, 6026.4]                                  ; Terrell et al. 2004
  loadct, 3
  
  ; MUST Run part zero if you redfine this
  Europa_Airglow         = [Na, K, O, S, SO, O2_plus]                    ; Get line list for individual emission components
  Undefine, Europa_Airglow_params
  Europa_Airglow_params  = {line:         Europa_Airglow, $
    brightness:     fltarr(N_elements(Europa_Airglow), N_elements(Europa_frames)), $
    err_brightness: fltarr(N_elements(Europa_Airglow), N_elements(Europa_frames)), $
    linewidth:      fltarr(N_elements(Europa_Airglow), N_elements(Europa_frames)), $
    linecenter:     fltarr(N_elements(Europa_Airglow), N_elements(Europa_frames))}

  Airglow_threshold  = 2.e2 ; we don't care much about faint telluric airglow, so threshold it here
  ;  Files = file_search('D:\DATA\Solar and Telluric Spectra\Telluric_airglow\UVES_Sky\*.tfits', count = n_files) ; Load a telluric emission line list
  ;  AG_WL = [] & AG_flux = [] & AG_FWHM = []
  ;  for i = 0, N_files-1 do begin
  ;    UVES    = MRDFITS(files[i], 1, header, /USE_COLNUM )
  ;    AG_WL   = [AG_WL,UVES.c1]
  ;    AG_flux = [AG_flux,UVES.c2]
  ;    AG_FWHM = [AG_FWHM,UVES.c3]
  ;  endfor
  ;  keep   = where(AG_flux gt Airglow_threshold, /Null)
  ;  Telluric_Airglow = AG_WL[Keep] & Telluric_Airglow_Flux = AG_Flux[Keep]

  ;--------------------------------------------------------------------Basic Reductions------------------------------------------------------------------
  if part eq 0 then begin

    big_array = fltarr(2139, 4096, N_elements(biases))

    for i = 0, N_elements(biases)-1 do begin
      big_array[*,*,i] = [mrdfits(Dir+'\hires' + biases[i] + '.fits', 3, header, /fscale), $
                          mrdfits(Dir+'\hires' + biases[i] + '.fits', 2, header, /fscale), $
                          mrdfits(Dir+'\hires' + biases[i] + '.fits', 1, header, /fscale)]
    endfor
    bias = median(big_array, dim = 3, /even)
    window, 2;, xs = 4096, ys = 2200
    cgimage, rotate(transpose(reform(bias)),7), minv=900, maxv=1000
    
;    tester = big_Array[*,*,2]
;    stop

    big_array = fltarr(2139, 4096, N_elements(flats))
    for i = 0, N_elements(flats)-1 do begin
      big_array[*,*,i] = [mrdfits(Dir+'\hires' + flats[i] + '.fits', 3, header, /fscale), $
                          mrdfits(Dir+'\hires' + flats[i] + '.fits', 2, header, /fscale), $
                          mrdfits(Dir+'\hires' + flats[i] + '.fits', 1, header, /fscale)]
    endfor
    junk = mrdfits(Dir+'\hires' + flats[0] + '.fits', 0, header, /fscale)  ; get the proper header
    flat = median(big_array, dim = 3, /even) - bias                  ; WE WILL NORMALIZE & USE THIS FLAT LATER, ORDER BY ORDER
    SXADDPAR, Header, 'BZERO', 0.0
    MWRFITS, rotate(transpose(FLAT),7), Dir+'\Processed\FLAT_BS.fits', header, /create, /silent ;super flat, bias subtracted but ***NOT NORMALIZED***
    flat = flat / mean(flat[where(flat gt 2000.)])                   ; normalize the flat, this seems to work better than order by order flat fielding (2000 is the whole aperture)
;             window, 2;, xs = 4096, ys = 2200
;             cgimage, bytscl(rotate(transpose(reform(flat)),7))
;             tv, bytscl(rotate(transpose(reform(flat)),7), 0, 2)


;tester = big_Array[*,*,2] - bias
;            stop

    big_array = fltarr(2139, 4096, N_elements(lamps))
    for i = 0, N_elements(lamps)-1 do begin
      big_array[*,*,i] = [mrdfits(Dir+'\hires' + lamps[i] + '.fits', 3, header, /fscale), $
                          mrdfits(Dir+'\hires' + lamps[i] + '.fits', 2, header, /fscale), $
                          mrdfits(Dir+'\hires' + lamps[i] + '.fits', 1, header, /fscale)]
    endfor
    junk    = mrdfits(Dir+'\hires' + lamps[0] + '.fits', 0, header, /fscale)
    ThAr    = median(big_array, dim = 3, /even) - bias
    SXADDPAR, Header, 'BZERO', 0.0
    MWRFITS, rotate(transpose(ThAr),7), Dir+'\Processed\ThAr.fits', header, /create, /silent

    Jupiter_array  = fltarr(2139, 4096, n_elements(Jupiter_frames))
    for i = 0, n_elements(Jupiter_frames)-1 do begin
      filename = '\hires' + Jupiter_frames[i] + '.fits'

      Jupiter_array[*,*,i] = [mrdfits(Dir+'\hires' + Jupiter_frames[i] + '.fits', 3, header, /fscale), $
                              mrdfits(Dir+'\hires' + Jupiter_frames[i] + '.fits', 2, header, /fscale), $
                              mrdfits(Dir+'\hires' + Jupiter_frames[i] + '.fits', 1, header, /fscale)]

      junk                 =  mrdfits(Dir+'\hires' + Jupiter_frames[0] + '.fits', 0, header, /fscale)
      new_filename         = STRMID(filename, 0, strpos(filename,'.fits'))

      Jupiter_array[*,*,i] = Jupiter_array[*,*,i] - bias
;      Jupiter_array[*,*,i] = Jupiter_array[*,*,i] / flat & PRINT, 'Do not apply this flat, it needs spectral *and* spatial normalization to unity'
      Jupiter_array[*,*,i] = Jupiter_array[*,*,i] / float(Sxpar(header, 'EXPTIME')) ;normalize to 1 second exposure time
      write_file = rotate(transpose(reform(Jupiter_array[*,*,i])),7)

      SXADDPAR, Header, 'BZERO', 0.0
      MWRFITS, write_file, Dir+'\Processed\' + new_filename + '.Cleaned.fits', header, /create, /silent
    endfor

    Star = [mrdfits(Dir+'\hires' + star_frames + '.fits', 3, header, /fscale), $
            mrdfits(Dir+'\hires' + star_frames + '.fits', 2, header, /fscale), $
            mrdfits(Dir+'\hires' + star_frames + '.fits', 1, header, /fscale)]
    junk =  mrdfits(Dir+'\hires' + star_frames + '.fits', 0, header, /fscale)

    Star = Star - bias ; Don't flat divide, since we'll want to use Star to find the trace in each order, and the misaligned flat would throw off the fit.
    Star = Star / float(Sxpar(header, 'EXPTIME')) ; normalize to 1 second expsosure time
    SXADDPAR, Header, 'BZERO', 0.0
    MWRFITS, rotate(transpose(Star),7), Dir+'\Processed\Star.Trace.fits', header, /create, /silent

    ; -------------------------------------------------------  Reduce Europa frames ----------------------------------------------------

    ;    READCOL,'D:\DATA\___Calibration___\Carl_centrifugal_equator.txt', torus_deg, torus_lat_out, skipline = 1, /Silent
    Europa_array               = fltarr(2139, 4096, n_elements(Europa_frames))
    ET_array                   = dblarr(N_elements(Europa_frames))
    exptime_array              = fltarr(N_elements(Europa_frames))
    torus_lat_array            = fltarr(N_elements(Europa_frames))
    Solar_Well_2_Europa_Dshift = fltarr(N_elements(Europa_frames))
    for i = 2, n_elements(Europa_frames)-1 do begin
      filename = '\hires' + Europa_frames[i] + '.fits'
      Europa_array[*,*,i] = [mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 3, header, /fscale), $
                             mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 2, header, /fscale), $
                             mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 1, header, /fscale)]

      junk                =  mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 0, header, /fscale)
      new_filename        =  STRMID(filename, 0, strpos(filename,'.fits'))
      
    test = rotate(transpose(reform(europa_array[*,*,i])),7)
;    
    window, 0, title='EUROPA BEFORE BIAS SUBTRACTION'
    cgimage, test[*,0:713]
;    cgplot, test[100,0:700], /ynozero
    

;    for i = 7, n_elements(Europa_frames)-1 do begin
      Europa_array[*,*,i] = Europa_array[*,*,i] - bias
;      Europa_array[*,*,i] = Europa_array[*,*,i] / flat & PRINT, 'Do not apply this flat, it needs spectral *and* spatial normalization to unity'
      Europa_array[*,*,i] = Europa_array[*,*,i] / float(Sxpar(header, 'EXPTIME')) ; normalize to 1 second expsosure time
      
      test = rotate(transpose(reform(europa_array[*,*,i])),7)
;      
      window, 1, title='EUROPA AFTER BIAS SUBTRACTION'
      cgimage, test[*,0:713]
;      cgplot, test[100,0:700], /ynozero
      
;      window, 3
;      cgimage, bytscl(test)
;stop
      ; find the instantaneous Earth-Europa Doppler Shift
      cspice_UTC2ET, sxpar(header, 'DATE_BEG'), ET
      ET_mid_exposure = ET + float(sxpar(header, 'EXPTIME'))/2.
      cspice_et2utc, ET_mid_exposure, 'C', 0, utcstr
      cspice_spkezr, 'Europa', ET_mid_exposure, 'J2000', 'LT+S', 'Earth', Europa_Earth_State, ltime
      theta  = cspice_vsep(Europa_Earth_state[0:2], Europa_Earth_state[3:5])
      Europa_wrt_Earth_Dopplershift = cos(theta) * norm(Europa_Earth_State[3:5])

      ; Find the Dopplershift between the solar wells and Europa's emission lines.
      cspice_spkezr, 'Sun', ET_mid_exposure - ltime, 'J2000', 'LT+S', 'Europa', Sun_Europa_State, ltime_up_leg
      theta  = cspice_vsep(Sun_Europa_State[0:2], Sun_Europa_State[3:5])
      Europa_wrt_Sun_Dopplershift = -1. * cos(theta) * norm(Sun_Europa_State[3:5]) ; km / s

      ;      cspice_subpnt, 'Near point: ellipsoid', 'Jupiter', ET_mid_exposure - ltime, 'IAU_Jupiter', 'None', Europa, Sub_Europa, trgepc, srfvec ;get sub solar point in cartesian IAU Jupiter coords
      ;      cspice_bodvrd, 'Jupiter', 'RADII', 3, radii ;get Jupiter shape constants
      ;      re = radii[0]
      ;      rp = radii[2]
      ;      f = (re-rp)/re
      ;      obspos = Sub_Europa - srfvec
      ;cspice_recpgr, 'Jupiter', obspos, re, f, Europa_SysIII, Europa_SysIII_LATITUDE, opgalt
      ;torus_lat_array[i] = interpol(torus_lat_out, reverse(torus_deg), Europa_SysIII*!radeg) ; Europa's latitude in the torus using Phil Phipp's arrays

      ET_array[i]                   = ET_mid_exposure
      exptime_array[i]              = float(sxpar(header, 'EXPTIME'))
      Solar_Well_2_Europa_Dshift[i] = Europa_wrt_Sun_Dopplershift

      ;SXADDPAR, header, 'Sys3_Lon',  Europa_SysIII*!radeg,          ' Sub-Europa System III Longitude'
      ;SXADDPAR, header, 'Sys3_Lat',  Europa_SysIII_LATITUDE*!radeg, ' Sub-Europa System III Latitude'
      ;SXADDPAR, header, 'Torus_Lat', torus_lat_array[i],          ' Torus Latitude WRT the Europa JRM09 Dipole Approx'
      SXADDPAR, header, 'Europa_DOP',  Europa_wrt_Earth_Dopplershift, ' Europa-Earth V_radial in km/s (mid exposure)'
      SXADDPAR, header, 'UTC_Mid', utcstr,                        ' UTC time (mid-exposure)'

      write_file = rotate(transpose(reform(Europa_array[*,*,i])),7)
      SXADDPAR, Header, 'BZERO', 0.0
      MWRFITS, write_file, Dir+'\Processed' + new_filename + '.Cleaned.fits', header, /create, /silent

    endfor
    Europa_Airglow_params = create_struct( 'torus_lat', torus_lat_array,           $
                                           'ET', ET_array,                         $
                                           'ExpTime', exptime_array,               $
                                           'Sun2Euro', Solar_Well_2_Europa_Dshift, $
                                           Europa_Airglow_params                   )
    print, 'XD Angle:', Sxpar(header, 'XDANGL'), ' Echelle Angle:', Sxpar(header, 'ECHANGL')
    save, Europa_Airglow_params, filename = Dir+'\Processed\Europa_Airglow_params.sav'
    stop
  endif ; Part 0 (basic bias/flat/exptime reductions)

  order_38 = {guess_coeffs:[9273.60,.0444120,-7.96656e-007], low_bound:1958,  hi_bound:2071,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_38'}
  order_39 = {guess_coeffs:[9035.56,.0435485,-7.96656e-007], low_bound:1852,  hi_bound:1962,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_39'}
  order_40 = {guess_coeffs:[8809.80,.0425047,-7.96656e-007], low_bound:1753,  hi_bound:1858,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_40'}
  order_41 = {guess_coeffs:[8594.70,.0416007,-7.96656e-007], low_bound:1658,  hi_bound:1764,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_41'}
  order_42 = {guess_coeffs:[8389.69,.0408341,-7.85080e-007], low_bound:1569,  hi_bound:1671,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_42'}
  order_43 = {guess_coeffs:[8194.70,.0398006,-7.54592e-007], low_bound:1483,  hi_bound:1585,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_43'}
  order_44 = {guess_coeffs:[8008.65,.0388470,-7.36656e-007], low_bound:1386,  hi_bound:1502,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_44'}
  order_45 = {guess_coeffs:[7830.85,.0379826,-7.06656e-007], low_bound:1308,  hi_bound:1411,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_45'}
  order_46 = {guess_coeffs:[7660.50,0.0370141,-6.36462e-007], low_bound:1234,  hi_bound:1335,  WL_range:[5576., 5578], aperture_limit:[4,39], name:'order_46'}                ; K order
  order_47 = {guess_coeffs:[7497.70,0.0362730,-6.36462e-007], low_bound:1163,  hi_bound:1262,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_47'}
  order_48 = {guess_coeffs:[7341.23,0.0357843,-6.91336e-007], low_bound:1094,  hi_bound:1193,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_48'}
  order_49 = {guess_coeffs:[7191.38,0.0348490,-6.36462e-007], low_bound:1029,  hi_bound:1125,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_49'}
  order_50 = {guess_coeffs:[7047.80,0.0341660,-6.36462e-007], low_bound:966,  hi_bound:1061,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_50'} ;dodgy 5th order fit, 3rd fine
  order_51 = {guess_coeffs:[6909.42,0.0336272,-6.36462e-007], low_bound:904,  hi_bound:999,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_51'}
  order_52 = {guess_coeffs:[6776.62,0.0330356,-6.45259e-007], low_bound:846,  hi_bound:940,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_52'}
  order_53 = {guess_coeffs:[6648.54,0.0325630,-6.61965e-007], low_bound:789,  hi_bound:885,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_53'}
  order_54 = {guess_coeffs:[6525.79,0.0313552,-4.81000e-007], low_bound:737,  hi_bound:828,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_54'}
  order_55 = {guess_coeffs:[6367.40,0.0310416,-5.14465e-007], low_bound:665,  hi_bound:775,  WL_range:[5576., 5578], aperture_limit:[4,41], name:'order_55'};dodgy!
  order_56 = {guess_coeffs:[6292.60,0.0305972,-5.90291e-007], low_bound:620,  hi_bound:699,  WL_range:[5576., 5578], aperture_limit:[5,40], name:'order_56'}
  order_57 = {guess_coeffs:[6182.21,0.0300263,-5.73491e-007], low_bound:564,  hi_bound:652,  WL_range:[5576., 5578], aperture_limit:[6,41], name:'order_57'}
  order_58 = {guess_coeffs:[6075.60,0.0293140,-5.14465e-007], low_bound:518,  hi_bound:607,  WL_range:[5576., 5578], aperture_limit:[7,42], name:'order_58'}
  order_59 = {guess_coeffs:[5972.71,0.0289892,-5.60534e-007], low_bound:470,  hi_bound:560,  WL_range:[5576., 5578], aperture_limit:[7,42], name:'order_59'}
  order_60 = {guess_coeffs:[5873.10,0.0285014,-5.54225e-007], low_bound:448,  hi_bound:501,  WL_range:[5576., 5578], aperture_limit:[7,42], name:'order_60'}                  ; Na order
;  order_60 = {guess_coeffs:[5878.99,0.0145895,-4.83663e-007], low_bound:448,  hi_bound:501,  WL_range:[5576., 5578], aperture_limit:[7,42], name:'order_60'}
;  order_60 = {guess_coeffs:[5873.26,0.0281397,-3.93779e-007], low_bound:448,  hi_bound:501,  WL_range:[5576., 5578], aperture_limit:[7,42], name:'order_60'}
  order_61 = {guess_coeffs:[5776.98,0.0279140,-5.14465e-007], low_bound:383,  hi_bound:473,  WL_range:[5576., 5578], aperture_limit:[3,42], name:'order_61'}
  order_62 = {guess_coeffs:[5683.90,0.0272618,-4.66680e-007], low_bound:340,  hi_bound:431,  WL_range:[5576., 5578], aperture_limit:[4,42], name:'order_62'}
  order_63 = {guess_coeffs:[5593.67,0.0269984,-5.01777e-007], low_bound:301,  hi_bound:389,  WL_range:[5576., 5578], aperture_limit:[5,42], name:'order_63'}
  order_64 = {guess_coeffs:[5506.12,0.0267284,-5.30678e-007], low_bound:263,  hi_bound:349,  WL_range:[5576., 5578], aperture_limit:[6,45], name:'order_64'}
  order_65 = {guess_coeffs:[5421.34,0.0263753,-5.34238e-007], low_bound:223,  hi_bound:308,  WL_range:[5575., 5580], aperture_limit:[8,43], name:'order_65'}
  order_66 = {guess_coeffs:[5339.50,0.0257361,-4.86579e-007], low_bound:188,  hi_bound:271,  WL_range:[5575., 5580], aperture_limit:[9,42], name:'order_66'}
  order_67 = {guess_coeffs:[5259.43,0.0256491,-5.34282e-007], low_bound:151,  hi_bound:235,  WL_range:[5575., 5580], aperture_limit:[10,40], name:'order_67'}
  order_68 = {guess_coeffs:[5182.29,0.0251111,-4.98049e-007], low_bound:118,  hi_bound:194,  WL_range:[5575., 5580], aperture_limit:[11,39], name:'order_68'}
  order_69 = {guess_coeffs:[5107.44,0.0245935,-4.72139e-007], low_bound:84,  hi_bound:158,  WL_range:[5575., 5580], aperture_limit:[12,38], name:'order_69'}
  order_70 = {guess_coeffs:[5034.47,0.0242993,-4.84337e-007], low_bound:51,  hi_bound:122,  WL_range:[5575., 5580], aperture_limit:[13,38], name:'order_70'}

  ;  orders   = [ order_38, order_39, order_40, order_41, order_42, order_43, order_44, order_45, order_46, $
  ;    order_47, order_48, order_49, order_50, order_51, order_52, order_53, order_54, order_55, order_56, order_57, order_58,$
  ;    order_59, order_60, order_61, order_62, order_63, order_64, order_65, order_66, order_67, order_68, order_69, order_70]
  ;
  ; ONLY reduce on K D, O 6300 and, Na D orders
  orders   = [order_46, order_56, order_60]


  ;x = findgen(22)+45
  ;y = orders.guess_coeffs[1]
  ;
  ;xinterp = findgen(30)+38
  ;
  ;
  ;coeffs = poly_fit(x, y, 2)
  ;yinterp = poly(xinterp, coeffs)
  ;
  ;cgplot, xinterp, yinterp, psym=1, /ynozero
  ;cgplot, x, orders.guess_coeffs[1], psym=2, color ='red', /overplot
  ;print, [transpose(xinterp), transpose(yinterp)]
  ;stop

  if part eq 1 then begin ; Straighten, extract and get the wavelength solutions

    slit_length_pix = 44 ; rough width of the aperture in pixels

    READCOL,'Z:\DATA\___Calibration___\thar_uves.dat', F='X,F', ThAr_WL, STRINGSKIP = '#', skipline = 800, numline = 1600

    Na_cube         = fltarr(4001, 44, 36) ; N_elements(Europa_frames))                                 ; the number 44 is based on the aperture upper and lower bounds
    K_cube          = fltarr(4001, 36, N_elements(Europa_frames))                                 ; the number 26 is based on the aperture upper and lower bounds
    O_cube          = fltarr(4001, 36, N_elements(Europa_frames))                                 ; the number 26 is based on the aperture upper and lower bounds

    Na_Jupiter_cube = fltarr(4001, 44, N_elements(Jupiter_frames))                                ; See how it's expecting more than 1 jupiter frame here? that's so you can compare each flux calibration later
    K_Jupiter_cube  = fltarr(4001, 36, N_elements(Jupiter_frames))
    O_Jupiter_cube  = fltarr(4001, 36, N_elements(Jupiter_frames))

    for h = 2, N_elements(orders) - 1 do begin
      order = orders[h]
      ;if order.name ne 'order_60' then continue

      ; Use a star to find the trace in within the orders of interest
      Star = mrdfits(Dir+'\Processed\Star.Trace.fits', 0, star_header) ; do not use multiple Ganymede frames here, it's important it's fixed.
      ThAr = mrdfits(Dir+'\Processed\ThAr.fits', 0, ThAr_header)
      flat = mrdfits(Dir+'\Processed\FLAT_BS.fits', 0, Flat_header)    ; not yet normalized to unity

      print, 'ECHANGL = ', sxpar(star_header, 'ECHANGL')
      print, ' XDANGL = ', sxpar(star_header, 'XDANGL' )

      ; Guess and check at a linear wavelength solution
      WL          = poly(findgen(4001), order.guess_coeffs)    ; Only roughly accurate
      WL_0        = WL
;      xr          = minmax(WL)
      xr = [5880., 5920.]
      order_lines = ThAr_WL[where( (ThAr_WL gt xr[0]) and (ThAr_WL lt xr[1]))]
      ID          = make_array(N_elements(order_lines), value = 1.e5)

      Frames = 'hires'+strcompress(Europa_frames, /rem)+'.Cleaned.fits'
      
; comparing on- and off-disk measurements. on-disk measurements are MUCH brighter here
;      img1 = mrdfits(Dir+'\Processed\' + Frames[3], 0, header)                            ; tests an EW on disk frame
;      img2 = mrdfits(Dir+'\Processed\' + Frames[6], 0, header)                            ; tests an EW 10 N frame
;      window, 0, title='EW ON DISK : COMPARE MAX VALUE'
;      cgplot, total(img1[500:1000,0:1000], 1, /nan)
;      cgtext, 0.8, 0.9, 'max: '+strcompress(max(img1[500:1000,0:1000])), color='red'
;      window, 1, title='EW ON DISK'
;      cgimage, img1[500:1000,0:1000]
;      window, 2, title='EW 10 N : COMPARE MAX VALUE'
;      cgplot, total(img2[500:1000,0:1000], 1, /nan)
;      cgtext, 0.8, 0.9, 'max: '+strcompress(max(img2[500:1000,0:1000])), color='red'
;      window, 3, title='EW 10 N'
;      cgimage, img2[500:1000,0:1000], minv=0, maxv=0.5
      
      
      if filt eq 'Na' then begin
        labels = strarr(38)
        if filt eq 'Na' then begin
          labels[2:5  ] = 'EW on disk'
          labels[6:7  ] = 'EW10 N'
          labels[8:9  ] = 'EW20 N'
          labels[10:14] = 'Juno flyby'
          labels[15:16] = 'EW10 S'
          labels[17:18] = 'EW20 S'
          labels[19:23] = 'NS on disk'                     ; file 147 is skipped because it's the sunlight spectrum
          labels[24:25] = 'NS10 W'
          labels[26:27] = 'NS20 W'
          labels[28:29] = 'NS10 E'
          labels[30:31] = 'NS20 E'
          labels[32:34] = 'NS on disk'
          labels[35:36] = 'EW on disk'
        endif
      endif
      statuses = []
      for frame = 0, 35 do begin ; files 2-37 are all of the Na filtered spectra! ;26 for worst cr frame
        Europa   = mrdfits(Dir+'\Processed\' + Frames[frame], 0, header)
        WL      = WL_0          ; start with the guessed WL

        ; -------------------------------------- Wavelength Solution in the Orders of Interest ------------------------------------------------------------------
        ;        Take_a_long_hard_look   = mrdfits('D:\DATA\Keck\Europa Na\'+Dir+'\Processed\' + Frames[0], 0, header)
        ;        window, 0, xs = 4096,ys = 2139
        ;        ;window, 0, xs = 4096,ys = 1050
        ;        tv, bytscl(Take_a_long_hard_look, 0, .3 )

        print, 'ECHANGL = ', sxpar(star_header, 'ECHANGL')
        print, 'XDANGL = ', sxpar(star_header, 'XDANGL')
        ;              rdpix, Take_a_long_hard_look
        ;              stop

        ; trace each order's footprint on the chip using the centroid of the Star Spectrum

        subframe = Star[0:4000, order.low_bound:order.hi_bound]
        ;subframe = Europa[0:4000, order.low_bound:order.hi_bound]

        s = size(subframe, /dim)
        trace = fltarr(s[0])
        for i = 0, s[0] - 1 do begin
          junk = max(subframe[i,*], loc, /nan)
          trace[i] = loc

          ;          ; Or, the much slower, but slightly better way to trace...
          ;          junk     = mpfitpeak(findgen(s[1]), subframe[i,*], a, STATUS = STATUS, estimates = [500., loc, 1., 0.], /positive)
          ;          trace[i] = a[1]
        endfor

        trace_old = trace
        replace_left = trace[0:2000]
        replace_right= trace[2001:*]
        replace_left[where(replace_left gt 41.)] = !values.F_nan
        replace_right[where( (replace_right lt 30.) or (replace_right gt 100.) )] = !values.F_nan
        trace[0:2000] = replace_left
        trace[2001:*] = replace_right
        keep = where(finite(trace))
        x    = findgen(s[0])
        keep = keep[where(x[keep] lt 4060., /null)] ;edge effects are rather dodgy so omit that part of the trace
        ; Manual over-ride, when this technique fails, especially at the ccd boundaries.
        ;  if ((date eq 'UT221124') and (order.name eq 'order_55')) then keep = keep[where(keep gt 1150., /null)] ;... for example

        coeffs = poly_fit(x[keep], trace[keep], 4, Yfit = Yerror)          ; trace continuum source's position with a 4th order polynomial

;        window, 1, xs = 1800, ys=800, title = 'WAVELENGTH SOLUTION FOR: ' + ORDER.NAME
;        cgplot, findgen(s[0]), trace_old, psym = 4, /ynozero, /yst, /xst
;        cgplot, findgen(s[0]), trace, psym = 4, /overplot, color = 'red'
;        cgplot, x, POLY( findgen(s[0]), coeffs), /overplot, color = 'blue' ; IF THIS FIT TO THE RED DATA POINTS IS BAD, EVERYTHING DOWNSTREAM FAILS
;        stop
        ; THIS IS A GOOD SPOT TO INSPECT THE TRACE FITTING, MOVING UPPER & LOWER BOUNDS

        ThAr_order_straight = fltarr(s[0], slit_length_pix)
        Flat_order_straight = fltarr(s[0], slit_length_pix)
        Star_order_straight = fltarr(s[0], slit_length_pix)
        Europa_order_straight = fltarr(s[0], slit_length_pix)
        
        for i = 0, s[0] - 1 do begin ; for column of wavelength, shift the spectra up & down to straighten them out.
          interpt_at               = float(order.low_bound) + findgen(slit_length_pix) - float(slit_length_pix)/2. + POLY( float(i), coeffs)

          Star_order_straight[i,*] = interpolate(Star[i,*], interpt_at, cubic = -0.5)
          ThAr_order_straight[i,*] = interpolate(ThAr[i,*], interpt_at, cubic = -0.5)
          flat_order_straight[i,*] = interpolate(flat[i,*], interpt_at, cubic = -0.5)
          Europa_order_straight[i,*] = interpolate(Europa[i,*], interpt_at, cubic = -0.5)
        endfor
        
        
;;;;;;;;;;;;;;;;; i'm running this test to see if on-disk spectra are already brighter than off-disk and they ARE.
;;;;;;;;;;;;;;;;; frame 3 (EW on disk) at 10 R_E west ~250 DN/s
;;;;;;;;;;;;;;;;; frame 24 (NS10W) on disk ~100 DN/s
;;;;;;;;;;;;;;;;; these should match up, but frame 3 at the correct location is already brighter than frame 24.
;;;;;;;;;;;;;;;;; tested this before and after CR sub too

;        plate_scale = 0.358    "/pix
;        suncol = where(Europa_order_straight[700,*] eq max(Europa_order_straight[700,*]))
;        ang_radius = 0.5
;        yr = [-(suncol*plate_scale/ang_radius), ((44.-suncol)*plate_scale/ang_radius)]
;        window, 11, title=labels[frame]
;        cgplot, total(Europa_order_straight[500:1000,*], 1), xr=[0,44], yr=[0,3.e4], xtickformat='(A1)', /ynozero
;        cgaxis, xaxis=0, xtit='Europa Radii', xr=yr, xstyle=1, xticklen=-0.05
;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        
        

        Flat_aperture = Flat_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;nomalize_with = mean(Flat_aperture, dim = 2, /NAN)
        flat_coeffs = poly_fit(findgen(s[0]), mean(Flat_aperture, dim = 2, /NAN), 3)
        nomalize_with = poly(findgen(s[0]), flat_coeffs)
        Flat_aperture = Flat_aperture / rebin(nomalize_with, s[0],order.aperture_limit[1] - order.aperture_limit[0] + 1) ; Normalize the flat field to unity, now its spectral AND spatially normalized
        aperture      = Europa_order_straight;[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;aperture      = aperture / Flat_aperture ; don't flatten twice!!!

        ; Evaluate the guess at the wavelength solution
        ThAr_order_straight = abs(ThAr_order_straight);[where(ThAr_order_straight lt 0., /null)] = !Values.F_Nan        I DON'T KNOW IF THIS IS ALLOWED. HACK HACK HACK. having NaN's was messing up the WL soln
        ThAr_measured = total(ThAr_order_straight,2);[*, order.aperture_limit[0]:order.aperture_limit[1]], 2)
        window, 2
        P = cglayout([1,2], ygap = 0.)
        cgplot, WL, ThAr_measured, /ylog, yr = [2.e3, 5.e5], /xstyle, pos = p[*,0], xtickformat = '(A1)', ytitle = 'ThAr Lamp Counts', xr=xr

        ; find the peaks and do a 3rd order wavelength solution in each order
        identwave_1 = [] & identpixl_1 = []
        search    = 50  ; pixel distance from the expected position based of the fits header to search for an arc line
        expected_pixel = round( interpol( findgen(N_elements(WL)), WL, order_lines ) )
        parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.,0.]}, 4)
        parinfo[2].fixed = 1
        parinfo[2].value = 4.0
        for i = 0, n_elements(order_lines) - 1 do begin
          if ( (expected_pixel[i]-search lt 0) or (expected_pixel[i]+search gt (s[0])-1) ) then continue
          y = ThAr_measured[expected_pixel[i]-search:expected_pixel[i]+search]
          if total(finite(y), /Nan) eq 0. then continue
          result = mpfitpeak(findgen(search*2. + 1), y, a, /POSITIVE, PARINFO = PARINFO, STATUS = STATUS, NFREE = 2, /nan)
          if (a eq !null) then continue
          if (not finite(total(a))) then continue ;else print, a
          if ( (status gt 0) and (a[0] gt 1.e4) and (a[1] gt 0.) and (a[2] gt 2.) and (a[2] lt 6.)) then begin
            identwave_1 = [identwave_1, order_lines[i]]
            identpixl_1 = [identpixl_1, a[1] - search + expected_pixel[i]]
          endif
        endfor
        coeff_1 = ROBUST_POLY_FIT(identpixl_1, identwave_1, 3, yfit_1, SIG)
        WL = poly(findgen(N_elements(WL)), coeff_1)

        ; iterate and try for 5th order
        identwave_2 = [] & identpixl_2 = []
        search    = 10  ; pixel distance from the expected position based of the fits header to search for an arc line
        expected_pixel = round( interpol( findgen(N_elements(WL)), WL, order_lines ) )
        for i = 0, n_elements(order_lines) - 1 do begin
          if ( (expected_pixel[i]-search lt 0) or (expected_pixel[i]+search gt s[0]) ) then continue
          result = mpfitpeak(findgen(search*2. + 1), float(ThAr_measured[expected_pixel[i]-search:expected_pixel[i]+search]), a, /POSITIVE, PARINFO = PARINFO, STATUS = STATUS, NFREE = 2)
          print, a
          statuses = [statuses, status]
          
          if not finite(total(a)) then continue
          if ( (status gt 0) and (a[0] gt 1.e3) and (a[1] gt 0.) and (a[2] gt 1.) and (a[2] lt 6.)) then begin
            identwave_2 = [identwave_2, order_lines[i]]
            identpixl_2 = [identpixl_2, a[1] - search + expected_pixel[i]]
          endif
        endfor
        
        coeff_2 = ROBUST_POLY_FIT(identpixl_2, identwave_2, 5, yfit_2, SIG)

        ;        ; Manual override, when this technique fails
        if order.name eq 'order_46' then coeff_2 = coeff_1 ; the potassium order's wavelength solution
        ;        if order.name eq 'order_44' then coeff_2 = order.GUESS_COEFFS
        ;        if order.name eq 'order_41' then coeff_2 = order.GUESS_COEFFS
        ;        if order.name eq 'order_40' then coeff_2 = order.GUESS_COEFFS
        ;        if order.name eq 'order_39' then coeff_2 = order.GUESS_COEFFS
        if order.name eq 'order_60' then coeff_2 = order.GUESS_COEFFS
        cgplot, order_lines, make_array(N_elements(order_lines), value = 1.e5), /overplot, psym=4                                                           ; plot cataloged lines
        cgplot, interpol(WL_0, findgen(N_elements(WL_0)), identpixl_1), make_array(N_elements(identpixl_1), value = 2.e5), /overplot, psym=4, color = 'red' ; plot which lines were found in iteration 1
        cgplot, interpol(WL_0, findgen(N_elements(WL_0)), identpixl_2), make_array(N_elements(identpixl_2), value = 3.e5), /overplot, psym=4, color = 'blue'; plot which lines were found in iteration 2

        cgplot, identwave_1, identpixl_1, pos = p[*,1], /noerase, XR = XR, /xstyle, psym=5, color = 'red'
        cgplot, identpixl_1, yfit_1, psym = 5, /overplot, color = 'red'
        cgplot, identwave_2, identpixl_2, psym=5, /overplot, color = 'blue'
        cgplot, identpixl_2, yfit_2, psym = 5, /overplot, color = 'blue'
        cgplot, poly(findgen(N_elements(WL)), order.guess_coeffs), findgen(N_elements(WL)), /overplot
        cgplot, poly(findgen(N_elements(WL)), coeff_1), findgen(N_elements(WL)), /overplot, color = 'red'
        cgplot, poly(findgen(N_elements(WL)), coeff_2), findgen(N_elements(WL)), /overplot, color = 'blue'

        Print, 'Fitting Wavelength Solution From'+ strcompress(N_elements(identpixl_2)) + ' ThAr lines.
        Print, 'Guess coefficients for [xstart, dispersion/pixel, 2nd order]', order.guess_coeffs
        Print, 'A better guess would have been:', ROBUST_POLY_FIT(identpixl_2, identwave_2, 2, junk, SIG)
        print, 'Poly fit coefficients:', coeff_2
        WL = poly(findgen(N_elements(WL)), coeff_2)                                                               ; this is the final wavelength solution (5th order)
        window, 3
        cgplot, WL, /ynozero
;        STOP ; <--- INSPECT WL SOLUTION HERE
        
;                             window, 4
;                             cgplot, WL, ThAr_measured, yr = [2.e3, 5.e5], /xstyle, xr=[5890,5900], ytitle = 'ThAr Lamp Counts'
;                             
;                             ; ThAr line at ~5891.5. use this to get HIRES resolution.
;                             
;                             test_WL = 5891.4
;                             range = [625, 675]
;                              indices  = where(WL eq WL_range)
;                             
;                             cgplot, WL[range[0]:range[1]], ThAr_Measured[range[0]:range[1]], /ylog, yr = [2.e3, 5.e5], /xstyle
;                             
;                             p  = [max(ThAr_Measured[range[0]:range[1]])-min(ThAr_Measured[range[0]:range[1]]), test_WL, 0.05, min(ThAr_Measured[range[0]:range[1]])]
;                             quickfunct  = { x:WL[range[0]:range[1]], y:ThAr_Measured[range[0]:range[1]], err:100.*sqrt(abs(ThAr_Measured[range[0]:range[1]])) }
;                             a  = mpfit('Gaussian_for_MPFIT', p, funct=quickfunct, STATUS = Did_it_work);, xtol=5D-9, ftol=1D-6, gtol=1D-9)
;                             
;                             cgplot, WL[range[0]:range[1]], gaussian(WL[range[0]:range[1]], a), /overplot, color='red'
;                             fwhm_val = a[2] * ( 2 * SQRT(2* ALOG(2)) )
;                             
;                             HIRES_res = test_WL / fwhm_val
;                             print, 'HIRES Res = ', HIRES_res
;                             
;                             YES IT WORKS!!!
;                             stop
        
        
        

        ; write the fits files and run the Cosmic Ray Corrections. Some orders overlap the CCD edge in the extraction so exclude these regions in the CR correction
        SXADDPAR, Header, 'BZERO', 0
        SXADDPAR, Header, 'BSCALE', 0
        MWRFITS, aperture, Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_' + Frames[frame], header, /create

        Case 1 of
          order.name eq 'order_43': statsec = '[0:2925,*]'
          order.name eq 'order_52': statsec = '[500:4000,*]'
          order.name eq 'order_53': statsec = '[0:2620,*]'
          order.name eq 'order_60': statsec = '[500:900,*]'
          else: junk = temporary(statsec)
        endcase

        gain = -1.0
        sigclip = 12.5
        la_cosmic, Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_' + Frames[frame],outsuff = "CR", sigclip = sigclip, statsec = statsec, gain = gain
        CR_result_1 = mrdfits(Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_hires'+strcompress(Europa_frames[frame], /rem)+'.CleanedCR.fits', 0, junk_header)
        CR_Mask     = mrdfits(Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_hires'+strcompress(Europa_frames[frame], /rem)+'.Cleaned-mask.fits', 0, junk_header)
        
;        D2Cen               = 609
;        LSF_fitting_ind2  = where( abs(wl - wl[D2cen]) lt 0.2, /NULL)
;        
        
;
        ; the LA Cosmic CR removal algorithm can sometimes introduce negative pixels, particularly at the edge, fix those
        n_sigma = 4. ; 1.*sigclip
        RESISTANT_Mean, CR_result_1, n_sigma, Mean_CR_result_1, Sigma_Mean_CR_result_1, Num_RejECTED
        Sig_CR_result_1 = ROBUST_SIGMA( CR_result_1 )
        CR_result_2 = CR_result_1
        Bad_Pixel_Indicies = where((CR_result_1 lt 0.) or (CR_result_1 gt (Mean_CR_result_1+50.*Sig_CR_result_1)), /NULL, complement = Good_Pixel_Indicies) ; hack
        s = size(aperture, /dim)
        bad_pixel_list = intarr(s[0],s[1])
        bad_pixel_list[Good_Pixel_Indicies] = 1
        fixpix, CR_result_1, bad_pixel_list, CR_result_2     ; replace Bad_Pixel_Indicies with avg of neighbor pixels

        ; the hot pixel correction really messes up the on-disk observations of europa, so i will only use CR_Result_1 for those frames (0, 1, 2, 11, 12, 13). definitely need
        ; hot pixel correction for all frames off the disk though.
        
;        imzap3, CR_result_1, CR_result_2
        
;        window, 10, title='Test CR sub'
;        cgimage, cr_result_1[500:1000,*],  minv=0.25*mean(aperture[500:1000,*]), maxv= 2.5*mean(aperture[500:1000,*])
;        window, 11, title='CR row'
;        cgplot, CR_result_1[500:1000,32]
;        
;        window, 12, title='Test hot pixel removal'
;        cgimage, cr_result_2[500:1000,*],  minv=0.25*mean(aperture[500:1000,*]), maxv= 2.5*mean(aperture[500:1000,*])
;        window, 13, title='CR row'
;        cgplot, CR_result_2[500:1000,32]
;        stop
        
        if order.name eq 'order_46' then begin
          K_cube[*,*,frame] = CR_result_2
          if (frame eq 0) or (frame eq 1) or (frame eq 2) or (frame eq 11) or (frame eq 12) or (frame eq 13) then begin
            K_cube[*,*,frame] = CR_result_1
          endif
        endif
        if order.name eq 'order_60' then begin
          Na_cube[*,*,frame] = CR_result_2
          if labels[frame] eq 'NS on disk' or labels[frame] eq 'EW on disk' or labels[frame] eq 'Juno flyby' then begin
            Na_cube[*,*,frame] = CR_result_1
          endif
          ;          window, 5, title=europa_frames[frame]
          ;          cgplot, CR_result_2, color='red'
          ;          cgplot, CR_result_1, /overplot
          ;          cgplot, Na_cube[*,*,frame], /overplot, color='blue', psym=12
          ;          stop                                                                     ; uncomment this stuff if you want to make sure the frame-by-frame is saving as CR_Result_1/2 is saving correctly
        endif
        if order.name eq 'order_56' then begin
          O_cube[*,*,frame] = CR_result_2
          if (frame eq 0) or (frame eq 1) or (frame eq 2) or (frame eq 11) or (frame eq 12) or (frame eq 13) then begin
            O_cube[*,*,frame] = CR_result_1
          endif
        endif
        
        

        ;CR_mask[Bad_Pixel_Indicies] = 1b

        ; Inspect Cosmic Ray Results
        loadct, 3
        window, 0, xs = 3400, ys = slit_length_pix, title = 'Bias and Flat Corrected Frame: '+ Europa_frames[frame]
        tv, bytscl(aperture, 0.5*mean(aperture[500:1000,*]), 1.5*mean(aperture[500:1000,*]))
        window, 4, xs = 3400, ys = slit_length_pix, ypos = 100, title = 'Cosmic Ray Corrected Frame '+ Europa_frames[frame]
        tv, bytscl(CR_result_1, 0.5*mean(aperture[500:1000,*]), 1.5*mean(aperture[500:1000,*]))
        window, 8, xs = 3400, ys = slit_length_pix, ypos = 200, title = 'Hot/Cold Pixel Filtered Cosmic Ray Corrected Frame '+ Europa_frames[frame]
        tv, bytscl(CR_result_2, 0.5*mean(aperture[500:1000,*]), 1.5*mean(aperture[500:1000,*]))
        window, 12, xs = 3400, ys = slit_length_pix, ypos = 300, title = 'Cosmic Ray Mask'
        tv, bytscl(CR_mask, 0, 1)
        
;        window, 2, xs=800, ys=800, title=labels[frame]
;        cgimage, aperture[450:1000,*], minv=-5, maxv=50
        if frame eq 0 then sunframe = aperture                              ; this saves one of the on-disk spectra BEFORE CR correction so that it can be used for solar subtraction later
        
        
        
;;;;;;;;;;;;;;;;; i'm running this test to see if on-disk spectra are already brighter than off-disk and they ARE.
;;;;;;;;;;;;;;;;; frame 3 (EW on disk) at 10 R_E west ~250 DN/s
;;;;;;;;;;;;;;;;; frame 24 (NS10W) on disk ~100 DN/s
;;;;;;;;;;;;;;;;; these should match up, but frame 3 at the correct location is already brighter than frame 24.
;;;;;;;;;;;;;;;;; tested this before and after CR sub too
        
;        plate_scale = 0.358
;        suncol = where(cr_result_1[700,*] eq max(cr_result_1[700,*]))
;        ang_radius = 0.5
;        yr = [-(suncol*plate_scale/ang_radius), ((s[1]-suncol)*plate_scale/ang_radius)]
;        window, 11, title=labels[frame]
;        cgplot, total(cr_result_1[500:1000,*], 1), xr=[0,44], yr=[0,1.e3], xtickformat='(A1)'
;        cgaxis, xaxis=0, xtit='Europa Radii', xr=yr, xstyle=1, xticklen=-0.05
;        stop
        
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
;        D2Cen               = 609
;        D1Cen               = 825
;
;        LSF_fitting_ind1  = where( abs(wl - wl[D1cen]) lt 0.2, /NULL)
;        LSF_fitting_ind2  = where( abs(wl - wl[D2cen]) lt 0.2, /NULL)
;        cloud_ind1        = findgen(42)+774.; where( abs(wl - wl[D1cen]) lt 1.5, /NULL)
;        cloud_ind2        = findgen(44)+553.; where( abs(wl - wl[D2cen]) lt 1.5, /NULL)
;        
;        dummy_cr                   = cr_result_1
;        dummy_cr[LSF_fitting_ind1] = !Values.F_Nan
;        dummy_cr[LSF_fitting_ind2] = !Values.F_Nan
;        dummy_cr[cloud_ind1      ] = !Values.F_Nan
;        dummy_cr[cloud_ind2      ] = !Values.F_Nan
;        
;        interp_to_mean = mean(dummy_cr, dim=2)
;        CR_result_2    = interpol(interp_to_mean, WL, WL, /NAN)
;        
;        window, 7
;        cgplot, interp_to_mean[500:1000]
;        cgplot, cr_result_2[500:1000], /overplot, color='red'
        
        
;        window, 9, title='Test CR sub'
;        cgimage, cr_result_1[500:1000,*], minv=0, maxv=0.2
;        stop
;        wait, 2
      endfor ; frame (Europa frames)

      ; Now do the same process for Jupiter
      Frames = '\hires'+Jupiter_frames+'.Cleaned.fits'
      for frame = 0, n_elements(Frames)-1 do begin
        Jupiter  = mrdfits(Dir+'\Processed\' + Frames[frame], 0, header)
        Jupiter_order_straight = fltarr(s[0], slit_length_pix)
        ;Flat_order_straight = fltarr(s[0], 30)
        for i = 0, s[0] - 1 do begin
          Jupiter_order_straight[i,*] = interpolate(Jupiter[i,*], order.low_bound + findgen(slit_length_pix) - slit_length_pix/2. + POLY( i, coeffs))
        endfor

        ; prepare the field field
        Flat_aperture = Flat_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;nomalize_with = mean(Flat_aperture, dim = 2, /NAN)
        coeffs = poly_fit(findgen(s[0]), mean(Flat_aperture, dim = 2, /NAN), 3)
        nomalize_with = poly(findgen(s[0]), coeffs)
        Flat_aperture = Flat_aperture / rebin(nomalize_with, s[0],order.aperture_limit[1] - order.aperture_limit[0] + 1) ; Normalize the flat field to unity, now its spectral AND spatially normalized

        if order.name eq 'order_46' then K_Jupiter_cube[*,*,frame] = Jupiter_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]] ;/ Flat_aperture ; flat field jupiter
        if order.name eq 'order_56' then O_Jupiter_cube[*,*,frame] = Jupiter_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]] ;/ Flat_aperture ; flat field jupiter
        if order.name eq 'order_60' then Na_Jupiter_cube[*,*,frame] = Jupiter_order_straight;[*, order.aperture_limit[0]:order.aperture_limit[1]]; / Flat_aperture ; flat field jupiter

      endfor ; frames (Jupiter frame number)

      if order.name eq 'order_46' then  K_WL = WL
      if order.name eq 'order_56' then  O_WL = WL
      if order.name eq 'order_60' then Na_WL = WL

      save, Na_cube, K_cube, O_cube, Na_Jupiter_cube, K_Jupiter_cube, O_Jupiter_cube, Na_WL, K_WL, O_WL, sunframe, filename = Dir+'\Processed\.sav'
    endfor ; h (order number)
    stop
  endif

  if part eq 1.5 then begin ; ==================== FLUX CALIBRATE ALL ORDERS INTO RAYLEIGHS PER ANGSTROM UNITS ===============================


    ; --------------------------------------------------- Determine Sensitivy for flux calibration ------------------------------------------

    ; Compare publications of Jupiter's spectral albedo at disk center

    ; absolute brightness: Digitized Plot from Woodman et al. 1979.
    READCOL,'Z:\DATA\___Calibration___\Jupiter_Reflectivity\Woodman_et_al_1979_plot1.txt', F='A,A', WL_1, Albedo_1, STRINGSKIP = '#', /Silent;wavelength in Angstroms I/F unitless
    READCOL,'Z:\DATA\___Calibration___\Jupiter_Reflectivity\Woodman_et_al_1979_plot2_new.txt', F='A,A', WL_2, Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
    Woodman_WL = float([WL_1, WL_2])             ; STITCH THESE TOGETHER
    Woodman_Albedo = Float([albedo_1, albedo_2]) ; STITCH THESE TOGETHER

    ; absolute brightness: from Karkoschka (1998) Icarus on the PDS as ID # ESO-J/S/N/U-SPECTROPHOTOMETER-4-V1.0
    READCOL,'Z:\DATA\___Calibration___\Jupiter_Reflectivity\Karkoschka_1995low.tab', F='X,A,X,A', Karkoschka_WL, Karkoschka_Albedo, STRINGSKIP = '#', /Silent ;wavelength in nm I/F unitless
    READCOL,'Z:\DATA\___Calibration___\Jupiter_Reflectivity\Karkoschka_1995high.tab', F='X,A,X,A', Karkoschka_WL_2, Karkoschka_Albedo_2, STRINGSKIP = '#', /Silent ;wavelength in Angstroms I/F unitless
    Karkoschka_WL = float([Karkoschka_WL, Karkoschka_WL_2])             ; STITCH THESE TOGETHER
    Karkoschka_Albedo = Float([Karkoschka_albedo, Karkoschka_albedo_2]) ; STITCH THESE TOGETHER
    Karkoschka_Albedo = Karkoschka_Albedo[sort(Karkoschka_WL)]
    Karkoschka_WL = Karkoschka_WL[sort(Karkoschka_WL)]

    ; compare the two
    cgplot, Woodman_WL / 10., Woodman_Albedo, color = 'blue', xstyle = 1., psym = 3, Xtitle = 'Wavelength (nm)', ytitle = 'I/F Reflectivity'
    cgplot, Karkoschka_WL, Karkoschka_Albedo*1.35, color = 'red', /overplot
    cgtext, 340, .1, 'EQUATOR AT CENTRAL MERIDIAN (Woodman et al. 1979)', color = 'blue'
    cgtext, 340, .16, 'FULL DISK scaled by 1.35 (Karkoschka 1998)', color = 'red'

    ; make an informed choice via scaling
    Karkoschka_wl = Karkoschka_wl * 10. ;nm to Angstroms
    Karkoschka_Albedo = Karkoschka_Albedo*1.35 ;Scale the "Full disk albedo" to the "Central Meridian Equatorial Absolute Reflectivity"

    ; Get the Solar Spectral Irradiance at 1 AU
    READCOL,'Z:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\Kurucz\Kurucz_2005_irradthuwl.dat', F='A,A', WL_nm, flux, STRINGSKIP = '#', /Silent ;flux is in W/m2/nm
    start  = where(WL_nm eq '299.100')
    WL_nm = float(WL_nm[start:*])
    flux = float(flux[start:*])

    ; change flux units from W/m^2/nm to photons/(cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
    conversion = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
    flux = flux * conversion      ; Cross-checked this result against Huebner et al. 1992.
    WL_A = temporary(WL_nm) * 10. ; Wavelength from nm into angstroms
    VACTOAIR, WL_A, WL_A_Air      ; Vacuum to air wavelength conversion
    WL_A = temporary(WL_A_Air)

    Jupiter_center_header = headfits(Dir+'\Processed\hires'+strcompress(Jupiter_frames[0], /rem)+'.Cleaned.fits')
    cspice_UTC2ET, sxpar(Jupiter_center_header, 'DATE-OBS') + ' '+ sxpar(Jupiter_center_header, 'UTC'), ET
    cspice_spkezr, 'Jupiter', ET, 'J2000', 'LT+S', 'Sun', Jupiter_Sun_State, ltime
    solar_distance = norm(Jupiter_Sun_State[0:2]) / 149597871.                                    ; converts to AU
    flux_at_jupiter = flux / solar_distance^2.

    ; Multiply incident solar irradiance x spectral albedo to determine the theoreitcal brightness of Jupiter at disk center
    Albedo = INTERPOL(Karkoschka_Albedo, Karkoschka_WL, WL_A)
    Rayleighs_per_angstrom = 4.*flux_at_jupiter*albedo / 1.e6
    if keyword_set(debug) then window, 5, Title = 'Instantaneous Rayleighs per Angstrom: Center of Jupiter''s Disk'
    if keyword_set(debug) then plot, WL_A, Rayleighs_per_angstrom, xr = [5885, 5900], charsize = 2 ;compare to 5.5 MR per angstrom by Brown & Schneider, 1981

    ; Adjust the expected absolute flux for Jupiter's Instantaneous Doppler shift

    ; Get the body-fixed non-inertial IAU_Earth Coordinates of Keck Observatory in Cartesian
    cspice_bodvrd, 'EARTH', 'RADII', 3, radii
    flat = (radii[0] - radii[2])/radii[0]
    OBSERVATORY_cs, 'Keck', obs_struct ;longitude is in degrees *West*, need to switch to East Longit for J2000 conversion via cspice_georec
    cspice_georec, (360. - obs_struct.longitude) * CSPICE_RPD(), obs_struct.latitude * CSPICE_RPD(), obs_struct.altitude / 1.e3, $
      radii[0], flat, obs_IAU_Earth
    origin = [obs_IAU_Earth, 0, 0, 0]

    ; Calculate radial velocity and adjust wavelength array for instantaneous Dopplershift
    cspice_spkezr, 'Jupiter', et, 'IAU_EARTH', 'lt+s', 'Earth', Jupiter_Earth_State, ltime
    state = Jupiter_Earth_State - origin
    Jupiter_wrt_Earth_Dopplershift = total(state[0:2]*state[3:5]/sqrt(total(state[0:2]^2)))
    WL_A = WL_A + WL_A*Jupiter_wrt_Earth_Dopplershift/cspice_clight()

    ; Find the angular radius of Europa
    cspice_bodvrd, 'Europa', 'RADII', 3, radii
    ang_radius = 206265.*tan(radii[0]/sqrt(total(state[0:2]^2)))
    
    ; --------------------------------- Flux Calibrate Using Distance and Absolute Spectral Reflectivity of Jupiter -------------------------------------

    ; D3 smooth_by = 3.6
    ; C3 smooth_by = 2.1

    smooth_by = 2.1

    FOR order = 2, N_elements(orders) - 1 do begin
      restore, dir + '\Processed\.sav'

      if orders[order].name eq 'order_46' then Jupiter_cube = K_Jupiter_cube
      if orders[order].name eq 'order_56' then Jupiter_cube = O_Jupiter_cube
      if orders[order].name eq 'order_60' then Jupiter_cube = Na_Jupiter_cube
      if orders[order].name eq 'order_46' then WL =  K_WL
      if orders[order].name eq 'order_56' then WL =  O_WL
      if orders[order].name eq 'order_60' then WL = Na_WL

      ; Determine the instrumental sensitivity from the expected versus measured flux at Jupiter Disk Center. This should be a smooth function
      expected_flux          = interpol(Rayleighs_per_angstrom, WL_A, WL)       ; move expected flux to the Jovian Doppler-shift, UNITS are R / A
      smoothed_expected_flux = GAUSS_SMOOTH(expected_flux, smooth_by, /edge_truncate) ; this smoothing looks about right for HIRES D3

      ; find the count rate in the middle 3 rows at slit center

      s = size(Jupiter_cube, /dim)
      Jupiter_DN_per_s_at_disk_center = median(Jupiter_cube[*, s[1]/2.-1:s[1]/2.+1], dim=2)

      ; The ThAr wavelength solution isn't always perfect match, which will throw the sensitivity fit.
      ; Find any subtle Doppler shifting and correct for it. Cross-correlate versus lag to align them to the nearest pixel
      lag         = indgen(11) - 5
      correlation = C_CORRELATE(Jupiter_DN_per_s_at_disk_center, smoothed_expected_flux, lag)
      junk        = max(correlation, corr_ind)
      shifted_expected_flux = shift(expected_flux, -lag[corr_ind])
      matched_expected_flux = GAUSS_SMOOTH(shifted_expected_flux, 3.6, /edge_truncate)

      ; Calculate sensitivity and fit the spline, avoiding spectral lines
      Sensitivity          = Jupiter_DN_per_s_at_disk_center / matched_expected_flux            ; Sensitivity in (DN / S) / (R / A)
      fit_Sensitivity      = medsmooth(Sensitivity, 100)

      ; Inspect the measured sensitivity. If the flat field worked, then this should not be curved with the echelle blaze, but rather a flatish line with some spectral artefacts.
      window, 6, xs=1024, ys=1024
      cgplot, WL, sensitivity, /xs, ytitle = '(DN / S) / (R / A)', xtitle = 'angstroms', title = 'measured (black) vs fit (red) sensitivity'
      cgplot, WL, fit_Sensitivity, color = 'red', /overplot

      ; Calibrate both Jupiter and Europa spectra INTO RAYLEIGHS PER ANGSTROM UNITS for this order
      _2D_sensitivity = rebin(fit_Sensitivity, s[0], s[1])
      Jupiter_cube = Jupiter_cube / _2D_sensitivity                                             ; jupiter flux in R / A
;      cgplot, mean(jupiter_cube, dim=2)                                                         ; should be about 5.5 MR / A
;stop
      for f = 0, N_elements(Na_cube[0,0,*])-1 do begin
;        window, 7, xs=4000, ys=44
;        cgimage, Na_cube[*,*,f], minv=0, maxv=30
;        stop
        Na_cube[*,*,f] = Na_cube[*,*,f] / _2D_sensitivity ; convert to R / A units
        K_cube[*,*,f]  =  K_cube[*,*,f] / _2D_sensitivity ; convert to R / A units
        O_cube[*,*,f]  =  O_cube[*,*,f] / _2D_sensitivity ; convert to R / A units
        print, 'you neglected differential airmass when scaling the sensitivity calculation from Jupiter to Europa!'
      endfor

      if orders[order].name eq 'order_46' then K_Jupiter_cube  = Jupiter_cube
      if orders[order].name eq 'order_56' then O_Jupiter_cube  = Jupiter_cube
      if orders[order].name eq 'order_60' then Na_Jupiter_cube = Jupiter_cube

      ; In a given row above be sure the answer matches the following notes...
      ;   -Trafton (1980) does the same thing I do: (page 118 here https://doi.org/10.1016/0019-1035(80)90249-3). A pi on both sides, but needs the factor of 4 to put it in Rayleighs
      ;   -Brown (1981) (doi:10.1086/158777) quotes 5.6 MR/A at 6300A. I get about 5.9 MR/A if Jupiter were its average distance of 5.2 AU
      ;   -Brown & Shemansky (1982) doi:10.1086/160515 quote 5.4 MR/A at 6724A. I get 5.7 MR/A if Jupiter were its average distance of 5.2 AU

      save, Na_cube, K_cube, O_cube, Na_Jupiter_cube, K_Jupiter_cube, O_Jupiter_cube, Na_WL, K_WL, O_WL, sunframe, ang_radius, _2D_sensitivity, filename = Dir+'\Processed\_Flux_Calibrated.sav'
      print, 'Flux-calibrated ', orders[order].name
    endfor
    stop
    PRINT, 'To Do: did not yet inspect flux calibration of K D and O 6300 orders!'
  endif

  if part eq 1.6 then begin
    restore, Dir+'\Processed\Europa_Airglow_params.sav'
    help, Europa_Airglow_params

    filenames = []
    for i = 0, n_elements(Europa_frames)-1 do begin
      filename = dir+'\Processed\Cosmic Rays\order_60_CR_hires' + Europa_frames[i] + '.CleanedCR.fits'
      header   =  headfits(filename)
      filenames= [filenames, filename]
    endfor
    

    SOLAR_SPECTRUM_FILE = 'Z:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\Coddington_2022\hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
    NCDF_LIST, SOLAR_SPECTRUM_FILE, /VARIABLES, /DIMENSIONS, /GATT, /VATT

    ID = NCDF_open(TEMPORARY(SOLAR_SPECTRUM_FILE))

    SSI = NCDF_VARID(id, 'SSI')
    NCDF_VARGET, id, SSI, Flux

    Vac_WL = NCDF_VARID(id, 'Vacuum Wavelength')
    NCDF_VARGET, id, Vac_WL, WL_nm

    NCDF_close, ID

    ; change flux units from W/m^2/nm to photons / (cm^2 s A)
    ; multiply by ((lambda / hc) / 1 W)  * (1 m^2 / 1e4 cm^2) * (1 nm / 10 A)
    conversion   = ((WL_nm*1.e-9)/(6.62606957e-34*299792458.D)) * (1./1.e4) * (1./10.)
    flux = flux * conversion                                                                  ; photons / (cm^2 s A)
    WL_A = temporary(WL_nm) * 10.

    junk = min(abs(WL_A - 5894.6), Na_ind)
    junk = min(abs(WL_A - 7680.), K_ind)


    ; ============ SUNLIGHT SUBTRACTION. I'm using Europa's on-disk observations as sunlight spectrum by assuming it's perfectly ⋆｡˚ ☁︎ ˚shiny｡⋆｡˚☽˚｡⋆ ========

    restore, Dir+'\Processed\_Flux_Calibrated.sav'
    for order_index = 2,2 do begin
      if order_index eq 2 then begin
        order = orders[order_index]
        cube  = Na_cube
        WL    = Na_WL
      endif
      if order_index eq 0 then begin
        order = orders[order_index]
        cube  = K_cube
        WL    = K_WL
      endif
      if order_index eq 1 then begin
        order = orders[order_index]
        cube  = O_cube
        WL    = O_WL
      endif

      if filt eq 'gg475' then begin
        EW0__   = 0.5*(REFORM(cube[*,*,0 ]+cube[*,*,1 ]+cube[*,*,2 ]))                                        ; EW orientation on-disk
        EW10_N  = 0.5*(REFORM(cube[*,*,3 ]+cube[*,*,4 ]))
        EW20_N  = 0.5*(REFORM(cube[*,*,5 ]+cube[*,*,6 ]))
        EW10_S  = 0.5*(REFORM(cube[*,*,7 ]+cube[*,*,8 ]))
        EW20_S  = 0.5*(REFORM(cube[*,*,9 ]+cube[*,*,10]))
        NS0__   = 0.5*(REFORM(cube[*,*,11]+cube[*,*,12]+cube[*,*,13]+cube[*,*,22] $
                             +cube[*,*,23]+cube[*,*,28]+cube[*,*,29]))                                        ; NS orientation on-disk
        NS10_W  = 0.5*(REFORM(cube[*,*,14]+cube[*,*,15]+cube[*,*,30]+cube[*,*,31] $
                             +cube[*,*,24]+cube[*,*,25]))
        NS20_W  = 0.5*(REFORM(cube[*,*,16]+cube[*,*,17]))
        NS10_E  = 0.5*(REFORM(cube[*,*,18]+cube[*,*,19]+cube[*,*,32]+cube[*,*,33] $
                             +cube[*,*,26]+cube[*,*,27]))
        NS20_E  = 0.5*(REFORM(cube[*,*,20]+cube[*,*,21]))
      endif


      ; we'll use NS0 only for now and that's the 12th index in the spectra we observed...
      sun2europa = Europa_Airglow_params.sun2euro[12]
      s = size(cube[*,*,0])

      if filt eq 'Na'    then orientations = fltarr(s[1], s[2], 11)                         ; there's definitely a better way to do this but whatever
      if filt eq 'gg475' then orientations = fltarr(s[1], s[2], 10)

      
      ; -------------------------------------------------- Get sunlight image from on-disk obs. -----------------------------------------------
      
; going to take an on-disk spectrum, get the brightest row and say that's sunlight. then i will rebin that one row into a 2D image and use that as my sunlight spectrum
      
      if order_index eq 2 then xr    = [5888.5, 5898.5] 
      junk  = min(abs(xr[0]- WL), index0)       
      junk  = min(abs(xr[1]- WL), index1)       
      
      sunspectrum = 3
      ondisk                      = fltarr(s[1], s[2])
      ondisk[index0:index1,*]     = REFORM(cube[index0:index1,*,sunspectrum]) ;  0.5*(REFORM(cube[*,*,20]))                            ; took a random one, this is EW on disk
      
      ondisk                      = sunframe
      
      
      ondisk1d                    = total(ondisk, 1, /nan)
      suncol                      = WHERE(ondisk1d EQ MAX(ondisk1d),count)
      
      sunspec = REBIN(ondisk[*,suncol], s[1], s[2])
      dummy_sunspec = sunspec
      
      
      sunmax = []
      newimg    = fltarr(s[1], s[2])
      sunimg    = fltarr(s[1], s[2])
      eurimg    = fltarr(s[1], s[2])
      scalesubb = fltarr(s[1], s[2])
      fakeio    = fltarr(s[1])
      no_io     = fltarr(s[1])

      ; ------------------------------------------- Find the Dispersion & Sunlight vs Exosphere Indicies --------------------------------------

      D2Cen               = 609
      D1Cen               = 825
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;WORK DONE ON 6/13. maybe epic fail. who knows.;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;     
;
; 
; 
;      sun_ind1          = where( abs(wl - wl[D1cen]) lt 0.8, /NULL)
;      ;      sun_ind1          = LSF_fitting_ind1[where( LSF_fitting_ind1 lt 1000., /NULL)]
;      sun_ind2          = where( abs(wl - wl[D2cen]) lt 0.8, /NULL)
;      ;      LSF_fitting_ind2  = LSF_fitting_ind2[where( LSF_fitting_ind2 lt 700. and LSF_fitting_ind2 gt 600.)]
;      cloud_ind1        = findgen(42)+774.; where( abs(wl - wl[D1cen]) lt 1.5, /NULL)
;      cloud_ind2        = findgen(44)+553.; where( abs(wl - wl[D2cen]) lt 1.5, /NULL)
;
;      windowwidth         = 30.
;      spec_1D             = total(ondisk, 2, /Nan)
;      result              = mpfitpeak(findgen(N_elements(sun_ind2)), spec_1D[sun_ind2], a, STATUS = STATUS)
;      D2_Solar            = sun_ind2[0] + a[1]
;
;      window, 2, xs=800, ys=500
;      cgplot, wl[sun_ind2], spec_1D[sun_ind2]
;      cgplot, wl[sun_ind2], gaussian(wl[sun_ind2], a), /overplot, color='red'                                         ; all plots in DN/s
;      stop
;      
;      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
      
      windowwidth         = 30.
      spec_1D             = total(ondisk, 2, /Nan)
      result              = mpfitpeak(findgen(windowwidth*2. + 1), spec_1D[D2cen - windowwidth:D2cen + windowwidth], a, STATUS = STATUS)
      D2_Solar            = D2cen - windowwidth + a[1]
      result              = mpfitpeak(findgen(windowwidth*2. + 1), spec_1D[D1cen - windowwidth:D1cen + windowwidth], a, STATUS = STATUS)
      D1_Solar            = D1cen - windowwidth + a[1]
      dispersion          = (5895.92424 - 5889.95095) / ( D1_Solar - D2_Solar )                                 ; A/pixel

      d2_sep              = sun2europa * 5889.95095 / (cspice_clight() * dispersion)                            ; # of pixels between solar wells and Europa's emission
      d1_sep              = sun2europa * 5895.92424 / (cspice_clight() * dispersion)                            ; # of pixels between solar wells and Europa's emission

      Europa_D2_pixel     = D2_Solar + D2_sep
      Europa_D1_pixel     = D1_Solar + D1_sep

      ; Get the pixel indices where the spectrum consists of just scattered sunlight
      center              = round((Europa_D1_pixel + Europa_D2_pixel) / 2.)
      fitindices          = where( (abs(center - 200 + findgen(400) - Europa_D1_pixel) gt 5) and $
                                   (abs(center - 200 + findgen(400) - Europa_D2_pixel) gt 5), /null)            ; Excludes the emission from Europa
      fitindices          = fitindices + center - 200
      
      labels = strarr(38);N_elements(cube[0,0,*]))
      if filt eq 'Na' then begin
        labels[2:5  ] = 'EW on disk'
        labels[6:7  ] = 'EW10 N'
        labels[8:9  ] = 'EW20 N'
        labels[10:14] = 'Juno flyby'
        labels[15:16] = 'EW10 S'
        labels[17:18] = 'EW20 S'
        labels[19:23] = 'NS on disk'                     ; file 147 is skipped because it's the sunlight spectrum
        labels[24:25] = 'NS10 W'
        labels[26:27] = 'NS20 W'
        labels[28:29] = 'NS10 E'
        labels[30:31] = 'NS20 E'
        labels[32:34] = 'NS on disk'
        labels[35:36] = 'EW on disk'
      endif

      P_returned       = fltarr(5,s[2])                  ; Three coefficients + MPFIT's "Status"
      P_guessed        = fltarr(4,s[2])                  ; Initial Guess that we throw at MPFIT
      gg475_new_images = fltarr(s[1], s[2], N_elements(labels)-1)
      Na_new_images    = fltarr(s[1], s[2], N_elements(labels))
      Nacolumns        = fltarr(s[2], N_elements(labels))
      int_Nacolumns    = fltarr(s[2], N_elements(labels))
      GGcolumns        = fltarr(s[2], N_elements(labels))
      
      LSF_fitting_ind1  = where( abs(wl - wl[D1cen]) lt 0.4, /NULL)
      LSF_fitting_ind1  = LSF_fitting_ind1[where( LSF_fitting_ind1 lt 1000., /NULL)]
      LSF_fitting_ind2  = where( abs(wl - wl[D2cen]) lt 0.4, /NULL)
      LSF_fitting_ind2  = LSF_fitting_ind2[where( LSF_fitting_ind2 lt 700. and LSF_fitting_ind2 gt 600.)]
      cloud_ind1        = findgen(42)+774.; where( abs(wl - wl[D1cen]) lt 1.5, /NULL)
      cloud_ind2        = findgen(44)+553.; where( abs(wl - wl[D2cen]) lt 1.5, /NULL)
      
      sunrow    = sunspec[*, suncol]
;      cgplot, sunrow[fitindices]
      
      
      guess_peaks = []
      all_EW      = []
      FOR orientation = 2, N_Elements(cube[0,0,*])-1 DO BEGIN ;cube[0,0,*]) - 1 DO BEGIN       ; N_Elements(orientations[0,0,*]) - 1 DO BEGIN
        if orientation eq sunspectrum then continue                            ; skips the NS on-disk frame i used for the fake sunlight spectrum
        
        mpfitD2emission = []
        mpfitD1emission = []
        int_D2_emission = []
        int_D1_emission = []
        sunsubbed       = fltarr(s[1],s[2])
        
        europa        = cube[*,*,orientation]
        dummy_europa  = europa
        
        dummy_sunspec[LSF_fitting_ind1,*] = !Values.F_Nan
        dummy_sunspec[LSF_fitting_ind2,*] = !Values.F_Nan
        dummy_sunspec[cloud_ind1      ,*] = !Values.F_Nan
        dummy_sunspec[cloud_ind2      ,*] = !Values.F_Nan

        dummy_europa[LSF_fitting_ind1,*] = !Values.F_Nan
        dummy_europa[LSF_fitting_ind2,*] = !Values.F_Nan
        dummy_europa[cloud_ind1      ,*] = !Values.F_Nan
        dummy_europa[cloud_ind2      ,*] = !Values.F_Nan
        
        
        FOR i = 0, s[2] - 1 DO BEGIN
          row         = europa[*,i]
          guess_scale = median(row[index0:index1] / sunrow[index0:index1])
          row_err     = sqrt(abs(row))
          
          ; Fit a y = A*Gauss_smooth(x,C) + B function to the spectrum, where x is the reference solar spectrum
          p0 = [guess_scale, max(sunrow[index0:index1]), 0.1, 0.0]                                           ; Guess at initial coefficients
          parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p0))
          
          parinfo.value         = p0
          parinfo[0].limited    = [1, 1]
          parinfo[0].limits     = [0.0, 2.]
          parinfo[1].limited    = [1, 1]
          parinfo[1].limits     = [0, 4.e6]
          parinfo[2].limited    = [1, 1]
          parinfo[2].limits     = [0.0, 10.]
          
          
          WEIGHTS = 1./(abs(findgen(s[1]) - Europa_D2_pixel))^.4 + 1./(abs(findgen(s[1]) - Europa_D1_pixel))^.4
          
          fa = {x:sunrow[index0:index1], y:row[index0:index1], err:sqrt(abs(row[index0:index1]))}
          p = mpfit('match_scattered_sunlight', p0, PERROR = err_P, functargs=fa, status=status, parinfo=parinfo, quiet=quiet)
          P_guessed[*,i]  = p0
          p_returned[*,i] = [p, status]
          scaled_sunlight = scale_fit_sunlight(P, sunrow)                   ; Puts it into y = A*shift(Gauss_smooth(x,D),C) + B form
;          
;                window, 2, xs=800, ys=500
;                cgplot, wl[index0:index1], sunrow[index0:index1]
;                cgplot, wl[index0:index1], scaled_sunlight[index0:index1], /overplot, color='red'                                         ; all plots in DN/s
;                stop
          
          sunimg[*,i] = scaled_sunlight
          sub         = guess_scale * sunrow                                ; If you JUST want the multiplicative correction (no scattered sunlight accounted for w mpfit)
          
          totsubtrd   = row  - scaled_sunlight                              ; WITH the fit sunlight...except am I even able to fit a gaussian to the whole sun spec?
          totsubtrd   = row  - sub                                          ; JUST multiplicative factor
;          if i eq suncol then continuum = scaled_sunlight                   ; this saves the non-sun subtracted continuum row so that i can reference it later
         
;          metrics=[]
;          qualitymetric = stddev(row, /nan) / abs(total(row, /nan))
;          metrics = [metrics, qualitymetric]

          if (labels[orientation] eq 'NS on disk') or (labels[orientation] eq 'EW on disk') $
            or (labels[orientation] eq 'Juno flyby') then begin
            if i gt 16 and i lt 28 then begin                                                         ; HACK, i'm eyeballing how many pixels i should block out based on solar subtraction
              ;                if qualitymetric lt 0 or qualitymetric gt 0.001 then begin             ; sets threshold and gets rid of sunlight over disk
              totsubtrd = !values.F_nan                                                               

              ;                endif
            endif
          endif
          sunsubbed[*,i] = totsubtrd
        endfor
        
        cgimage, sunsubbed[index0:index1,*]
        stop
;        loadct, 3
;        window, 0, title=labels[orientation]+' sun subtracted'
;        cgimage, sunsubbed[500:1000,*]
;        
;        window, 2, title=labels[orientation]
;        cgimage, europa[500:1000,*]
;        
;        window, 1, title='totaled along spatial dimension'
;        cgplot, total(sunsubbed[500:1000,*], 2, /nan)
;        
;        window, 3, title='totaled along spatial dimension'
;        cgplot, total(europa[500:1000,*], 2, /nan)
;        print, orientation, '  ', europa_frames[orientation]
;        stop
        

        
; ------------------------------------------------- io subtraction ------------------------------------------------------------
; this is me trying to subtract io cloud. first, i block out the europa emission & continuum and then interpolate over these blocked out regions.
        dummy_sunsubbed = sunsubbed
        dummy_sunsubbed[LSF_fitting_ind1,*] = !Values.F_NaN
        dummy_sunsubbed[LSF_fitting_ind2,*] = !Values.F_NaN                                                   ; HACK HACK HACK. WL solution is weird.
        
        meanimg = median(dummy_sunsubbed, dim=2)
        fakeio[index0:index1] = interpol(meanimg[index0:index1], WL[index0:index1], WL[index0:index1], /NAN)
;        window, 6
;        cgplot, fakeio[index0:index1,*]
;       stop
        
        for i = 0, s[2] - 1 DO BEGIN
          dummy_row = dummy_sunsubbed[*,i]
          
          if total(dummy_row) eq !Values.F_Nan then continue
          if max(dummy_row) eq 0. then scaled_row = dummy_row
          if max(dummy_row) ne 0. then scaled_row = dummy_row * (max(meanimg) / max(dummy_row))        ; scales row to the spatially resolved mean of sunlight-subtracted image

;          
;;          eurimg[index0,i] = 1.e6                                                                     ; use this to track which row you're looking at

;          window, 4, title=labels[orientation]
;          cgimage, eurimg[index0:index1,*];, $
            ;minv=-3.5*mean(eurimg[index0:index1,*], /nan), maxv=3.5*mean(eurimg[index0:index1,*], /nan)
          
          no_io[index0:index1] = dummy_row[index0:index1] - fakeio[index0:index1]
          
;          window, 3, xs=1000, ys=400, title=labels[orientation]+' io sub row '+string(i)
;          cgplot, WL[index0:index1], dummy_row[index0:index1], yr= [-100,100];[min(no_io[index0:index1]),max(dummy_row[index0:index1])]             ; black = sun-subbed
;          cgplot, WL[index0:index1], scaled_row[index0:index1], /overplot, color='green'                          ; green = sun-subbed row, scaled
;          cgplot, wl[index0:index1], fakeio[index0:index1], /overplot, color='red', psym=12                       ; red   = fake io, where we interpolated over D-lines

;;                  no_io[LSF_fitting_ind1] = totsubtrd[LSF_fitting_ind1]
;;                  no_io[LSF_fitting_ind2] = totsubtrd[LSF_fitting_ind2]
                  noIo = scaled_row - fakeio
                  
          if dummy_row[0] eq !values.F_NAN then totsubtrd = !values.F_NAN
          if dummy_row[0] eq !values.F_nan then continue
          
;          cgplot, WL[index0:index1], noio[index0:index1], /overplot, symsize=3, color='teal'                      ; teal  = scaled row - fake io row
;          cgplot, WL[index0:index1], no_io, /overplot, color='pink'                                ; pink  = black (sun-subbed, NOT spatially resolved) row - fake io row
          
          newimg[*,i] = no_Io ; noIo                                                                              ; choosing between if i should use scaled subtraction for io cloud or not... right now the scaled works better
        
          
          if (labels[orientation] eq 'NS on disk') or (labels[orientation] eq 'EW on disk') $
            or (labels[orientation] eq 'Juno flyby') then begin
              
              newimg[LSF_fitting_ind1, *] = sunsubbed[LSF_fitting_ind1, *]
              newimg[LSF_fitting_ind2, *] = sunsubbed[LSF_fitting_ind2, *]
             
;              newimg[LSF_fitting_ind1, 0:16] = sunsubbed[LSF_fitting_ind1, 0:16]
;              newimg[LSF_fitting_ind1, 28:*] = sunsubbed[LSF_fitting_ind1, 28:*]
;              newimg[LSF_fitting_ind2, 0:16] = sunsubbed[LSF_fitting_ind2, 0:16]
;              newimg[LSF_fitting_ind2, 28:*] = sunsubbed[LSF_fitting_ind2, 28:*]
              
                                                            ; quick calculation to derive Na exosphere temperature from scale height. i'm using this sunspectrum to find
                                                            ; the fall off rate of Na to get scale height.
                                              
;                                                            window, 0
;                                                            cgimage, newimg[LSF_fitting_ind2,*], minv=0, maxv=1.e6
;                                              
;                                                            sun_cutoff = 16.                                   ; pretty arbitrary, but chose closest row that is not contaminated by sunlight
;                                                            NaOnDisk = max(newimg[LSF_fitting_ind2,sun_cutoff]);median(ondisk[LSF_fitting_ind2,sun_cutoff])
;                                              
;                                                            scale_height_val = NaOnDisk / 2.7182818284590452
;                                                            NaAlongSlit = []
;                                                            for i = 0, 20 do begin;'s[2]-1 do begin
;                                                              Narowbyrow = median(newimg[LSF_fitting_ind2,i])
;                                                              NaAlongSlit = [NaAlongSlit, Narowbyrow]
;                                                            endfor
;                                                            differences  = abs(NaAlongSlit - scale_height_val)
;                                                            scale_height = where(differences eq min(differences))
;                                              
;                                                            window, 1, title="minimize this array to find fall off rate"
;                                                            cgplot, differences, /ynozero
;                                                            window, 3
;                                                            cgplot, Naalongslit, psym=2, symsize=2
;                                                            cgplot, fltarr(44)+scale_height_val, /overplot, color='red'
;                                              
;                                                            plate_scale = 0.358
;                                                            scale_height_radius = abs(sun_cutoff-scale_height)*plate_scale/ang_radius
;                                                            ;                                scale_height_radius = 5.0
;                                                            scale_height_m = scale_height_radius * 1560. * 1.e3
;                                              
;                                                            Namass = 3.817E-26 ; kg
;                                                            europamass = 4.8e22 ; kg, googled it
;                                                            dist_bt_row16_and_europacent = abs(suncol - sun_cutoff)*plate_scale/ang_radius
;                                                            ;                                europag = !const.G * europamass / (dist_bt_row16_and_europacent)^2.
;                                                            europag = 1.315 / scale_height_radius^2    ; m/s^2
;                                              
;                                                            Na_T = scale_height_m * Namass * europag / !const.k
;                                                            print, "Scale height of Na is",scale_height_radius,"Re, corresponding to T =",Na_T,"K"
;                                              
;                                                            stop
              
              endif
              
              
              
;          endif else begin
            newimg[LSF_fitting_ind1, *] = sunsubbed[LSF_fitting_ind1, *]
            newimg[LSF_fitting_ind2, *] = sunsubbed[LSF_fitting_ind2, *]
;          endelse
          
        ENDFOR ; io subtraction
        
        
        window, 0, title=labels[orientation]+' Io subtracted'
        cgimage, newimg[index0:index1,*]
        
        window, 1, title=labels[orientation]+' Io subtracted'
        cgplot, wl[index0:index1], total(newimg[index0:index1,*], 2, /nan)
        stop
;        yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]
;        window, 2
;        cgplot, total(newimg[500:1000,*], 1, /nan), xr=[0,44];, xtickformat='(A1)'
;        cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
;        stop
          
; ---------------------------------------------------- Na D1 line --------------------------------------------------------
        
        plota0 = []
        plota1 = []
        plota2 = []
        plota3 = []
        statuses=[]
        
        FOR i = 0, s[2]-1 DO BEGIN
          totsubtrd        =  newimg[*,i]
          guess_peak       =  max(totsubtrd[LSF_fitting_ind1], loc)                                     ; initial guesses [amplitude, peak centroid, hwhm, vertical shift]
          guess_low        =  min(totsubtrd[LSF_fitting_ind1], minloc)
          loc              =  loc + LSF_fitting_ind1[0]
          
          if strmid(labels[orientation],2,1) eq '1' then $
            p = [abs(guess_peak-(median(totsubtrd[index0:index1]))), 5896.25, 0.065, abs(median(totsubtrd[index0:index1]))]
          if strmid(labels[orientation],2,1) eq '2' then $
            p = [abs(guess_peak-(median(totsubtrd[index0:index1]))), 5896.25, 0.080, abs(median(totsubtrd[index0:index1]))]
          
          parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p))
          parinfo.value         = p
          parinfo[0].limited    = [1, 1]                                                              ; 
          parinfo[0].limits     = [abs(guess_low), guess_peak]                                             ;
;          parinfo[1].limited    = [1, 1]                                                              ;
;          parinfo[1].limits     = [loc-0.2, loc+0.2]                                                  ;
          parinfo[1].fixed      = 1                                                                   ;
          parinfo[1].value      = p[1]
;          parinfo[2].limited    = [1, 1]                                                              ;
;          parinfo[2].limits     = [0., 0.075]                                                            ;
          if strmid(labels[orientation],2,1) eq '1' then begin
            parinfo[2].fixed      = 1                                                                   ;
            parinfo[2].value      = 0.065
          endif
          if strmid(labels[orientation],2,1) eq '2' then begin
            parinfo[2].fixed      = 1                                                                   ;
            parinfo[2].value      = 0.09
          endif
          parinfo[3].limited    = [1, 1]                                                              ;
          parinfo[3].limits     = [0.,5000.]                                                            ;
          
          y1                    = totsubtrd[LSF_fitting_ind1]
          D1fa                  = { x:wl[LSF_fitting_ind1], y:y1, err:10.*sqrt(abs(y1)) }
          
          a                     = mpfit('Gaussian_for_MPFIT', p, parinfo=parinfo, funct=D1fa, STATUS = Did_it_work, xtol=5D-9, ftol=1D-6, gtol=1D-9)

          if labels[orientation] eq 'EW on disk' or $
             labels[orientation] eq 'NS on disk' or $
             labels[orientation] eq 'Juno flyby' then begin
              p                     = [abs(guess_peak-(median(totsubtrd[index0:index1]))), 5896.25, 0.06, abs(median(totsubtrd[index0:index1]))]
              parinfo               = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p))
              parinfo.value         = p
              parinfo[0].limited    = [1, 1]                                                              ;
              parinfo[0].limits     = [abs(guess_low), guess_peak]                                        ;
              parinfo[1].limited    = [1, 1]                                                              ;
              parinfo[1].limits     = [WL[loc]-0.25, WL[loc]+0.25]                                                                 ;
;              parinfo[1].fixed      = 1
;              parinfo[1].value      = p[1]
              parinfo[2].limited    = [1, 1]                                                              ;
              parinfo[2].limits     = [0., 0.2]
              parinfo[3].limited    = [1, 1]                                                              ;
              parinfo[3].limits     = [0.,guess_peak]
              a                     = mpfit('Gaussian_for_MPFIT', p, parinfo=parinfo, funct=D1fa, STATUS = Did_it_work, xtol=1D-9, ftol=1D-2)
          endif
          
          
          print, 'STATUS: ',did_it_work
          
          D1_area               = a[0] * a[2] * sqrt(2.*!pi)
          D1_int                = INT_TABULATED(D1fa.x, y1)
          mpfitD1emission       = [mpfitD1emission, D1_area]
          int_D1_emission       = [int_D1_emission, D1_int]
        ; Inspect:
                          window, 2, xs=800, ys=500, title=labels[orientation]+' D1, row '+strcompress(i)
                          cgplot, D1fa.x, D1fa.y, err_yhigh = D1fa.err, err_ylow = D1fa.err
                          cgplot, D1fa.x, gaussian(D1fa.x, a), /overplot, color='red'                                         ; all plots in DN/s
                          stop
                
           
; ---------------------------------------------------- Na D2 line --------------------------------------------------------

          guess_peak       = max(totsubtrd[LSF_fitting_ind2], loc)                                     ; initial guesses [amplitude, peak centroid, hwhm, vertical shift]
          guess_low        = min(totsubtrd[LSF_fitting_ind2], minloc)
          loc              = loc + LSF_fitting_ind2[0]
          if strmid(labels[orientation],2,1) eq '1' then $ 
            p = [abs(guess_peak-(median(totsubtrd[index0:index1]))), 5890.25, 0.068, abs(median(totsubtrd[index0:index1]))]
          if strmid(labels[orientation],2,1) eq '2' then $  
            p = [abs(guess_peak-(median(totsubtrd[index0:index1]))), 5890.25, 0.078, abs(median(totsubtrd[index0:index1]))]

          parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p))
          parinfo.value         = p
          parinfo[0].limited    = [1, 1]                                                              ; 
          parinfo[0].limits     = [abs(guess_low), guess_peak]                                             ;
          parinfo[1].limited    = [1, 1]                                                          ;
          parinfo[1].limits     = [WL[loc]-0.25, WL[loc]+0.25]                                               ;
;          parinfo[1].fixed      = 1                                                                   ;
;          parinfo[1].value      = 5890.25
          parinfo[2].limited    = [1, 1]                                                              ;
          parinfo[2].limits     = [0,0.1]                                                            ;
;          parinfo[2].fixed      = 1                                                                   ;
;          parinfo[2].value      = 0.08
          parinfo[3].limited    = [1, 1]                                                              ;
          parinfo[3].limits     = [0.,5000.]                                                            ;

          y2                    = totsubtrd[LSF_fitting_ind2]
          D2fa                  = { x:wl[LSF_fitting_ind2], y:y2, err:10.*sqrt(abs(y2)) }
          a                     = mpfit('Gaussian_for_MPFIT', p, parinfo=parinfo, funct=D2fa, STATUS = Did_it_work, xtol=1D-9, ftol=1D-2);, gtol=1D10)
          
          if labels[orientation] eq 'EW on disk' or $
             labels[orientation] eq 'NS on disk' or $
             labels[orientation] eq 'Juno flyby' then begin
             
               p                     = [abs(guess_peak-(median(totsubtrd[index0:index1]))), 5890.25, 0.06, abs(median(totsubtrd[index0:index1]))]
               parinfo               = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p))
               parinfo.value         = p
               parinfo[0].limited    = [1, 1]                                                              ;
               parinfo[0].limits     = [abs(guess_low), guess_peak];
               parinfo[1].limited    = [1, 1]                                                          ;
               parinfo[1].limits     = [WL[loc]-0.25, WL[loc]+0.25]
;               parinfo[1].fixed      = 1                                                                   ;
;               parinfo[1].value      = 5890.25
               parinfo[2].limited    = [1, 1]                                                              ;
               parinfo[2].limits     = [0,0.2]
               parinfo[3].limited    = [1, 1]                                                              ;
               parinfo[3].limits     = [0.,3000.]
               a                     = mpfit('Gaussian_for_MPFIT', p, parinfo=parinfo, funct=D2fa, STATUS = Did_it_work, xtol=5D-8, ftol=1D-2)
          endif
           
          print, 'STATUS: ' , did_it_work
          D2_area               = a[0] * a[2] * sqrt(2.*!pi)
          D2_int                = TSUM(D2fa.x, y2)
          mpfitD2emission       = [mpfitD2emission, D2_area]
          int_D2_emission       = [int_D2_emission, D2_int]
;          ; Inspect:
;                  window, 3, xs=800, ys=500, title=labels[orientation]+' D2'
;                  cgplot, D2fa.x, D2fa.y, err_yhigh = D2fa.err, err_ylow = D2fa.err
;                  cgplot, D2fa.x, gaussian(D2fa.x, a), /overplot, color='red'                                         ; all plots in DN/s
;                  stop
                  
                  
          plota0 = [plota0, a[0]]
          plota1 = [plota1, a[1]]
          plota2 = [plota2, a[2]]
          plota3 = [plota3, a[3]]
          statuses = [statuses, did_it_work]
          
        ENDFOR ; each row of ONE orientation
        
        
;        window, 0, title='amplitude a0'
;        cgplot, plota0
;        window, 1, title='center wavelength a1'
;        cgplot, plota1, /ynozero
;        window, 2, title='gaussian width a2'
;        cgplot, plota2, /ynozero
;        window, 3, title='dc offset a3'
;        cgplot, plota3, /ynozero
;        plate_scale = 0.358
;        yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]
;        window, 5, title='comparing gauss fits to data'+labels[orientation]
;        cgplot, mpfitD2emission, /ynozero, title='total fits over all rows', xr=[0,n_elements(mpfitd1emission)],xtickformat='(A1)'       ; units?
;        cgplot, int_D2_emission, /overplot, color='red'                           ; units [R]
;        cgaxis, xaxis=0, xtit='Europa Radii', xr=yr, xstyle=1, xticklen=-0.05
;        window, 6, title='status updates'
;        cgplot, statuses, title='STATUS PER ROW', yr=[0,8], psym=5
;;        
;        stop
        
        plate_scale = 0.358
        yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]
        
              window, 2, title='D1 mpfit'
              cgplot, mpfitD1emission, xr=[0,n_elements(mpfitD1emission)], xtickformat='(A1)', title='black = fitted, red = int'
              cgplot, int_D1_emission, /overplot, color='red'
              cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
              window, 3, title='D2 mpfit'
              cgplot, mpfitD2emission, xr=[0,n_elements(mpfitD2emission)], xtickformat='(A1)', title='black = fitted, red = int'
              cgplot, int_D2_emission, /overplot, color='red'
              cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
              stop
        angstrom_per_pixel = mean(deriv(WL))
        ;
        ;; below, i calculate units of rayleighs to match leblanc (2005) plots.
        ;
        nogaussfit = newimg;[LSF_fitting_ind1, *] + newimg[LSF_fitting_ind2, *]
        collapse   = mean(nogaussfit, dim=1)
        ;
        ;very rough calculations here
        ind = min(abs(wl[500:900] - 5890.2),loc)
        area_profile = mpfitD1emission + mpfitD2emission
        int_profile  = int_D1_emission + int_D2_emission
        
        
              window, 3
              cgplot, area_profile, title = labels[orientation], ytitle = 'Rayleighs', xr=[0,n_elements(area_profile)], xtickformat='(A1)';, yr=[0,max(collapse)], xticklen=0
              cgplot, int_profile, /overplot, color='red'
              cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
              cglegend, colors=['black', 'red'], titles=['Gaussian fitted', 'Sum under data'], length=0.01, symsize=0.1, /Box, Location=[0.15, 0.70], charsize=1.0, /Background, vspace=1
;              save, area_profile, yr, filename= dir+'\Processed\test.sav'
;             stop
              
;              window, 4
;              cgplot, int_profile, title = labels[orientation], ytitle = 'Rayleighs', xr=[0,n_elements(int_profile)], xtickformat='(A1)';, yr=[0,max(collapse)], xticklen=0
;              cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
              
        ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% apples2apples leblanc compare %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        P = cglayout([2,2], ygap = 0., oxmargin = [14, 2], oymargin = [9, 5], xgap = 0.)
        axis_format = {XTicklen:-.01, yticklen:-0.01 }

        cgPS_Open, filename = Dir+'\Figures\falloff_Rayleighs_'+order.name+'_'+labels[orientation]+'.eps', /ENCAPSULATED, xsize = 7.5, ysize = 6
        !P.font=1
        loadct, 3
        device, SET_FONT = 'Helvetica Bold', /TT_FONT

        title = 'HIRES 2022-09-29 : '+labels[orientation]

        cgplot, int_profile, title = labels[orientation], ytitle = 'Emission D1+D2 (Rayleighs)', xr=[0,n_elements(area_profile)], xtickformat='(A1)', xticks=1, xminor=1;, xticklen=-0.1
        cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr*2., xstyle=1, xticklen=-0.02

        cgps_Close
;        stop
        ;---------------------------------------------------- calculate the g-value -------------------------------------------------------
        ; first, will use SPICE to get europa's heliocentric velocity... need observations time stamp so reload headers

        dates = []

        ; need to match up the time stamps to orientations to get the specific g-values
        if filt eq 'gg475' then begin
          if orientation eq 0  then file = 49             ; hires0177
          if orientation eq 1  then file = 56             ; hires0184
          if orientation eq 2  then file = 57             ; hires0185
          if orientation eq 3  then file = 52             ; hires0180
          if orientation eq 4  then file = 53             ; hires0181
          if orientation eq 5  then file = 38             ; hires0166
          if orientation eq 6  then file = 40             ; hires0168
          if orientation eq 7  then file = 42             ; hires0170
          if orientation eq 8  then file = 44             ; hires0172
          if orientation eq 9  then file = 46             ; hires0174
        endif

        filename      = order.name+'_CR_hires' + Europa_frames[orientation] + '.Cleaned.fits'
        europa_header = headfits(dir+'\Processed\Cosmic Rays\'+filename)
        date          = strcompress(sxpar(europa_header, 'DATE'),  /remove_all)

        cspice_str2et, date, et
        dates         = [dates, et]

        target        = 'Europa'
        coframe       = 'J2000'
        abcorr        = 'LT'                                                          ; I found that using LT instead of LT+S gives a more accurate answer based on the NASA JPL Horizons System
        observer      = 'Earth'
        observatory   = 'keck'

        cspice_spkezr, target, et, coframe, abcorr, observer, state, earth_europa_ltime

        timediff      = et - earth_europa_ltime

        cspice_spkezr, target, timediff, coframe, abcorr, 'sun', state, sun_europa_ltime
        obspos        = state[0:2]                                                         ; position of europa in cartesian coordinates
        obsvel        = state[3:5]                                                         ; velocity of europa in cartesian coordinates
        theta         = cspice_vsep(obspos, obsvel)
        sun_europa    = cos(theta) * norm(obsvel)
        distance      = 3.e5 * sun_europa_ltime

        GVALUE, 'Na-D', sun_europa * (10^3), distance / 1.496e8, WL_A[Na_ind-9000:Na_ind+9000], Flux[Na_ind-9000:Na_ind+9000], g_Na
        print, 'g-value for Na D1+D2 using horizons', g_na
        GVALUE,  'K-D', sun_europa * (10^3), distance / 1.496e8, WL_A[k_ind-9000:k_ind+9000], Flux[k_ind-9000:k_ind+9000], g_K
        print,  'g-value for K D1+D2 using horizons', g_K

        column                                               = area_profile * 10.e6 / g_Na
        int_column                                           = int_profile  * 10.e6 / g_Na
               
        if filt eq 'Na' then Nacolumns[*, orientation]       = column
        if filt eq 'Na' then int_Nacolumns[*, orientation]   = int_column
        if filt eq 'gg475' then GGcolumns[*, orientation]    = column
        ;      scale_height = WHERE(column[15:*] eq 3.5E10) - WHERE(column[15:*] eq 3.5E10 / 2.718281828)
        scale_height = 3.                                         ; 3 europa radii scale height
        scale_height = scale_height * 1.5608E8                    ; converts to cm ?
        volume_dens  = column / scale_height
        

        ; Compare raw vs sun-subtracted spectra.
;        xr    = [wl[500],wl[900]]    ;[5887.5, 5897.9]     ;
;        junk  = min(abs(xr[0]- WL), index0)
;        junk  = min(abs(xr[1]- WL), index1)
        WL_xr = [index0, index1]
        
        print, 'average column density for ', labels[orientation], ' is ', mean(column, /nan)
        
        ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3-panel plots w col. dens. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
        P = cglayout([2,2], ygap = 0., oxmargin = [14, 2], oymargin = [9, 5], xgap = 0.)
        axis_format = {XTicklen:-.05, yticklen:-0.01 }

        cgPS_Open, filename = Dir+'\Figures\'+europa_frames[orientation]+' '+labels[orientation]+' 3panel.eps', /ENCAPSULATED, xsize = 7.5, ysize = 6
        !P.font=1
        loadct, 3
        device, SET_FONT = 'Helvetica Bold', /TT_FONT

        title = 'HIRES 2022-09-29 : '+labels[orientation]

        cgplot, wl[500:900], total(newimg[index0:index1, *], 2, /NAN)/1.e4, /xs, xr = [wl[500], wl[900]], pos = p[*,0], xtickformat = '(A1)', $
          title=title, yticklen=-0.02 , ytitle = cgsymbol('times')+'10!U4!N Rayleighs / ' + cgsymbol('Angstrom')

        cgimage, newimg[index0:index1, *], /axes, xr = xr, pos = p[*,2], yr = yr, /noerase, $
          xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', AXKEYWORDS = axis_format
        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01

        cgplot, fltarr(N_elements(newimg[index0:index1, *])), findgen(N_elements(newimg[index0:index1, *])), $
          /ynozero, /noerase, pos = p[*,1], ys= 5, xtickformat="(A1)", ytickformat="(A1)", xstyle=8
;        cgaxis, xr = [min(volume_dens), max(volume_dens)], xaxis=1, xminor=1, charsize=1.3, color='maroon'
;        cgtext, 0.60, 0.64, 'Estimated Volume Density Along Juno Trajectory (atoms / cm!U3!N)', color='maroon', /normal, charsize=1

        cgplot, column/1.e10, findgen(N_elements(newimg[index0:index1, *])), /ynozero, /noerase, pos = p[*,3], ys= 5, $
          xtitle = 'Na Column Density ('+cgsymbol('times')+'10!U10!N atoms / cm!U2!N)', xstyle=9;, xr=[0,400]
                cgaxis, xaxis = 1, xr = [min(area_profile)/1.e3, max(area_profile)/1.e3], charsize=1.5
;                cgaxis, xr = [min(area_profile)/1000., max(area_profile)/1000.], charsize=1, xticklen=-0.03, xminor=1
        cgtext, 0.65, 0.59, 'D1 + D2 ('+cgsymbol('times')+'10!U3!N Rayleighs)', color='black', /normal, charsize=1.5

        cgcolorbar, POSITION=[0.563, 0.541, 0.575, 0.912], range = [0.0,max(column/1.e10)], charsize=1.2, /top, /vertical, /right
        cgtext, 0.63, 0.67, 'Na Col. Dens.', color='black', /normal, charsize=1.5, orientation=90.
        cgtext, 0.66, 0.65, '('+cgsymbol('times')+'10!U10!N atoms / cm!U2!N)', color='black', /normal, charsize=1.3, orientation=90.

        cgps_Close
        
;        ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sun-sub comparison figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
;        ;P  = cglayout([2,2], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
;        plate_scale = 0.358
;        yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]
;
;        cgPS_Open, filename = 'Z:\DATA\Keck\Europa Na\HIRES_20220928\Figures\'+order.name+'_'+labels[orientation]+'_suncomparison.eps', $
;          /ENCAPSULATED, xsize = 10, ysize = 8
;        !P.font=1
;        loadct, 3
;        device, SET_FONT = 'Helvetica Bold', /TT_FONT
;
;        title = 'HIRES 2022-09-29 : '+labels[orientation]
;
;        ; non-sunlight subtracted on the left hand side
;        cgplot, wl, europa[*,suncol], /xs, xr = xr, pos = p[*,0], xtickformat = '(A1)', $
;          title=labels[orientation]+' Raw', ytitle = 'Rayleighs / ' + cgsymbol('Angstrom')
;
;        cgplot, wl, continuum, /overplot, color='red'
;
;        cgtext, .29, .58, 'Above disk', color = 'black', /normal
;        cgtext, .29, .55, 'Fit Reflectance', color = 'red', /normal
;
;        axis_format = {XTicklen:-.01, yticklen:-0.01 }
;        cgimage, orientations[index0:index1,*,orientation], minv=0.75*mean(orientations[*,*,orientation]), maxv=3.0*mean(orientations[*,*,orientation]), $
;          /axes, xr = xr, pos = p[*,2], /noerase, $
;          ytit = 'Europa Radii', xtitle = 'Angstroms', AXKEYWORDS = axis_format, yr = yr
;
;        cgcolorbar, POSITION=[p[0,0], 0.045, p[2,0], 0.065], range = minmax(orientations[*,*,orientation])
;
;        ; sunlight subtracted on the right hand side
;        cgplot, wl, total(newimg, 2, /NAN), /xs, xr = xr, pos = p[*,1], xtickformat = '(A1)', $
;          ytickformat = '(A1)', /noerase, title=labels[orientation]+' Sun-Subtracted'
;        cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
;
;        axis_format = {XTicklen:-.01, yticklen:-0.01, ystyle:5}
;        cgimage, newimg[index0:index1, *], minv=-100, maxv=100, /axes, xr = xr, pos = p[*,3], yr = yr, /noerase, $
;          xtitle = 'Angstroms', AXKEYWORDS = axis_format
;        cgaxis, yaxis = 1, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
;
;        cgcolorbar, POSITION=[p[2,0], 0.045, 0.979, 0.065], range = [-5.e3,5.e3]
;        cgps_Close
;
;        ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5 panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
;        P  = cglayout([3,2], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
;        plate_scale = 0.358
;        yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]
;
;        cgPS_Open, filename = dir+'\Figures\test_figures\'+order.name+'_'+labels[orientation]+'_5panel_process.eps', $
;          /ENCAPSULATED, xsize = 20, ysize = 10
;        !P.font=1
;        loadct, 3
;        device, SET_FONT = 'Helvetica Bold', /TT_FONT
;
;        title = 'HIRES 2022-09-29 : '+labels[orientation]+' Subtraction Process'
;
;        ;   non-sunlight subtracted on the upper left hand side
;        cgimage, europa[index0:index1, *], pos = p[*,0], title=labels[orientation]+' Raw', ytitle = 'Rayleighs / '+cgsymbol('Angstrom')
;        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
;        cgtext, .2, .58, 'Raw', color = 'white', /normal
;
;        ;   sunlight subtracted on the lower left hand side
;        cgimage, eurimg[index0:index1, *], pos = p[*,3], /axes, /noerase, AXKEYWORDS = axis_format, xr = xr, xtitle ='Angstroms'
;        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
;        cgtext, .2, .15, 'Sun Sub', color = 'white', /normal
;
;        ;   sun and io subtracted on lower middle
;        cgimage, newimg[index0:index1, *], /axes, xr = xr, pos = p[*,4], yr = yr, /noerase, $
;          xtitle = 'Angstroms', AXKEYWORDS = axis_format
;        cgtext, .5, .15, 'Sun and Io Sub', color = 'white', /normal
;
;        ;   1D spectra of sun- and io-subbed on the upper middle
;        cgplot, wl, total(newimg[index0:index1, *], 2, /NAN), /xs, xr = xr, pos = p[*,1], xtickformat = '(A1)', $
;          ytickformat = '(A1)', /noerase, title=labels[orientation]+' Sun-Subtracted'
;        cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
;
;        ;   1D spectra of sun- and io-subbed on the lower right hand side
;        cgplot, column/1.e10, findgen(N_elements(newimg[index0:index1, *])), /ynozero, /noerase, pos = p[*,5], ys= 5, $
;          xtitle = 'Na Column Density ('+cgsymbol('times')+'10!U10!N atoms / cm!U2!N)';, xr=[0,400]
;        cgaxis, xaxis = 1, xtit = 'D1 + D2 (Rayleighs)', xr = [min(area_profile), max(area_profile)]
;
;        cgPS_Close
;
;        ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6 panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
        P  = cglayout([3,2], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
        plate_scale = 0.358
        yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]

        cgPS_Open, filename = dir+'\Figures\test_figures\'+europa_frames[orientation]+' '+labels[orientation]+'_6panel_compare.eps', $
          /ENCAPSULATED, xsize = 20, ysize = 10
        !P.font=1
        loadct, 3
        device, SET_FONT = 'Helvetica Bold', /TT_FONT

        title = 'HIRES 2022-09-29 : '+labels[orientation]+' Subtraction Comparisons'

        ;   raw on the left hand side
        cgimage, europa[index0:index1, *], pos = p[*,0]
        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01, charsize=1.8
        cgtext, 0.25, 0.85, 'Basic redux', align = 0.5, /normal, charsize=2, color='white'
        ;   raw 1D
        cgplot, wl[index0:index1], total(europa[index0:index1, *], 2, /NAN)/1.e6, pos = p[*,3], xr = xr, /xs, ytickformat = '(A1)', xtickformat = '(A1)', /noerase
;        cgplot, wl[index0:index1], sunrow[index0:index1]/max(total(europa[index0:index1, *], 2, /NAN)), /overplot, color='orange'
;        cglegend, colors=['black', 'orange'], psym=[0,0], titles=['Raw spectrum', 'Fake sun spectrum'], length=0.01, symsize=0.1, /Box, Location=[0.11, 0.50], charsize=1.0, /Background, vspace=1
        cgaxis, yaxis = 0, ytitle = '10!U6!N Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1, charsize=1.8
        cgaxis, xaxis = 0, xtitle = 'Wavelength (' + cgsymbol('Angstrom') + ')', xstyle=1, charsize=1.8
;        cgcolorbar, POSITION=[0.105, 0.905, 0.369, 0.92], range = [min(europa[index0:index1, *])/1.e5,max(europa[index0:index1, *])/1.e5], charsize=1.8, /top, title='Na Column Density (10!U5!N Rayleighs / '+ cgsymbol('Angstrom')+')'

        ;   sunlight subtracted in the middle
        cgimage, sunsubbed[index0:index1, *], pos = p[*,1], /noerase, title='Sun-Sub'
        ;   sun sub 1D
        cgplot, wl[index0:index1], total(sunsubbed[index0:index1, *], 2, /NAN)/1.e4, pos = p[*,4], xr = xr, /xs, ytickformat = '(A1)', xtickformat = '(A1)', /noerase
        cgplot, wl[index0:index1], fltarr(index1-index0), /overplot, color='red', linestyle=2
        cgtext, 0.51, 0.85, 'Sun-subtracted', align = 0.5, /normal, charsize=2, color='white'
        cgaxis, xaxis = 0, xtitle = 'Wavelength (' + cgsymbol('Angstrom') + ')', xstyle=1, charsize=1.8
        cgcolorbar, POSITION=[0.369, 0.905, 0.630, 0.92], range = [min(sunsubbed[index0:index1, *])/1.e3,max(sunsubbed[index0:index1, *])/1.e3], charsize=1.8, /top, title='Na Column Density (10!U3!N Rayleighs / '+ cgsymbol('Angstrom')+')'
        
        ;   sun and io subtracted on the right hand side
        cgimage, newimg[index0:index1, *], pos = p[*,2], yr = yr, /noerase, title='Sun-Io-Sub'
        ;   sun io sub 1D
        cgplot, wl[index0:index1], total(newimg[index0:index1, *], 2, /NAN)/1.e4, pos = p[*,5], ytickformat = '(A1)', xtickformat = '(A1)', /noerase, xr = xr, /xs;, yr=[min(total(eurimg, 2)), max(total(eurimg, 2))+1.e3]
        cgplot, wl[index0:index1], fltarr(index1-index0), /overplot, color='red', linestyle=2
        cgaxis, yaxis = 1, ytitle = '10!U4!N Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1, charsize=1.8
        cgtext, 0.78, 0.85, 'Io-subtracted', align = 0.5, /normal, charsize=2, color='white'
        cgaxis, xaxis = 0, xtitle = 'Wavelength (' + cgsymbol('Angstrom') + ')', xstyle=1, charsize=1.8
        cgcolorbar, POSITION=[0.630, 0.905, 0.895, 0.92], range = [min(newimg[index0:index1, *])/1.e3,max(newimg[index0:index1, *])/1.e3], charsize=1.8, /top, title='Na Column Density (10!U5!N Rayleighs / '+ cgsymbol('Angstrom')+')'
        
        cgPS_Close
        
        if europa_frames[orientation] eq '0140' then begin
          P = cglayout([2,2], ygap = 0., oxmargin = [14, 2], oymargin = [9, 5], xgap = 0.)
          cgPS_Open, filename = dir+'\Figures\test_figures\juno_flyby_spectrum.eps', /ENCAPSULATED, xsize=7.5, ysize=6
          !P.font=1
          loadct, 3
          device, SET_FONT = 'Helvetica Bold', /TT_FONT

          title = 'HIRES 2022-09-29 : '+labels[orientation]

          
          cgimage, newimg[index0:index1, *], /axes, xr = xr, pos = p[*,2], yr = yr, /noerase, $
             xtitle = 'Wavelength ('+cgsymbol('Angstrom')+')', AXKEYWORDS = axis_format
          cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
          cgcolorbar, POSITION=[0.563, 0.165, 0.575, 0.537], range = [0.0,max(column/1.e10)], charsize=1.2, /top, /vertical, /right
          cgtext, 0.62, 0.17, 'Na Col. Dens. ('+cgsymbol('times')+'10!U10!N atoms / cm!U2!N)', color='black', /normal, charsize=1, orientation=90.
          
          cgPS_Close
        endif 
        
;        if labels[orientation] eq 'EW on disk' then save, /all, filename = Dir+'\Processed\'+europa_frames[orientation]+'_EWondisk_rayleighs.sav'
;        if labels[orientation] eq 'NS on disk' then save, /all, filename = Dir+'\Processed\'+europa_frames[orientation]+'_NSondisk_rayleighs.sav'
;        if labels[orientation] eq 'EW10 N' then save, /all, filename = Dir+'\Processed\EW10N_rayleighs.sav'
        
        if filt eq 'gg475' then gg475_new_images[*,*,orientation] = newimg
        if filt eq 'Na'    then    Na_new_images[*,*,orientation] = newimg
      endfor ; orientation
    endfor ; order
    save, /all, filename = Dir+'\Processed\sun_subbed_europa.sav'
    stop
  endif                               ; part 1.6


  if part eq 2 then begin                                                         ; mapping brightnesses around europa in grid structure
    RESTORE, Dir+'\Processed\sun_subbed_europa.sav'

    guiders          = FILE_SEARCH(dir+'\MAGIQ files', 'hiresslit*'+'*.fits')
    spectra          = FILE_SEARCH(dir+'\Processed\Cosmic Rays', 'order_60_*'+'*CR.fits')
    guider_times     = []
    end_times_et     = []
    start_times_et   = []
    exp_times        = []

    ; first, get the time stamps for the guider frames, convert to et
    FOR i = 0, N_elements(guiders)-1 DO BEGIN
      guider_header  = headfits(guiders[i])
      cspice_UTC2ET,   sxpar(guider_header, 'DATE-OBS') + 'T'+ sxpar(guider_header, 'UTC'), guider_ET
      guider_times   = [guider_times, guider_ET]

      ;    guider_frame   = mrdfits(guiders[i], 0, guider_header, /fscale)
      ;    window, 4, xs=512, ys=512, title=strmid(guiders[i], 50)
      ;    cgimage, guider_frame, minv=0.75*mean(guider_frame), maxv=1.5*mean(guider_frame)
      ;    stop

    ENDFOR
    ; now i'm going to unpack one just so that i can get image dimensions and put them into an array later
    guider_img = mrdfits(guiders[0], 0, guider_header, /fscale)
    s          = size(guider_img)


    ; now, get time stamps for spectra, convert to et AND get exposure times
    FOR j = 0, N_elements(spectra)-1 DO BEGIN
      spectra_header = headfits(spectra[j])
      cspice_UTC2ET,   sxpar(spectra_header, 'DATE'), spectra_ET
      exposure       = sxpar(spectra_header, 'EXPTIME')
      exp_times      = [exp_times, exposure]
      end_times_et   = [end_times_et, spectra_ET]
      start_times_et = [start_times_et, spectra_ET - exposure]
    ENDFOR

    mid_exposure_times                = (start_times_et + end_times_et ) / 2.

    magiq_per_frame                   = fltarr(s[1], s[2], n_elements(guiders)); + !Values.F_NaN
    dummy_cube                        = fltarr(s[1], s[2], n_elements(guiders))
    keepindex                         = []

    NSshort                           = [239., 242.]
    EWshort                           = [236., 244.]
    NSlongs                           = [191., 356.]
    EWlongs                           = [191., 356.]

    ; after rotation
    NSshort                           = [268., 275.]
    EWshort                           = [236., 243.]
    NSlongs                           = [158., 321.]
    EWlongs                           = [156., 319.]

    slit                              = fltarr(NSshort[1]-NSshort[0]+1., NSlongs[1]-NSlongs[0]+1.)
    slitsize                          = size(slit)


    ;;;; this is for Na filtered data
    if filt eq 'Na' then begin
      FOR run = 0, 1 DO BEGIN                                                                   ; TWO runs: one with column = area under D-line gaussians, second with column = integral of the D-line gaussians
        magiq_per_frame                   = fltarr(s[1], s[2], n_elements(guiders)); + !Values.F_NaN
        dummy_cube                        = fltarr(s[1], s[2], n_elements(guiders))
        keepindex                         = []
        FOR file = 3, 37-1 DO BEGIN  ; juno frames: 11,37-1 do begin                                                           ; these are the Na filter files
          if file eq sunspectrum then continue
  
          NSslit_locations                  = fltarr(s[1], s[2]) 
          EWslit_locations                  = fltarr(s[1], s[2]) 
          junoslit_shift                    = fltarr(s[1], s[2]) 
  
          slit_location_Juno                = ROT(NSslit_locations, 316)
          slit_location_344                 = ROT(NSslit_locations, 350)
  
          h = europa_frames[file]
          
  
;          if (h ge 130) and (h le 132)   then frame = 0  ; $
;;            or (h ge 162) and (h le 164) then frame = 0                                         ;  EW on disk                                                  orientations[*,*,0 ] = EW0__
;          if (h ge 133) and (h le 134)   then frame = 1                                         ;  EW 10 N                                                     orientations[*,*,1 ] = EW10_N
;          if (h ge 135) and (h le 136)   then frame = 2                                         ;  EW 20 N                                                     orientations[*,*,2 ] = EW20_N
;          if (h ge 137) and (h le 141)   then frame = 10                                        ;  juno                                                        orientations[*,*,3 ] = EW10_S
;          if (h ge 142) and (h le 143)   then frame = 3                                         ;  EW 10 S                                                     orientations[*,*,4 ] = EW20_S
;          if (h ge 144) and (h le 145)   then frame = 4                                         ;  EW 20 S                                                     orientations[*,*,5 ] = NS0__
;          if (h ge 146) and (h le 150)   then frame = 5  ; $                                                                                                   orientations[*,*,6 ] = NS10_W
;;            or (h ge 159) and (h le 161) then frame = 5                                         ;  NS on disk                                                  orientations[*,*,7 ] = NS20_W
;          if (h ge 151) and (h le 152)   then frame = 6                                         ;  NS 10 W                                                     orientations[*,*,8 ] = NS10_E
;          if (h ge 153) and (h le 154)   then frame = 7                                         ;  NS 20 W                                                     orientations[*,*,9 ] = NS20_E
;          if (h ge 155) and (h le 156)   then frame = 8                                         ;  NS 10 E                                                     if filt eq 'Na' then orientations[*,*,10] = Juno
;          if (h ge 157) and (h le 158)   then frame = 9                                         ;  NS 20 E
;  
;          orientation  = Na_new_images[*,*,frame]

          
          
          
          if run eq 0 then column       = Nacolumns[*, file]
          if run eq 1 then column       = int_Nacolumns[*, file]
          
          slitfiller1d = CONGRID(column, slitsize[2], slitsize[3], /interp)
          slitfiller2d = REBIN(slitfiller1d, slitsize[2], slitsize[1])
          if column[0] eq 0. then slitfiller2d = make_array(slitsize[2], slitsize[1], value=!Values.F_Nan)
          
              loadct, 3
              window, 0
              cgplot, total(slitfiller2d, 2)
  
;              slitfiller2d  = fltarr(slitsize[2], slitsize[1]) +!Values.F_NaN                       ; comment this line out if you want to include col dens map
             
          NSslit_locations[NSshort[0]:NSshort[1],NSlongs[0]:NSlongs[1]] = TRANSPOSE(slitfiller2d)
          EWslit_locations[EWlongs[0]:EWlongs[1],EWshort[0]:EWshort[1]] = slitfiller2d
  
          window, 2, xs=512, ys=512, title=h+' '+labels[file]
          cgimage, EWslit_locations
;          stop
          junoslit_shift[NSshort[0]:NSshort[1],NSlongs[0]:NSlongs[1]]     = TRANSPOSE(slitfiller2d)
          slit_location_Juno                = ROT(junoslit_shift, 316, 1, mean(NSshort), mean(NSlongs), /pivot)
          slit_location_344                 = ROT(NSslit_locations, 350)
          
          
          print, labels[file]
  
          cspice_ET2UTC, start_times_et[file], "ISOC", 2, start_times
          date = sxpar(headfits(spectra[file]), 'DATE')
          print, spectra[file] + ' ' + sxpar(headfits(spectra[file]), 'DATE-OBS') + ' START: '+  STRMID(start_times, 11) + '      END: '+ STRMID(date, 11)
          print, 'has MAGIQ files..............................'
          cspice_str2et, date, et
  
          FOR k = 0, N_elements(guider_times)-1 DO BEGIN
            IF guider_times[k] LT end_times_et[file] AND $
              guider_times[k] GT start_times_et[file] then begin
              print, strmid(guiders[k],50)+ ' ' + sxpar(headfits(guiders[k]), 'DATE-OBS') + ' '+ sxpar(headfits(guiders[k]), 'UTC')+'   ', k
  
              cspice_ET2UTC, guider_times[k], "ISOC", 2, checking
              print, 'CHECK TIME: ' + STRMID(checking, 11)
  
              if k eq 158 then continue                                                              ; idk what this frame is but it's not 10 S like the log says
  
              magiq_per_frame[*,*,k] = mrdfits(guiders[k], 0, guider_header, /fscale)              ; saves all the magiq files per spectral observation
              print, round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.)
  
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 334 then begin          ; NS Oriented
                magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
                magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 1)
                dummy_cube[*,*,k]      = rotate(dummy_cube[*,*,k], 5)
                dummy_cube[*,*,k]      = rotate(dummy_cube[*,*,k], 1)
              endif
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 244 then begin          ; EW Oriented
                magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
                dummy_cube[*,*,k]      = rotate(dummy_cube[*,*,k], 5)
              endif
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 245 then begin          ; EW Oriented
                magiq_per_frame[*,*,k] = rotate(magiq_per_frame[*,*,k], 5)
                dummy_cube[*,*,k]      = rotate(dummy_cube[*,*,k], 5)
              endif
  
              window, 0, xs=512, ys=512, title=STRMID(spectra[file], 62, 9)+' --> '+strmid(guiders[k], 50, 17)
              cgimage, magiq_per_frame[*,*,k], minv=0.75*median(magiq_per_frame[*,*,k]), maxv=1.5*median(magiq_per_frame[*,*,k])
  
  
              ; let's do centroids now
  
              maxes = []
              ylocs = []
  
              for cols = 0, s[2]-1 do begin
                colsd = magiq_per_frame[*,cols,k]
                maxes = [maxes, max(colsd)]
              endfor
  
              centroid_loc  = max(maxes, yloc)
              xloc = WHERE(magiq_per_frame[*,yloc,k] eq max(magiq_per_frame[*,yloc,k]))
  
              CNTRD, magiq_per_frame[*,*,k], xloc, yloc, xcen, ycen, 100
  
              if xcen[-1] eq -1. or ycen[-1] eq -1. then begin
                magiq_per_frame[*,*,k] = fltarr(s[1], s[2])
                continue
              endif
  
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 245 then magiq_per_frame[*,*,k] = EWslit_locations   + magiq_per_frame[*,*,k]         ; for some EW frames, we used 244.5 and some were 244
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 244 then magiq_per_frame[*,*,k] = EWslit_locations   + magiq_per_frame[*,*,k]
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 334 then magiq_per_frame[*,*,k] = NSslit_locations   + magiq_per_frame[*,*,k]
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq  44 then magiq_per_frame[*,*,k] = slit_location_Juno + magiq_per_frame[*,*,k]         ; juno flyby byebye
  
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 245 then dummy_cube[*,*,k] = EWslit_locations   + dummy_cube[*,*,k]         ; for some EW frames, we used 244.5 and some were 244
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 244 then dummy_cube[*,*,k] = EWslit_locations   + dummy_cube[*,*,k]
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq 334 then dummy_cube[*,*,k] = NSslit_locations   + dummy_cube[*,*,k]
              if round(sxpar(headfits(spectra[file]), 'ROTPOSN') - 90.) eq  44 then dummy_cube[*,*,k] = slit_location_Juno + dummy_cube[*,*,k]         ; juno flyby byebye
  
              shift_in_x = s[1]/2. - xcen
              shift_in_y = s[2]/2. - ycen
  
              magiq_per_frame[*,*,k] = shift(magiq_per_frame[*,*,k], shift_in_x, shift_in_y)
              dummy_cube[*,*,k]      = shift(dummy_cube[*,*,k], shift_in_x, shift_in_y)
  
              if median(magiq_per_frame[*,*,k]) ne !Values.F_NaN then keepindex = [keepindex, k]
  
              window, 3, xs=512, ys=512, title='centroid ; '+labels[file]+ '  ' + strmid(spectra[file], 77,4)
              cgimage, magiq_per_frame[*,*,k], minv=0.75*median(magiq_per_frame[*,*,k]), maxv=1.5*median(magiq_per_frame[*,*,k])
              
            ENDIF             ; guider frame & spectra matching
          ENDFOR              ; guider frames
  
          ;      layered = median(magiq_per_frame, dim=3, /even)
          layered = total(magiq_per_frame, 3)
  
          loadct, 3
          window, 1, xs=512, ys=512, title='TOTALED MAGIQ FRAMES FOR THIS SPECTRUM'
          cgimage, layered, minv=0, maxv=5.e4 ; minv=0.75*median(layered), maxv=1.5*median(layered)
          cgcolorbar, POSITION=[0.7, 0.35, 0.72, 0.65], range = [min(layered), max(layered)], /vertical, /right, color='white';, Format='(F0.2)'
        endfor ; h is the number of sodium filtered spectra
  
        dummy_cube = dummy_cube[*,*,keepindex]
  
        zeros   = WHERE(dummy_cube eq 0.)
        dummy_cube[zeros] = !Values.F_Nan
  
        justmap = mean(dummy_cube, dim=3, /nan)
        yr = [-8.*(slitsize[1]*plate_scale/ang_radius), 8.*(slitsize[1]*plate_scale/ang_radius)]      ; hack
        
        window, 2, xs=512, ys=512, title='JUST MAP'
        cgimage, justmap, minv=0, maxv=3.e10
        cgcolorbar, POSITION=[0.7, 0.35, 0.72, 0.65], range = [min(justmap), max(justmap)], /vertical, /right, color='white'
        cgaxis, 40., 40., color='white', /data, yaxis=0, yr = yr
        cgaxis, 40., 40., color='white', /data, xaxis=0, xr = yr
        print, '    '
        
        if run eq 0 then begin
          cgPS_Open, filename = dir+'\Figures\Na_guider_map_Juno.eps', /ENCAPSULATED, xsize = 10, ysize = 10
          !P.font=1
          loadct, 3
          device, SET_FONT = 'Helvetica Bold', /TT_FONT
    
          title = 'HIRES 2022-09-29 Na filter Slit Orientations Around Europa'
          ;      cgimage, layered, minv=0, maxv=5.e4
          cgimage, justmap, minv=0, maxv=3.5e10;, minv=0.75*median(layered), maxv=1.5*median(layered)
          cgcolorbar, POSITION=[0.7, 0.35, 0.72, 0.65], range = [min(justmap), max(justmap)]/1.e10, /vertical, /right, color='white', $
            title='Na Column Density (x10!U10!N)', charsize=1.8
;          cgaxis, 0.2, 0.3, color='white', /data, yaxis=0, yr = yr
;          cgaxis, 0.4, 0.2, color='white', /data, xaxis=0, xr = yr
    
    
          cgPS_Close
          
          save, justmap, filename = Dir+'\Processed\Na_map_junoflyby.sav'
          stop
        endif
        
        if run eq 1 then begin
          cgPS_Open, filename = dir+'\Figures\Na_guider_map_Juno_INTEGRATED.eps', $
            /ENCAPSULATED, xsize = 10, ysize = 10
          !P.font=1
          loadct, 3
          device, SET_FONT = 'Helvetica Bold', /TT_FONT

          title = 'HIRES 2022-09-29 Na filter Slit Orientations Around Europa'
          ;      cgimage, layered, minv=0, maxv=5.e4
          cgimage, justmap, minv=0, maxv=3.5e10;, minv=0.75*median(layered), maxv=1.5*median(layered)
          cgcolorbar, POSITION=[0.7, 0.35, 0.72, 0.65], range = [min(justmap), max(justmap)]/1.e10, /vertical, /right, color='white', $
              title='Na Column Density (x10!U10!N)', charsize=1.8
          cgaxis, 40., 100., color='white', /data, yaxis=0, yr = yr
          cgaxis, 250., 40., color='white', /data, xaxis=0, xr = yr


          cgPS_Close
        endif
        stop
      ENDFOR                        ; both runs (area under gaussian curve fitted to D lines vs. int_tabulated area under D lines (NO GAUSS FIT))
    endif    ; Na filter data
    stop



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
            stop
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
            stop
          ENDIF
        ENDFOR
        goodframes = []
        layered = TOTAL(magiq_per_frame, 3)
        window, 1, xs=512, ys=512, title='TOTALED MAGIQ FRAMES FOR THIS SPECTRUM'
        cgimage, layered, minv=0.75*median(layered), maxv=1.5*median(layered)
        print, '    '
        stop
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



    stop

  endif
  ;
  ;
  ;        for order = 0, n_elements(orders)-1 do begin
  ;
  ;          if order eq 2 then begin ; Na order_60
  ;
  ;
  ;
  ;          endif
  ;
  ;          ;
  ;          ;            if order eq 0 then begin ; K order_46
  ;          ;
  ;          ;              axis_format = {XTicklen:-.01, yticklen:-0.005 }
  ;          ;
  ;          ;
  ;          ;              P = cglayout([1,2]);, ygap = 0., oxmargin = [12, .9], oymargin = [9, 5])
  ;          ;              cgPS_Open, filename = 'C:\Users\elovett\EuropaResearch\Europa_Flyby\3 Panel Plots\'+orders_to_calibrate[order]+'_'+labels[orientation]+'.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5
  ;          ;              !P.font=1
  ;          ;              device, SET_FONT = 'Helvetica Bold', /TT_FONT
  ;          ;              cgplot, WL, newimg[*,suncol], /xs, xr = [wl[0], wl[s[2]]], pos = p[*,0], xtickformat = '(A1)', $
  ;          ;                 ytitle = 'Rayleighs / ' + cgsymbol('Angstrom')
  ;          ;
  ;          ;              cgimage, newimg, minv=0.5*mean(newimg), maxv=1.5*mean(newimg), /axes, xr = [wl[0], wl[1200]], pos = p[*,1], yr = [-(14.*.5/ang_radius), (14.*.5/ang_radius)]-.7, /noerase, $
  ;          ;                ytit = 'Europa Radii', xtitle = 'Angstroms', AXKEYWORDS = axis_format
  ;          ;              cgps_Close
  ;          ;            endif
  ;
  ;
  ;        ENDFOR ; each order within one orientation
  ;
  ;
  ;
  ;
  ;    loadct, 0
  ;
  ;    WINDOW, 1, XS = 3500, YS = 1000
  ;    cgimage, _10_RsubE_east_NS, minval = 0, maxval = 1., /axes, xr = [wl[0], wl[1200]]
  ;
  ;    S = SIZE(_10_RsubE_east_NS, /DIM)
  ;
  ;    wl = double(WL)
  ;
  ;    _10_RsubE_east_NS = double(_10_RsubE_east_NS)
  ;
  ;    window, 2
  ;
  ;    K_Fit_Params  = {brightness:     fltarr(2, s[1]), $
  ;                     err_brightness: fltarr(2, s[1]), $
  ;                     linewidth:      fltarr(2, s[1]), $
  ;                     linecenter:     fltarr(2, s[1])}
  ;
  ;    ; K D2
  ;      for i = 0, s[1] - 1 do begin
  ;        LSF_fitting_ind  = where( abs(wl- 7665.33) lt 1., /NULL)
  ;
  ;        scale = mean(_10_RsubE_east_NS[LSF_fitting_ind, i] / sunlight)
  ;        y = _10_RsubE_east_NS[LSF_fitting_ind, i] - scale*sunlight
  ;        fa               = { x:wl[LSF_fitting_ind], y:y, err:sqrt(y/30.) }
  ;        a                = mpfit('Gaussian_for_MPFIT', [.3, 7665.33, 0.1], funct=fa, STATUS = Did_it_work, /quiet)
  ;        a = abs(a)
  ;
  ;        ; Inspect:
  ;          cgplot, fa.x, fa.y, err_yhigh = fa.err, err_ylow = fa.err
  ;          cgplot, fa.x, gaussian(fa.x, a), /overplot
  ;
  ;          K_Fit_Params.brightness[0, i]     = A[0]*A[2]*SQRT(2*!DPI)
  ;          K_Fit_Params.err_brightness[0, i] = SQRT(K_Fit_Params.brightness[0, i] )
  ;          K_Fit_Params.linewidth[0, i]      = A[2]
  ;          K_Fit_Params.linecenter[0, i]     = A[1]
  ;         ; stop
  ;      endfor
  ;
  ;
  ;    ; K D1
  ;      for i = 0, s[1] - 1 do begin
  ;        LSF_fitting_ind  = where( abs(wl- 7699.43) lt 1., /NULL)
  ;        fa               = { x:wl[LSF_fitting_ind], y:_10_RsubE_east_NS[LSF_fitting_ind, i] - mean(_10_RsubE_east_NS[*, i]), err:sqrt(_10_RsubE_east_NS[LSF_fitting_ind, i])/30. }
  ;        a                = mpfit('Gaussian_for_MPFIT', [.3, 7699.43, 0.1], funct=fa, STATUS = Did_it_work, /quiet)
  ;        a = abs(a)
  ;
  ;        ; Inspect:
  ;        cgplot, fa.x, fa.y, err_yhigh = fa.err, err_ylow = fa.err
  ;        cgplot, fa.x, gaussian(fa.x, a), /overplot
  ;
  ;        K_Fit_Params.brightness[1, i]     = A[0]*A[2]*SQRT(2*!DPI)
  ;        K_Fit_Params.err_brightness[1, i] = SQRT(K_Fit_Params.brightness[0, i] )
  ;        K_Fit_Params.linewidth[1, i]      = A[2]
  ;        K_Fit_Params.linecenter[1, i]     = A[1]
  ;      endfor
  ;
  ;    window, 3
  ;    cgplot, K_Fit_Params.brightness[0, *]
  ;    cgplot, K_Fit_Params.brightness[1, *], /overplot
  ;    cgplot, .05*mean( _10_RsubE_east_NS , dim = 1), /overplot

  stop
  ;    ENDFOR ; orientation

  if part eq 2 then begin ; mapping observations in grid structure, finding distance of slit from europa's disk

    frames = Dir+'\hires'+strcompress(Europa_frames, /rem)+'.fits'

    ets = []
    amass = []
    distances = []

    for frame = 0, n_elements(frames) - 1 do begin
      array = fltarr(2139, 4096)
      array = [mrdfits(frames[frame], 3, header, /fscale, /silent), $
        mrdfits(frames[frame], 2, header, /fscale, /silent), $
        mrdfits(frames[frame], 1, header, /fscale, /silent)]
      junk  =  mrdfits(frames[frame], 0, header, /fscale, /silent)

      filter  = strcompress(sxpar(header, 'FIL2NAME'), /remove_all)
      if filter ne 'clear' then continue
      date    = strcompress(sxpar(header, 'DATE'),  /remove_all)
      keckra  = ten(strcompress(sxpar(header, 'RA'),  /remove_all) )
      keckdec = ten(strcompress(sxpar(header, 'DEC'),  /remove_all))

      cspice_str2et, date, et

      ; there is no position angle in headers. another way to check?
      ; https://koa.ipac.caltech.edu/UserGuide/WCS/wcs.html says PA = ROTPOSN + 90 to get NORTH

      posang = float(strcompress(sxpar(header, 'ROTPOSN'),/remove_all)) - 90.     ; off by 10 degrees... maybe our log is wrong. idk why its rotposn - 90.

      ; now use SPICE to get predicted location of europa

      target      = 'Europa'
      coframe     = 'J2000'
      abcorr      = 'LT'                                                          ; I found that using LT instead of LT+S gives a more accurate answer based on the NASA JPL Horizons System
      observer    = 'Earth'
      observatory = 'keck'

      cspice_spkezr, target, et, coframe, abcorr, observer, state, ltime
      obspos = state[0:2]                                                         ; position of europa in cartesian coordinates
      obsvel = state[3:5]                                                         ; velocity of europa in cartesian coordinates

      distance = 3.e5 * ltime  ; d (km) = rate (km/s) * time (s)
      radius = 1560.80    ; km

      halfang_rad = atan( radius / distance )   ; radians
      halfang_arc = halfang_rad * 206265.       ; arcsecs
      halfang_deg = halfang_arc / 3600.         ; degrees
      diameter    = 2. * halfang_deg            ; angular diameter in degrees
      europarad   = diameter / 2.

      N_times  = size(et, /N_El)

      ; Get the body-fixed non-inertial IAU_Earth coordinates of the observatory in cartesian
      cspice_bodvrd, observer, 'RADII', 3, radii
      flat = (radii[0] - radii[2])/radii[0]
      OBSERVATORY_cs, observatory, obs_struct                                      ; longitude is in degrees *West*, need to switch to East Longit for J2000 conversion via cspice_georec
      cspice_georec, (360. - obs_struct.longitude) * CSPICE_RPD(), $
        obs_struct.latitude * CSPICE_RPD(), obs_struct.altitude / 1.e3, radii[0], flat, obs_IAU_Earth

      ; Convert to IAU_Earth to J2000
      cspice_pxform, 'IAU_EARTH', coframe, et, Earth_to_J2000_transform_matrix
      obs_J2000 = MAKE_ARRAY(3, N_times, /DOUBLE)
      FOR i=0, N_times-1 DO obs_J2000[*,i] = TRANSPOSE(Earth_to_J2000_transform_matrix[*,*,i]) # obs_IAU_Earth

      cspice_spkpos, target, et, coframe, abcorr, observer, Europa_Earth_vector, ltime
      cspice_recrad, (Europa_Earth_vector - obs_J2000), Europa_distance_WRT_obs, predra, preddec

      predra = ten(predra * !radeg ) / 15.
      preddec= ten(preddec * !radeg)

      diffra   = predra - keckra
      diffdec  = preddec - keckdec
      differ   = ((cos(preddec * !pi / 180.)*diffra)^2 + diffdec^2)^0.5  * 3600.   ; distance from europa's disk in arcseconds
      ;print, frames[frame], differ
      ;if differ gt 20 then continue

      dif_radii = differ / halfang_arc                                             ; converts distance to europa radii

      distances = [distances, dif_radii]

      ets = [ets, et]
      ;amass = [amass, strcompress(sxpar(header, 'AIRMASS'),/remove_all)]

      ;    if frames[frame] eq dir+'hires0127.fits' then stop
      ;    if frames[frame] eq dir+'hires0136.fits' then stop
      ;    if frames[frame] eq dir+'hires0145.fits' then stop
      ;    if frames[frame] eq dir+'hires0155.fits' then stop
      ;    if frames[frame] eq dir+'hires0164.fits' then stop

      if posang lt 245 and posang gt 243.99 then begin
        if dif_radii lt 2 then print, frames[frame], '   On-disk ', dif_radii
        if dif_radii gt 2 and diffdec lt 0 then print, frames[frame], '   N ', dif_radii
        if dif_radii gt 2 and diffdec gt 0 then print, frames[frame], '   S ', dif_radii
      endif

      if posang lt 335 and posang gt 333.99 then begin
        if dif_radii lt 2 then print, frames[frame], '   On-disk ', dif_radii
        if dif_radii gt 2 and diffra lt 0 then print, frames[frame], '   E ', dif_radii
        if dif_radii gt 2 and diffra gt 0 then print, frames[frame], '   W ', dif_radii
      endif
      ;print, files[file], date, dif_radii

    endfor ;files

    ;    differences = (shouldbe - distances)
    ;  window, 1, xs=900, ys=900
    ;  cgplot, ets, amass, yr=[min(amass),max(amass)], xr=[min(ets),max(ets)], $
    ;    title='Airmass Over the Night of 9/29/2022', xtitle='Ephemeris Time [s]', ytitle='Airmass', color='red', psym=-46, symcolor= 'crimson', symsize=2, linestyle=0
    ;
    ;  window, 2, xs=900, ys=900
    ;  cgplot, ets, differences, yr=[min(differences),max(differences)], xr=[min(ets),max(ets)], $
    ;    title='Offness', xtitle='Ephemeris Time [s]', ytitle='Should-Be - Actual Distance', color='blue', psym=-46, symcolor= 'sky blue', symsize=2, linestyle=0
    ;  cgplot, /overplot, ets, fltarr(N_elements(ets)), linestyle=2, symthick=2

    ;    window, 1, xs=900, ys=900
    ;    cgplot, string(indgen(n_elements(europa_frames))+127, format='(I4.4)'), amass, yr=[min(amass),max(amass)], xr=[127,164], $
    ;      title='Airmass Over the Night of 9/29/2022', xtitle='File #', ytitle='Airmass', color='red', psym=-46, symcolor= 'crimson', symsize=2, linestyle=0

    ;    window, 2, xs=900, ys=900
    ;    cgplot, string(indgen(n_elements(europa_frames))+127, format='(I4.4)'), differences, yr=[min(differences),max(differences)], xr=[127,164], $
    ;      title='Offness', xtitle='File #', ytitle='Should-Be - Actual Distance', color='blue', psym=-46, symcolor= 'sky blue', symsize=2, linestyle=0
    ;    cgplot, string(indgen(n_elements(europa_frames))+127, format='(I4.4)'), fltarr(38), linestyle=2, symthick=2, /overplot


    stop
  endif ; part 2 end

end