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
  fitted = p[0]*exp(-1.*z^2) + abs(p[3])
  return, (y - fitted)/err
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
  return, abs(y - fit)/err
end ; matched with "scale_fit_sunlight" below

FUNCTION scale_fit_sunlight, p, x
  return, P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]
end

PRO Keck_Europa_Flyby, part = part, dir = dir

  case dir of
    'Z:\DATA\Keck\Europa Na\HIRES_20220928': begin
      ;Europa           = '1003517' ; aka C/2017 K2 (PANSTARRS)
      biases          = string(indgen(10)+4, format='(I4.4)')
      Flats           = string(indgen(10)+15, format='(I4.4)')
      Lamps           = string(indgen(5)+25, format='(I4.4)')
      Star_Frames     = string(indgen(2)+78, format='(I4.4)')
      Europa_frames   = string(indgen(5)+137, format='(I4.4)')                        ; these are the frames for the Juno flyby, Na filter only though...
      Europa_frames   = [Europa_frames, string(indgen(31)+165, format='(I4.4)')]      ;this is for K and Na data; for just Na, use string(indgen(38)+127, format='(I4.4)')
      Jupiter_frames  = string(116, format='(I4.4)')                 ; Jupiter disk center post eclipse
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
    ;         window, 0, xs = 4096, ys = 2200
    ;         tv, bytscl(rotate(transpose(reform(flat)),7), 0, 2)
    ;        stop

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
      Jupiter_array[*,*,i] = Jupiter_array[*,*,i] / flat & PRINT, 'Do not apply this flat, it needs spectral *and* spatial normalization to unity'
      Jupiter_array[*,*,i] = Jupiter_array[*,*,i] / float(Sxpar(header, 'EXPTIME')) ;normalize to 1 second exposure time
      write_file = rotate(transpose(reform(Jupiter_array[*,*,i])),7)
      
      SXADDPAR, Header, 'BZERO', 0.0
      MWRFITS, write_file, Dir+'\Processed\' + new_filename + '.Cleaned.fits', header, /create, /silent
    endfor

    Star = [mrdfits(Dir+'\hires' + star_frames[0] + '.fits', 3, header, /fscale), $
            mrdfits(Dir+'\hires' + star_frames[0] + '.fits', 2, header, /fscale), $
            mrdfits(Dir+'\hires' + star_frames[0] + '.fits', 1, header, /fscale)]
    junk =  mrdfits(Dir+'\hires' + star_frames[0] + '.fits', 0, header, /fscale)

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
    for i = 0, n_elements(Europa_frames)-1 do begin
      filename = '\hires' + Europa_frames[i] + '.fits'
      Europa_array[*,*,i] = [mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 3, header, /fscale), $
                             mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 2, header, /fscale), $
                             mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 1, header, /fscale)]

      junk                =  mrdfits(Dir+'\hires' + Europa_frames[i] + '.fits', 0, header, /fscale)
      new_filename        =  STRMID(filename, 0, strpos(filename,'.fits'))

      Europa_array[*,*,i] = Europa_array[*,*,i] - bias
      Europa_array[*,*,i] = Europa_array[*,*,i] / flat & PRINT, 'Do not apply this flat, it needs spectral *and* spatial normalization to unity'
      Europa_array[*,*,i] = Europa_array[*,*,i] / float(Sxpar(header, 'EXPTIME')) ; normalize to 1 second expsosure time
      
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
    Europa_Airglow_params = create_struct( 'torus_lat', torus_lat_array, $
                                           'ET', ET_array, $
                                           'ExpTime', exptime_array, $
                                           'Sun2Euro', Solar_Well_2_Europa_Dshift, $
                                           Europa_Airglow_params )
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
  order_60 = {guess_coeffs:[5873.10,0.0285314,-5.54225e-007], low_bound:448,  hi_bound:501,  WL_range:[5576., 5578], aperture_limit:[7,42], name:'order_60'}                  ; Na order
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
    
    Na_cube          = fltarr(4001, 36, N_elements(Europa_frames))                                 ; the number 26 is based on the aperture upper and lower bounds
     K_cube          = fltarr(4001, 36, N_elements(Europa_frames))                                 ; the number 26 is based on the aperture upper and lower bounds
     O_cube          = fltarr(4001, 36, N_elements(Europa_frames))                                 ; the number 26 is based on the aperture upper and lower bounds
     
    Na_Jupiter_cube = fltarr(4001, 36, N_elements(Jupiter_frames))                                ; See how it's expecting more than 1 jupiter frame here? that's so you can compare each flux calibration later
     K_Jupiter_cube = fltarr(4001, 36, N_elements(Jupiter_frames))
     O_Jupiter_cube = fltarr(4001, 36, N_elements(Jupiter_frames))
    
    for h = 0, N_elements(orders) - 1 do begin

      order = orders[2]
      
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
      xr          = minmax(WL)
      order_lines = ThAr_WL[where( (ThAr_WL gt xr[0]) and (ThAr_WL lt xr[1]))]
      ID          = make_array(N_elements(order_lines), value = 1.e5)

      Frames = 'hires'+strcompress(Europa_frames, /rem)+'.Cleaned.fits'
      for frame = 0, n_elements(Frames)-1 do begin

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

        window, 1, xs = 1800, ys=800, title = 'WAVELENGTH SOLUTION FOR: ' + ORDER.NAME
        cgplot, findgen(s[0]), trace_old, psym = 4, /ynozero, /yst, /xst 
        cgplot, findgen(s[0]), trace, psym = 4, /overplot, color = 'red'
        cgplot, x, POLY( findgen(s[0]), coeffs), /overplot, color = 'blue' ; IF THIS FIT TO THE RED DATA POINTS IS BAD, EVERYTHING DOWNSTREAM FAILS
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

;        E2 = Europa
;        for i = 0, s[0] - 1 do begin ; for column of wavelength, shift the spectra up & down tgo straighten them out.
;          shift_by               = POLY( float(i), coeffs) - POLY( 0., coeffs)
; 
;          sz = size(Europa[i,*])
;          E2[i,*] = interpolate(Europa[i,*], findgen(sz(1)), findgen(sz(2))+shift_by, /grid)
;        endfor
;        window, 0, xs = 4001, ys =2139, title = 'straightened '
;        tv, bytscl(E2, 0, 2000)
;        stop

;
;          interpt_at          = POLY(findgen(s[0]), coeffs)
;stop
;          Star_order_straight = interpolate(Star, findgen(s[0]), interpt_at, cubic = -0.5, /grid)
;          ThAr_order_straight = interpolate(ThAr, findgen(s[0]), interpt_at, cubic = -0.5, /grid)
;          flat_order_straight = interpolate(flat, findgen(s[0]), interpt_at, cubic = -0.5, /grid)
;          Europa_order_straight = interpolate(Europa, findgen(s[0]), interpt_at, cubic = -0.5, /grid)

;  a = Europa[*,order.low_bound:order.hi_bound]
;  S_a = size(a)
;          XI = x[keep]
;          YI = trace[keep]
;          ; Enter the Xos and Yos:
;          XO = x[keep]
;          YO = replicate(50., N_elements(keep)) 
;          ; Run POLYWARP to obtain a Kx and Ky:
;          POLYWARP, XI, YI, XO, YO, 3, KX, KY
;          ; Create a warped image based on Kx and Ky with POLY_2D:
;          B = POLY_2D(A, KX, KY)
;window, 0, xs = s[0], ys = s[1]
;tv, bytscl(a, 0, 2100)
;window, 1, xs = s[0], ys = s[1]
;tv, bytscl(b, 0, 2100)
;
;stop

;        window, 1
;        plot, total(Europa_order_straight[*,10:30],2)
;        ; Uh oh, I can see a the pixelated stair case in the resulting values on any give row, need a better way to straighten!! 
;        window, 0, xs = 4001, ys = 10*slit_length_pix, title = 'straightened '
;        tv, bytscl(alog(rebin(Europa_order_straight, 4001, 10*slit_length_pix))) 
;        stop
        
        
        Flat_aperture = Flat_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;nomalize_with = mean(Flat_aperture, dim = 2, /NAN)
        flat_coeffs = poly_fit(findgen(s[0]), mean(Flat_aperture, dim = 2, /NAN), 3)
        nomalize_with = poly(findgen(s[0]), flat_coeffs)
        Flat_aperture = Flat_aperture / rebin(nomalize_with, s[0],order.aperture_limit[1] - order.aperture_limit[0] + 1) ; Normalize the flat field to unity, now its spectral AND spatially normalized
        aperture      = Europa_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]]
        ;aperture      = aperture / Flat_aperture ; don't flatten twice!!!

        ; Evaluate the guess at the wavelength solution
        ThAr_order_straight[where(ThAr_order_straight lt 0., /null)] = !Values.F_Nan
        ThAr_measured = total(ThAr_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]], 2)
        P = cglayout([1,2], ygap = 0.)
        cgplot, WL, ThAr_measured, /ylog, yr = [2.e3, 5.e5], /xstyle, pos = p[*,0], xtickformat = '(A1)', ytitle = 'ThAr Lamp Counts'
        
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
       ;STOP ; <--- INSPECT WL SOLUTION HERE

; write the fits files and run the Cosmic Ray Corrections. Some orders overlap the CCD edge in the extraction so exclude these regions in the CR correction
        SXADDPAR, Header, 'BZERO', 0
        SXADDPAR, Header, 'BSCALE', 0
        MWRFITS, aperture, Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_' + Frames[frame], header, /create
        
        Case 1 of
          order.name eq 'order_43': statsec = '[0:2925,*]'
          order.name eq 'order_52': statsec = '[500:4000,*]'
          order.name eq 'order_53': statsec = '[0:2620,*]'
          else: junk = temporary(statsec)
        endcase
        
        gain = 1.0
        sigclip = 8.5
        la_cosmic, Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_' + Frames[frame],outsuff = "CR", sigclip = sigclip, statsec = statsec, gain = gain
        CR_result_1 = mrdfits(Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_hires'+strcompress(Europa_frames[frame], /rem)+'.Cleaned.fits', 0, junk_header)
        CR_Mask     = mrdfits(Dir+'\Processed\Cosmic Rays\'+order.name+'_CR_hires'+strcompress(Europa_frames[frame], /rem)+'.Cleaned-mask.fits', 0, junk_header)
        
; the LA Cosmic CR removal algorithm can sometimes introduce negative pixels, particularly at the edge, fix those
        n_sigma = 1.*sigclip
        RESISTANT_Mean, CR_result_1, n_sigma, Mean_CR_result_1, Sigma_Mean_CR_result_1, Num_RejECTED
        Sig_CR_result_1 = ROBUST_SIGMA( CR_result_1 )
        CR_result_2 = CR_result_1
        Bad_Pixel_Indicies = where((CR_result_1 lt 0.) or (CR_result_1 gt (Mean_CR_result_1+n_sigma*Sig_CR_result_1)), /NULL, complement = Good_Pixel_Indicies)
        s = size(aperture, /dim)
        bad_pixel_list = intarr(s[0],s[1])
        bad_pixel_list[Good_Pixel_Indicies] = 1
        fixpix, CR_result_1, bad_pixel_list, CR_result_2     ; replace Bad_Pixel_Indicies with avg of neighbor pixels
        
; the hot pixel correction really messes up the on-disk observations of europa, so i will only use CR_Result_1 for those frames (0, 1, 2, 11, 12, 13). definitely need
; hot pixel correction for all frames off the disk though.
        
        if order.name eq 'order_46' then begin
          K_cube[*,*,frame] = CR_result_2
          if (frame eq 0) or (frame eq 1) or (frame eq 2) or (frame eq 11) or (frame eq 12) or (frame eq 13) then begin
            K_cube[*,*,frame] = CR_result_1
          endif
        endif
        if order.name eq 'order_60' then begin
          Na_cube[*,*,frame] = CR_result_2
          if (frame eq 0) or (frame eq 1) or (frame eq 2) or (frame eq 11) or (frame eq 12) or (frame eq 13) then begin
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
          tv, bytscl(aperture);, 0, 100) ;0.5*mean(aperture), 1.5*mean(aperture))
          window, 4, xs = 3400, ys = slit_length_pix, ypos = 100, title = 'Cosmic Ray Corrected Frame '+ Europa_frames[frame]
          tv, bytscl(CR_result_1);, 0, 100) ;0.5*mean(aperture), 1.5*mean(aperture))
          window, 8, xs = 3400, ys = slit_length_pix, ypos = 200, title = 'Hot/Cold Pixel Filtered Cosmic Ray Corrected Frame '+ Europa_frames[frame]
          tv, bytscl(CR_result_2);, 0, 100) ;0.5*mean(aperture), 1.5*mean(aperture))
          window, 12, xs = 3400, ys = slit_length_pix, ypos = 300, title = 'Cosmic Ray Mask'
          tv, bytscl(CR_mask, 0, 1)
          ;wait, 2
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

        if order.name eq 'order_46' then K_Jupiter_cube[*,*,frame] = Jupiter_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]] / Flat_aperture ; flat field jupiter
        if order.name eq 'order_56' then O_Jupiter_cube[*,*,frame] = Jupiter_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]] / Flat_aperture ; flat field jupiter
        if order.name eq 'order_60' then Na_Jupiter_cube[*,*,frame] = Jupiter_order_straight[*, order.aperture_limit[0]:order.aperture_limit[1]] / Flat_aperture ; flat field jupiter
        
      endfor ; frames (Jupiter frame number)
      
      if order.name eq 'order_46' then  K_WL = WL
      if order.name eq 'order_56' then  O_WL = WL
      if order.name eq 'order_60' then Na_WL = WL
      
      save, Na_cube, K_cube, O_cube, Na_Jupiter_cube, K_Jupiter_cube, O_Jupiter_cube, Na_WL, K_WL, O_WL, filename = Dir+'\Processed\.sav'
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
       cgplot, mean(jupiter_cube, dim=2)                                                         ; should be about 5.5 MR / A
       
       for f = 0, N_elements(Europa_frames)-4 do begin
         Na_cube[*,*,f] = Na_cube[*,*,f] / _2D_sensitivity ; convert to R / A units
          K_cube[*,*,f] =  K_cube[*,*,f] / _2D_sensitivity ; convert to R / A units
          O_cube[*,*,f] =  O_cube[*,*,f] / _2D_sensitivity ; convert to R / A units
         print, 'you neglected differential airmass when scaling the sensitivity calculation from Jupiter to Europa!'
       endfor
       
       if orders[order].name eq 'order_46' then K_Jupiter_cube  = Jupiter_cube
       if orders[order].name eq 'order_56' then O_Jupiter_cube  = Jupiter_cube
       if orders[order].name eq 'order_60' then Na_Jupiter_cube = Jupiter_cube
       
     ; In a given row above be sure the answer matches the following notes...
     ;   -Trafton (1980) does the same thing I do: (page 118 here https://doi.org/10.1016/0019-1035(80)90249-3). A pi on both sides, but needs the factor of 4 to put it in Rayleighs
     ;   -Brown (1981) (doi:10.1086/158777) quotes 5.6 MR/A at 6300A. I get about 5.9 MR/A if Jupiter were its average distance of 5.2 AU
     ;   -Brown & Shemansky (1982) doi:10.1086/160515 quote 5.4 MR/A at 6724A. I get 5.7 MR/A if Jupiter were its average distance of 5.2 AU  
              
     save, Na_cube, K_cube, O_cube, Na_Jupiter_cube, K_Jupiter_cube, O_Jupiter_cube, Na_WL, K_WL, O_WL, ang_radius, _2D_sensitivity, filename = Dir+'\Processed\_Flux_Calibrated.sav'
     print, 'Flux-calibrated ', orders[order].name
   endfor  
  stop
  PRINT, 'To Do: did not yet inspect flux calibration of K D and O 6300 orders!' 
  endif    
    
  if part eq 1.6 then begin
    restore, Dir+'\Processed\Europa_Airglow_params.sav'
    help, Europa_Airglow_params
    
    for i = 0, n_elements(Europa_frames)-1 do begin
      filename = dir+'\hires' + Europa_frames[i] + '.fits'
      header   =  headfits(filename)
      print, sxpar(header, 'ROTPOSN') - 90.
    endfor

    
    SOLAR_SPECTRUM_FILE = 'Z:\DATA\___Calibration___\Solar_and_Stellar_Calibration_Spectra\Coddington_2022\hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
    NCDF_LIST, SOLAR_SPECTRUM_FILE, /VARIABLES, /DIMENSIONS, /GATT, /VATT

    ID = NCDF_open(SOLAR_SPECTRUM_FILE)

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
      
      Juno    = 0.5*(REFORM(cube[*,*,0]+cube[*,*,1]+cube[*,*,2]+cube[*,*,3]+cube[*,*,4]))
      EW0__   = 0.5*(REFORM(cube[*,*,5 ]+cube[*,*,6 ]+cube[*,*,7]+cube[*,*,27]+cube[*,*,28]+cube[*,*,29]))             ; EW orientation on-disk
      EW10_N  = 0.5*(REFORM(cube[*,*,8 ]+cube[*,*,9 ]))
      EW20_N  = 0.5*(REFORM(cube[*,*,10]+cube[*,*,11]))
      EW10_S  = 0.5*(REFORM(cube[*,*,12]+cube[*,*,13]))
      EW20_S  = 0.5*(REFORM(cube[*,*,14]+cube[*,*,15]))
      NS0__   = 0.5*(REFORM(cube[*,*,16]+cube[*,*,17]+cube[*,*,18]+cube[*,*,34]+cube[*,*,35]))             ; NS orientation on-disk
      NS10_W  = 0.5*(REFORM(cube[*,*,19]+cube[*,*,20]+cube[*,*,30]+cube[*,*,31]))
      NS20_W  = 0.5*(REFORM(cube[*,*,21]+cube[*,*,22]))
      NS10_E  = 0.5*(REFORM(cube[*,*,23]+cube[*,*,24]+cube[*,*,32]+cube[*,*,33]))
      NS20_E  = 0.5*(REFORM(cube[*,*,25]+cube[*,*,26]))
      
      ; we'll use NS0 only for now and that's the 12th index in the spectra we observed...
        sun2europa = Europa_Airglow_params.sun2euro[12]
        s = size(NS0__)
      
      orientations = fltarr(s[1], s[2], 11)                         ; there's definitely a better way to do this but whatever
      orientations[*,*,0] = NS0__ 
      orientations[*,*,1] = NS10_E
      orientations[*,*,2] = NS20_E
      orientations[*,*,3] = NS10_W
      orientations[*,*,4] = NS20_W
      orientations[*,*,5] = EW0__ 
      orientations[*,*,6] = EW10_N
      orientations[*,*,7] = EW20_N
      orientations[*,*,8] = EW10_S
      orientations[*,*,9] = EW20_S
      orientations[*,*,10]= Juno
      s[2] = s[2] - 4.                                              ; HACK! just wanted to crop out the weird order for the Na order; don't know if it'll be the same for K and O
      
      labels = ['NS0__','NS10_E','NS20_E','NS10_W','NS20_W','EW0__','EW10_N','EW20_N','EW10_S','EW20_S','JUNO flyby']
    
    sunmax = []
    newimg    = fltarr(s[1], s[2])
    sunimg    = fltarr(s[1], s[2])
    eurimg    = fltarr(s[1], s[2])
    
    
; ------------------------------------------- Find the Dispersion & Sunlight vs Exosphere Indicies --------------------------------------

    D2Cen               = 609
    D1Cen               = 825
    
    windowwidth         = 30.
    spec_1D             = total(NS0__, 2, /Nan)
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
      
    P_returned  = fltarr(4,s[2])                  ; Three coefficients + MPFIT's "Status"
    P_guessed   = Fltarr(3,s[2])                  ; Initial Guess that we throw at MPFIT
    
    
    FOR orientation = 10, N_Elements(orientations[0,0,*]) - 1 DO BEGIN
      mpfitD2emission = []
      mpfitD1emission = []
      dont_include    = []
      FOR i = 0, s[2] - 1 DO BEGIN

    ; generate an initial guess for multipliciative scaling
        europa      = orientations[*,*,orientation]
        row         = europa[*,i]
        
    ; Finds the disk center spectrum in each column and puts that row into a 1D array
        FOR c = 0, s[1] - 1 DO BEGIN
          column = europa[c,*]
          suncol = WHERE(column[0:20] EQ MAX(column[0:20]),count)
          suncol = suncol[0]
          sunmax = [sunmax, suncol[0]]
        ENDFOR

        pixelsofsun = 1;3
        sunlight    = europa[*, suncol - pixelsofsun : suncol + pixelsofsun]
        ;TV, bytscl(sunlight)
        sunlight_1d = TOTAL(sunlight, 2, /Nan)
        
        guess_scale = median(row[fitindices] / sunlight_1D[fitindices])
        row_err     = sqrt(abs(row))
        
    ; Fit a y = A*Gauss_smooth(x,C) + B function to the spectrum, where x is the reference solar spectrum
        p0 = [guess_scale, 0.0, 0.5]                                           ; Guess at initial coefficients
        parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p0))
        
        parinfo.value         = p0
        parinfo[0].limited    = [1, 1]
        parinfo[0].limits     = [0.0, 2.]
        parinfo[1].limited    = [1, 1]
        parinfo[1].limits     = [-20.0, 20.]
        parinfo[2].limited    = [1, 1]
        parinfo[2].limits     = [0.0, 20.]

        WEIGHTS = 1./(abs(findgen(s[1]) - Europa_D2_pixel))^.4 + 1./(abs(findgen(s[1]) - Europa_D1_pixel))^.4
        
        fa = {x:sunlight_1d[fitindices], y:row[fitindices], err:sqrt(abs(row[fitindices]))}
        p = mpfit('match_scattered_sunlight', p0, PERROR = err_P, functargs=fa, status=status, parinfo=parinfo, quiet=quiet)
        P_guessed[*,i]  = p0
        p_returned[*,i] = [p, status]
        scaled_sunlight = scale_fit_sunlight(P, sunlight_1d)             ; Puts it into y = A*shift(Gauss_smooth(x,D),C) + B form
        
        sunimg[*,i] = scaled_sunlight
        sub         = guess_scale * sunlight_1d                          ; If you JUST want the multiplicative correction (no scattered sunlight accounted for w mpfit)
        totsubtrd   = row  - scaled_sunlight                             ; Change back to row - sub to get just the multiplicative factor
        if i eq suncol then continuum = scaled_sunlight                  ; this saves the non-sun subtracted continuum row so that i can reference it later
        
; for each row, fit a gaussian to the D1 and D2 lines. find area under curve --> emission!
        windowwidth         = 30.
        
; ---------------------------------------------------- Na D1 line --------------------------------------------------------
        LSF_fitting_ind1  = where( abs(wl - wl[D1cen]) lt 0.2, /NULL)
        LSF_fitting_ind2  = where( abs(wl - wl[D2cen]) lt 0.2, /NULL)
        cloud_ind1        = where( abs(wl - wl[D1cen]) lt 1.5 , /NULL)
        cloud_ind2        = where( abs(wl - wl[D2cen]) lt 1.5, /NULL)

;;;;;;; this is me trying to subtract io cloud...not working well. first, i block out the europa emission & continuum and then interpolate over these blocked out regions.
;;;;;;; i then subtract this "fake" io background from the og image and add the europa emission back to the row. this doesn't work for all files.
  
        qualitymetric = 0
        qualitymetric = stddev(totsubtrd[fitindices]) / total(row)
        if qualitymetric lt 5.0E-05 then begin
          totsubtrd = !values.F_nan;fltarr(4001)                                       ; sets threshold and gets rid of sunlight over disk
          dont_include = [dont_include, i]
          ;continue                                                      ; include the continue if you do not want to include blocked out continuum
        endif
        eurimg[*,i]  = totsubtrd
      endfor
      meanimg = median(eurimg, dim=2)
      meanimg[LSF_fitting_ind1] = !values.F_nan                                                 ; first, i'm masking out europa's emission (continuum is already masked out)
      meanimg[LSF_fitting_ind2] = !values.F_nan
      fakeio = interpol(meanimg, WL, WL, /NAN)
      
      for i = 0, s[2] - 1 DO BEGIN
        dummy_row = eurimg[*,i]
;        dummy_row[LSF_fitting_ind1] = !values.F_nan                                                 ; first, i'm masking out europa's emission (continuum is already masked out)
;        dummy_row[LSF_fitting_ind2] = !values.F_nan
        if max(dummy_row) eq 0. then scaled_row = dummy_row
        if max(dummy_row) ne 0. then scaled_row = dummy_row * (max(meanimg) / max(dummy_row))                                     ; scales row to the spatially resolved mean of sunlight-subtracted image
        
        if order_index eq 2 then xr    = [5888.1, 5898.5]
        junk  = min(abs(xr[0]- WL), index0)
        junk  = min(abs(xr[1]- WL), index1)
        
        no_io = dummy_row - fakeio
        
        window, 5, title='D1 io sub'
        cgplot, WL[index0:index1], dummy_row[index0:index1], yr= [-3000.,10000.]
        cgplot, wl[index0:index1], fakeio[index0:index1], /overplot, color='red', psym=12
      
;        no_io[LSF_fitting_ind1] = totsubtrd[LSF_fitting_ind1]
;        no_io[LSF_fitting_ind2] = totsubtrd[LSF_fitting_ind2]
;        noIo = scaled_row - fakeio
        totsubtrd = no_Io
        if dummy_row[0] eq !values.F_NAN then totsubtrd = !values.F_NAN
        if dummy_row[0] eq !values.F_nan then continue
        cgplot, WL[index0:index1], totsubtrd[index0:index1], /overplot, color='blue'
        ;cgplot, WL[index0:index1], scaled_row[index0:index1], /overplot, color='green'
        newimg[*,i] = totsubtrd
        
;; now that the io cloud is subtracted, fit the emission
        
        guess_peak       = max(totsubtrd[LSF_fitting_ind1], loc)                                     ; initial guesses [amplitude, peak centroid, hwhm]
        guess_low        = min(totsubtrd[LSF_fitting_ind1], minloc)
        p                = [abs(guess_peak-(guess_low)), wl[LSF_fitting_ind1[loc]], 0.07, median(totsubtrd[index0:index1])]
        
        parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p))
        parinfo.value         = p
        parinfo[0].limited    = [1, 1]
        parinfo[0].limits     = [0.0, guess_peak]
        parinfo[1].limited    = [1, 1]
        parinfo[1].limits     = [minloc, loc]
        parinfo[2].limited    = [1, 0]
        parinfo[2].limits     = [0.0, 0.5]
        parinfo[3].limited    = [1, 1]
        parinfo[3].limits     = [abs(guess_low), abs(guess_low)+1000.]
        
        y1                    = totsubtrd[LSF_fitting_ind1]
        D1fa                  = { x:wl[LSF_fitting_ind1], y:y1, err:1./sqrt(abs(y1)) }
        a                     = mpfit('Gaussian_for_MPFIT', p, funct=D1fa, parinfo=parinfo, STATUS = Did_it_work)
        D1_area               = a[0] * a[2] * sqrt(2.*!pi)
        mpfitD1emission       = [mpfitD1emission, D1_area]
        
        ; Inspect:
;        window, 3, title='D1'
;        cgplot, D1fa.x, D1fa.y, err_yhigh = D1fa.err, err_ylow = D1fa.err
;        cgplot, D1fa.x, gaussian(D1fa.x, a), /overplot, color='red'
        
; ---------------------------------------------------- Na D2 line -------------------------------------------------------- 
        
        guess_peak       = max(totsubtrd[LSF_fitting_ind2], loc)                                     ; initial guesses [amplitude, peak centroid, hwhm, vertical shift]
        guess_low        = min(totsubtrd[LSF_fitting_ind2], minloc)
        p                = [guess_peak-(guess_low), wl[LSF_fitting_ind2[loc]], 0.07, median(totsubtrd[index0:index1])] 
        
        parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p))
        parinfo.value         = p
        parinfo[0].limited    = [1, 1]
        parinfo[0].limits     = [0.0, guess_peak]
        parinfo[1].limited    = [1, 1]
        parinfo[1].limits     = [0, loc]
        parinfo[2].limited    = [1, 0]
        parinfo[2].limits     = [0.0, 0.5]
        parinfo[3].limited    = [1, 1]
        parinfo[3].limits     = [guess_low, guess_low+1000.]
        
        y2                    = totsubtrd[LSF_fitting_ind2]
        D2fa                  = { x:wl[LSF_fitting_ind2], y:y2, err:sqrt(y2*100.) }
        a                     = mpfit('Gaussian_for_MPFIT', p, funct=D2fa, parinfo=parinfo, STATUS = Did_it_work)
        D2_area               = a[0] * a[2] * sqrt(2.*!pi)
        mpfitD2emission       = [mpfitD2emission, D2_area]

        ; Inspect:
;        window, 2, title='D2'
;        cgplot, D2fa.x, D2fa.y, err_yhigh = D2fa.err, err_ylow = D2fa.err
;        cgplot, D2fa.x, gaussian(D2fa.x, a), /overplot, color='red'                                         ; all plots in DN/s
      ENDFOR ; each row of ONE orientation
      
      window, 0
      loadct, 3
      cgimage, newimg[index0:index1, *], minv=-1e4, maxv=1e4
      window, 1
      cgplot, wl[index0:index1], total(newimg[index0:index1, *], 2, /nan)
      
      plate_scale = 0.358
      yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]
;      if newimg[index0:index1, *] ne 
;      window, 1, title='D1 mpfit'
;      cgplot, mpfitD1emission, xr=[0,n_elements(mpfitD1emission)], xtickformat='(A1)'
;      cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
;      window, 2, title='D2 mpfit'
;      cgplot, mpfitD2emission, xr=[0,n_elements(mpfitD2emission)], xtickformat='(A1)'
;      cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
;      
      angstrom_per_pixel = mean(deriv(WL))
;              
;; below, i calculate units of rayleighs to match leblanc (2005) plots. 
;          
      nogaussfit = newimg;[LSF_fitting_ind1, *] + newimg[LSF_fitting_ind2, *]
      collapse   = mean(nogaussfit, dim=1)
;      
      ;very rough calculations here
      ind = min(abs(wl[500:900] - 5890.2),loc)
      profile = mpfitD1emission + mpfitD2emission
;      
      window, 3
      cgplot, profile, title = labels[orientation], ytitle = 'Rayleighs', xr=[0,n_elements(profile)], xtickformat='(A1)', yr=[0,max(collapse)], xticklen=0
      cgplot, collapse, /overplot, color='red'
      cgaxis, xaxis = 0, xtit = 'Europa Radii', xr = yr, xstyle=1, xticklen=-0.01
      
;---------------------------------------------------- calculate the g-value -------------------------------------------------------
; first, will use SPICE to get europa's heliocentric velocity... need observations time stamp so reload headers

      dates = []
      
      ; need to match up the time stamps to orientations to get the specific g-values
      
      if orientation eq 0 then file = 11
      if orientation eq 1 then file = 18
      if orientation eq 2 then file = 20  
      if orientation eq 3 then file = 14  
      if orientation eq 4 then file = 16  
      if orientation eq 5 then file = 0   
      if orientation eq 6 then file = 3   
      if orientation eq 7 then file = 5   
      if orientation eq 8 then file = 7   
      if orientation eq 9 then file = 9   
      
      filename      = order.name+'_CR_hires' + Europa_frames[file] + '.Cleaned.fits'
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
      
      column = profile * 9. * angstrom_per_pixel * 1.5 * 10.e6 / g_Na
      
      
      ; Compare raw vs sun-subtracted spectra.
      xr    = [5888.1, 5898.5]
      junk  = min(abs(xr[0]- WL), index0)
      junk  = min(abs(xr[1]- WL), index1)
      WL_xr = [index0, index1]
      
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3-panel plots w col. dens. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      P = cglayout([2,2], ygap = 0., oxmargin = [14, 2], oymargin = [9, 5], xgap = 0.)
      axis_format = {XTicklen:-.01, yticklen:-0.01 }
      
      cgPS_Open, filename = Dir+'\Figures\3panel_'+order.name+'_'+labels[orientation]+'.eps', /ENCAPSULATED, xsize = 7.5, ysize = 5
        !P.font=1
        loadct, 3
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        
        title = 'HIRES 2022-09-29 : '+labels[orientation]
  
        cgplot, wl[500:900], total(newimg[index0:index1, *], 2, /NAN), /xs, xr = [wl[500], wl[900]], pos = p[*,0], xtickformat = '(A1)', $
          title=title, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom')
  
        cgimage, newimg[index0:index1, *], minv=-5.e3, maxv=5.e3, /axes, xr = xr, pos = p[*,2], yr = yr, /noerase, $
          xtitle = 'Angstroms', AXKEYWORDS = axis_format
        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01

        cgcolorbar, POSITION=[0.57, 0.55, 0.59, 0.95], range = [0.0,max(column/1.e10)], /vertical, /right
  
        cgplot, column/1.e10, findgen(N_elements(newimg[index0:index1, *])), /ynozero, /noerase, pos = p[*,3], ys= 5, $
          xtitle = 'Na Column Density ('+cgsymbol('times')+'10!U10!N atoms / cm!U2!N)';, xr=[0,400]
        cgaxis, xaxis = 1, xtit = 'D1 + D2 (Rayleighs)', xr = [min(profile), max(profile)]

      cgps_Close
      
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sun-sub comparison figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      
      ;P  = cglayout([2,2], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
      plate_scale = 0.358
      yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]

      cgPS_Open, filename = 'Z:\DATA\Keck\Europa Na\HIRES_20220928\Figures\'+order.name+'_'+labels[orientation]+'_suncomparison.eps', $
        /ENCAPSULATED, xsize = 10, ysize = 8
        !P.font=1
        loadct, 3
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
  
        title = 'HIRES 2022-09-29 : '+labels[orientation]
  
        ; non-sunlight subtracted on the left hand side
        cgplot, wl, europa[*,suncol], /xs, xr = xr, pos = p[*,0], xtickformat = '(A1)', $
          title=labels[orientation]+' Raw', ytitle = 'Rayleighs / ' + cgsymbol('Angstrom')
  
        cgplot, wl, continuum, /overplot, color='red'
  
        cgtext, .29, .58, 'Above disk', color = 'black', /normal
        cgtext, .29, .55, 'Fit Reflectance', color = 'red', /normal
  
        axis_format = {XTicklen:-.01, yticklen:-0.01 }
        cgimage, orientations[index0:index1,*,orientation], minv=0.75*mean(orientations[*,*,orientation]), maxv=3.0*mean(orientations[*,*,orientation]), $
          /axes, xr = xr, pos = p[*,2], /noerase, $
          ytit = 'Europa Radii', xtitle = 'Angstroms', AXKEYWORDS = axis_format, yr = yr
  
        cgcolorbar, POSITION=[p[0,0], 0.045, p[2,0], 0.065], range = minmax(orientations[*,*,orientation])
  
        ; sunlight subtracted on the right hand side
        cgplot, wl, total(newimg, 2, /NAN), /xs, xr = xr, pos = p[*,1], xtickformat = '(A1)', $
          ytickformat = '(A1)', /noerase, title=labels[orientation]+' Sun-Subtracted'
        cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
  
        axis_format = {XTicklen:-.01, yticklen:-0.01, ystyle:5}
        cgimage, newimg[index0:index1, *], minv=-100, maxv=100, /axes, xr = xr, pos = p[*,3], yr = yr, /noerase, $
          xtitle = 'Angstroms', AXKEYWORDS = axis_format
        cgaxis, yaxis = 1, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
  
        cgcolorbar, POSITION=[p[2,0], 0.045, 0.979, 0.065], range = [-5.e3,5.e3]
      cgps_Close
      
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5 panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      P  = cglayout([3,2], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
      plate_scale = 0.358
      yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]

      cgPS_Open, filename = dir+'\Figures\test_figures\'+order.name+'_'+labels[orientation]+'_5panel_process.eps', $
        /ENCAPSULATED, xsize = 20, ysize = 10
        !P.font=1
        loadct, 3
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        
        title = 'HIRES 2022-09-29 : '+labels[orientation]+' Subtraction Process'
        
        ;   non-sunlight subtracted on the upper left hand side
        cgimage, europa[index0:index1, *], pos = p[*,0], title=labels[orientation]+' Raw', ytitle = 'Rayleighs / '+cgsymbol('Angstrom')
        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
        cgtext, .2, .58, 'Raw', color = 'white', /normal
        
        ;   sunlight subtracted on the lower left hand side
        cgimage, eurimg[index0:index1, *], pos = p[*,3], /axes, /noerase, AXKEYWORDS = axis_format, xr = xr, xtitle ='Angstroms'
        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
        cgtext, .2, .15, 'Sun Sub', color = 'white', /normal
        
        ;   sun and io subtracted on lower middle
        cgimage, newimg[index0:index1, *], /axes, xr = xr, pos = p[*,4], yr = yr, /noerase, $
          xtitle = 'Angstroms', AXKEYWORDS = axis_format
        cgtext, .5, .15, 'Sun and Io Sub', color = 'white', /normal
        
        ;   1D spectra of sun- and io-subbed on the upper middle
        cgplot, wl, total(newimg[index0:index1, *], 2, /NAN), /xs, xr = xr, pos = p[*,1], xtickformat = '(A1)', $
          ytickformat = '(A1)', /noerase, title=labels[orientation]+' Sun-Subtracted'
        cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
        
        ;   1D spectra of sun- and io-subbed on the lower right hand side
        cgplot, column/1.e10, findgen(N_elements(newimg[index0:index1, *])), /ynozero, /noerase, pos = p[*,5], ys= 5, $
          xtitle = 'Na Column Density ('+cgsymbol('times')+'10!U10!N atoms / cm!U2!N)';, xr=[0,400]
        cgaxis, xaxis = 1, xtit = 'D1 + D2 (Rayleighs)', xr = [min(profile), max(profile)]
        
      cgPS_Close

; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6 panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      P  = cglayout([3,2], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
      plate_scale = 0.358
      yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]

      cgPS_Open, filename = dir+'\Figures\test_figures\'+order.name+'_'+labels[orientation]+'_6panel_compare.eps', $
        /ENCAPSULATED, xsize = 20, ysize = 10
        !P.font=1
        loadct, 3
        device, SET_FONT = 'Helvetica Bold', /TT_FONT
        
        title = 'HIRES 2022-09-29 : '+labels[orientation]+' Subtraction Comparisons'
        
        ;   raw on the left hand side
        cgimage, europa[index0:index1, 0:32], pos = p[*,0], title=labels[orientation]+' Raw', ytitle = 'Rayleighs / '+cgsymbol('Angstrom')
        cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
        ;   raw 1D
        cgplot, wl[index0:index1], total(europa[index0:index1, *], 2, /NAN), pos = p[*,3], xr = xr, /xs, ytickformat = '(A1)', /noerase
        cgaxis, yaxis = 0, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
        
        ;   sunlight subtracted in the middle
        cgimage, eurimg[index0:index1, *], pos = p[*,1], /noerase, title='Sun-Sub'
        ;   sun sub 1D
        cgplot, wl[index0:index1], total(eurimg[index0:index1, *], 2, /NAN), pos = p[*,4], xr = xr, /xs, ytickformat = '(A1)', /noerase
        
        ;   sun and io subtracted on the right hand side
        cgimage, newimg[index0:index1, *], pos = p[*,2], yr = yr, /noerase, title='Sun-Io-Sub'
        ;   sun io sub 1D
        cgplot, wl[index0:index1], total(newimg[index0:index1, *], 2, /NAN), pos = p[*,5], ytickformat = '(A1)', /noerase, xr = xr, /xs;, yr=[min(total(eurimg, 2)), max(total(eurimg, 2))+1.e3]
        cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1, xtitle = 'Angstroms', xaxis = 0
        
      cgPS_Close
      
    endfor ; orientation
    stop
  endfor ; order
  save, /all, filename = Dir+'\Processed\sun_subbed_europa.sav'
        stop
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
  endif
  
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