pro quick_gvalue_calc

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
  
  
  dir = 'Z:\DATA\Keck\Europa Na\HIRES_20220928'
  
  RESTORE, Dir+'\Processed\sun_subbed_europa.sav'
  
  stop
  date          = '1999-12-28T05:30:00'

  cspice_str2et, date, et

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

  column = profile * 9. * angstrom_per_pixel * 1.5 * 10.e6 / g_Na
  
  window, 0
  cgplot, column/1.e10, findgen(N_elements(newimg[index0:index1, *]))
;  cgaxis, xaxis = 1, xtit = 'D1 + D2 (Rayleighs)', xr = [min(profile), max(profile)]

  stop


end