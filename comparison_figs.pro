pro comparison_figs


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
  
  ; ----------------------------------------------------------------------------------------------------------------
  
  dir = 'Z:\DATA\Keck\Europa Na\HIRES_20220928'
  RESTORE, dir+'\Processed\calcg.sav'
  
  EWondisk_brown = read_csv('Z:\DATA\Keck\Europa Na\HIRES_20220928\leblanc2005EW_ondisk.csv')
  NSondisk_brown = read_csv('Z:\DATA\Keck\Europa Na\HIRES_20220928\leblanc2005NS_ondisk.csv')
  EW10N_brown = read_csv('C:\Users\elovett\EuropaResearch\Paper 1\leblanc2005_EW10N.csv')
; calculate g-value
  browndate =  '2022-09-29T10:45:37'

  cspice_str2et, browndate, et

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

  EWcolumn = EWondisk_brown.field2; * 10.e6 / g_Na
  NScolumn = NSondisk_brown.field2; * 10.e6 / g_Na
  EW10Ncol = EW10N_brown.field2; * 10.e6 / g_Na
  RESTORE, Dir+'\Processed\0132_EWondisk_rayleighs.sav'
  EW1 = area_profile
  RESTORE, dir+'\Processed\0163_EWondisk_rayleighs.sav'
  EW2 = area_profile
  RESTORE, dir+'\Processed\0164_EWondisk_rayleighs.sav'
  EW3 = area_profile
  
  RESTORE, dir+'\Processed\0146_NSondisk_rayleighs.sav'
  NS1 = area_profile
  RESTORE, dir+'\Processed\0147_NSondisk_rayleighs.sav'
  NS2 = area_profile
;  RESTORE, dir+'\Processed\0159_NSondisk_rayleighs.sav'
;  NS3 = area_profile
;  RESTORE, dir+'\Processed\EW10N_rayleighs.sav'
;  EW10N_area_profile = area_profile
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; EW PLOTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;  window, 0, title='EW Comparison'
;  cgplot, EWondisk_brown.field1, EWcolumn, psym=2, title='Comparison to LeBlanc et al. (2005)',ytitle='D1+D2 Emission (Rayleighs)', xtitle='Distance from Europa (Europa radii)', color='green', xr=[-20,20], /ynozero
;  cgplot, EWondisk_brown.field1, EW_area_profile, /overplot, color='black', psym=4
;  cglegend, colors=['green', 'black'],psym=[2,1], titles=['LeBlanc et al. (2005)', 'Keck/HIRES'], length=0.01, symsize=0.1, /Box, Location=[0.15, 0.70], charsize=1.0, /Background, vspace=1
  
  xaxis = findgen(n_elements(area_profile), start=yr[0], increment=(yr[1]-yr[0])/n_elements(area_profile))
  
  P  = cglayout([1,1], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
        
  cgPS_Open, filename = dir+'\Figures\leblanc05_compareEW.eps', $
    /ENCAPSULATED, xsize = 8, ysize = 6
      !P.font=1
      loadct, 3
      device, SET_FONT = 'Palatino Linotype', /TT_FONT
    
      title = 'Comparing Keck/HIRES 2022/09/29 to LeBlanc (2005) EW'
      
      cgplot, EWondisk_brown.field1, EWcolumn, psym=9, title='EW Comparison to LeBlanc et al. (2005)',ytitle='D1+D2 Emission (Rayleighs)', xtitle='Distance from Europa (Europa radii)', $
        color='black', /ynozero, thick=3.0, xr=[min(EWondisk_brown.field1), max(EWondisk_brown.field1)]
      cgplot, xaxis, EW1, /overplot, color='hot pink', psym=16, symsize=1.5, /ylog;, ERR_xLOW = error_bars, ERR_xHigh = error_bars
      cgplot, xaxis, EW2, /overplot, color='hot pink', psym=16, symsize=1.5, /ylog
      cgplot, xaxis, EW3, /overplot, color='hot pink', psym=16, symsize=1.5, /ylog
    ;  cgplot, EWondisk_brown.field1, EW2, /overplot, color='black', psym=16, symsize=1.5
    ;  cgplot, EWondisk_brown.field1, EW3, /overplot, color='black', psym=16, symsize=1.5
      cglegend, colors=['black', 'hot pink'],psym=[9,16], titles=['LeBlanc et al. (2005)', 'Lovett et al. (in prep)'], length=0, symsize=1.5, /Box, Location=[0.15, 0.80], charsize=2.0, /Background, vspace=2
  
  cgPS_Close
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; NS PLOTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  
;  window, 1, title='NS Comparison'
;  cgplot, NSondisk_brown.field1, NScolumn, psym=2, title='Comparison to LeBlanc et al. (2005)',ytitle='D1+D2 Emission (Rayleighs)', xtitle='Distance from Europa (Europa radii)', color='green', xr=[-20,20], /ynozero
;  cgplot, EWondisk_brown.field1, NS_area_profile, /overplot, color='black', psym=4
;  cglegend, colors=['green', 'black'],psym=[2,1], titles=['LeBlanc et al. (2005)', 'Keck/HIRES'], length=0.5, symsize=0.5, /Box, Location=[0.15, 0.70], charsize=2.0, /Background, vspace=1
  
  
  P  = cglayout([1,1], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
        
  cgPS_Open, filename = dir+'\Figures\leblanc05_compareNS.eps', $
    /ENCAPSULATED, xsize = 10, ysize = 6
      !P.font=1
      loadct, 3
      device, SET_FONT = 'Palatino Linotype', /TT_FONT
    
      title = 'Comparing Keck/HIRES 2022/09/29 to LeBlanc (2005) NS'
      
      cgplot, NSondisk_brown.field1, NScolumn, psym=9, title='NS Comparison to LeBlanc et al. (2005)',ytitle='D1+D2 Emission (Rayleighs)', xtitle='Distance from Europa (Europa radii)', $
        color='black', /ynozero, thick=3.0, xr=[min(EWondisk_brown.field1), max(EWondisk_brown.field1)]
      cgplot, xaxis, NS1, /overplot, color='hot pink', psym=16, symsize=1.5, /ylog;, ERR_xLOW = error_bars, ERR_xHigh = error_bars
      cgplot, xaxis, NS2, /overplot, color='hot pink', psym=16, symsize=1.5, /ylog
    ;  cgplot, EWondisk_brown.field1, EW2, /overplot, color='black', psym=16, symsize=1.5
    ;  cgplot, EWondisk_brown.field1, EW3, /overplot, color='black', psym=16, symsize=1.5
      cglegend, colors=['black', 'hot pink'],psym=[9,16], titles=['LeBlanc et al. (2005)', 'Lovett et al. (in prep)'], length=0, symsize=1.5, /Box, Location=[0.15, 0.80], charsize=2.0, /Background, vspace=2
  
  cgPS_Close
  stop
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; EW10 N PLOTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  window, 0, title='EW10N Comparison'
  cgplot, EWondisk_brown.field1, EW10Ncol, psym=2, title='Comparison to LeBlanc et al. (2005)',ytitle='D1+D2 Emission (Rayleighs)', xtitle='Distance from Europa (Europa radii)', color='green', xr=[-20,20], /ynozero
  cgplot, EWondisk_brown.field1, EW10N_area_profile, /overplot, color='black', psym=4
  cglegend, colors=['green', 'black'],psym=[2,1], titles=['LeBlanc et al. (2005)', 'Keck/HIRES'], length=0.01, symsize=0.1, /Box, Location=[0.15, 0.70], charsize=1.0, /Background, vspace=1



  P  = cglayout([1,1], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
  plate_scale = 0.358
  yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]

  cgPS_Open, filename = dir+'\Figures\leblanc05_compareEW10N.eps', $
    /ENCAPSULATED, xsize = 10, ysize = 10
  !P.font=1
  loadct, 3
  device, SET_FONT = 'Helvetica Bold', /TT_FONT

  title = 'Comparing Keck/HIRES 2022/09/29 to LeBlanc (2005) EW10N'

  cgplot, EWondisk_brown.field1, EW10Ncol, psym=9, title='EW10N Comparison to LeBlanc et al. (2005)',ytitle='D1+D2 Emission (Rayleighs)', xtitle='Distance from Europa (Europa radii)', color='green', xr=[-20,20], /ynozero
  cgplot, EWondisk_brown.field1, EW10N_area_profile, /overplot, color='black', psym=16, symsize=1.5
  cglegend, colors=['green', 'black'],psym=[9,16], titles=['LeBlanc et al. (2005)', 'Keck/HIRES'], length=0, symsize=1.5, /Box, Location=[0.15, 0.80], charsize=2.0, /Background, vspace=2

  cgPS_Close
  
  

stop

end