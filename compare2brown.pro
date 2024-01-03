pro compare2brown


  RESTORE, 'Z:\DATA\Keck\Europa Na\HIRES_20220928\Processed\sun_subbed_europa.sav'
  
  brownpath = 'C:\Users\elovett\EuropaResearch\Europa_Flyby\brown_data.csv'
  
  READCOL, brownpath, ewx, ewy, nsx, nsy, skipline=1
  
  cgPS_Open, filename = Dir+'\Figures\EW_brown_compare.eps', /ENCAPSULATED, xsize = 7.5, ysize = 6
  !P.font=1
  loadct, 3
  device, SET_FONT = 'Helvetica Bold', /TT_FONT

  title = 'EW Brightness Profile'

  cgplot, ewx, ewy, psym=9, color='maroon', xtitle='Europa Radii', ytitle='Emission D1 + D2 (Rayleighs)', xr=[-25,25], title='Europa Na EW On-Disk', thick=2
  cgplot, findgen(n_elements(profile))-(n_elements(profile)/2.)+2, EWprofile, /overplot, xtickformat='(A1)', xticks=1, xminor=1, color='hot pink', thick=2                                                                        
  cglegend, SymColors=['maroon','hot pink'], PSyms=[9, 0], linestyles=1, symsize=1, Titles=["LeBlanc (2005) Brown's Data","Lovett+ (in prep)"],$
    Length=0.0, /Box, Location=[0.65, 0.7], charsize=1.2, /Background, vspace=1  
  cgps_Close
  
  
  
  cgPS_Open, filename = Dir+'\Figures\NS_brown_compare.eps', /ENCAPSULATED, xsize = 7.5, ysize = 6
    !P.font=1
    loadct, 3
    device, SET_FONT = 'Helvetica Bold', /TT_FONT
    
    title = 'NS Brightness Profile'
    cgplot, nsx, nsy, psym=9, color='maroon', xtitle='Europa Radii', ytitle='Emission D1 + D2 (Rayleighs)', xr=[-25,25], title='Europa Na NS On-Disk', thick=2
    cgplot, findgen(n_elements(profile))-(n_elements(profile)/2.)+2, NSprofile, /overplot, xtickformat='(A1)', xticks=1, xminor=1, color='hot pink', thick=2
    cgps_Close
  
  stop



end