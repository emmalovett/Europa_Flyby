pro juno_flyby_EWframes, dir=dir

  restore, Dir+'\Processed\sun_subbed_europa.sav'
  
;  
;  for order_index = 2,2 do begin
;    if order_index eq 2 then begin
;      order = orders[order_index]
;      cube  = Na_cube
;      WL    = Na_WL
;    endif
;    if order_index eq 0 then begin
;      order = orders[order_index]
;      cube  = K_cube
;      WL    = K_WL
;    endif
;    if order_index eq 1 then begin
;      order = orders[order_index]
;      cube  = O_cube
;      WL    = O_WL
;    endif
  
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5 panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    cgplot, wl, total(newimg, 2), /xs, xr = xr, pos = p[*,1], xtickformat = '(A1)', $
      ytickformat = '(A1)', /noerase, title=labels[orientation]+' Sun-Subtracted'
    cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
    
;   1D spectra of sun- and io-subbed on the lower right hand side
    cgplot, column/1.e10, findgen(N_elements(newimg[index0:index1, *])), /ynozero, /noerase, pos = p[*,5], ys= 5, $
      xtitle = 'Na Column Density ('+cgsymbol('times')+'10!U10!N atoms / cm!U2!N)';, xr=[0,400]
    cgaxis, xaxis = 1, xtit = 'D1 + D2 (Rayleighs)', xr = [min(profile), max(profile)]
    
  cgPS_Close
  
  
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6 panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    cgimage, europa[index0:index1, *], pos = p[*,0], title=labels[orientation]+' Raw', ytitle = 'Rayleighs / '+cgsymbol('Angstrom')
    cgaxis, yaxis = 0, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01
;   raw 1D
    cgplot, wl, total(europa, 2), pos = p[*,3], xr = xr, /xs, ytickformat = '(A1)', /noerase
    cgaxis, yaxis = 0, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
    
;   sunlight subtracted in the middle
    cgimage, eurimg[index0:index1, *], pos = p[*,1], /noerase, title='Sun-Sub'
;   sun sub 1D
    cgplot, wl, total(eurimg, 2), pos = p[*,4], xr = xr, /xs, ytickformat = '(A1)', /noerase, yr=[min(total(eurimg, 2)), max(total(eurimg, 2))+1.e3]
    ;cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
    
;   sun and io subtracted on the right hand side
    cgimage, newimg[index0:index1, *], pos = p[*,2], yr = yr, /noerase, title='Sun-Io-Sub'
;   sun io sub 1D
    cgplot, wl, total(newimg, 2), pos = p[*,5], ytickformat = '(A1)', /noerase, xr = xr, /xs, yr=[min(total(eurimg, 2)), max(total(eurimg, 2))+1.e3]
    cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1, xtitle = 'Angstroms', xaxis = 0
    
  cgPS_Close
  
 ; endfor ; each order



end