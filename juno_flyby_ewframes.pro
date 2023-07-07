pro juno_flyby_EWframes, dir=dir

  restore, Dir+'\Processed\sun_subbed_europa.sav'
  
  
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
  
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sun-sub comparison figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  P  = cglayout([3,2], ygap = 0., oxmargin = [10,10], oymargin = [8,5], xgap = 0.)
  plate_scale = 0.358
  yr = [-(suncol*plate_scale/ang_radius), ((s[2]-suncol)*plate_scale/ang_radius)]
  
  cgPS_Open, filename = dir+'\Figures\'+order.name+'_'+labels[orientation]+'_subtraction_process.eps', $
    /ENCAPSULATED, xsize = 20, ysize = 8
  !P.font=1
  loadct, 3
  device, SET_FONT = 'Helvetica Bold', /TT_FONT

  title = 'HIRES 2022-09-29 : '+labels[orientation]

; non-sunlight subtracted on the upper left hand side
  cgimage, europa[index0:index1, *], pos = p[*,0], title=labels[orientation]+' Raw', ytitle = 'Rayleighs / ' + cgsymbol('Angstrom')
    
; sunlight subtracted on the lower left hand side
  cgimage, eurimg[index0:index1, *], /xs, xr = xr, pos = p[2,*], /axes, /noerase, title=labels[orientation]+' Sun-Subtracted'
  cgaxis, yaxis = 0, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1
  stop
    

  cgplot, wl, continuum, /overplot, color='red'

  cgtext, .29, .58, 'Above disk', color = 'black', /normal
  cgtext, .29, .55, 'Fit Reflectance', color = 'red', /normal

  axis_format = {XTicklen:-.01, yticklen:-0.01 }
  cgimage, orientations[index0:index1,*,orientation], minv=0.75*mean(orientations[*,*,orientation]), maxv=3.0*mean(orientations[*,*,orientation]), $
    /axes, xr = xr, pos = p[*,2], /noerase, $
    ytit = 'Europa Radii', xtitle = 'Angstroms', AXKEYWORDS = axis_format, yr = yr

  cgcolorbar, POSITION=[p[0,0], 0.045, p[2,0], 0.065], range = minmax(orientations[*,*,orientation])

  ; sunlight subtracted on the right hand side
  cgplot, wl, total(newimg, 2), /xs, xr = xr, pos = p[*,1], xtickformat = '(A1)', $
    ytickformat = '(A1)', /noerase, title=labels[orientation]+' Sun-Subtracted'
  cgaxis, yaxis = 1, ytitle = 'Rayleighs / ' + cgsymbol('Angstrom'), ystyle=1

  axis_format = {XTicklen:-.01, yticklen:-0.01, ystyle:5}
  cgimage, newimg[index0:index1, *], minv=-100, maxv=100, /axes, xr = xr, pos = p[*,3], yr = yr, /noerase, $
    xtitle = 'Angstroms', AXKEYWORDS = axis_format
  cgaxis, yaxis = 1, ytit = 'Europa Radii', yr = yr, ystyle=1, yticklen=-0.01

  cgcolorbar, POSITION=[p[2,0], 0.045, 0.979, 0.065], range = [-5.e3,5.e3]
  cgps_Close














  endfor ; each order



end