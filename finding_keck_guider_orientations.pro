pro finding_keck_guider_orientations




dir = 'Z:\DATA\Keck\Europa Na\HIRES_20220928\MAGIQ files'
dir = 'Z:\DATA\Keck\Io Eclipse HIRES\08-09-2023\Guider'
ewfile  = dir+'\hiresSlit000_0066.fits'
nsfile  = dir+'\hiresSlit000_0075.fits'
guiders = FILE_SEARCH(dir, '*.fits')

for i = 21, N_ELEMENTS(guiders)-1 DO BEGIN ;64
  frame = mrdfits(guiders[i], 0, header, /fscale) 
  frame = rotate(frame, 5)
  frame = rotate(frame, 2)
  
  window, 6, xs=512, ys=512, title=STRMID(guiders[i], 50)
  cgimage, frame, minv=500, maxv=1600
  print, sxpar(header, 'ROTPOSN') - 90.
  
  
  stop
endfor
stop
nsframe = mrdfits(nsfile, 0, nsheader, /fscale) 
ewframe = mrdfits(ewfile, 0, ewheader, /fscale)



window, 0, xs=512, ys=512, title=STRMID(nsfile, 50)+'  NS'
cgimage, nsframe, minv=0.75*mean(nsframe), maxv=1.5*mean(nsframe)

window, 1, xs=512, ys=512, title=STRMID(ewfile, 50)+'  EW'
cgimage, ewframe, minv=0.75*mean(ewframe), maxv=1.5*mean(ewframe)


; transpose first. remember, io is below jupiter's equator


transposed_NS = ROTATE(nsframe, 5)
window, 2,  xs=512, ys=512, title=STRMID(nsfile, 50)+'  TRANSPOSED '
cgimage, transposed_NS, minv=0.75*mean(nsframe), maxv=1.5*mean(nsframe) 

transposed_EW = ROTATE(ewframe, 5)
window, 3,  xs=512, ys=512, title=STRMID(ewfile, 50)+'  TRANSPOSED '
cgimage, transposed_EW, minv=0.75*mean(ewframe), maxv=1.5*mean(ewframe)

; now rotate.

rotated_NS = ROTATE(transposed_NS, 1)
window, 4,  xs=512, ys=512, title=STRMID(nsfile, 50)+'  ROTATED '       ; <------------- USE THIS FOR ORIENTATION YAAA
cgimage, rotated_NS, minv=0.75*mean(nsframe), maxv=1.5*mean(nsframe)


rotated_EW = transposed_EW
window, 5,  xs=512, ys=512, title=STRMID(ewfile, 50)+'  ROTATED '       ; <------------- USE THIS FOR ORIENTATION YAAA
cgimage, rotated_EW, minv=0.75*mean(ewframe), maxv=1.5*mean(ewframe)


;;;; why are these different???

testpath1 = dir+'\hiresSlit000_0205.fits'
testpath2 = dir+'\hiresSlit000_0206.fits'
testimg1  = mrdfits(testpath1, 0, header1, /fscale)
testimg2  = mrdfits(testpath2, 0, header2, /fscale)

print, sxpar(header1, 'ROTPOSN') - 90.
print, sxpar(header2, 'ROTPOSN') - 90.


testimg1 = rotate(testimg1, 5)
testimg1 = rotate(testimg1, 1)
                  
testimg2 = rotate(testimg2, 5)
testimg2 = rotate(testimg2, 1)

window, 6,  xs=512, ys=512
tv, bytscl(testimg1)

window, 7,  xs=512, ys=512
tv, bytscl(testimg2)
stop



;;; now let's test this rotation on other files. first, an EW Io frame (hires99)

;path = 'Z:\DATA\Keck\Europa Na\HIRES_20220928\just RAW data\hires0200.fits'
;
;ioheader = headfits(path)
;print, sxpar(headfits(path), 'ROTPOSN')-90.
;cspice_UTC2ET,   sxpar(ioheader, 'DATE-OBS') + ' '+ sxpar(ioheader, 'UTC'), io_ET
;start_time = io_et - sxpar(ioheader, 'EXPTIME')
;
;
;guiders          = FILE_SEARCH(dir, 'hiresslit*'+'*.fits')
;guider_times     = []
;
;FOR i = 0, N_elements(guiders)-1 DO BEGIN
;  guider_header  = headfits(guiders[i])
;  cspice_UTC2ET,   sxpar(guider_header, 'DATE-OBS') + ' '+ sxpar(guider_header, 'UTC'), guider_ET
;  IF guider_et LT io_et AND guider_et GT START_time then begin
;    print, strmid(guiders[i],50)+ ' ' + sxpar(headfits(guiders[i]), 'DATE-OBS') + ' '+ sxpar(headfits(guiders[i]), 'UTC')+'   ', i
;  ENDIF
;ENDFOR
;
;window, 6,  xs=512, ys=512, title=STRMID(path, 50)
;cgimage, ioimg, minv=0.75*mean(ioimg), maxv=1.5*mean(ioimg)
;
;window, 7,  xs=512, ys=512, title=STRMID(path, 50)+'  ROTATED'
;cgimage, ROTATE(ioimg, 5), minv=0.75*mean(ioimg), maxv=1.5*mean(ioimg)
;stop
nsjupfile = dir+'\hiresSlit000_0275.fits'             ; 100 for NS Jupiter, 275 for EW Jupiter
jupimg  = mrdfits(nsjupfile, 0, header, /fscale)

window, 6,  xs=512, ys=512, title=STRMID(nsjupfile, 50)
tv, bytscl(jupimg)


jupimg  = rotate(jupimg, 5)
;jupimg  = rotate(jupimg, 1)

window, 7,  xs=512, ys=512, title=STRMID(nsjupfile, 50)+'  ROTATED'
tv, bytscl(jupimg)

; test for ew jupiter?



 

stop










end