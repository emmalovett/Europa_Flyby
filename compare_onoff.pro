pro compare_onoff
  
  dir             = 'Z:\DATA\Keck\Europa Na\HIRES_20220928'
  Europa_frames   = string(indgen(37)+127, format='(I4.4)')
  
    img1                = [mrdfits(Dir+'\hires' + Europa_frames[3] + '.fits', 3, header, /fscale), $
                           mrdfits(Dir+'\hires' + Europa_frames[3] + '.fits', 2, header, /fscale), $
                           mrdfits(Dir+'\hires' + Europa_frames[3] + '.fits', 1, header, /fscale)]
   img2                 = [mrdfits(Dir+'\hires' + Europa_frames[6] + '.fits', 3, header, /fscale), $
                           mrdfits(Dir+'\hires' + Europa_frames[6] + '.fits', 2, header, /fscale), $
                           mrdfits(Dir+'\hires' + Europa_frames[6] + '.fits', 1, header, /fscale)]                       
    img1 = rotate(transpose(img1),7)
    img2 = rotate(transpose(img2),7)
    
    
    window, 0, title='EW ON DISK : COMPARE MAX VALUE'
    cgplot, total(img1[500:1000,450:480], 1, /nan)
    cgtext, 0.8, 0.9, 'max: '+strcompress(max(img1[500:1000,450:480])), color='red'
    window, 1, title='EW ON DISK'
    cgimage, img1[500:1000,450:480], minv=900, maxv=1100
    window, 2, title='EW 10 N : COMPARE MAX VALUE'
    cgplot, total(img2[500:1000,450:480], 1, /nan)
    cgtext, 0.8, 0.9, 'max: '+strcompress(max(img2[500:1000,450:480])), color='red'
    window, 3, title='EW 10 N'
    cgimage, img2[500:1000,450:480], minv=900, maxv=1100
    
    stop
  
  
end