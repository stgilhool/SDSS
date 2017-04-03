; Look at specific lines in APOGEE spectra
pro line_viewer, wl_ctr, width

common viewer_info, wl_grid, wl_grid_over, spectra, ln_num
pmultisav = !p.multi
!p.multi=0
; Find the (oversampled) pixel where wl is 
; nearest to the supplied central feature wl
;diff_wl = abs(wl_grid_over - wl_ctr)
diff_wl = abs(wl_grid - wl_ctr)
ctr_idx = where(diff_wl eq min(diff_wl))

; Make vectors of width width, centered about ctr
dwin = width/2
pix_idx = lindgen(width)-dwin+ctr_idx[0]
winspec = spectra[ctr_idx-dwin:ctr_idx+dwin, *]
winwl = wl_grid[ctr_idx-dwin:ctr_idx+dwin]


; Optionally mask pixels


; Get some info about the window
nspec = n_elements(spectra[0,*]) 
winmed = median(winspec, dimension = 2)
winmean = mean(winspec, dimension = 2, /nan)
winsig = stddev(winspec, dimension=2, /nan)

; Plot

; for the vertical line at wl_ctr
;diff_wl_over = abs(wl_grid_over - wl_ctr)
;ctr_idx_over = where(diff_wl_over eq min(diff_wl_over))
;vline_x = wl_grid_over[ctr_idx_over]

; for the plot range
yplotmin = min(winmed)-5*mean(winsig)
yplotmax = max(winmed)+5*mean(winsig)

; make the plots
;plot, winwl/1d4, winmed, ps = 6, symsize = 0.5, yr = [yplotmin,
;yplotmax], /xs, $
plot, pix_idx, winmed, ps = 6, symsize = 0.5, yr = [yplotmin, yplotmax], /xs, $
  title = "Feature number: " + strtrim(ln_num,2) + " at wl = " + strtrim(wl_ctr,2), $
  charsize = 1.5
;oplot, winwl/1d4, winmed, ps = 6, symsize = 0.5, color=99999
oplot, pix_idx, winmed, ps = 6, symsize = 0.5, color=99999

;oploterror, winwl/1d4, winmean, winsig, /nohat
oploterror, pix_idx, winmean, winsig, /nohat
;oplot, winwl/1d4, winmean, color=200
oplot, pix_idx, winmean, color=200
oplot, [!x.crange[0], !x.crange[1]], [1d0, 1d0], linestyle = 2
;oplot, [wl_ctr/1d4, wl_ctr/1d4], [!y.crange[0], !y.crange[1]],
;linestyle = 2
oplot, [ctr_idx, ctr_idx], [!y.crange[0], !y.crange[1]], linestyle = 2
;for spec = 0, nspec-1 do begin
;    oplot, winwl, winspec[*,spec], ps=3
;endfor
dxstop
!p.multi = pmultisav

end



pro view_lines

common viewer_info, wl_grid, wl_grid_over, spectra, ln_num

; Constants
wl0 = 4.179d0
wl1 = 6d-6
npix = 8575L
oversamp = 7L
width = 35L

; Read in line list
lines = readin_lines()
nlines = n_elements(lines) 

; Read in spectra
apg = mrdfits('/home/stgilhool/APOGEE/absk_corr/absk_corr_init.fits',1)
spectra = apg.spec
nspec = n_elements(apg) 

; WLs
logwl = dindgen(npix)*wl1 + wl0
wl_grid = 10d0^logwl
; Oversample the wavelength grid
pix_vec = lindgen(npix)
pix_vec_over = (lindgen(npix*oversamp)-(oversamp/2))/oversamp
logwl_over = interpol(logwl, pix_vec, pix_vec_over)
wl_grid_over = 10d0^logwl_over

; View the lines
ti_lines_idx = where(lines.species eq 'Ti', ticnt)

if ticnt gt 0 then ti_lines = lines[ti_lines_idx] $ 
else message, 'Error finding Ti lines'

lines_in = ti_lines
nlines_in = n_elements(ti_lines) 

window, 0, xs = 1600, ys = 900

for ln_num = 0, nlines_in-1 do begin


    wl_ctr = lines_in[ln_num].wl_ctr
    
    line_viewer, wl_ctr, width

endfor

stop


end
