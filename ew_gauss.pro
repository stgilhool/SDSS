pro show_spectra, pixnum

common block1, spectra, wl_grid, gfitco_all

;HARDCODE
npix_window = 11
npix_halfwin = npix_window/2
nspec = 25

window_x = dindgen(npix_window) - npix_halfwin + pixnum

spec = spectra[window_x, *]



for snum = 0, nspec - 1 do begin
    
    fitcos = gfitco_all[long(pixnum),snum,*]
    fitcos = reform(fitcos)
    fit = gaussian(window_x, fitcos)

    print, fitcos


    plot, window_x, spec[*,snum], /ys, /xs, yr=[0.5,1.1]
    oplot, window_x, fit, color=200
    
    dxstop

endfor


end


pro ew_gauss

common block1, spectra, wl_grid, gfitco_all

;restore, 'ew_gausstest.sav';
;restore, 'ew_gausstest2.sav'
restore, 'ew_gausstest3.sav'

; for each pixel and spectrum number, check if a gaussian amplitude
; was returned
fin_mtx = finite(gfitco_all[*,*,0])

; now sum up the columns
fin_arr = total(fin_mtx, 2)

; check status
;gdstatus_idx = where(status_fit gt 0, ngdstat)

;gdstatus = status_fit eq 1

gdcont = (gfitco_all[*,*,3] ge 0.8) and (gfitco_all[*,*,3] le 1.2) and (gfitco_all[*,*,3] ne 1.0d0)

;gdall = gdstatus and gdcont
gdall = gdcont

;gdstat_col = total(gdstatus, 2)
gdstat_col = total(gdall, 2)



; let's require that 25/25 have valid results
;gdpix_idx = where(fin_arr ge 24, ngd)
gdpix_idx = where(gdstat_col ge 10, ngd)
gdpix_idx = [463,1053,1132,1297,1348,1473,1531,1546,2177,2190,2244,2463,2657,2772,2785]
ngd = n_elements(gdpix_idx) 
; Get median, mean and stddev of chi_fit columns
med_chi = median(chi_fit, dimension=2)
mean_chi = mean(chi_fit, dimension=2, /nan)
sig_chi = stddev(chi_fit, dimension=2, /nan)

med_gfitco = median(gfitco_all, dimension=2)
mean_gfitco = mean(gfitco_all, dimension=2, /nan)
sig_gfitco = stddev(gfitco_all, dimension=2, /nan)


if ngd gt 0 then begin
    plot, gdpix_idx, med_chi[gdpix_idx], ps=6

    stop
endif else message, "requirement is too strict"

stop

end
