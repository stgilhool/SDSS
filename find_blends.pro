function corr_calc, v1, v2

v1 = reform(v1)
v1mean = mean(v1, /nan)
v1diff = v1 - v1mean

v2 = reform(v2)
v2mean = mean(v2, /nan)
v2diff = v2 - v2mean

r = total(v1diff * v2diff, /double)/(sqrt(total(v1diff^2, /double)) * sqrt(total(v2diff^2, /double)))

return, r

end



pro find_blends

restore, 'checklines.dat'

corr_x0b = [720L, 2000L] + bmin

corr_x0g = [300, 400, 510, 1090, 1170, 1290, 1350, 1390, 1580, $
            1610, 1790, 1810, 2080, 2110, 2130, 2200] + gmin

corr_x0r = [12, 100, 120, 150, 890, 910, 950, 1050, 1080] + rmin

x0 = [corr_x0b, corr_x0g, corr_x0r]


n_lines = n_elements(x0)

pcorr_arr = dblarr(n_lines)
bestpcorr_arr = dblarr(n_lines)

bestx0  = dblarr(n_lines)

maxshift=18
!p.multi=[0,1,2]
foreach pix0, x0, xidx do begin
    
    pcorr_i = dblarr(maxshift)
    
    for i = 0, maxshift-1 do begin

        fp = pix0 - (npix_win/2) + (i-(maxshift/2))


        lp = pix0 + (npix_win/2) + (i-(maxshift/2))

        ew = total((1d0-spectra[fp:lp,*]), 1,  /double)


        pcorr = corr_calc(ew, feh_afk)

        pcorr_i[i] = pcorr
    endfor

    bestshift_idx = where(pcorr_i eq max(pcorr_i))
    bestx0[xidx] = pix0 + bestshift_idx - (maxshift/2)
    
    fp = pix0 - (npix_win/2)
    
    
    lp = pix0 + (npix_win/2)
    
    ew = total((1d0-spectra[fp:lp,*]), 1,  /double)
    
    
    pcorr = corr_calc(ew, feh_afk)

    pcorr_arr[xidx] = pcorr
    bestpcorr_arr[xidx] = pcorr_i[bestshift_idx]


    plot, lindgen(npix_win) + fp, spectra[fp:lp, 43], /xs, yr = [0.6,1.1]
    
    fp = bestx0[xidx] - (npix_win/2)

    lp = bestx0[xidx] + (npix_win/2)
    
    plot, lindgen(npix_win) + fp, spectra[fp:lp, 43], /xs, yr = [0.6,1.1]
    ;dxstop

endforeach

print, pcorr_arr
print, bestpcorr_arr
print, ''
print, bestpcorr_arr - pcorr_arr

print, ''
print, x0
print, bestx0

print, ''
print, bestx0 - x0  

xx = rebin(lindgen(nx), nx, n)

window, 0, xsize = 1000, ysize=1000

terminate = 0
iter = 0
while terminate eq 0 do begin
    
    xrmin = iter*50 + bmin
    xrmax = xrmin + 100
    plot, xx, spectra, ps=3, yr=[0.3,1.3], xr=[xrmin,xrmax]
    dxstop
    
    iter++
    if (xrmax + 50) gt nx then terminate = 1 
    
endwhile


stop

end
