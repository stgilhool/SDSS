pro show_grid_results_2, fnum

 f='/home/stgilhool/APOGEE/vsini_results/newtest/rfile'+strtrim(fnum,2)+'.fits'
r = mrdfits(f,0)
;pfeh = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0',
;'-0.5', '-0.0', $
pfeh = ['-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
        '+0.3', '+0.5']
feh = float(pfeh)

vstart = 5d0
vstep = 5d0
nvel = 10
vsini = dindgen(nvel)*vstep+vstart
nvel = 8
vsini = [4d0, 6d0, 10d0, 16d0, 24d0, 34d0, 46d0, 60d0, 76d0]


nteff=16
teff_vec = lindgen(nteff)*100+2400

lmin = dblarr(nteff)
lmin_idx = lonarr(nteff,2)

z_max = max(r)*1.05
z_min = min(r)*0.95



for t=0,nteff-1 do begin 
    p = reform(r[t,*,*]) 
    teff=teff_vec[t]
    cgsurface, p, feh, vsini, xtitle='FeH', ytitle='vsini', zr=[z_min,z_max], title='Teff = '+strtrim(teff,2), charsize=2.0 
    lmin[t] = min(p, lminidx, /nan)
    lmin_idx[t,*]=array_indices(p, lminidx)
endfor

globmin = min(lmin, globidx, /nan)
best_teff = teff_vec[globidx]

best_feh_idx = lmin_idx[globidx,0]
best_feh = feh[best_feh_idx]

best_vel_idx = lmin_idx[globidx,1]
best_vel = vsini[best_vel_idx]

print, "Best fit is: "
print, "Teff  = " + strtrim(best_teff,2)
print, "FeH   = " + strtrim(best_feh,2)
print, "vsini = " + strtrim(best_vel,2)
print, "chi2  = " + strtrim(globmin,2)
plot, vsini, r[globidx,best_feh_idx,*], xtitle='vsini', ytitle='chi2', title='Teff = '+strtrim(best_teff,2)+' | FeH = '+strtrim(best_feh,2), charsize=2.0

stop 


 f='/home/stgilhool/APOGEE/vsini_results/mptest/rfile'+strtrim(fnum,2)+'.fits'
r = mrdfits(f,0)
z_max = max(r)*1.05
z_min = min(r)*0.95



for t=0,nteff-1 do begin 
    p = reform(r[t,*,*]) 
    teff=teff_vec[t]
    cgsurface, p, feh, vsini, xtitle='FeH', ytitle='vsini', zr=[z_min,z_max], title='Teff = '+strtrim(teff,2), charsize=2.0, /nan 
    lmin[t] = min(p, lminidx, /nan)
    lmin_idx[t,*]=array_indices(p, lminidx)
endfor

globmin = min(lmin, globidx, /nan)
best_teff = teff_vec[globidx]

best_feh_idx = lmin_idx[globidx,0]
best_feh = feh[best_feh_idx]

best_vel_idx = lmin_idx[globidx,1]
best_vel = vsini[best_vel_idx]

print, "Best fit is: "
print, "Teff  = " + strtrim(best_teff,2)
print, "FeH   = " + strtrim(best_feh,2)
print, "vsini = " + strtrim(best_vel,2)
print, "chi2  = " + strtrim(globmin,2)
plot, vsini, r[globidx,best_feh_idx,*], xtitle='vsini', ytitle='chi2', title='Teff = '+strtrim(best_teff,2)+' | FeH = '+strtrim(best_feh,2), charsize=2.0

stop

end
