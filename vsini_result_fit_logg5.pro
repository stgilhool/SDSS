pro vsini_result_fit_logg5

nspec = 1350

; Grab APOGEE RESULTS
apg = readin_apo(hdu=1, nfiles=nspec)

;outdata = []

for fnum = 0, nspec-1 do begin 
    
    apghead = *(apg[fnum].header)
    teff_apg = fxpar(apghead, 'RVTEFF')
    feh_apg = fxpar(apghead, 'RVFEH')
    


    f='/home/stgilhool/APOGEE/vsini_results/logg5/rfile'+strtrim(fnum,2)+'.fits'

    if file_test(f) then r = mrdfits(f,0) else begin
        ; outstr = {TEFF_VEC:!values.d_nan, $
;                   FEH_VEC:!values.d_nan, $
;                   VSINI_VEC:!values.d_nan, $
;                   VSINI_OVER_VEC:!values.d_nan, $
;                   TEFF_GRID:!values.d_nan, $
;                   TEFF_POLY:!values.d_nan, $
;                   TEFF_SPLINE:!values.d_nan, $
;                   TEFF_APG:!values.d_nan, $
;                   TEFF_BEST:!values.d_nan, $
;                   FEH_GRID:!values.d_nan, $
;                   FEH_POLY:!values.d_nan, $
;                   FEH_SPLINE:!values.d_nan, $
;                   FEH_APG:!values.d_nan, $
;                   FEH_BEST:!values.d_nan, $
;                   VSINI_GRID:!values.d_nan, $
;                   VSINI_POLY:!values.d_nan, $
;                   VSINI_SPLINE:!values.d_nan, $
;                   VSINI_BEST:!values.d_nan, $
;                   CHI2_GRID:!values.d_nan, $
;                   CHI2_POLY:!values.d_nan, $
;                   CHI2_SPLINE:!values.d_nan, $
;                   CHI2_BEST:!values.d_nan, $
;                   CHIVEL_GRID:!values.d_nan, $
;                   CHIVEL_POLY:!values.d_nan, $
;                   CHIVEL_SPLINE:!values.d_nan, $
;                   MODEL_BEST:!values.d_nan, $
;                   TEFF_NE_FLAG:!values.d_nan, $
;                   FEH_NE_FLAG:!values.d_nan, $
;                   TEFF_APG_FLAG:!values.d_nan, $
;                   FEH_APG_FLAG:!values.d_nan $
;                  }
        ; THIS IS WRONG
        outdata[fnum] = outstr
        continue
    endelse
    ;pfeh = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0',
    ;'-0.5', '-0.0', $
    ;pfeh = ['-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
    pfeh = ['-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
            '+0.3', '+0.5']
    feh = float(pfeh)
    nfeh = n_elements(feh) 
    
    ;vstart = 5d0
    ;vstep = 5d0
    ;vsini = dindgen(nvel)*vstep+vstart
    vsini = [4d0, 6d0, 10d0, 16d0, 24d0, 34d0, 46d0, 60d0, 76d0]
    vsini_over = dindgen(264)/4d0 + 4d0
    nvel = n_elements(vsini)
    nvel_over = n_elements(vsini_over)
    
    nteff=16
    teff_vec = lindgen(nteff)*100+2400
    
    lmin = dblarr(nteff)
    lminover = dblarr(nteff)
    lminspline = dblarr(nteff)
    lmin_idx = lonarr(nteff,2)
    lminover_idx = lonarr(nteff,2)
    lminspline_idx = lonarr(nteff,2)
    
    r_over = dblarr(nteff, nfeh,nvel_over) 
    r_spline = dblarr(nteff, nfeh,nvel_over) 
    err_over = dblarr(nteff, nfeh,nvel_over) 
    err_spline = dblarr(nteff, nfeh,nvel_over) 
    
    
    z_max = max(r)*1.05
    z_min = min(r)*0.95
    
    
    
    for t=0,nteff-1 do begin 
        p = reform(r[t,*,*]) 
        teff=teff_vec[t]
        for met = 0, nfeh-1 do begin
            vel_fit_coeff = poly_fit(vsini, p[met,*], 3, /double, chisq = vel_fit_chi) 
            vel_fit = poly(vsini_over, vel_fit_coeff)
            vel_spline = interpol(p[met,*], vsini, vsini_over, /spline)
            vel_linterp = interpol(p[met,*], vsini, vsini_over)
            r_over[t,met,*]= vel_fit
            err_over[t,met,*]=vel_linterp-vel_fit
            r_spline[t,met,*]= vel_spline
            err_spline[t,met,*]=vel_linterp-vel_spline
        endfor
        p_over = reform(r_over[t,*,*])
        p_spline = reform(r_spline[t,*,*])
        
        ;cgsurface, p_over, feh, vsini, xtitle='FeH', ytitle='vsini', zr=[z_min,z_max], title='Teff = '+strtrim(teff,2), charsize=2.0 
        lminover[t] = min(p_over, lminoveridx, /nan)
        lminspline[t] = min(p_spline, lminsplineidx, /nan)
        lmin[t] = min(p, lminidx, /nan)
        
        lmin_idx[t,*]=array_indices(p, lminidx)
        lminover_idx[t,*]=array_indices(p_over, lminoveridx)
        lminspline_idx[t,*]=array_indices(p_spline, lminsplineidx)
    endfor
    
    
    globmin = min(lmin, globidx, /nan)
    globmin_over = min(lminover, globoveridx, /nan)
    globmin_spline = min(lminspline, globsplineidx, /nan)
    
    best_teff = teff_vec[globidx]
    bestover_teff = teff_vec[globoveridx]
    bestspline_teff = teff_vec[globsplineidx]
    
    best_feh_idx = lmin_idx[globidx,0]
    best_feh = feh[best_feh_idx]
    bestover_feh_idx = lminover_idx[globoveridx,0]
    bestover_feh = feh[bestover_feh_idx]
    bestspline_feh_idx = lminspline_idx[globsplineidx,0]
    bestspline_feh = feh[bestspline_feh_idx]
    
    best_vel_idx = lmin_idx[globidx,1]
    best_vel = vsini[best_vel_idx]
    bestover_vel_idx = lminover_idx[globoveridx,1]
    bestover_vel = vsini_over[bestover_vel_idx]
    bestspline_vel_idx = lminspline_idx[globsplineidx,1]
    bestspline_vel = vsini_over[bestspline_vel_idx]
    
    best_globmin = min([globmin, globmin_over, globmin_spline])
    
    ;Determine which method gives the best result
    case 1 of
        
        best_globmin eq globmin: begin
            best_best_teff = best_teff
            best_best_feh = best_feh
            best_best_vel = best_vel
            best_best_chi2 = globmin
            model_best = "GRID"
        end

        best_globmin eq globmin_over and best_globmin eq globmin_spline: begin
            best_best_teff = bestover_teff
            best_best_feh = bestover_feh
            best_best_vel = bestover_vel
            best_best_chi2 = globmin_over
            model_best = "BOTH"
        end

        best_globmin eq globmin_over: begin
            best_best_teff = bestover_teff
            best_best_feh = bestover_feh
            best_best_vel = bestover_vel
            best_best_chi2 = globmin_over
            model_best = "POLY"
        end

        best_globmin eq globmin_spline: begin
            best_best_teff = bestspline_teff
            best_best_feh = bestspline_feh
            best_best_vel = bestspline_vel
            best_best_chi2 = globmin_spline
            model_best = "SPLINE"
        end

        else: begin
            best_best_teff = best_teff
            best_best_feh = best_feh
            best_best_vel = best_vel
            best_best_chi2 = globmin
            model_best = "NONE"
        end
    endcase

    ;Check to see if we get different FEH's and TEFF's
    if best_teff ne best_best_teff or bestover_teff ne best_best_teff $
      or bestspline_teff ne best_best_teff then teff_ne_flag = 1 $
      else teff_ne_flag = 0
    if best_feh ne best_best_feh or bestover_feh ne best_best_feh $
      or bestspline_feh ne best_best_feh then feh_ne_flag = 1 $
      else feh_ne_flag = 0

    ;Check if APOGEE result is consistent
    if abs(best_best_teff-teff_apg) ge 200 then teff_apg_flag = 1 $
      else teff_apg_flag = 0
    if abs(best_best_feh-feh_apg) ge 1 then feh_apg_flag = 1 $
      else feh_apg_flag = 0
    
    ;Check smoothness of chi2 surface somehow
        


    outstr = {TEFF_VEC:teff_vec, $
              FEH_VEC:feh, $
              VSINI_VEC:vsini, $
              VSINI_OVER_VEC:vsini_over, $
              TEFF_GRID:best_teff, $
              TEFF_POLY:bestover_teff, $
              TEFF_SPLINE:bestspline_teff, $
              TEFF_APG:teff_apg, $
              TEFF_BEST:best_best_teff, $
              FEH_GRID:best_feh, $
              FEH_POLY:bestover_feh, $
              FEH_SPLINE:bestspline_feh, $
              FEH_APG:feh_apg, $
              FEH_BEST:best_best_feh, $
              VSINI_GRID:best_vel, $
              VSINI_POLY:bestover_vel, $
              VSINI_SPLINE:bestspline_vel, $
              VSINI_BEST:best_best_vel, $
              CHI2_GRID:globmin, $
              CHI2_POLY:globmin_over, $
              CHI2_SPLINE:globmin_spline, $
              CHI2_BEST:best_globmin, $
              CHIVEL_GRID:reform(r[globidx,best_feh_idx,*]), $
              CHIVEL_POLY:reform(r_over[globoveridx,bestover_feh_idx,*]), $
              CHIVEL_SPLINE:reform(r_spline[globsplineidx,bestspline_feh_idx,*]), $
              MODEL_BEST:model_best, $
              TEFF_NE_FLAG:teff_ne_flag, $
              FEH_NE_FLAG:feh_ne_flag, $
              TEFF_APG_FLAG:teff_apg_flag, $
              FEH_APG_FLAG:feh_apg_flag $
             }
    
    if fnum eq 0 then outdata = replicate(outstr, nspec)
    outdata[fnum] = outstr
    ;outdata = [outdata, outstr]
    
    ;print, "Best fit is: "
    ;print, "Teff  = " + strtrim(best_teff,2)
    ;print, "FeH   = " + strtrim(best_feh,2)
    ;print, "vsini = " + strtrim(best_vel,2)
    ;print, "chi2  = " + strtrim(globmin,2)
    ;print, "--------------------------------------"
    ;print, "Teff  = " + strtrim(bestover_teff,2)
    ;print, "FeH   = " + strtrim(bestover_feh,2)
    ;print, "vsini = " + strtrim(bestover_vel,2)
    ;print, "chi2  = " + strtrim(globmin_over,2)
    ;plot, vsini, r[globidx,best_feh_idx,*], xtitle='vsini', ytitle='chi2', title='Teff = '+strtrim(best_teff,2)+' | FeH = '+strtrim(best_feh,2), charsize=2.0
    ;oplot, vsini_over, r_over[globoveridx,bestover_feh_idx,*], color=200
    ;stop 
   
endfor
;Write results
outfile = '/home/stgilhool/APOGEE/vsini_results/logg5/rfile_master.fits'

mwrfits, outdata, outfile, /create

filenum = lindgen(n_elements(outdata)) 
plot, filenum, outdata.vsini_best, xs = 2, ps=6, title = 'Vsini for APOGEE M dwarfs', xtitle = 'File Number', ytitle = 'vsini (km/s)', yr = [0,65], xr=[0,1350]

stop

end
