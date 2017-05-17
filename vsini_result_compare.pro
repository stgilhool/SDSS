pro vsini_result_compare

nspec = 1350

; Grab APOGEE RESULTS
apg = readin_apo(hdu=1, nfiles=nspec)

; Aspcap
asp = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits',1)

;outdata = []

for fnum = 0, nspec-1 do begin 
    
    apghead = *(apg[fnum].header)
    teff_apg = fxpar(apghead, 'RVTEFF')
    feh_apg = fxpar(apghead, 'RVFEH')
    


    f='/home/stgilhool/APOGEE/vsini_results/logg45_with_fit_full/rfile'+strtrim(fnum,2)+'.fits'
    f5='/home/stgilhool/APOGEE/vsini_results/logg5_with_fit_full/rfile'+strtrim(fnum,2)+'.fits'
    

    if file_test(f) then begin
        r = mrdfits(f,0) 
        rbest = mrdfits(f,1)
        r5 = mrdfits(f5,0) 
        rbest5 = mrdfits(f5,1)
        
    endif else begin
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
        outbest[fnum] = beststr
        nbest_entries = n_elements(beststr) 
        outbest[fnum].valid_flag = 0
        continue
    endelse
    ;pfeh = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0',
    ;'-0.5', '-0.0', $
    ;pfeh = ['-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
    pfeh = ['-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
            '+0.3', '+0.5']
    feh = float(pfeh)
    nfeh = n_elements(feh) 
    
    ; Initialize
    vsini = [4d0, 6d0, 10d0, 16d0, 24d0, 34d0, 46d0, 60d0, 76d0]
    vsini_over = dindgen(264)/4d0 + 4d0
    nvel = n_elements(vsini)
    nvel_over = n_elements(vsini_over)
    
    nteff=16
    teff_vec = lindgen(nteff)*100+2400
    
    ; More Initialize
    lmin = dblarr(nteff)
    lminspline = dblarr(nteff)
    lmin_idx = lonarr(nteff,2)
    lminspline_idx = lonarr(nteff,2)
    
    r_spline = dblarr(nteff, nfeh,nvel_over) 
    ; logg 5
    lmin5 = dblarr(nteff)
    lminspline5 = dblarr(nteff)
    lmin_idx5 = lonarr(nteff,2)
    lminspline_idx5 = lonarr(nteff,2)
    
    r_spline5 = dblarr(nteff, nfeh,nvel_over) 
        
        
    for t=0,nteff-1 do begin 
        ; get results at a given temp
        p = reform(r[t,*,*]) 
        p5 = reform(r5[t,*,*]) 
        teff=teff_vec[t]
        for met = 0, nfeh-1 do begin
            ; get chi2 vector at given temp and feh
            vel_spline = interpol(p[met,*], vsini, vsini_over, /spline)
            r_spline[t,met,*]= vel_spline
            
            vel_spline5 = interpol(p5[met,*], vsini, vsini_over, /spline)
            r_spline5[t,met,*]= vel_spline5
            
        endfor
        
        p_spline = reform(r_spline[t,*,*])
        p_spline5 = reform(r_spline5[t,*,*])
        
        ; Get minima
        lminspline[t] = min(p_spline, lminsplineidx, /nan)
        lmin[t] = min(p, lminidx, /nan)
        ; logg 5
        lminspline5[t] = min(p_spline5, lminsplineidx5, /nan)
        lmin5[t] = min(p5, lminidx5, /nan)
        
        ; Get 2D array idx's
        lmin_idx[t,*]=array_indices(p, lminidx)
        lminspline_idx[t,*]=array_indices(p_spline, lminsplineidx)

        lmin_idx5[t,*]=array_indices(p5, lminidx5)
        lminspline_idx5[t,*]=array_indices(p_spline5, lminsplineidx5)
    endfor
    
    ; Find the overall minimum logg=4.5
    globmin = min(lmin, globidx, /nan)
    globmin_spline = min(lminspline, globsplineidx, /nan)
    
    best_teff = teff_vec[globidx]
    bestspline_teff = teff_vec[globsplineidx]
    
    best_feh_idx = lmin_idx[globidx,0]
    best_feh = feh[best_feh_idx]
    bestspline_feh_idx = lminspline_idx[globsplineidx,0]
    bestspline_feh = feh[bestspline_feh_idx]
    
    best_vel_idx = lmin_idx[globidx,1]
    best_vel = vsini[best_vel_idx]
    bestspline_vel_idx = lminspline_idx[globsplineidx,1]
    bestspline_vel = vsini_over[bestspline_vel_idx]
    
    ; Find the overall minimum logg=5
    globmin5 = min(lmin5, globidx5, /nan)
    globmin_spline5 = min(lminspline5, globsplineidx5, /nan)
    
    best_teff5 = teff_vec[globidx5]
    bestspline_teff5 = teff_vec[globsplineidx5]
    
    best_feh_idx5 = lmin_idx5[globidx5,0]
    best_feh5 = feh[best_feh_idx5]
    bestspline_feh_idx5 = lminspline_idx5[globsplineidx5,0]
    bestspline_feh5 = feh[bestspline_feh_idx5]
    
    best_vel_idx5 = lmin_idx5[globidx5,1]
    best_vel5 = vsini[best_vel_idx5]
    bestspline_vel_idx5 = lminspline_idx5[globsplineidx5,1]
    bestspline_vel5 = vsini_over[bestspline_vel_idx5]

    ; Find the Vsini and Chi2 at the best point in other logg thang
    ; logg5 at logg4.5 minimum
    chivel_545_vec = r5[globsplineidx, bestspline_feh_idx, *]
    chivel_545_vec = reform(chivel_545_vec)
    chivel_545_splinevec = interpol(chivel_545_vec, vsini, vsini_over, /spline)

    
    if fnum eq 5 then begin
        help, chivel_545_splinevec
        help, chivel_545_vec
        stop
    endif

    best_chi_spline545 = min(chivel_545_splinevec)

    min_545_splineidx = where(chivel_545_splinevec eq best_chi_spline545)
        
    ;best_vel_545 = vsini[min_545_idx]
    best_vel_spline545 = vsini_over[min_545_splineidx]
    ;best_chi_545 = min(chivel_545_vec)
    ;best_chi_spline545 = min(chivel_545_splinevec)


    ; logg4.5 at logg5 minimum
    chivel_455_vec = r[globsplineidx5, bestspline_feh_idx5, *]
    chivel_455_vec = reform(chivel_455_vec)
    chivel_455_splinevec = interpol(chivel_455_vec, vsini, vsini_over, /spline)

    best_chi_spline455 = min(chivel_455_splinevec)

    min_455_splineidx = where(chivel_455_splinevec eq best_chi_spline455)
        
    ;best_vel_545 = vsini[min_545_idx]
    best_vel_spline455 = vsini_over[min_455_splineidx]
    ;best_chi_545 = min(chivel_545_vec)
    ;best_chi_spline545 = min(chivel_545_splinevec)

    if fnum eq 5 then begin
        help, chivel_455_splinevec
        help, chivel_455_vec
        stop
    endif

    
    
 ;;; ALWAYS USE SPLINE!  POLY IS UNRELIABLE
    ;best_globmin = min([globmin, globmin_over, globmin_spline])
    ;best_globmin = min([globmin, globmin_spline])
    ;best_globmin5 = min([globmin5, globmin_spline5])

    ;Determine which method gives the best result

    best_best_teff5 = bestspline_teff5
    best_best_feh5 = bestspline_feh5
    best_best_vel5 = bestspline_vel5
    best_best_chi25 = globmin_spline5
    model_best = "SPLINE"

    best_best_teff = bestspline_teff
    best_best_feh = bestspline_feh
    best_best_vel = bestspline_vel
    best_best_chi2 = globmin_spline
    model_best = "SPLINE"

;     case 1 of
        
;         best_globmin eq globmin: begin
;             best_best_teff = best_teff
;             best_best_feh = best_feh
;             best_best_vel = best_vel
;             best_best_chi2 = globmin
;             model_best = "GRID"
;         end

;         best_globmin eq globmin_spline: begin
;             best_best_teff = bestspline_teff
;             best_best_feh = bestspline_feh
;             best_best_vel = bestspline_vel
;             best_best_chi2 = globmin_spline
;             model_best = "SPLINE"
;         end

;         else: begin
;             best_best_teff = best_teff
;             best_best_feh = best_feh
;             best_best_vel = best_vel
;             best_best_chi2 = globmin
;             model_best = "NONE"
;         end
;     endcase

;     case 1 of
        
;         best_globmin5 eq globmin5: begin
;             best_best_teff5 = best_teff5
;             best_best_feh5 = best_feh5
;             best_best_vel5 = best_vel5
;             best_best_chi25 = globmin5
;             model_best = "GRID"
;         end

;         best_globmin5 eq globmin_spline5: begin
;             best_best_teff5 = bestspline_teff5
;             best_best_feh5 = bestspline_feh5
;             best_best_vel5 = bestspline_vel5
;             best_best_chi25 = globmin_spline5
;             model_best = "SPLINE"
;         end

;         else: begin
;             best_best_teff5 = best_teff5
;             best_best_feh5 = best_feh5
;             best_best_vel5 = best_vel5
;             best_best_chi25 = globmin5
;             model_best = "NONE"
;         end
;     endcase

    ;Check to see if we get different FEH's and TEFF's
;     if best_teff ne best_best_teff or bestover_teff ne best_best_teff $
;       or bestspline_teff ne best_best_teff then teff_ne_flag = 1 $
;       else teff_ne_flag = 0
;     if best_feh ne best_best_feh or bestover_feh ne best_best_feh $
;       or bestspline_feh ne best_best_feh then feh_ne_flag = 1 $
;       else feh_ne_flag = 0

;     ;Check if APOGEE result is consistent
;     if abs(best_best_teff-teff_apg) ge 200 then teff_apg_flag = 1 $
;       else teff_apg_flag = 0
;     if abs(best_best_feh-feh_apg) ge 1 then feh_apg_flag = 1 $
;       else feh_apg_flag = 0
    
    ;Check smoothness of chi2 surface somehow
        
print, chivel_455_vec eq reform(r[globsplineidx5,bestspline_feh_idx5,*])

    outstr = {TEFF_VEC:teff_vec, $
              FEH_VEC:feh, $
              VSINI_VEC:vsini, $
              VSINI_OVER_VEC:vsini_over, $
              TEFF_APG:teff_apg, $
              TEFF_ASP:asp[fnum].teff, $
              TEFF_BEST_45:best_best_teff, $
              TEFF_BEST_5:best_best_teff5, $
              FEH_APG:feh_apg, $
              FEH_BEST_45:best_best_feh, $
              FEH_BEST_5:best_best_feh5, $
              VSINI_ASP:asp[fnum].vsini, $
              VSINI_BEST_45:best_best_vel, $
              VSINI_BEST_5:best_best_vel5, $
              VSINI_BEST_45_5:best_vel_spline455, $
              VSINI_BEST_5_45:best_vel_spline545, $
              CHI2_BEST_45:best_best_chi2, $
              CHI2_BEST_5:best_best_chi25, $
              CHI2_BEST_45_5:best_chi_spline455, $
              CHI2_BEST_5_45:best_chi_spline545, $
              CHIVEL_GRID_45:reform(r[globsplineidx,bestspline_feh_idx,*]), $
              CHIVEL_GRID_5:reform(r5[globsplineidx5,bestspline_feh_idx5,*]), $
              CHIVEL_GRID_45_5:reform(r[globsplineidx5,bestspline_feh_idx5,*]), $
              CHIVEL_GRID_5_45:reform(r5[globsplineidx,bestspline_feh_idx,*]) $
             }
    
    ;;; BEST STUFF    
    ;;; Get the valid entries
    valid_idx = where(rbest.valid_flag eq 1, nvalid)

    ; Save the structure that corresponds to the best fit
    if nvalid gt 0 then begin
        rb = rbest[valid_idx]
        cbest_idx = where(rb.chi2_best eq min(rb.chi2_best, /nan), ncbest)
        if ncbest ne 1 then message, "problem finding best chi2"
        if rb[cbest_idx].chi2_best eq 0 then message, "Error: best chi2 = 0"
        beststr = rb[cbest_idx]
    endif else message, "No valid flags"




    if fnum eq 0 then begin
        outdata = replicate(outstr, nspec)
        outbest = replicate(beststr, nspec)
    endif
    outdata[fnum] = outstr
    outbest[fnum] = beststr
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
   
    
    plot, outstr.vsini_vec, outstr.chivel_grid_45
    oplot, [beststr.vsini_best], [beststr.chi2_best], ps=6, color=200

    ;dxstop
    
endfor
;Write results
outfile = '/home/stgilhool/APOGEE/vsini_results/logg45_logg5_compare/rfile_master.fits'

mwrfits, outdata, outfile, /create
mwrfits, outbest, outfile



stop

end
