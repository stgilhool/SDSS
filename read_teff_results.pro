pro read_teff_results

; LOOP

; read in a results file
rfilepath = '/home/stgilhool/APOGEE/teff_results/'
rfilename = 'teff_fit_69.fits'
rfile = rfilepath + rfilename

rstr = mrdfits(rfile, 1)

; make indices that depend on teff, feh, and logg
nteff = n_elements(rstr.model_teff)
nfeh = n_elements(rstr.model_feh) 
nlogg = n_elements(rstr.model_logg) 

; right now, logg = 5.0 only

; use the indices to plot chi2 vs teff
teff_mod = rstr.model_teff
teff_col = transpose(rstr.model_teff)
teff_col_arr = rebin(teff_col, nfeh, nteff)
teff_vec = reform(teff_col_arr, nfeh*nteff)

feh_arr = rebin(rstr.model_feh, nfeh, nteff)
feh_vec = reform(feh_arr, nfeh*nteff)
 
; find minimum for each model and global minimum
;!p.multi = [0,1,nfeh]
!p.multi = 0
chi2_arr=[]
chi2_over_arr = []
teff_min_vec = []
chi2_min_vec = []
oversamp = 7L
teff_over_grid = dindgen(nteff*oversamp)*(100d0/oversamp) + min(teff_mod)
for feh_idx = 0, nfeh-1 do begin
    feh_i = rstr.model_feh[feh_idx]
    feh_idx_vec = where(feh_vec eq feh_i, nfehmatch)
    if nfehmatch gt 0 then begin
        chi2_vec = rstr.chi2[feh_idx_vec]
        chi2_arr = [[chi2_arr],[chi2_vec]]
        chi2_over = interpol(chi2_vec, teff_mod, teff_over_grid, /spline)
        chi2_over_arr = [[chi2_over_arr], [chi2_over]]
    endif else message, "Error in making the teff vector for feh= "+strtrim(feh_i,2)
;    if feh_idx eq 0 then plot, teff_mod, chi2_vec, linestyle = feh_idx+1 $
 ;     else oplot, teff_mod, chi2_vec, linestyle = feh_idx + 1
    if feh_idx eq 0 then plot, teff_over_grid, chi2_over, linestyle = feh_idx+1 $
    else oplot, teff_over_grid, chi2_over, linestyle = feh_idx + 1
    teff_min_idx = where(chi2_over eq min(chi2_over))
    teff_min = teff_over_grid[teff_min_idx]
    teff_min_vec = [teff_min_vec, teff_min]
    chi2_min = chi2_over[teff_min_idx]
    chi2_min_vec = [chi2_min_vec, chi2_min]
endfor
help, teff_min_vec
print, teff_min_vec
help, chi2_min_vec
print, chi2_min_vec
stop
; plot model and report goodness of fit indicators (chi2, status

; rewrite structure


end
