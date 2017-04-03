pro teff_results_2d_interp, color

if n_elements(color) eq 0 then color = 'blue'


; LOOP
; read in sdss info
sdss_info_file = '/home/stgilhool/APOGEE/master_info5.fits'
sdss_info = mrdfits(sdss_info_file,1)

; read in a results file
rfilepath = '/home/stgilhool/APOGEE/teff_results/'

ridx_vec = [51,69,197,211,326,356,362,380,622,740,825,826,982,1060,1061,1062,1076,1216,1218,1269]
nridx = n_elements(ridx_vec) 

for ridx = 0, nridx-1 do begin
    
    ; Read in results file
    rfile_idx = ridx_vec[ridx]
    rfilename = 'teff_fit_'+color+'_'+strtrim(rfile_idx,2)+'.fits'
    rfile = rfilepath + rfilename

    rstr = mrdfits(rfile, 1)
    
    ; Get comparison results from irtf
    comp_info = sdss_info[rfile_idx]
    teff_apg  = comp_info.teff_apg
    teff_vfit = comp_info.teff_vfit
    teff_irtf = comp_info.teff_irtf
    feh_irtf  = comp_info.feh_irtf
    efeh_irtf = comp_info.efeh_irtf
    feh_vfit  = comp_info.feh_vfit
    feh_rgrss = comp_info.feh_rgrss

    ; Manage the grid elements
    nteff = n_elements(rstr.model_teff)
    nfeh = n_elements(rstr.model_feh) 
    nlogg = n_elements(rstr.model_logg) 
    
    ; right now, logg = 5.0 only
    
    ; Get the data for the plot
    teff_vec = rstr.model_teff
    
    feh_vec = rstr.model_feh
    
    chi2_vec = rstr.chi2
    
    chi2_arr = reform(chi2_vec, nfeh, nteff)
    
    ; Interpolate in 2d
    xgrid = maken(0d0, 3d0, 14)
    ygrid = maken(0d0, 14d0, 57)
    
    feh_grid = interpol(feh_vec, dindgen(4), xgrid)
    teff_grid = interpol(teff_vec, dindgen(15), ygrid)
    
    chi2_over_arr = interpolate(chi2_arr, xgrid, ygrid, /grid, $
                                missing=!values.d_nan)
    
    minchi = min(chi2_over_arr)
    minchi_idx = where(chi2_over_arr eq minchi)
    minchi_rowcol = array_indices(chi2_over_arr, minchi_idx)
    minchi_col = minchi_rowcol[0]
    minchi_row = minchi_rowcol[1]
    
    best_feh = feh_grid[minchi_col]
    best_teff = teff_grid[minchi_row]

    
    ;cgsurf, chi2_over_arr, feh_grid, teff_grid
    ;s =surface(chi2_over_arr, feh_grid, teff_grid)
    
    ;set_plot, 'ps'
    ;device, filename = 'teff_contour.ps'
    ;device, /color, bits=8
    ;device, xs = 13, ys= 8, /inches
    ;loadct, 13
    
    ;contour, chi2_over_arr, feh_grid, teff_grid, max_value = 15000,
    ;levels=[11000,12000,13000], /xs, /ys
    for i = 0, 9 do begin
        levels= lindgen(i+1)*2d-2*minchi + (minchi)
        ;xr = [min(feh_grid), max(feh_grid)]
        xr = [min(feh_grid), 0.5]
        yr = [min(teff_grid), max(teff_grid)]
        contour, chi2_over_arr, feh_grid, teff_grid, max_value = 2*minchi, levels=levels, /xs, /ys, /downhill, xr = xr, yr=yr
        
        plot, [best_feh],[best_teff], ps=4, xs=5, ys=5, /noerase, xr=xr, yr=yr
        plot, [feh_irtf],[teff_irtf], ps=6, color=200, xs=5, ys=5, /noerase, xr=xr, yr=yr
        print, levels
        dxstop
    endfor
    
    print, min(chi2_over_arr)
    print, where(chi2_over_arr eq min(chi2_over_arr))


endfor




; chi2_arr=[]
; chi2_over_arr = []
; teff_min_vec = []
; chi2_min_vec = []
; oversamp = 7L
; ;??????
; teff_over_grid = dindgen(nteff*oversamp)*(100d0/oversamp) + min(teff_mod)
; teff_over_idx = dindgen(nteff*oversamp)/oversamp
; feh_over_grid = dindgen(9)*1d-1 - 0.5
; feh_over_idx = dindgen(9)

; for feh_idx = 0, nfeh-1 do begin
;     feh_i = rstr.model_feh[feh_idx]
;     feh_idx_vec = where(feh_vec eq feh_i, nfehmatch)
;     if nfehmatch gt 0 then begin
;         chi2_vec = rstr.chi2[feh_idx_vec]
;         chi2_arr = [[chi2_arr],[chi2_vec]]
;         chi2_over = interpol(chi2_vec, teff_mod, teff_over_grid, /spline)
;         chi2_over_arr = [[chi2_over_arr], [chi2_over]]
;     endif else message, "Error in making the teff vector for feh= "+strtrim(feh_i,2)
; ;    if feh_idx eq 0 then plot, teff_mod, chi2_vec, linestyle = feh_idx+1 $
;  ;     else oplot, teff_mod, chi2_vec, linestyle = feh_idx + 1
;     if feh_idx eq 0 then plot, teff_over_grid, chi2_over, linestyle = feh_idx+1 $
;     else oplot, teff_over_grid, chi2_over, linestyle = feh_idx + 1
;     teff_min_idx = where(chi2_over eq min(chi2_over))
;     teff_min = teff_over_grid[teff_min_idx]
;     teff_min_vec = [teff_min_vec, teff_min]
;     chi2_min = chi2_over[teff_min_idx]
;     chi2_min_vec = [chi2_min_vec, chi2_min]
; endfor
; help, teff_min_vec
; print, teff_min_vec
; help, chi2_min_vec
; print, chi2_min_vec
; stop
; ; plot model and report goodness of fit indicators (chi2, status

; ; rewrite structure


end
