pro teff_results_interp_all, color

if n_elements(color) eq 0 then color = 'blue'


; LOOP
; read in sdss info
sdss_info_file = '/home/stgilhool/APOGEE/master_info5.fits'
sdss_info = mrdfits(sdss_info_file,1)

; read in empirical results file
mannfile = '/home/stgilhool/APOGEE/teff_results/teff_mann_empirical.fits'
mannstr = mrdfits(mannfile,1)
rel1_idx = mannstr.rel1_idx
rel2_idx = mannstr.rel2_idx
teff1_mann_vec = mannstr.teff_rel1
teff2_mann_vec = mannstr.teff_rel2


; read in a results file
rfilepath = '/home/stgilhool/APOGEE/teff_results/'

ridx_vec = [51,69,197,211,326,356,362,380,622,740,825,826,982,1060,1061,1062,1076,1216,1218,1269]
nridx = n_elements(ridx_vec) 

;com1 = set_intersection(ridx_vec, rel1_idx)
;com2 = set_intersection(ridx_vec, rel2_idx)

;print, com1
;print, com2
;stop

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

    ; Get comparison from empirical results from mann
    mann1_idx = where(rel1_idx eq rfile_idx, nmann1)
    mann2_idx = where(rel2_idx eq rfile_idx, nmann2)
    if nmann1 gt 0 then teff1_mann = teff1_mann_vec[mann1_idx] $
    else teff1_mann = !values.d_nan
    if nmann2 gt 0 then teff2_mann = teff2_mann_vec[mann2_idx] $
    else teff2_mann = !values.d_nan

    ; Manage the grid elements
    nteff = n_elements(rstr.model_teff)
    nfeh = n_elements(rstr.model_feh) 
    nlogg = n_elements(rstr.model_logg) 
    
    ; right now, logg = 5.0 only
    
    ; Get the data for the plot
    teff_vec = rstr.model_teff
    
    feh_vec = rstr.model_feh
    
    ; Interpolate in 2d
    xgrid = maken(0d0, 3d0, 14)
    ygrid = maken(0d0, 14d0, 57)
    
    feh_grid = interpol(feh_vec, dindgen(4), xgrid)
    teff_grid = interpol(teff_vec, dindgen(15), ygrid)
    
    ; Now loop over each color
    chi2_mtx = []
    chi2_over_mtx = []
    best_feh_arr = []
    best_teff_arr = []
    minchi_arr = []
    minchi_arr_template = []
    best_feh_arr_template = []
    best_teff_arr_template = []
    minchi_idx_arr_template = []
    dof_arr = []

    for color_idx = 0,2 do begin
        
        case color_idx of
            ;blue
            0: begin
                color='blue'
                rfilename = 'teff_fit_'+color+'_'+strtrim(rfile_idx,2)+'.fits'
                
            end
            ;green
            1: begin
                color='green'
                rfilename = 'teff_fit_'+color+'_'+strtrim(rfile_idx,2)+'.fits'
            end
            ;red
            2: begin
                color='red'
                rfilename = 'teff_fit_'+color+'_'+strtrim(rfile_idx,2)+'.fits'
            end

        endcase
        
        ;reread file
        rfile = rfilepath + rfilename
        rstr = mrdfits(rfile, 1)
        
        ;get results
        chi2_vec = rstr.chi2
        
        chi2_arr = reform(chi2_vec, nfeh, nteff)
        
        chi2_mtx = [[[chi2_mtx]],[[chi2_arr]]]

       
        
        ; Find which template corresponds to the minima with no
        ; interpolation
        minchi_template = min(chi2_arr, /nan)
        minchi_idx_template = where(chi2_arr eq minchi_template)
        minchi_rowcol_template = array_indices(chi2_arr, $
                                               minchi_idx_template)
        minchi_col_template = minchi_rowcol_template[0]
        minchi_row_template = minchi_rowcol_template[1]
        
        best_feh_template = feh_vec[minchi_col_template]
        best_teff_template = teff_vec[minchi_row_template]
        
        ;determine number of dof
        
        err = rstr.sdss_err
        npix_chip = n_elements(err) 
        
        mask_idx = where(err ge 0.5, nmask)
        
        nparam = n_elements(rstr.result_arr[*,0]) 
        dof = npix_chip - nmask - nparam
        dof_arr = [dof_arr, dof]
        ;chi2_dof = minchi_arr_template[c]/dof

        ;save template results
        minchi_arr_template = [minchi_arr_template, minchi_template]
        minchi_idx_arr_template = [minchi_idx_arr_template, minchi_idx_template]
        best_feh_arr_template = [best_feh_arr_template, best_feh_template]
        best_teff_arr_template = [best_teff_arr_template, best_teff_template]


        ; Interpolate to make a more finely gridded chi2 surface
        chi2_over_arr = interpolate(chi2_arr, xgrid, ygrid, /grid, $
                                missing=!values.d_nan, /cubic)
        
        chi2_over_mtx = [[[chi2_over_mtx]],[[chi2_over_arr]]]

        minchi = min(chi2_over_arr, /nan)
        
        minchi_arr = [minchi_arr, minchi]

        minchi_idx = where(chi2_over_arr eq minchi)
        minchi_rowcol = array_indices(chi2_over_arr, minchi_idx)
        minchi_col = minchi_rowcol[0]
        minchi_row = minchi_rowcol[1]
    
        best_feh = feh_grid[minchi_col]
        best_teff = teff_grid[minchi_row]

        best_feh_arr = [best_feh_arr,best_feh]
        best_teff_arr = [best_teff_arr,best_teff]
        
    
        ;cgsurf, chi2_over_arr, feh_grid, teff_grid
        ;s =surface(chi2_over_arr, feh_grid, teff_grid)
        
        ;set_plot, 'ps'
        ;device, filename = 'teff_contour.ps'
        ;device, /color, bits=8
        ;device, xs = 13, ys= 8, /inches
        ;loadct, 13
        
        ;contour, chi2_over_arr, feh_grid, teff_grid, max_value = 15000,
        ;levels=[11000,12000,13000], /xs, /ys
    endfor
    
    
    

    for c = 0,2 do begin
        case c of
            ;blue
            0: begin
                color='blue'
                rfilename = 'teff_fit_'+color+'_'+strtrim(rfile_idx,2)+'.fits'
                
            end
            ;green
            1: begin
                color='green'
                rfilename = 'teff_fit_'+color+'_'+strtrim(rfile_idx,2)+'.fits'
            end
            ;red
            2: begin
                color='red'
                rfilename = 'teff_fit_'+color+'_'+strtrim(rfile_idx,2)+'.fits'
            end

        endcase
        
        ;reread file
        rfile = rfilepath + rfilename
        rstr = mrdfits(rfile, 1)
          
        wl = rstr.sdss_wl
        data = rstr.sdss_spec
        npix_chip = n_elements(data) 
        err = rstr.sdss_err
        model = rstr.model_arr[*,minchi_idx_arr_template[c]]
        
        ; make mask and apply to best-fit model
        mask_idx = where(err ge 0.5, nmask)
        mask = replicate(0, npix_chip)
        if nmask gt 0 then begin
            mask[mask_idx] = 1
            model[mask_idx] = !values.d_nan
        endif
        nparam = n_elements(rstr.result_arr[*,0]) 
        ;dof = npix_chip - nmask - nparam
        dof = dof_arr[c]
        chi2_dof = minchi_arr_template[c]/dof
    
        ; data to display along with fit
        
        bfit_feh_st = strtrim(best_feh_arr_template[c], 2)
        bfit_teff_st = strtrim(best_teff_arr_template[c], 2)
        chi2_dof_st = strtrim(chi2_dof, 2)
        
        
        window, 0, xs = 1700, ys=400
        plot, wl, data, /xs, yr=[0.5,1.1],/ys, title = 'No: '+strtrim(rfile_idx,2)+' | Chip: '+color+' | Best Teff: '+bfit_teff_st+' | Best Fe/H: '+ $
          bfit_feh_st+' | Chi2/DOF: '+chi2_dof_st
        oplot, wl, model, color=200
        
        if c eq 0 then begin
            window, 2
            plot, wl, data, xr = [15363, 15383], yr=[0.5,1.1], /ys
            oplot, wl, model, color=200
            oplot, [15364.4,15364.4], [0,2], linestyle=2
            oplot, [15381.0,15381.0], [0,2], linestyle=2
        endif
        ;;;;; Show contours
        window, 1, xs = 1300, ys=600
        for i = 0, 5 do begin
            ;levels= lindgen(i+1)*2d-2*minchi + (minchi)
            minchi = minchi_arr[c]/dof
            ;levels= lindgen(i+1)*2d-2*minchi + (minchi)
            levels= lindgen(i+1) + (minchi)
            ;xr = [min(feh_grid), max(feh_grid)]
            xr = [min(feh_grid)-0.5, 0.5]
            yr = [min(teff_grid)-100, max(teff_grid)+100]
            chi2_display = chi2_over_mtx[*,*,c]/dof
            ;contour, chi2_over_arr, feh_grid, teff_grid, max_value =
            ;2*minchi, levels=levels, /xs, /ys, /downhill, xr = xr,
            ;yr=yr
            mc0=strtrim(minchi_arr[0]/dof_arr[0],2)
            mc1=strtrim(minchi_arr[1]/dof_arr[1],2)
            mc2=strtrim(minchi_arr[2]/dof_arr[2],2)
            mcstring = strtrim(minchi_arr/dof_arr,2)
            chi_string = strjoin(mcstring,', ')

            bt0=strtrim(best_teff_arr[0],2)
            bt1=strtrim(best_teff_arr[1],2)
            bt2=strtrim(best_teff_arr[2],2)
            btstring = strtrim(best_teff_arr,2)
            teff_string = strjoin(btstring,', ')

            bf0=strtrim(best_feh_arr[0],2)
            bf1=strtrim(best_feh_arr[1],2)
            bf2=strtrim(best_feh_arr[2],2)
            bfstring = strtrim(best_feh_arr,2)
            feh_string = strjoin(bfstring,', ')

            contour, chi2_display, feh_grid, teff_grid, max_value = 2*minchi, levels=levels, /xs, /ys, /downhill, xr = xr, yr=yr, title = color+' | Chi2 = '+chi_string+' | Teff: '+teff_string+' | FeH: '+feh_string+' | Teff_IRTF: '+strtrim(teff_irtf,2)+' | FeH_IRTF: '+strtrim(feh_irtf,2)
            
            plot, [best_feh_arr[0]],[best_teff_arr[0]], ps=4, xs=5, ys=5, /noerase, xr=xr, yr=yr
            plot, [best_feh_arr[1]],[best_teff_arr[1]], ps=5, xs=5, ys=5, /noerase, xr=xr, yr=yr
            plot, [best_feh_arr[2]],[best_teff_arr[2]], ps=1, xs=5, ys=5, /noerase, xr=xr, yr=yr
            plot, [feh_irtf],[teff_irtf], ps=6, color=200, xs=5, ys=5, /noerase, xr=xr, yr=yr
            plot, [feh_irtf],[teff1_mann], ps=4, color=200, xs=5, ys=5, /noerase, xr=xr, yr=yr
            plot, [feh_irtf],[teff2_mann], ps=5, color=200, xs=5, ys=5, /noerase, xr=xr, yr=yr
            print, levels
            dxstop
        endfor
    endfor
    
    print, min(chi2_over_arr, /nan)

    print, where(chi2_over_arr eq min(chi2_over_arr, /nan))


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
