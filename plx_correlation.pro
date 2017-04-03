; Find features in APG spectra that correlate with parallax
function plx_fit, p, ew=ew, plx=plx, err=err, vis=vis

common outputs, model, res

plx_col = reform(plx, 1, n_elements(plx))

model = ew ## p

model = reform(model)

data = plx

res = data - model

dev = res/err

chi2 = sqrt(total(dev^2, /double))

sorti = sort(data)
sortd = data[sorti]
sortm = model[sorti]
sortr = res[sorti]

if vis eq 1 then begin
    ;plot, data, ps = 6, xtitle = "File number", ytitle = "Metalicity", title = strtrim(chi2,2), /xs
    ;oplot, model, ps = 6, color = 200
    ;plot, res, ps = 6, title = "Residuals", /xs
    plot, sortd, ps = 6, xtitle = "File number", ytitle = "Parallax", title = strtrim(chi2,2), /xs
    oplot, sortm, ps = 6, color = 200
    plot, sortr, ps = 6, title = "Residuals with RMSE: "+strtrim(stddev(res),2), /xs
endif

return, dev

end


pro plx_correlation

; Options
npix_window = 11L 
gd_window_threshold = 60
max_eq_widths = 30
vis = 1

; Initialize 
npix = 8575L
pix_vec = lindgen(npix)

; Read headers to get wl info
head_str = readin_apo(hdu=0, nfiles=1)
head = *(head_str.header)
wl0 = fxpar(head, 'CRVAL1')
wl1 = fxpar(head, 'CDELT1')
log_wl = poly(pix_vec, [wl0,wl1])
wl_grid = 10d0^log_wl
; Assuming the wl_value is the pixel's central wavelength,
; figure out the delta_wl for each pixel
log_wl_left_pixel_bound  = log_wl - (0.5d0*wl1)
log_wl_right_pixel_bound = log_wl + (0.5d0*wl1)
delta_wl = 10d0^(log_wl_right_pixel_bound) - $
           10d0^(log_wl_left_pixel_bound)


; Read result of plx_correlation_init
info = mrdfits('/home/stgilhool/APOGEE/parallax/plx_corr_init.fits',1)

; Take only spectra with DTM PLX's
plx_idx = where(finite(info.plx_dtm), nplx)
plx_vec = info[plx_idx].plx_dtm
e_plx_vec = info[plx_idx].eplx_dtm

; Smooth the spectra
spectra = info[plx_idx].spec
errors  = info[plx_idx].err
masks   = info[plx_idx].mask

; Get equivalent width
depth = 1d0 - spectra
eq_width_per_pix = depth*rebin(delta_wl, npix, nplx)

;smoothed_depths = smooth(depth,[npix_window,1], /nan)
; Convolve with a tophat to sum up each window
tophat = replicate(1d0, npix_window)
eq_width_nomask = convol(eq_width_per_pix, tophat)
; Now replace masked pixels with nans, convolve again, and compare
; (because I want to only select windows with no masked pixels)
flg_idx = where(masks ne 0, flgcnt)
eq_width_per_pix_masked = eq_width_per_pix
eq_width_per_pix_masked[flg_idx] = !values.d_nan
eq_width_masked = convol(eq_width_per_pix_masked, tophat, /nan)
; Compare masked to unmasked
eq_width_diff = eq_width_nomask-eq_width_masked
gd_eq_wid_idx = where(eq_width_diff eq 0, gdeqwidcnt, complement=bd_eq_wid_idx)
if gdeqwidcnt eq 0 then message, "Need some good windows"
; Save eq_widths from windows with no masked pixels
eq_width = dblarr(npix,nplx)
eq_width[gd_eq_wid_idx] = eq_width_nomask[gd_eq_wid_idx]
eq_width[bd_eq_wid_idx] = !values.d_nan

; Find number of clean spectra per window and make a cut
gd_spec_per_window = total(finite(eq_width),2)
gd_spec_window_idx = where(gd_spec_per_window ge gd_window_threshold, gdwincnt, $
                          complement = bd_spec_window_idx)

; Calculate correlations
eq_width_final = eq_width
eq_width_final[bd_spec_window_idx,*] = !values.d_nan

mean_eq_width = mean(eq_width_final, dimension=2, /nan)
diff_eq_width = eq_width_final - rebin(mean_eq_width, npix, nplx)

plx_arr = rebin(reform(plx_vec,1,nplx),npix,nplx)
mean_plx = mean(plx_arr, dimension=2)
diff_plx = plx_arr - rebin(mean_plx, npix, nplx)

corr_coeff = total(diff_eq_width*diff_plx, 2, /double, /nan)/ $
  (sqrt(total(diff_eq_width^2, 2, /double, /nan)) * $
   sqrt(total(diff_plx^2, 2, /double, /nan)))

; Find peak correlations
best_window_pos = []
ccoeff = corr_coeff

plot, pix_vec, corr_coeff, ps=6

repeat begin
    ; Find max correlation position
    max_corr = max(ccoeff, /nan)
    window_pos = where(ccoeff eq max_corr)
    best_window_pos = [best_window_pos, window_pos]
    ; Mask the region around that position
    mask_pos = window_pos[0] + lindgen(npix_window) - (npix_window/2)
    ccoeff[mask_pos] = 0d0

    oplot, mask_pos, corr_coeff[mask_pos], ps = 6, color=200
    oplot, window_pos, corr_coeff[window_pos], ps = 6, color=99999
    oplot, mask_pos, ccoeff[mask_pos], ps=6

endrep until n_elements(best_window_pos) eq max_eq_widths

; Do Multiple Linear Regression
eq_width_fit = eq_width_final
eq_width_fit = eq_width_final[best_window_pos,*]
; Add column for constant term
eq_width_fit = [reform(replicate(1d0,nplx),1,nplx),eq_width_fit]
n_eq_width_coeffs = n_elements(eq_width_fit[*,0]) 
common outputs, model, res

functargs = {ew:eq_width_fit, $
             plx:plx_vec, $
             err:e_plx_vec, $
             vis:vis $
            }
    
parinfo = replicate({value:0.1d0},n_eq_width_coeffs+1)


r = mpfit('plx_fit', parinfo=parinfo, functargs=functargs, status=status, dof=dof, bestnorm=chi2)


print, "result:"
print, r



stop

end
