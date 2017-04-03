; Find features in APG spectra that correlate with parallax
function absk_fit, p, ew=ew, absk=absk, err=err, vis=vis

common outputs, model, res

absk_col = reform(absk, 1, n_elements(absk))

model = ew ## p

model = reform(model)

data = absk

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
    plot, sortd, ps = 6, xtitle = "File number", ytitle = "Absolute K Magnitude", title = strtrim(chi2,2), /xs
    oplot, sortm, ps = 6, color = 200
    plot, sortr, ps = 6, title = "Residuals with RMSE: "+strtrim(stddev(res),2), /xs
    wait, 0.001
endif

return, dev

end


pro absk_correlation_test2

; Options
npix_window = 11L 
gd_window_threshold = 60
max_eq_widths = 1
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
info = mrdfits('/home/stgilhool/APOGEE/absk_corr/absk_corr_init.fits',1)

; Take only spectra with DTM PLX's
;plx_idx = where(finite(info.plx_dtm), nplx)
plx_idx = where(finite(info.feh_irtf) + finite(info.plx_dtm) $
                + finite(info.vmag_apass) eq 3, nplx)
plx_vec = info[plx_idx].plx_dtm/1d3
feh_vec = info[plx_idx].feh_irtf
e_plx_vec = info[plx_idx].eplx_dtm/1d3
e_feh_vec = info[plx_idx].efeh_irtf
k_vec = info[plx_idx].kmag
v_vec = info[plx_idx].vmag_apass
; Get their absolute k magnitudes
absk_vec = k_vec + 5d0 + 5d0*alog10(plx_vec)

eabsk_vec = 5d0*0.434*(e_plx_vec/plx_vec)

; V-K
vk_vec = v_vec - k_vec

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
;eq_width_per_pix_masked[flg_idx] = !values.d_nan
eq_width_masked = convol(eq_width_per_pix_masked, tophat, /nan)
; Compare masked to unmasked
eq_width_diff = eq_width_nomask-eq_width_masked
gd_eq_wid_idx = where(eq_width_diff eq 0, gdeqwidcnt, complement=bd_eq_wid_idx)
if gdeqwidcnt eq 0 then message, "Need some good windows"
; Save eq_widths from windows with no masked pixels
eq_width = dblarr(npix,nplx)
eq_width[gd_eq_wid_idx] = eq_width_nomask[gd_eq_wid_idx]
;eq_width[bd_eq_wid_idx] = !values.d_nan

; Find number of clean spectra per window and make a cut
gd_spec_per_window = total(finite(eq_width),2)
gd_spec_window_idx = where(gd_spec_per_window ge gd_window_threshold, gdwincnt, $
                          complement = bd_spec_window_idx)

; Calculate correlations
eq_width_final = eq_width
;eq_width_final[bd_spec_window_idx,*] = !values.d_nan

mean_eq_width = mean(eq_width_final, dimension=2, /nan)
diff_eq_width = eq_width_final - rebin(mean_eq_width, npix, nplx)

absk_arr = rebin(reform(absk_vec,1,nplx),npix,nplx)
mean_absk = mean(absk_arr, dimension=2)
diff_absk = absk_arr - rebin(mean_absk, npix, nplx)

corr_coeff = total(diff_eq_width*diff_absk, 2, /double, /nan)/ $
  (sqrt(total(diff_eq_width^2, 2, /double, /nan)) * $
   sqrt(total(diff_absk^2, 2, /double, /nan)))

; Find peak correlations
best_window_pos = []
ccoeff = corr_coeff

;plot, pix_vec, corr_coeff, ps=6

repeat begin
    ; Find max correlation position
    max_corr = max(ccoeff, /nan)
    window_pos = where(ccoeff eq max_corr)
    window_pos = [6859L]
    max_corr = ccoeff[window_pos]
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

; Add column for feh, v-k
vk_col = reform(vk_vec,1,nplx)
feh_col = reform(feh_vec, 1,nplx)
eq_width_fit = [eq_width_fit, feh_col, vk_col]

; Get rid of NaNs
nan_idx = where(finite(eq_width_fit) eq 0, nancnt)
if nancnt gt 0 then begin
    eq_width_fit[nan_idx] = 0d0 
    stop

endif else print, "NICE, no NANs"
n_eq_width_coeffs = n_elements(eq_width_fit[*,0]) 
common outputs, model, res

functargs = {ew:eq_width_fit, $
             absk:absk_vec, $
             err:eabsk_vec, $
             vis:vis $
            }
    
parinfo = replicate({value:0.1d0},n_eq_width_coeffs)
parinfo[0].value=20d0

!p.multi = [0,1,2]


r = mpfit('absk_fit', parinfo=parinfo, functargs=functargs, status=status, dof=dof, bestnorm=chi2)


print, "result:"
print, r

stop
;stop
;;; Now, take that result, and make absolute k guesses for all spec
; Smooth the spectra
nspec   = n_elements(info) 
spectra = info.spec
errors  = info.err
masks   = info.mask

; Get equivalent width
depth = 1d0 - spectra
eq_width_per_pix = depth*rebin(delta_wl, npix, nspec)

;smoothed_depths = smooth(depth,[npix_window,1], /nan)
; Convolve with a tophat to sum up each window

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
eq_width = dblarr(npix,nspec)
eq_width[gd_eq_wid_idx] = eq_width_nomask[gd_eq_wid_idx]
eq_width[bd_eq_wid_idx] = !values.d_nan

; Make array of just eq_widths for our particular windows
eq_width_fit = eq_width[best_window_pos,*]
; Add column for constant term
eq_width_fit = [reform(replicate(1d0,nspec),1,nspec),eq_width_fit]
; Get rid of NaNs
nan_idx = where(finite(eq_width_fit) eq 0, nancnt)
if nancnt gt 0 then eq_width_fit[nan_idx] = 0d0

; Now make the absk's
absk_model = eq_width_fit ## r

absk_model = reform(absk_model)

;;;PLOT result
!p.multi=0
; V-K
vmag_apass = info.vmag_apass
kmag = info.kmag
abs_mag = absk_model

v_k = vmag_apass - abs_mag

; Get rid of nans
fin_idx = where(finite(vmag_apass) eq 1, fincnt)
ra = info[fin_idx].ra
dec = info[fin_idx].dec
abs_mag = abs_mag[fin_idx]
v_k = v_k[fin_idx]
eabs_mag = 0.25d0 * v_k ;assume 25% error
;eabs_mag = eabs_mag[fin_idx]
;plx_out = plx[fin_idx]

; Plot
window, 0, xs=1500, ys=900
;plot, v_k, abs_mag, ps=6, yr = [max(abs_mag)+1, min(abs_mag)-1], xr =
;[3,9], $
;plot, v_k, abs_mag, ps=6, yr = [max(abs_mag),min(abs_mag)], $; yr =
;[12,3], xr = [3,9], 
plot, v_k, abs_mag, ps=6, yr = [16,0], xr = [0,15], $
  xtit = 'V-K', ytit = 'M_k', /xs, /ys
;errplot, v_k, abs_mag+eabs_mag, abs_mag-eabs_mag

; Compare with Cullen
cinfo = mrdfits('/home/stgilhool/APOGEE/cullenparallax.fits',1)
cra = cinfo.ra[fin_idx]
cde = cinfo.dec[fin_idx]
cpa = cinfo.parallax[fin_idx]
cvmag = cinfo.vmag[fin_idx]
cmk = cinfo.mk[fin_idx]
svmag = vmag_apass[fin_idx]

stop

end
