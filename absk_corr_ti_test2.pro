; Find features in APG spectra that correlate with parallax
function ew_gauss_calc, ti_info, spectrum, error, mask

common ew_info, delta_wl, wl_grid
common temp, specnum, plx_idx, both_idx
; vectors and indices
lcont_range = ti_info.left_cont_range
rcont_range = ti_info.right_cont_range
npix = rcont_range[1]-lcont_range[0]+1
npix_lcont = lcont_range[1]-lcont_range[0]+1
npix_rcont = rcont_range[1]-rcont_range[0]+1
npix_feat = npix - npix_lcont - npix_rcont

winpix = lindgen(npix) + lcont_range[0]
winspec = spectrum[winpix]
winmask = mask[winpix]
winwl = wl_grid[winpix]
win_dwl = delta_wl[winpix]

idx_vec = lindgen(npix)
lcont_idx = idx_vec[0:npix_lcont-1]
rcont_idx = idx_vec[-1*npix_rcont:-1]
feat_idx = idx_vec[lcont_idx[1]+1:rcont_idx[0]-1]

wl_ctr = ti_info.wl_ctr
pix_ctr = ti_info.pix_ctr

; Mask
mask_idx = where(winmask ne 0, mskcnt)
winspec_copy = winspec

if mskcnt gt 0 then begin
    winspec[mask_idx] = !values.d_nan
endif

; Reject window if there are no pixels in either continuum area
if (total(finite(winspec[lcont_idx])) eq 0) or $
  (total(finite(winspec[rcont_idx])) eq 0) then begin
    ew = !values.d_nan
    return, ew
endif 

; Reject window if more than 1/3 of the feature pixels are masked
if total(finite(winspec[feat_idx])) le (2/3.)*npix_feat then begin
    ew = !values.d_nan
    return, ew
endif 



; Fit a Gaussian
amplitude = winspec[where(winpix eq pix_ctr)] - 1d0
centroid = wl_ctr
fwhm = 3d0* mean(win_dwl)
cont0 = 1d0
cont1 = 0d0

estimates = [amplitude, centroid, fwhm, cont0, cont1]
profile = mpfitpeak(winwl, winspec, coeffs, nterms=5, estimates=estimates, /nan)

cont_fit = poly(winwl, [coeffs[3],coeffs[4]])
ew_per_pix = ((cont_fit - profile)/cont_fit) * win_dwl
ew_per_pix_compare = ((cont_fit - winspec_copy)/cont_fit) * win_dwl

ew = total(ew_per_pix, /double)
ew_compare = total(ew_per_pix_compare, /double)

plot, winpix, winspec_copy, title = "EW(nomask) = "+strtrim(ew_compare,2)+ $
  "EW(mask) = " + strtrim(ew,2) + " Diff = " + strtrim(abs(ew-ew_compare),2), $
  charsize = 1.5
oplot, winpix, profile, color=200
oplot, winpix, cont_fit, color = 99999

;in_plx_test = where(both_idx eq specnum, inplxcnt)
in_plx_test = where(plx_idx eq specnum, inplxcnt)
if inplxcnt eq 1 then begin
    flag = ''
    read, flag, prompt = "Keep this window? y/n"
    
    if flag eq 'n' then ew = !values.d_nan
endif

;if ew gt 0.539 and ew lt 0.54 then begin
;     plot, winpix, winspec_copy, title = "EW(nomask) = "+strtrim(ew_compare,2)+ $
;       "EW(mask) = " + strtrim(ew,2) + " Diff = " + strtrim(abs(ew-ew_compare),2), $
;       charsize = 1.5
;     oplot, winpix, profile, color=200
    
;     stop
;    ew = !values.d_nan
;    return, ew
;endif


return, ew

end

function ew_cont_fit, lcont, rcont, spectrum, error, mask
;;;Fit a low-order polynomial to the continuum regions 
npix_window = rcont[1]-lcont[0]+1

winpix = lindgen(npix_window) + lcont[0]

winy = replicate(!values.d_nan, npix_window)

; Populate continuum areas
npix_lcont = lcont[1]-lcont[0]+1
npix_rcont = rcont[1]-rcont[0]+1
lcont_idx = lindgen(npix_lcont) + lcont[0]
rcont_idx = lindgen(npix_rcont) + rcont[0]

new_npix = npix_window - npix_lcont - npix_rcont + 2
new_winpix = lindgen(new_npix) + lcont[1]

newx0 = lcont[1]
newx1 = rcont[0]
newy0 = mean(spectrum[lcont_idx], /nan)
newy1 = mean(spectrum[rcont_idx], /nan)

slope = (newy1-newy0)/(newx1-newx0)
b = newy1 - slope*newx1

new_winy = slope * new_winpix + b

winy[0:npix_lcont-1] = spectrum[lcont_idx]
winy[-1*npix_rcont:-1] = spectrum[rcont_idx]

;Fit a new continuum
;cont_coeffs = robust_poly_fit(winpix, winy, 1, cont_fit, /double, numit=0)
cont_coeffs = poly_fit(new_winpix, new_winy, 1, /double, yfit = cont_fit)

return, cont_coeffs

end

function ew_calc, continuum_coefficients, feature_range, window_range, spectrum, $
                  error, mask

common ew_info, delta_wl, wl_grid

npix_window = window_range[1]-window_range[0]+1
winpix = lindgen(npix_window) + window_range[0]
winspec = spectrum[winpix]
win_dwl = delta_wl[winpix]

;featpix = winpix[feature_range[0]:feature_range[1]]

; Calculate EW
cont_fit = poly(winpix, continuum_coefficients)
ew_per_pix = ((cont_fit-winspec)/cont_fit)*win_dwl
ew = total(ew_per_pix, /double)

return, ew

end

function ti_lines, spectra, errors, masks
;;;Make an array of EW's for ti lines that I've defined

common ew_info, delta_wl, wl_grid
common temp, specnum, plx_idx, both_idx
;TI EW info
n_ti_lines = 6
nspec = n_elements(spectra[0,*]) 
ti_str = {ti_num:0L, $
          wl_ctr:0d0, $
          pix_ctr:0L, $
          feature_range:[0L,0L], $
          left_cont_range:[0L,0L], $
          right_cont_range:[0L,0L], $
          ew:dblarr(nspec) $
         }

ti_info = replicate(ti_str, n_ti_lines)

;Populate TI_INFO
ti_num = [3,7,10,15,22,28]
wl_ctr = [15339.038d0, $
          15548.009d0, $
          15719.872d0, $
          16211.129d0, $
          16344.993d0, $
          16405.992d0]
pix_ctr = [1133L,2112L,2908L,5136L,5686L,6000L]

left_cont_range = [[1125L,1127L], $
                   [2107L,2108L], $
                   [2897L,2899L], $
                   [5129L,5131L], $
                   [5681L,5683L], $
                   [5995L,5997L]]

right_cont_range = [[1138L,1139L], $
                    [2120L,2123L], $
                    [2920L,2922L], $
                    [5146L,5148L], $
                    [5691L,5693L], $
                    [6004L,6006L]]

; Save vectors
ti_info.ti_num = ti_num
ti_info.wl_ctr = wl_ctr
ti_info.pix_ctr = pix_ctr
      
;REMOVE
n_ti_lines = 1
start_line = 0
; Save arrays
for linenum = start_line, start_line+n_ti_lines-1 do begin
    
    lcont = left_cont_range[*,linenum]
    ti_info[linenum].left_cont_range = lcont
    
    rcont = right_cont_range[*,linenum]
    ti_info[linenum].right_cont_range = rcont    

    feature_left = lcont[1]+1
    feature_right = rcont[0]-1
    feature_range = [feature_left,feature_right]
    window_left = lcont[0]
    window_right = rcont[1]
    window_range = [window_left,window_right]
    ti_info[linenum].feature_range = feature_range
    
    ; Fit continuum
    for specnum = 0, nspec-1 do begin
        ; continuum_coeffs = ew_cont_fit(lcont, rcont, $
;                                        spectra[*,specnum], errors[*,specnum], $
;                                        masks[*,specnum])
;         ; Calculate EW
;         ew_i = ew_calc(continuum_coeffs, feature_range, window_range, $
;                        spectra[*,specnum], errors[*,specnum], $
;                        masks[*,specnum])
;         ti_info[linenum].ew[specnum] = ew_i

        ew_i = ew_gauss_calc(ti_info[linenum], spectra[*,specnum], errors[*,specnum], $
                             masks[*,specnum])
        ti_info[linenum].ew[specnum] = ew_i
    endfor
endfor

return, ti_info


end

    
    


function absk_fit, p, ew=ew, absk=absk, err=err, output=output, vis=vis

common fit_outputs, model, res

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

if output eq 'chi2' then return, dev

if output eq 'model' then return, model

end


pro absk_corr_ti_test2

common ew_info, delta_wl, wl_grid
common temp, specnum, plx_idx, both_idx
; Options
npix_window = 11L 
gd_window_threshold = 60
max_eq_widths = 75
vis = 1

; Initialize 
npix = 8575L
pix_vec = dindgen(npix)

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


; Read in spectra
spectra = info.spec
errors  = info.err
masks   = info.mask
nspec = n_elements(info)


; na1_wl = 16378.343d0
; na1_diff = abs(wl_grid-na1_wl)
; na1_pix = where(na1_diff eq min(na1_diff))

; na2_wl = 16393.327d0
; na2_diff = abs(wl_grid-na2_wl)
; na2_pix = where(na2_diff eq min(na2_diff))

; ; Plot
; winsz = 11L
; win1_vec = spectra[na1_pix-(winsz/2):na1_pix+(winsz/2),*]
; win2_vec = spectra[na2_pix-(winsz/2):na2_pix+(winsz/2),*]

; wl1_vec = wl_grid[na1_pix-(winsz/2):na1_pix+(winsz/2)]
; wl2_vec = wl_grid[na2_pix-(winsz/2):na2_pix+(winsz/2)]

; !p.multi=[0,1,2]
; window, 0, xs=1600, ys=1000
; plot, wl1_vec, win1_vec[*,0], yr = [0.5,1.2]
; for i = 1, nspec-1 do begin
;     oplot, wl1_vec, win1_vec[*,i];, ps=3
; endfor
; oplot, wl1_vec, median(win1_vec, dimension=2), color=200

; plot, wl2_vec, win2_vec[*,0],yr = [0.5,1.2]
; for i = 1, nspec-1 do begin
;     oplot, wl2_vec, win2_vec[*,i];, ps=3
; endfor
; oplot, wl2_vec, median(win2_vec, dimension=2), color=200

;stop
; Read in DTM plx's
plx_idx = where(finite(info.plx_dtm), nplx)
plx_vec = info[plx_idx].plx_dtm/1d3
e_plx_vec = info[plx_idx].eplx_dtm/1d3
k_vec = info[plx_idx].kmag
; Get their absolute k magnitudes
absk_vec = k_vec + 5d0 + 5d0*alog10(plx_vec)
eabsk_vec = 5d0*0.434*(e_plx_vec/plx_vec)

; Get Fe/H as well
feh_idx = where(finite(info.feh_irtf), nfeh)
feh_vec = info[feh_idx].feh_irtf
efec_vec = info[feh_idx].efeh_irtf

both_idx = set_intersection(plx_idx, feh_idx)
fehboth_vec = info[both_idx].feh_irtf
plxboth_vec = info[both_idx].plx_dtm/1d3
e_plxboth_vec = info[both_idx].eplx_dtm/1d3
kboth_vec = info[both_idx].kmag
; Get their absolute k magnitudes
abskboth_vec = kboth_vec + 5d0 + 5d0*alog10(plxboth_vec)
eabskboth_vec = 5d0*0.434*(e_plxboth_vec/plxboth_vec)

; Take only spectra with DTM plx's
; eq_width_plx = eq_width[plx_idx,*]

; Find number of clean spectra per window and make a cut
; gd_spec_per_window = total(finite(eq_width_plx),2)
; ; Only keep regions where there are more than gd_window_threshold EW's
; gd_spec_window_idx = where(gd_spec_per_window ge gd_window_threshold, gdwincnt, $
;                           complement = bd_spec_window_idx)


; eq_width_final[bd_spec_window_idx,*] = !values.d_nan

; mean_eq_width = mean(eq_width_final, dimension=2, /nan)
; diff_eq_width = eq_width_final - rebin(mean_eq_width, npix, nplx)

; absk_arr = rebin(reform(absk_vec,1,nplx),npix,nplx)
; mean_absk = mean(absk_arr, dimension=2)
; diff_absk = absk_arr - rebin(mean_absk, npix, nplx)

; corr_coeff = total(diff_eq_width*diff_absk, 2, /double, /nan)/ $
;   (sqrt(total(diff_eq_width^2, 2, /double, /nan)) * $
;    sqrt(total(diff_absk^2, 2, /double, /nan)))

; ; Find peak correlations
; best_window_pos = []
; ccoeff = corr_coeff

; plot, pix_vec, corr_coeff, ps=6

; repeat begin
;     ; Find max correlation position
;     max_corr = max(ccoeff, /nan)
;     window_pos = where(ccoeff eq max_corr)
;     best_window_pos = [best_window_pos, window_pos]
;     ; Mask the region around that position
;     mask_pos = window_pos[0] + lindgen(npix_window) - (npix_window/2)
;     ccoeff[mask_pos] = 0d0

;     oplot, mask_pos, corr_coeff[mask_pos], ps = 6, color=200
;     oplot, window_pos, corr_coeff[window_pos], ps = 6, color=99999
;     oplot, mask_pos, ccoeff[mask_pos], ps=6

; endrep until n_elements(best_window_pos) eq max_eq_widths

; ; Do Multiple Linear Regression
; eq_width_fit = eq_width_final
; eq_width_fit = eq_width_final[best_window_pos,*]
; ; Add column for constant term
; eq_width_fit = [reform(replicate(1d0,nplx),1,nplx),eq_width_fit]
; ; Get rid of NaNs
; nan_idx = where(finite(eq_width_fit) eq 0, nancnt)
; if nancnt gt 0 then eq_width_fit[nan_idx] = 0d0
; n_eq_width_coeffs = n_elements(eq_width_fit[*,0]) 

ti_info = ti_lines(spectra, errors, masks)


eq_width_fit = ti_info.ew[plx_idx]
eq_width_fit = transpose(eq_width_fit)

;eq_width_both = ti_info.ew[both_idx]
;eq_width_both = transpose(eq_width_both)


;nti = n_elements(ti_info) 
nti = 1
 eq_width_fit = ti_info[0].ew[plx_idx]
 fit_idx = where(finite(eq_width_fit), nfit)
;eq_width_both = ti_info[0].ew[both_idx]
;fit_idx = where(finite(eq_width_both), nfit)

if nfit eq 0 then message, "no EW's were fit"
; eq_width = eq_width_both[fit_idx]
; eq_width = reform(eq_width, 1, nfit)
eq_width = eq_width_fit[fit_idx]
eq_width = reform(eq_width, 1, nfit)

 absk_vec = absk_vec[fit_idx]
 eabsk_vec = eabsk_vec[fit_idx]
; absk_vec = abskboth_vec[fit_idx]
; eabsk_vec = eabskboth_vec[fit_idx]
feh_vec = fehboth_vec[fit_idx]
feh_vec = reform(feh_vec, 1, nfit)
; add col for constant term
; and add feh

;eq_width = [reform(replicate(1d0, nfit),1,nfit), eq_width, feh_vec]
eq_width = [reform(replicate(1d0, nfit),1,nfit), eq_width]

n_eq_width_coeffs = n_elements(eq_width[*,0]) 

common fit_outputs, model, res

functargs = {ew:eq_width, $
             absk:absk_vec, $
             err:eabsk_vec, $
             output:'chi2', $
             vis:vis $
            }
    
parinfo = replicate({value:0.1d0},n_eq_width_coeffs)
parinfo[0].value=5d0

!p.multi = [0,1,2]


r = mpfit('absk_fit', parinfo=parinfo, functargs=functargs, $
          status=status, dof=dof, bestnorm=chi2)

print, "result:"
print, r

stop

absk_model = absk_fit(r, ew=functargs.ew, absk=functargs.absk, err=functargs.err, $
                      output='model', /vis)

!p.multi=0
window, 0, xs = 1500, ys=800
plot, reform(eq_width[1,*]), absk_model, ps=6, $
  title="EW vs. M_K(fit) | Red pts are data", $
  xtit='M_K(fit)', ytit='EW', charsize = 1.5
oplot, reform(eq_width[1,*]), absk_vec, ps=4, color=200

eqs = sort(eq_width[1,*])
eq_sort = eq_width[1,eqs]
absk_model_sort = absk_model[eqs]
absk_vec_sort = absk_vec[eqs]
res = absk_vec_sort-absk_model_sort
scatter = stddev(res)
print, "Scatter is: "+strtrim(scatter,2)


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
