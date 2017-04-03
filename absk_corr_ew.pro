; Determine the correlation between features
; absolute K-band magnitude, binned by metallicity


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


function ew_calc_rough, window_range, spectrum, $
                        error, mask



npix_window = window_range[1]-window_range[0]+1
winpix = lindgen(npix_window) + window_range[0]
winspec = spectrum[winpix]

; Calculate EW

; mask pixels roughly
winspec_masked = winspec
winmask = mask[winpix]
mask_idx = where(winmask ne 0, nmask)
nmask=0 ;FIX
if nmask gt 0 then begin
    winspec_masked_copy = winspec_masked
    winspec_masked_copy[mask_idx] = !values.d_nan
    ;interpolate
    
    winspec_masked = interpol(winspec_masked_copy, winpix, winpix, /spline)
    plot, winspec, /xs, yr=[0.5,1.1], /ys
    oplot, winspec_masked_copy, color=200
    oplot, winspec_masked, color=99999
    dxstop
endif


ew = total(1d0-winspec_masked, /double)

return, ew

end

pro determine_which_spectra_to_use, testidx
; Determine which spectra to use as a calibration set
infofile = '/home/stgilhool/APOGEE/master_info5.fits'

; Readin the information structure that holds
; various info on each SDSS spectrum
info = mrdfits(infofile, 1)

; Criteria are:
;  -Has Dittman parallax (absK)
;  -Is bright
;  -Preferrably has irtf Fe/H
; Note: there are 25 stars with both DTM and IRTF, and HMAG le 9.5
testidx = where(finite(info.feh_irtf) and $
                finite(info.plx_dtm) and info.hmag le 9.5, ntest)



end



pro absk_corr_ew, feature0


; HARDCODE
;feature0 = 300
npix_window = 11
feature1 = feature0+npix_window

feature_range = [feature0, feature1]

feature_idx = lindgen(feature1-feature0)+feature0

testidx = []
; get indices of spectra to use
determine_which_spectra_to_use, testidx
; dummy indices for testidx
nspec = n_elements(testidx) 
testidxidx = lindgen(nspec) 

; Load those spectra
specinfo = mrdfits('/home/stgilhool/APOGEE/absk_corr/absk_corr_init.fits',1)
specinfo = specinfo[testidx]

spectra_global = specinfo.spec
error_global = specinfo.err
mask_global = specinfo.mask

; Take only the bitmasks we want
; initialize the new mask to be the same size, but all 0's
mask = replicate(0, size(mask_global,/dim))
mask_bits = [0,12]
mask_dec = long(total(2L^mask_bits))
; this does a bitwise and, which will be zero if the mask_bits are not
; flagged for a given pixel, and will be non-zero if one or more of
; mask_bits is set
mask_arr = mask_global and mask_dec
mask_idx = where(mask_arr ne 0,nmasked)
; if any pixels are flagged, set their value to 1
if nmasked gt 0 then begin
    mask[mask_idx] = 1
endif

; Now, make another mask that excludes regions with too many masked
; pixels
mask_exclude = replicate(0, size(mask, /dim))
mask_smooth = double(mask)
; I want to exclude any pixel that falls in a npix_window region, with
; (npix_window/2) or more bad pixels
; Actually, changing to 3 pixel window with 2 bad pix
mask_window = 3
badpix_per_window_thresh = 2
;smooth_width = [npix_window, 0]
smooth_width = [mask_window, 0]
nbadpix_per_window = smooth(mask_smooth, smooth_width, /edge_mirror)*smooth_width[0]
exclude_idx = where(nbadpix_per_window ge badpix_per_window_thresh, nexclude)
if nexclude gt 0 then begin
    mask_exclude[exclude_idx] = 1
endif

interp_idx = where((mask eq 1) and (mask_exclude eq 0), ninterp)
interp_mask = replicate(0, size(mask, /dim))
if ninterp gt 0 then begin
    interp_mask[interp_idx] = 1
endif

spectra=spectra_global
spectra[exclude_idx] = !values.d_nan


;diff_test = mask_exclude - mask
;print, total(diff_test)
;stop

;x0 = 500
;x1 = 1000
;plot, mask[x0:x1, 0], ps = 6, yr=[-0.5, 1.5], /ys

;oplot, mask_exclude[x0:x1, 0], ps=6, color=200
;stop


; and get wl info
npix = 8575L
pix_vec = dindgen(npix)

head_str = readin_apo(hdu=0, nfiles=1)
head = *(head_str.header)
wl0 = fxpar(head, 'CRVAL1')
wl1 = fxpar(head, 'CDELT1')
log_wl = poly(pix_vec, [wl0,wl1])
wl_grid = 10d0^log_wl
;

;plot, spectra[*,0], /xs, yr=[0,2]
;stop
plot, wl_grid[feature_idx], median(spectra[feature_idx,*], dimension=2), /xs, yr=[0,2],/ys
;stop


; Calculate an EW
ew_vec = []
for i = 0, nspec-1 do begin

    ;interp over just bad pixels
    

    ew_i = ew_calc_rough(feature_range, spectra_global[*,i], error_global[*,i], mask[*,i])
    ew_vec = [ew_vec, ew_i]

endfor

; Get absolute mag
plx_vec = specinfo.plx_dtm/1d3
e_plx_vec = specinfo.eplx_dtm/1d3
k_vec = specinfo.kmag
; Get their absolute k magnitudes
absk_vec = k_vec + 5d0 + 5d0*alog10(plx_vec)
eabsk_vec = 5d0*0.434*(e_plx_vec/plx_vec)

;;;;;;;;;;;;;;;;


; Bin the metallicities
feh = specinfo.feh_irtf
efeh = specinfo.efeh_irtf

; Look at the binned metallicities
;;; A reasonable first attempt is:
;;; Fe/H < -0.1
;;; -0.1 < Fe/H < 0.0
;;;  0.0 < Fe/H < 0.1
;;;  0.0 < Fe/H < 0.1
;;; Fe/H > 0.1

bsize = 0.025

fehhist = histogram(feh, binsize = bsize, omin=om)
histx = dindgen(n_elements(fehhist))*bsize + om

plot, histx, fehhist, ps=10
;stop 

feh_bounds = [-0.1, 0.1]

fehsort_idx = sort(feh)
fehsort = feh[fehsort_idx]

bin_placement = value_locate(feh_bounds, fehsort)

; this depends on the bounds being placed within the range of feh's
nbins = n_elements(feh_bounds) + 1

binidx = ptrarr(nbins, /allocate_heap)

color_vec = [999, 99999, 9999999, 99999999, 9999999999, 99999999999, 999999999999] 

; Plot the correlations
plot, ew_vec, absk_vec, ps = 6, /yno

; Find the entries for each bin, fit a line, and overplot
for bini = 0, nbins-1 do begin
    color_bin = color_vec[bini]
    bin_val = bini-1
    bin_idx = where(bin_placement eq bin_val, nvals)
    if nvals eq 0 then message, "feh_bounds not placed properly"
    ; find the indices for feh that belong in the bini-th bin
    actual_idx = fehsort_idx[bin_idx]
    ; and store them in the pointer array
    *(binidx[bini]) = actual_idx

    ; now fit a line
    ew_bin = ew_vec[actual_idx]
    
    absk_bin = absk_vec[actual_idx]
    eabsk_bin = eabsk_vec[actual_idx]

    fit_bin_coeffs = poly_fit(ew_bin, absk_bin, 1, yfit = fit_bin, $
                              /double, measure_errors=eabsk_bin, $
                              yband=yband, yerror=yerror, sigma=sigma_bin)
    oplot, ew_bin, absk_bin, ps=6, color=color_bin
    oplot, ew_bin, fit_bin, color=color_bin
    oplot, ew_bin, fit_bin-yerror, linest=2, color=color_bin
    oplot, ew_bin, fit_bin+yerror, linest=2, color=color_bin
    
    print, ew_bin
    stop
endfor

; bin0idx = where(feh lt -0.1, n0)
; ;bin1idx = where(feh ge -0.1 and feh lt 0.0, n1)
; ;bin2idx = where(feh ge 0.0 and feh lt 0.1, n2)
; bin2idx = where(feh ge -0.1 and feh lt 0.1, n2)

; bin3idx = where(feh ge 0.1, n3)





; Plot the correlations
; plot, ew_vec, absk_vec, ps = 6, /yno
; oplot, ew_vec[bin0idx], absk_vec[bin0idx], ps = 6
; ;oplot, ew_vec[bin1idx], absk_vec[bin1idx], ps = 6, color=200
; oplot, ew_vec[bin2idx], absk_vec[bin2idx], ps = 6, color=99999
; oplot, ew_vec[bin3idx], absk_vec[bin3idx], ps = 6, color=200
stop

end
