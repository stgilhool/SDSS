; Determine the correlation between features
; absolute K-band magnitude, binned by metallicity

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


function ew_cont_fit, lcont, rcont, spectrum, error, mask
common ewinfo, wl_grid, wl_grid_over

if n_elements(rcont) eq 1 then rcont = [rcont, rcont]
if n_elements(lcont) eq 1 then lcont = [lcont, lcont]

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






pro ew_mc, left_continuum, right_continuum, spectra, error, mask, ew_mean, ew_sig, ew_arr
; Do an MC simulation to determine the distribution of EW's we get
; given a reasonable error in our definition of the continuum

;Let's start with just taking the value at left_continuum and
;right_continuum as the continuum at those points, and then fitting a
;line between them.  We'll allow the choice of left and right
;continuum to vary by +/- 2 pixels

common ewinfo, wl_grid, wl_grid_over

; Oversample spectra
oversamp = 11
npix = n_elements(wl_grid) 
nspec = n_elements(spectra[0,*]) 
xvec = dindgen(npix)
xvec_over = (dindgen(npix*oversamp)-(oversamp/2))/oversamp

wl_grid_over = interpol(wl_grid, xvec, xvec_over)

spec_over = []
for snum = 0, nspec - 1 do begin
    
    spec_over_temp = interpol(spectra[*,snum], wl_grid, wl_grid_over)
    spec_over = [[spec_over], [spec_over_temp]]
endfor

; Determine the delta_wl at each point
wl_grid_plus = interpol(wl_grid, xvec, xvec+1)
wl_grid_plus_over = interpol(wl_grid, xvec, xvec_over+1)
delta_wl = wl_grid_plus - wl_grid
delta_wl_over = wl_grid_plus_over - wl_grid_over


; Create big vectors of left and right continuum points 
ninstance = 500
;delta_cont_float = randomu(seed, ninstance)*5d0 - 2d0
;delta_cont_int = round(delta_cont_float)
;delta_cont_int_left = round(randomu(seed1, ninstance)*6d0 - 2.5d0)
;delta_cont_int_right = round(randomu(seed2, ninstance)*6d0 - 2.5d0)
delta_cont_float_left = randomu(seed1, ninstance)*2d0 - 1d0
delta_cont_float_right = randomu(seed2, ninstance)*2d0 - 1d0


;lcont_vec = delta_cont_int_left + left_continuum
;rcont_vec = delta_cont_int_right + right_continuum
lcont_vec = delta_cont_float_left + left_continuum
lcont_vec_over = value_locate(xvec_over, lcont_vec)
rcont_vec = delta_cont_float_right + right_continuum
rcont_vec_over = value_locate(xvec_over, rcont_vec)

display_left = left_continuum - 3
display_right =  right_continuum + 3
npix_display = display_right - display_left + 1
display_x = dindgen(npix_display) + display_left

; Loop through and fit a continuum line to each one
ew_arr = dblarr(ninstance, nspec)
ew_mean = dblarr(nspec)
ew_sig = dblarr(nspec)
for snum = 0, nspec - 1 do begin
    spec_i = spectra[*,snum]
    spec_over_i = spec_over[*,snum]
    err_i = error[*,snum]
    mask_i = mask[*,snum]

    display_spec = spec_i[display_x]
   

    ; Loop through each trial
    for inum = 0, ninstance - 1 do begin
        ; Fit continuum
        ;contco_i = ew_cont_fit(lcont_vec[inum], rcont_vec[inum],
        ;spec_i, err_i, mask_i)
        contco_i = ew_cont_fit(lcont_vec_over[inum], rcont_vec_over[inum], spec_over_i, err_i, mask_i)
        
        ; Subtract spectrum from continuum
        ;npix_feat = rcont_vec[inum]-lcont_vec[inum]+1
        npix_feat = rcont_vec_over[inum]-lcont_vec_over[inum]+1
        ;feat_x = dindgen(npix_feat) + lcont_vec[inum]
        feat_x = dindgen(npix_feat) + lcont_vec_over[inum]
        cont_i = poly(feat_x, contco_i)
        feat_spec = spec_over_i[feat_x]
        feat_ew = cont_i - feat_spec
        
        dwl = delta_wl_over[feat_x]

        feat_x_dec = xvec_over[feat_x]
        
        ; Fit gaussian OR
        
        ; Numerically integrate
        ew_i = total(feat_ew*dwl, /double)
        plot, display_x, display_spec, /xs, yr=[0.3, 1.1], /ys, ps=6
        oplot, feat_x_dec, feat_spec, color=200, ps=6
        oplot, feat_x_dec, cont_i, color=99999
        case snum of
            22: begin
                print, "G err"
                wait, 0.01
            end
            16: begin
                print, "G err"
                wait, 0.01
            end
            2: begin
                print, "G err"
                wait, 0.01
            end
            14: begin
                print, "B err"
                wait, 0.01
            end
            18: begin
                print, "B err"
                wait, 0.01
            end
            4: begin
                print, "B err"
                wait, 0.01
            end
            else: begin
            end
        endcase
        
        ;wait, 0.01
        ew_arr[inum, snum] = ew_i
    endfor
; Get the distribution of EW's
; View the distribution
; Fit a gaussian to the distribution
ew_mean[snum] = mean(ew_arr[*,snum])
ew_sig[snum] = stddev(ew_arr[*,snum])

endfor

end


pro show_ew_dist, snum, bsize
common distinfo, ew_arr, ew_mean, ew_sig

ew_vec = ew_arr[*,snum]

hist = histogram(ew_vec, binsize = bsize, omin = om)
histx = dindgen(n_elements(hist))*bsize + om

plot, histx, hist, ps=10, title = 'Mean: '+strtrim(ew_mean[snum],2)+ $
  ' | StdDev: '+strtrim(ew_sig[snum],2)

stop

end















function ew_calc_errtest, window_range, spectrum, $
                        error, mask

common ewinfo, wl_grid, wl_grid_over

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

pro determine_which_spectra_to_use_all, testidx
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
                finite(info.plx_dtm), ntest)



end

pro absk_newtry


; OVERALL: fit a gaussian at each pixel, with amplitude and sigma as
; free parameters.  Then determine the points at which the fit is best
; for the spectra in question.  This should tell us where the good
; features are



;;; Load the set of most promising spectra

; get indices of spectra to use
determine_which_spectra_to_use, testidx_cal
determine_which_spectra_to_use_all, testidx



; dummy indices for testidx
nspec = n_elements(testidx)
nspec_cal = n_elements(testidx_cal)+19
testidxidx = lindgen(nspec) 

; find the indices of the calibration set in testidx
cal_idx = value_locate(testidx, testidx_cal)  
cal_idx = lindgen(nspec)

; Load those spectra
specinfo = mrdfits('/home/stgilhool/APOGEE/absk_corr/absk_corr_init.fits',1)
specinfo = specinfo[testidx]

spectra_global = specinfo.spec
error_global = specinfo.err
mask_global = specinfo.mask

common ew_info, delta_wl, wl_grid


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
delta_wl = wl_grid-shift(wl_grid,1)


;;; Apply masks and interpolate over isolated pixels

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




;;; Determine which regions to analyze (let's say 20/25 spectra must
;;; have a clean span of 11 pixels about the centroid)
;;; What's written is a bit different, but should be good enough

; Take a median spectrum for display, and finding out which pixels are
; masked in which spectrum
mask_exclude_sum = total(mask_exclude, 2)
; exclude this pixel if it's masked in 20 or more spectra
exclude_col = where(mask_exclude_sum ge 20, nexc)

medspec_unmasked = median(spectra, dimension=2)
medspec = median(spectra, dimension=2)

if nexc gt 0 then begin
    mask_exclude_sum[exclude_col] = 1
    medspec[exclude_col] = !values.d_nan
endif


plot, medspec_unmasked, /xs, /ys, yr=[0.5,1.1], linest=2
oplot, medspec, color=200

;;; Fit Gaussians
nterms=5
amp_fit = dblarr(npix, nspec)
ctr_fit = dblarr(npix, nspec)
sig_fit = dblarr(npix, nspec)
cnt0_fit = dblarr(npix, nspec)
cnt1_fit = dblarr(npix, nspec)
gfitco_all = dblarr(npix, nspec, nterms)
status_fit = dblarr(npix, nspec)
npegged_fit = dblarr(npix, nspec)

ew_fit = dblarr(npix, nspec)

chi_fit = dblarr(npix, nspec)

;HARDCODE TAKE OUT
exclude_idx = lindgen(npix)
;keep_idx =
;[463,1053,1132,1297,1348,1473,1531,1546,2177,2190,2244,2463,2657,2772,2785]
;keep_idx = [1843,2319,2474,4047,5125,5894,6790]
;keep_idx = [943,1133,1843,2112,2189,2319,2474,4047,5125,5894,6790]
keep_idx =[943,1133,1843,2112,2189,2319,2474,2908,3039,3062,4047,5125,5894,6790]
;keep_idx = [943,1133,1843,2112,2189,2319,2474,2908,3039,3062,4047,5125,5894]
nkeep = n_elements(keep_idx) 
exclude_idx[keep_idx] = 0

for pixnum = 0, npix-1 do begin
    ; Determine if we want to fit a Gaussian on this pixel
    mask_check = where(exclude_idx eq pixnum, skip_pix)
    if skip_pix then begin
        amp_fit[pixnum,*] = !values.d_nan
        ctr_fit[pixnum,*] = !values.d_nan
        sig_fit[pixnum,*] = !values.d_nan
        cnt0_fit[pixnum,*] = !values.d_nan
        cnt1_fit[pixnum,*] = !values.d_nan
        gfitco_all[pixnum,*,*] = !values.d_nan
        status_fit[pixnum,*] = -1
        npegged_fit[pixnum,*] = !values.d_nan
        chi_fit[pixnum,*] = !values.d_nan
        continue
    endif else begin
        ; Take chunk of spectrum
        npix_window = 9 ;HARDCODE
        npix_halfwin = npix_window/2
        oversamp = 7
        
        win_x = lindgen(npix_window) + pixnum - npix_halfwin
        win_x_over = (dindgen(npix_window*oversamp)-(oversamp/2))/oversamp + pixnum - npix_halfwin
        
        
        parinfo_str = {limited:[0,0], limits:[0d0,0d0], fixed:0}
        parinfo = replicate(parinfo_str, nterms)
        
        parinfo[0].limited = [1,1]
        parinfo[0].limits  = [-5d0, 0d0]
        parinfo[1].limited = [1,1]
        parinfo[1].limits  = [double(pixnum)-1.5d0, double(pixnum)+1.5d0]
        parinfo[2].limited = [1,1]
        parinfo[2].limits  = [0.01d0, 1.5d0]
	parinfo[3].limited = [1,1]
        parinfo[3].limits  = [0.85d0, 1.15d0]

        ;parinfo[1].fixed = 1 ; HARDCODE
        clock = tic()
        for snum = 0, nspec-1 do begin
            spec_chunk = spectra[win_x, snum]
            
            fin_spec_chunk = finite(spec_chunk)
            spec_chunk_nan = where(fin_spec_chunk eq 0, nnan)
            
             ;Disqualify the chunk if there are 3 or more bad pixels
             if nnan ge 3 then begin
                 amp_fit[pixnum,snum] = !values.d_nan
                 ctr_fit[pixnum,snum] = !values.d_nan
                 sig_fit[pixnum,snum] = !values.d_nan
                 cnt0_fit[pixnum,snum] = !values.d_nan
                 cnt1_fit[pixnum,snum] = !values.d_nan
                 gfitco_all[pixnum,snum,*] = !values.d_nan
                 status_fit[pixnum,snum] = -1
                 npegged_fit[pixnum,snum] = !values.d_nan
                 chi_fit[pixnum,snum] = !values.d_nan
                 continue
             endif else begin
                
                
                spec_chunk_over = interpol(spec_chunk, win_x, win_x_over)
                
                estimates = dblarr(nterms)
                estimates[0] = spec_chunk[npix_halfwin] - 1d0
                estimates[1] = double(pixnum)
                estimates[2] = 1d0
                estimates[3] = 1d0
                estimates[4] = 1d-2
                
                
                
                gfit = mpfitpeak(win_x_over, spec_chunk_over, gfitco, nterms=nterms, $
                                 estimates=estimates, parinfo=parinfo, chisq=chi2, $
                                 npegged=npegged, status=status, /nan)
               
                

                amp_fit[pixnum,snum] = gfitco[0]
                ctr_fit[pixnum,snum] = gfitco[1]
                sig_fit[pixnum,snum] = gfitco[2]
                cnt0_fit[pixnum,snum] = gfitco[3]
                cnt1_fit[pixnum,snum] = gfitco[4]
                gfitco_all[pixnum,snum,*] = gfitco
                status_fit[pixnum,snum] = status
                ;npegged_fit[pixnum,snum] = npegged
                npegged_fit[pixnum,snum] = !values.d_nan
            ;;; Collect chi2's
                chi_fit[pixnum,snum] = chi2
                
                cont = poly(win_x_over, gfitco[3:4])
                
                
                ew_sum = ew_calc_rough(minmax(win_x), spectra_global[*,snum], $
                                       error_global[*,snum], mask[*,snum])

                ;ew_sum = ew_calc([gfitco[3],gfitco[4]],minmax(win_x), $
                                 ;minmax(win_x), spectra_global[*,snum], $
                                 ;error_global[*,snum], mask[*,snum])
                ew = total(cont-gfit, /double)

                if finite(ew) eq 0 then ew = 0d0
                ;ew_fit[pixnum,snum] = ew
                ew_fit[pixnum,snum] = ew_sum
                
                plot, win_x_over, spec_chunk_over, /xs, /ys, yr = [0.5, 1.1]
                oplot, win_x_over, gfit, color=200
                print, ew               
;                dxstop
            endelse  
        endfor
        print, strtrim(pixnum,2)
        toc, clock
    endelse
endfor


; Get median, mean and stddev of chi_fit columns
med_chi = median(chi_fit, dimension=2)
mean_chi = mean(chi_fit, dimension=2, /nan)
sig_chi = stddev(chi_fit, dimension=2, /nan)

med_gfitco = median(gfitco_all, dimension=2)
mean_gfitco = mean(gfitco_all, dimension=2, /nan)
sig_gfitco = stddev(gfitco_all, dimension=2, /nan)

plot, med_chi, ps=6
;stop


for pixel = 0, nkeep-1 do begin

pixnum = keep_idx[pixel]

ew_vec = ew_fit[pixnum, *]



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
;oploterror, ew_vec, absk_vec, ew_sig, eabsk_vec
oploterror, ew_vec, absk_vec, eabsk_vec


; Find the entries for each bin, fit a line, and overplot
; for bini = 0, nbins-1 do begin
;     color_bin = color_vec[bini]
;     bin_val = bini-1
;     bin_idx = where(bin_placement eq bin_val, nvals)
;     if nvals eq 0 then message, "feh_bounds not placed properly"
;     ; find the indices for feh that belong in the bini-th bin
;     actual_idx = fehsort_idx[bin_idx]
;     ; and store them in the pointer array
;     *(binidx[bini]) = actual_idx

;     ; now fit a line
;     ew_bin = ew_vec[actual_idx]
    
;     absk_bin = absk_vec[actual_idx]
;     eabsk_bin = eabsk_vec[actual_idx]
;     fit_bin_coeffs = poly_fit(ew_bin, absk_bin, 1, yfit = fit_bin, $
;                               /double, measure_errors=eabsk_bin, $
;                               yband=yband, yerror=yerror, sigma=sigma_bin)
;     oplot, ew_bin, absk_bin, ps=6, color=color_bin
;     oplot, ew_bin, fit_bin, color=color_bin
;     oplot, ew_bin, fit_bin-yerror, linest=2, color=color_bin
;     oplot, ew_bin, fit_bin+yerror, linest=2, color=color_bin
;     color2 = ['red', 'yellow', 'cyan']
;     ;oploterror, ew_bin, absk_bin, ew_sig[actual_idx], eabsk_bin,
;     ;errcolor=color2[bini]
;     oploterror, ew_bin, absk_bin, eabsk_bin, errcolor=color2[bini]
;     print, ew_bin
;     ;stop
; endfor

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
dxstop
endfor

vis = 1

;;;;; Changing below to use only the calibration set

; I want to include the faintest star in the calibration set
;faintest_idx = where(absk_vec eq max(absk_vec), nfaint)
;if nfaint eq 1 then cal_idx=[cal_idx, faintest_idx] else message,"faintest error"

ew_arr = ew_fit[keep_idx,*]

ew_mean = mean(ew_arr, dimension=2,/nan)

ew_arr_cal = ew_arr[*, cal_idx]

feh_col = reform(feh, 1, nspec)
feh_col_cal = reform(feh[cal_idx], 1, nspec_cal)

eq_width = [reform(replicate(1d0, nspec),1,nspec), ew_arr, feh_col]
eq_width_cal = [reform(replicate(1d0, nspec_cal),1,nspec_cal), $
                ew_arr_cal, feh_col_cal]

absk_vec_cal = absk_vec[cal_idx]
eabsk_vec_cal = eabsk_vec[cal_idx]



;n_eq_width_coeffs = n_elements(eq_width[*,0]) 
n_eq_width_coeffs = n_elements(eq_width_cal[*,0]) 

common fit_outputs, model, res

functargs = {ew:eq_width_cal, $
             absk:absk_vec_cal, $
             err:eabsk_vec_cal, $
             output:'chi2', $
             vis:vis $
            }
    
parinfo = replicate({value:0.1d0},n_eq_width_coeffs)
parinfo[0].value=5d0

!p.multi = [0,1,2]


r = mpfit('absk_fit', parinfo=parinfo, functargs=functargs, status=status, dof=dof, bestnorm=chi2)

print, "result:"
print, r

;absk_model = absk_fit(r, ew=functargs.ew, absk=functargs.absk,
;err=functargs.err, $
absk_model = absk_fit(r, ew=eq_width, absk=absk_vec, err=eabsk_vec, $
                      output='model', /vis)

absk_model_cal = absk_model[cal_idx]
!p.multi=0
window, 0, xs = 1500, ys=800

aksort = sort(absk_vec)
aksort_cal = sort(absk_vec_cal)

rmse = stddev(absk_model-absk_vec, /nan)
rmse_cal = stddev(absk_model_cal-absk_vec_cal, /nan)

filenum = lindgen(nspec)
filenum_cal = value_locate(absk_vec[aksort], absk_vec_cal[aksort_cal])

; set_plot, 'ps'
;  device, filename = 'absk_corr_fit.eps'
;  device, /color, bits=8
;  device, xs = 13, ys= 8, /inches
;  loadct, 13



plot, filenum, absk_vec[aksort], ps = 6, title=$
  "Agreement between fit and data | RMSE(all) = "$
  +strtrim(rmse,2) + " mags | RMSE(cal) = "+ strtrim(rmse_cal,2) + " mags",$
  xtit = "Spectrum", ytit = "M_K", yr = [(min(absk_vec)<min(absk_model,/nan))-1, (max(absk_vec)>max(absk_model,/nan))+1]
oploterror, filenum, absk_vec[aksort], eabsk_vec[aksort], ps=6
oplot, filenum, absk_vec[aksort], ps=6
oploterror, filenum_cal, absk_vec_cal[aksort_cal], eabsk_vec_cal[aksort_cal], ps = 6, color=200
oplot, filenum_cal, absk_vec_cal[aksort_cal], ps=6, color=200
oplot, absk_model[aksort], ps=6, color=400

;device, /close_file

;print, "plot printed"






;plot, reform(eq_width[1,*]), absk_model, ps=6, $
;  title="EW vs. M_K(fit) | Red pts are data", $
;  xtit='M_K(fit)', ytit='EW', charsize = 1.5
;oplot, reform(eq_width[1,*]), absk_vec, ps=4, color=200


stop



plot, lindgen(n_elements(mean_chi)), mean_chi, ps=6, yr=[-0.1,0.1]
oploterror, lindgen(n_elements(mean_chi)), mean_chi, sig_chi

stop
;;; Develop some criteria for what constitutes a good region
    ;;; Display the features and their fits


stop


;common ewinfo, wl_grid, wl_grid_over
;common distinfo, ew_arr, ew_mean, ew_sig
; HARDCODE
;feature0 = 300
;npix_window = 75
; feature1 = feature0+npix_window

; feature_range = [feature0, feature1]

; feature_idx = lindgen(feature1-feature0)+feature0
feature_idx = lindgen(right_continuum-left_continuum+1)+left_continuum

; testidx = []






;diff_test = mask_exclude - mask
;print, total(diff_test)
;stop

;x0 = 500
;x1 = 1000
;plot, mask[x0:x1, 0], ps = 6, yr=[-0.5, 1.5], /ys

;oplot, mask_exclude[x0:x1, 0], ps=6, color=200
;stop



    

;plot, spectra[*,0], /xs, yr=[0,2]
;stop
;plot, wl_grid[feature_idx],
;median(spectra[feature_idx,*],dimension=2), /xs, yr=[0,2],/ys

plot, feature_idx, medspec_all, $
  /xs, yr=[0,2],/ys, linest = 2
oplot, feature_idx, medspec, color=200
;plot, feature_idx, median(spectra[feature_idx,*], dimension=2), /xs, yr=[0.3,1.1],/ys
;stop


; Calculate an EW
;Do an MC approach at defining the continuum pixels, and the value
;of the continuum at that point
;ew_dist = ew_mc(left_continuum, right_continuum, spectra,
;error_global, mask)
ew_mc, left_continuum, right_continuum, spectra, error_global, mask, ew_mean, ew_sig, ew_arr



;    ew_i = ew_calc_rough(feature_range, spectra_global[*,i], error_global[*,i], mask[*,i])
;    ew_vec = [ew_vec, ew_i]
ew_vec = ew_mean

;endfor

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
oploterror, ew_vec, absk_vec, ew_sig, eabsk_vec

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
    color2 = ['red', 'yellow', 'cyan']
    oploterror, ew_bin, absk_bin, ew_sig[actual_idx], eabsk_bin, errcolor=color2[bini]
    print, ew_bin
    ;stop
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
