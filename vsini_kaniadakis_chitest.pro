; Also calculates the expected number of non-detections and
; probability of getting none

; Changed binsize to 100 K and range tops out at 3900K

; Added cut on VSINI_BAD and VSINI_WARN flags

; Added a fitted line to the result


; Implementing maximum likelihood estimation for censored
; data. Previous version just fit a distribution using mpfit
; Added penalties to keep parameters within appropriate limits
; Added a goto statement which re-fits if the amoeba fit fails
; (that is, no nan's in the results)

; Reproducing Cullen's vsini dist parameter vs. Teff analysis

; Cullen's description of his analysis.  
; Hi Steve,

; OK, here is the first plot from my analysis.

; Basically, what I've done here is take the DR13 stars that meet the following criteria:

; flag bit 19 set
; 2600<teff<4000 (so, excluding the -99999 guys)
; vsini > 0.

; This gives me 719 guys.

; Then, for each bin of 200K (starting at 2600) take the vsini values and fit a three parameter Kaniadakis function (see the paper from Carvalho et al. 2009 - equations 4 and 8). The parameters are the kappa, the sigma, and a center of the distribution. I did this in a maximum likelihood sense. For the ones below 5 km/s, I integrated the Kaniadakis function and included a liklihood for all of those non-detections combined.

; I bootstrapped the errors on the parameters by running the same fit 100 times, but adding a noise  of 2 km/s for the > 5km/s guys, and shuffling the 5 km/s guys around. I checked to make sure that the bootstrap didn't push a guy with, say, 7 km/s below 5 km/s. I think that this is the correct way to handle this, but I am happy to hear arguments otherwise.

; The plot shows a clear correlation between the center of the distribution (V_0sini) and teff. I have overplotted approximate spectral types. Also, here is a table of info about temps, spectral types, etc that I put together from a range of sources.



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Fit the Kaniadakis function to the vsini distribution
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function kaniadakis_fit, params

common amoeba_info, vsini_boot, vsini_xvec, vsini_xvec_midpt, d_vsini, nnondet, thresh, vsini_nondet, vis, chi2, oversamp, chi2dof
; Get params
vsini_0 = params[0]
kappa = params[1]
sigma = params[2]

; Make pdf of kaniadakis distribution
y = vsini_xvec
yp = y-vsini_0

sqrt_term = sqrt(1d0 + (kappa^2)*(yp^4)/(sigma^4))
other_term = kappa * (-1d0 * (yp^2)/(sigma^2))
combined_term = (sqrt_term + other_term)^(1/kappa)

kan_dist = y * combined_term

; Normalize the kan_dist pdf
pdf = kan_dist/total(kan_dist*d_vsini)

; Get the CDF
cdf = total(pdf*d_vsini, /cumulative)

; Integrate (downsample) the pdf for chi2 calculation purposes
integrated_pdf = downsample_tophat(pdf, oversamp)

; Get the likelihood for a single non-detection
censored_prob = double(interpol(cdf, vsini_xvec, thresh))

; Get the probability of each detection, given the distribution
bigl_vec = interpol(pdf, vsini_xvec, vsini_boot)

; Combine the detection probs and non-detection probs
if nnondet gt 0 then begin
    ;vsini_nondet1 = randomu(seed, nnondet)*thresh
    ;bigl_nondet = interpol(pdf, vsini_xvec, vsini_nondet1)
    ;bigl_nondet = interpol(cdf, vsini_xvec, vsini_nondet)
    bigl_nondet = replicate(censored_prob, nnondet)
    bigl_vec = [bigl_nondet, bigl_vec]
    vsini_all = [vsini_nondet,vsini_boot]
    nspec = n_elements(vsini_all) 
    chi2_vec_0 = ((nnondet-(nspec*censored_prob))^2d0)/(nspec*censored_prob)
endif else begin
    vsini_all = vsini_boot
    chi2_vec_0 = 0d0
endelse

integrated_pdf_remainder = integrated_pdf[fix(thresh):-1]
expected_vec = integrated_pdf_remainder * nspec

actual_counts = histogram(vsini_boot, binsize = 1, min=fix(thresh), nbins=n_elements(expected_vec))

;chi2_vec_remainder = ((vsini_boot - expected_vec)^2d0)/expected_vec
chi2_vec_remainder = ((actual_counts - expected_vec)^2d0)/expected_vec

chi2_vec = [chi2_vec_0, chi2_vec_remainder]
;help, actual_counts
;help, expected_vec
;help, chi2_vec
;dxstop

; Calculate chi2
chi2 = total(chi2_vec, /double)

dof = n_elements(chi2_vec)-3

chi2dof = chi2/double(dof)


if vis then begin
    
    ;histx = dindgen(max(vsini_xvec))
    ;histy = histogram(vsini_all, min=0d0)
    histx = vsini_xvec_midpt
    histy = histogram(vsini_all, min=0d0, nbins=n_elements(vsini_xvec_midpt) )
    plot, histx, histy, ps=10
    oplot, vsini_xvec, pdf*n_elements(vsini_all), color=200
    plot, histx, total(histy, /cumulative), ps = 10
    oplot, vsini_xvec, cdf*n_elements(vsini_all), color=200
    wait, 0.0001
endif

; Take log likelihood
logl_vec = alog(bigl_vec)
; Sum the vector of log likelihoods
logl = total(logl_vec)

if vsini_0 lt 0 then v_penalty = abs(logl)*0.1d0*exp(-1d0*vsini_0) $
else v_penalty=0
if kappa lt 0 then klow_penalty = abs(logl)*0.1d0*exp(-1d0*kappa) $
else klow_penalty=0
if kappa gt 1 then khigh_penalty = abs(logl)*0.1d0*exp(kappa-1d0) $
else khigh_penalty=0
if sigma lt 0 then s_penalty = abs(logl)*0.1d0*exp(-1d0*sigma) $
else s_penalty=0

return, -1d0*logl+v_penalty+klow_penalty+khigh_penalty+s_penalty

end


; Calculate the expected number of non-detections, and the probability
; of 0 non-detections
function kaniadakis_analysis, params

common amoeba_info, vsini_boot, vsini_xvec, vsini_xvec_midpt, d_vsini, nnondet, thresh, vsini_nondet, vis, chi2, oversamp, chi2dof
; Get params
vsini_0 = params[0]
kappa = params[1]
sigma = params[2]

; Make pdf of kaniadakis distribution
y = vsini_xvec
yp = y-vsini_0

sqrt_term = sqrt(1d0 + (kappa^2)*(yp^4)/(sigma^4))
other_term = kappa * (-1d0 * (yp^2)/(sigma^2))
combined_term = (sqrt_term + other_term)^(1/kappa)

kan_dist = y * combined_term

; Normalize the kan_dist pdf
pdf = kan_dist/total(kan_dist*d_vsini)

; Get the CDF
cdf = total(pdf*d_vsini, /cumulative)


; Get the probability for a single non-detection
censored_prob = double(interpol(cdf, vsini_xvec, thresh))
; prob of a single detection
detection_prob = 1d0 - censored_prob

; Probability for 0 non-detections
nstars = n_elements(vsini_boot) + nnondet
all_det_prob = detection_prob^nstars

; Expect number of non-detections
n_nondet_expect = censored_prob*nstars


; Also do the same for number of points below minimum of data or
; thresh, whichever is highest
if nnondet gt 0 then begin
    min_data = thresh
    min_nondet_prob = double(interpol(cdf, vsini_xvec, min_data))
    min_det_prob = 1d0 - min_nondet_prob
endif else begin
    min_data = min(vsini_boot)
    min_nondet_prob = double(interpol(cdf, vsini_xvec, min_data))
    min_det_prob = 1d0 - min_nondet_prob
endelse

; Prob to have no data below the minimum
none_below_prob = min_det_prob^nstars
; Expected number of data below the minimum
n_below_expect = min_nondet_prob * nstars



return, [all_det_prob, n_nondet_expect, none_below_prob, n_below_expect]

end


function kaniadakis_out, params

common amoeba_info, vsini_boot, vsini_xvec, vsini_xvec_midpt, d_vsini, nnondet, thresh, vsini_nondet, vis, chi2, oversamp, chi2dof
; Get params
vsini_0 = params[0]
kappa = params[1]
sigma = params[2]

; Make pdf of kaniadakis distribution
y = vsini_xvec
yp = y-vsini_0

sqrt_term = sqrt(1d0 + (kappa^2)*(yp^4)/(sigma^4))
other_term = kappa * (-1d0 * (yp^2)/(sigma^2))
combined_term = (sqrt_term + other_term)^(1/kappa)

kan_dist = y * combined_term

; Normalize the kan_dist pdf
pdf = kan_dist/total(kan_dist*d_vsini)

; Get the CDF
cdf = total(pdf*d_vsini, /cumulative)

return, [[pdf],[cdf]]

end











pro vsini_kaniadakis_chitest, vis=vis

skip_to_plot = 1

if skip_to_plot then goto, plot_the_data

if n_elements(vis) eq 0 then vis=1


; Read in the ASPCAP data
aspcap_file = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'
aspcap_info = mrdfits(aspcap_file,1)

; Make list of Teff b/w 2600 and 4000, and vsini gt 0
teff_all = aspcap_info.teff
teff_err_all = aspcap_info.teff_err
vsini_all = aspcap_info.vsini

aspcapflag = aspcap_info.aspcapflag

;cut the aspcap VSINI_WARN and VSINI_BAD guys
vsini_warn_bit = 14
vsini_bad_bit = 30
star_bad_bit = 23
;flag_bits = [vsini_warn_bit, vsini_bad_bit]
flag_bits = [vsini_warn_bit, star_bad_bit, vsini_bad_bit]
;flag_bits = [star_bad_bit]
flag_dec = long(total(2L^flag_bits))

flag_idx = where((aspcapflag and flag_dec) ne 0, nflg, complement=unflag_idx, ncomplement=nunflg)

sel_idx = where(teff_all gt 3500 and teff_all lt 4000 and vsini_all gt 0 and ((aspcapflag and flag_dec) eq 0), nsel)
;sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0, nsel)

print, nsel
dxstop


;sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0, nsel)
; take just our selection
sel_info = aspcap_info[sel_idx]
teff = sel_info.teff
teff_err = sel_info.teff_err
vsini = sel_info.vsini

maxvsini = max(vsini)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; Fit the Kaniadakis function 100 times
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

param_arr = dblarr(3, 100)
status_arr = lonarr(100)
chi2_arr = dblarr(100)
chi2_dof_arr = dblarr(100)
nondet_prob_arr = dblarr(100)
nondet_expect_arr = dblarr(100)

mindet_prob_arr = dblarr(100)
below_expect_arr = dblarr(100)
min_boot = dblarr(100)
nstars_per_bin = nsel



common amoeba_info, vsini_boot, vsini_xvec, vsini_xvec_midpt, d_vsini, nnondet, thresh, vsini_nondet, vis1, chi2, oversamp, chi2dof
vis1=vis

thresh = 5d0

; make a vsini x vector for the pdf later
minvsini = 0d0
vsini_range = maxvsini+5d0-minvsini
oversamp=11
;d_vsini = 0.1d0
d_vsini = 1d0/oversamp
npts_xvec = fix(vsini_range/d_vsini)


;vsini_xvec = maken(minvsini, maxvsini, npts_xvec)
vsini_xvec = dindgen(vsini_range*oversamp)/oversamp + (1/(2d0*oversamp)) 

vsini_xvec_midpt= downsample_tophat(vsini_xvec, oversamp)

vsini_base = vsini
    
nvel = n_elements(vsini_base) 

nstars_per_bin = nvel

; Make vectors of vsini_detected and vsini_nondetected
det_idx = where(vsini_base gt thresh, ndet, complement=nondet_idx, $
                ncomplement=nnondet)
    
nondet_data_arr = nnondet
vsini_det = vsini_base[det_idx]
if nnondet gt 0 then vsini_nondet = vsini_base[nondet_idx]

; Loop 100 times
for iter = 0, 99 do begin
    refit:
    ; Add error (of 2 km/s)
    chi2=[]
    std_error = 2d0
    vsini_err = randomn(seed, ndet)*std_error
    
    vsini_boot = vsini_det + vsini_err
    
    ; Make sure the detections stay in the detection range
    repeat begin
        too_low_idx = where(vsini_boot lt thresh, ntoo_low)
        if ntoo_low gt 0 then begin
            new_err = randomn(seed, ntoo_low)
            vsini_boot[too_low_idx] = vsini_det[too_low_idx]+new_err
        endif
    endrep until ntoo_low eq 0
    
    
    vsini_guess = median(vsini_boot)
    kappa_guess = 0.5d0
    sigma_guess = 4d0
    
    vsini_scale = stddev(vsini_boot)
    kappa_scale = 0.5d0
    sigma_scale = 3d0
    
    guess = [vsini_guess, kappa_guess, sigma_guess]
    scale = [vsini_scale, kappa_scale, sigma_scale]
    
    ; Fit
    !p.multi = [0,1,2]
    result = amoeba3(1d-10, function_name = 'kaniadakis_fit', p0=guess, $
                     scale=scale, nmax=1000)
    !p.multi = 0
    if n_elements(result) gt 1 then begin
        param_arr[*,iter] = result
        status_arr[iter] = 1
        
        prob_analysis = kaniadakis_analysis(result)
        nondet_prob_arr[iter] = prob_analysis[0]
        nondet_expect_arr[iter] = prob_analysis[1]
        mindet_prob_arr[iter] = prob_analysis[2]
        below_expect_arr[iter] = prob_analysis[3]
        min_boot[iter] = min(vsini_boot)
        chi2_arr[iter] = chi2
        chi2_dof_arr[iter] = chi2dof
        
    endif else if n_elements(result) eq 1 then begin
        goto, refit
        ;param_arr[*,bin,iter] = replicate(!values.d_nan,3)
        ;status_arr[bin,iter] = -1
        
    endif
    print, result
endfor


save, /variables, filename='vsini_kan_chitest.sav'

plot_the_data:

if skip_to_plot then restore, 'vsini_kan_chitest.sav'

; For some reason, the usual ps approach didn't work
  set_plot, 'ps'
  device, filename = 'vsini_kan_chitest.eps', /encapsulated
 device, /color, bits=8
  device, xs = 13, ys= 8, /inches
  loadct, 13

mean_param_arr = mean(param_arr, dimension=2)

final_dist_arr = kaniadakis_out(param_arr)

final_pdf = final_dist_arr[*,0]
final_cdf = final_dist_arr[*,1]


data_hist = histogram(vsini, min=0d0, nbins=n_elements(vsini_xvec_midpt) )



;!p.multi=[0,1,2]
;window,0,xs=1200,ys=900
plot, vsini_xvec_midpt, data_hist, ps=6, xr=[4d0, 7d1], yr=[1d-2, 50d0], /ys, /xs, /xlog, xtitle = "Vsini", ytit = "Number of Detections per bin", $
  title = "Comparison of Kaniadakis Distribution to Data", charsize = 2
oplot, vsini_xvec, final_pdf*nsel
;plot, vsini_xvec_midpt, total(data_hist,/cumulative), ps=6, xr=[4d0, 7d1], /ylog,yr=[5d1, 3d2], /ys, /xs
;oplot, vsini_xvec,final_cdf*nsel

print, mean(chi2_arr)
print, minmax(chi2_arr)

print, mean(chi2_dof_arr)
print, minmax(chi2_dof_arr)


device, /close_file

stop

set_plot, 'ps'
device, filename = 'vsini_kaniadakis_chitest.eps', /encapsulated
device, /color, bits=8
device, xs = 13, ys= 8, /inches
loadct, 13



device, /close_file














print, teff_bin_vec_ctr
print, vsini_mean
print, nondet_data_arr
print, mean(nondet_expect_arr, dimension=2)
print, mean(nondet_prob_arr, dimension=2)
print, mean(below_expect_arr, dimension=2)
print, mean(mindet_prob_arr, dimension=2)


nondet_expected = mean(nondet_expect_arr, dimension=2)
nondet_probability = mean(nondet_prob_arr, dimension=2)
belowmin_expected = mean(below_expect_arr, dimension=2)
nonebelow_probability = mean(mindet_prob_arr, dimension=2)
min_vsini_arr = mean(min_boot, dimension=2)

big_arr = [[teff_bin_vec_0], [teff_bin_vec_1], [vsini_mean], [vsini_stddev], [nstars_per_bin],[nondet_data_arr], [min_vsini_arr], [nondet_expected], [nondet_probability], [belowmin_expected], [nonebelow_probability]]

print, transpose(big_arr)

stop

outstr = {teff_bin_ctr:teff_bin_vec_ctr, $
          vsini_mean:vsini_mean, $
          vsini_stddev:vsini_stddev, $
          nstars_per_bin:nstars_per_bin, $
          nnondet:nondet_data_arr, $
          min_vsini:min_vsini_arr, $
          nondet_expected:nondet_expected, $
          nondet_probability:nondet_probability, $
          belowmin_expected:belowmin_expected, $
          nonebelow_probability:nonebelow_probability}

mwrfits, outstr, 'vsini_kaniadakis_chitest.fits', /create

stop
end
















