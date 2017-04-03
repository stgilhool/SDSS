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

common amoeba_info, vsini_boot, vsini_xvec, d_vsini, nnondet, thresh, vsini_nondet, vis
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
    vsini_all = [vsini_boot,vsini_nondet]
endif else vsini_all = vsini_boot

if vis then begin
    
    histx = dindgen(max(vsini_xvec))
    histy = histogram(vsini_all, min=0d0)
    plot, histx, histy, ps=10
    oplot, vsini_xvec, pdf*n_elements(vsini_all), color=200
    plot, histx, total(histy, /cumulative), ps = 10
    oplot, vsini_xvec, cdf*n_elements(vsini_all), color=200
    wait, 0.001
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











pro vsini_kaniadakis_3, vis=vis

skip_to_plot = 0

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

;sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0 and ((aspcapflag and flag_dec) eq 0), nsel)
sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0, nsel)

print, nsel
dxstop


;sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0, nsel)
; take just our selection
sel_info = aspcap_info[sel_idx]
teff = sel_info.teff
teff_err = sel_info.teff_err
vsini = sel_info.vsini

maxvsini = max(vsini)

; Bin by Teff
binsize = 200
min_temp = 2600

teff_hist = histogram(teff, binsize = 200, min = min_temp, reverse_indices=ri)

nbins = n_elements(teff_hist)

; teff bin vectors (low bound, high bound and center)
teff_bin_vec_0 = lindgen(nbins)*binsize + min_temp
teff_bin_vec_1 = (lindgen(nbins)+1)*binsize + min_temp
teff_bin_vec_ctr = (lindgen(nbins)+0.5)*binsize + min_temp

bin_num_vec = lindgen(nbins)

; Initialize hash to store data
data = hash(bin_num_vec, bin_num_vec)                   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Use a histogram to sort and bin data into a structure  ;;;
;;; (adapt this into a stand-alone function)               ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Loop through reverse indices and get the values and indices for each
; bin
for bin = 0, nbins-1 do begin
    
    ri0 = ri[bin]
    ri1 = ri[bin+1]

    if ri1 gt ri0 then begin
        ; get all the data
        bin_idx = ri[ri0:ri1-1]
        bin_aspcap_idx = sel_idx[bin_idx]
        bin_teff = teff[bin_idx]        
        bin_teff_err = teff_err[bin_idx]        
        bin_vsini = vsini[bin_idx]
        ; store all the data
	bin_str = {bin_num:bin, $
                   bin_teff:[teff_bin_vec_0[bin],teff_bin_vec_1[bin]], $
                   bin_idx:bin_idx, $
                   aspcap_idx:bin_aspcap_idx, $
                   teff:bin_teff, $
                   teff_err:bin_teff_err, $
                   vsini:bin_vsini}
    endif else begin
        ; store the null values
	bin_str = {bin_num:bin, $
                   bin_teff:[teff_bin_vec_0[bin],teff_bin_vec_1[bin]], $
                   bin_idx:!values.f_nan, $
                   aspcap_idx:!values.f_nan, $
                   teff:!values.f_nan, $
                   teff_err:!values.f_nan, $
                   vsini:!values.f_nan}
    endelse

    ; Replace the placeholder value in the hash with bin_str
    data[bin] = bin_str

endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; Loop through each bin and fit the Kaniadakis function 100 times
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

param_arr = dblarr(3, nbins, 100)
status_arr = lonarr(nbins, 100)
chi2_arr = dblarr(nbins, 100)

common amoeba_info, vsini_boot, vsini_xvec, d_vsini, nnondet, thresh, vsini_nondet, vis1
vis1=vis

thresh = 5d0

; make a vsini x vector for the pdf later
minvsini = 0d0
vsini_range = maxvsini+5d0-minvsini
d_vsini = 0.1d0
npts_xvec = fix(vsini_range/d_vsini)

vsini_xvec = maken(minvsini, maxvsini, npts_xvec)

for bin = 0, nbins-1 do begin

    ; Retrieve bin info
    bin_info = data[bin]
    
    vsini_base = bin_info.vsini
    
    nvel = n_elements(vsini_base) 


    ; Make vectors of vsini_detected and vsini_nondetected
    det_idx = where(vsini_base gt thresh, ndet, complement=nondet_idx, $
                    ncomplement=nnondet)
    
    vsini_det = vsini_base[det_idx]
    if nnondet gt 0 then vsini_nondet = vsini_base[nondet_idx]

    ; Loop 100 times
    for iter = 0, 99 do begin
        refit:
        ; Add error (of 2 km/s)
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
            param_arr[*,bin,iter] = result
            status_arr[bin,iter] = 1
        endif else if n_elements(result) eq 1 then begin
            goto, refit
            ;param_arr[*,bin,iter] = replicate(!values.d_nan,3)
            ;status_arr[bin,iter] = -1
            
        endif
        print, result
    endfor
    print, "Teff bin: " + strtrim(bin,2) + " finished."

endfor

save, /variables, filename='vsini_kan_3.sav'

plot_the_data:

restore, 'vsini_kan_3.sav'

; For some reason, the usual ps approach didn't work
; set_plot, 'ps'
; device, filename = 'vsini_teff_kan.ps'

; device, xs = 13, ys= 8, /inches
; loadct, 13


; Plot results
window,0,xs=1200,ys=900
; Get mean and stddev of vsini_0
vsini_arr = param_arr[0,*,*]
vsini_arr = reform(vsini_arr)
vsini_mean = mean(vsini_arr, dimension=2, /nan)
vsini_stddev = stddev(vsini_arr, dimension=2, /nan)

; fit a line
teff_fit_xvec = maken(2600d0,4000d0,3000)

teff_fit_coeff = poly_fit(teff_bin_vec_ctr, vsini_mean, 1, $
                          measure_errors=vsini_stddev, yfit=teff_fit)
teff_fit_yvec = poly(teff_fit_xvec, teff_fit_coeff)
co1_str = string(round(teff_fit_coeff[1]*1d3)/1d3, format='(D6.3)')
co0_str = string(round(teff_fit_coeff[0]*1d1)/1d1, format='(D4.1)')
;eq_string = 'Vsini = '+strtrim(teff_fit_coeff[1],2)+'*Teff + '+ $
  ;strtrim(teff_fit_coeff[0],2) + ' km/s'
; eq_string = 'Vsini = '+strtrim(co1_str,2)+'*Teff + '+ $
;   strtrim(co0_str,2) + ' km/s'
eq_string = 'Vsini = '+co1_str+'*Teff + '+ co0_str + ' km/s'

; cgplot, teff_bin_vec_ctr, vsini_mean, ps=6, xtit = "T_eff", ytit = "Vsini_0", $
;   tit=eq_string, charsize=2.5, /window
; cgplot, teff_fit_xvec, teff_fit_yvec, /overplot, /addcmd
; cgerrplot, teff_bin_vec_ctr, vsini_mean-vsini_stddev, vsini_mean+vsini_stddev, /addcmd
; cgcontrol, output='cgkan_test.eps', IM_width=1000

 plot, teff_bin_vec_ctr, vsini_mean, ps=6, xtit = "Teff", ytit = "Vsini_0", $
   tit=eq_string, charsize=2.5
 oplot, teff_fit_xvec, teff_fit_yvec
 oploterr, teff_bin_vec_ctr, vsini_mean, vsini_stddev

outstr = {teff_bin_ctr:teff_bin_vec_ctr, $
          vsini_mean:vsini_mean, $
          vsini_stddev:vsini_stddev}

mwrfits, outstr, 'vsini_temp_data.fits', /create

stop
end
















