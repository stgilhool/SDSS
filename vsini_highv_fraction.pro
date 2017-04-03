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







pro vsini_highv_fraction, fast_rot_thresh, vis=vis

if n_elements(vis) eq 0 then vis = 1

COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

; Readin the aspcap files
aspfile = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'
aspcap_info = mrdfits(aspfile,1)

aspcap_idx_all = lindgen(n_elements(aspcap_info)) 

; Make list of Teff b/w 2600 and 4000, and vsini gt 0
teff_all = aspcap_info.teff
teff_err_all = aspcap_info.teff_err
vsini_all = aspcap_info.vsini
snr_all = aspcap_info.snr

; Axe bad spectra
aspcapflag = aspcap_info.aspcapflag

;cut the aspcap VSINI_WARN and VSINI_BAD guys
vsini_warn_bit = 14
vsini_bad_bit = 30
star_bad_bit = 23

flag_bits = [vsini_warn_bit, star_bad_bit, vsini_bad_bit]

flag_dec = long(total(2L^flag_bits))

sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0 $
                and (aspcapflag and flag_dec) eq 0, nsel)

asp_info = aspcap_info[sel_idx]
aspcap_idx = aspcap_idx_all[sel_idx]
teff = teff_all[sel_idx]
teff_err = teff_err_all[sel_idx]
vsini = vsini_all[sel_idx]
snr = snr_all[sel_idx]
model = asp_info.aspcap_class

;;;;;;;;;;;;;;;;;;;
; Load my vfit results
sgfile = '/home/stgilhool/APOGEE/master_info5.fits'
sgstr = mrdfits(sgfile,1)
steff_all = sgstr.teff_vfit
svsini_all = sgstr.vsini_vfit

teff_sg = steff_all[sel_idx]
vsini_sg = svsini_all[sel_idx]

thresh = 5d0

nondet_asp = where(vsini le thresh and vsini_sg gt thresh, nnonasp)
nondet_sg = where(vsini_sg le thresh and vsini gt thresh, nnonsg)
nondet_both = where(vsini le thresh and vsini_sg le thresh, nnonboth)
det_both = where(vsini gt thresh and vsini_sg gt thresh, ndetboth)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bin by Teff
binsize = 100
min_temp = 2600

 
teff_hist = histogram(teff, binsize = binsize, min = min_temp, reverse_indices=ri)

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
        bin_vsini_sg = vsini_sg[bin_idx]
        ; store all the data
	bin_str = {bin_num:bin, $
                   bin_teff:[teff_bin_vec_0[bin],teff_bin_vec_1[bin]], $
                   bin_idx:bin_idx, $
                   aspcap_idx:bin_aspcap_idx, $
                   teff:bin_teff, $
                   teff_err:bin_teff_err, $
                   vsini:bin_vsini, $
                  vsini_sg:bin_vsini_sg}
    endif else begin
        ; store the null values
	bin_str = {bin_num:bin, $
                   bin_teff:[teff_bin_vec_0[bin],teff_bin_vec_1[bin]], $
                   bin_idx:!values.f_nan, $
                   aspcap_idx:!values.f_nan, $
                   teff:!values.f_nan, $
                   teff_err:!values.f_nan, $
                   vsini:!values.f_nan, $
                   vsini_sg:!values.f_nan}
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

common amoeba_info, vsini_boot, vsini_xvec, d_vsini, nnondet, thresh1, vsini_nondet, vis1
vis1=vis

thresh1 = thresh

; make a vsini x vector for the pdf later
minvsini = 0d0
maxvsini = max(vsini_sg)

vsini_range = maxvsini+5d0-minvsini
d_vsini = 0.1d0
npts_xvec = fix(vsini_range/d_vsini)

vsini_xvec = maken(minvsini, maxvsini, npts_xvec)


;fast_rot_thresh = 8
fastfrac_sg = dblarr(3, nbins)
fastfrac_asp = dblarr(3, nbins)

;;;; Make the plot showing the fraction of fast rotators per bin
for bin = 0, nbins-1 do begin

    ; recover bin data
    dstr = data[bin]

    vsg = dstr.vsini_sg
    vasp = dstr.vsini

    ; get number of fast rotators
    idx_sg = where(vsg ge fast_rot_thresh, nfast_sg)
    idx_asp = where(vasp ge fast_rot_thresh, nfast_asp)
    
    ntot_sg = n_elements(vsg) 
    ntot_asp = n_elements(vasp) 

    frac_sg = double(nfast_sg)/double(ntot_sg)
    frac_asp = double(nfast_asp)/double(ntot_asp)

    ; save
    fastfrac_sg[0,bin] = frac_sg
    fastfrac_asp[0,bin] = frac_asp
    fastfrac_sg[1,bin] = nfast_sg
    fastfrac_asp[1,bin] = nfast_asp
    fastfrac_sg[2,bin] = ntot_sg
    fastfrac_asp[2,bin] = ntot_asp

endfor

sg_frac = fastfrac_sg[0,*]
asp_frac = fastfrac_asp[0,*]

; Get confidence intervals
npts = n_elements(sg_frac) 
sg_ul_99 = dblarr(npts)
sg_ll_99 = dblarr(npts)
asp_ul_99 = dblarr(npts)
asp_ll_99 = dblarr(npts)
sg_ul_68 = dblarr(npts)
sg_ll_68 = dblarr(npts)
asp_ul_68 = dblarr(npts)
asp_ll_68 = dblarr(npts) 

; For each Teff bin, given the total number of stars, and the
; number of fast rotators, compute the 1, and 5 sigma limits
for dpt = 0, npts-1 do begin

    ; Totals
    bin_tot_sg = fastfrac_sg[2,dpt]
    bin_tot_asp = fastfrac_asp[2,dpt]
    ; Fast guys
    bin_fast_sg = fastfrac_sg[1,dpt]
    bin_fast_asp = fastfrac_asp[1,dpt]

    ; Limits
    limits_99_sg = binomial_stat_errors(bin_tot_sg, bin_fast_sg, 99d0)
    limits_68_sg = binomial_stat_errors(bin_tot_sg, bin_fast_sg, 68d0)
    limits_99_asp = binomial_stat_errors(bin_tot_asp, bin_fast_asp, 99d0)
    limits_68_asp = binomial_stat_errors(bin_tot_asp, bin_fast_asp, 68d0)

    ; Save limits
    sg_ul_99[dpt] = limits_99_sg[1]
    sg_ll_99[dpt] = limits_99_sg[0]
    sg_ul_68[dpt] = limits_68_sg[1]
    sg_ll_68[dpt] = limits_68_sg[0]
    asp_ul_99[dpt] = limits_99_asp[1]
    asp_ll_99[dpt] = limits_99_asp[0]
    asp_ul_68[dpt] = limits_68_asp[1]
    asp_ll_68[dpt] = limits_68_asp[0]


endfor


; cgDisplay, 600, 500, title = "Fractions of 'Fast' Rotators"

; cgPlot, teff_bin_vec_ctr, sg_frac, /NoData, YRange=[-0.5,1.5], YStyle=1

; cgColorFill, 




plot, teff_bin_vec_ctr, sg_frac, ps=6, yr=[0,1]
oplot, teff_bin_vec_ctr, asp_frac, ps=1, color=200

;oplot, teff_bin_vec_ctr, sg_ul_68, linest=2
;oplot, teff_bin_vec_ctr, sg_ll_68, linest=2

 oplot, teff_bin_vec_ctr, sg_ul_99, linest=3
 oplot, teff_bin_vec_ctr, sg_ll_99, linest=3

 oplot, teff_bin_vec_ctr, asp_ul_99, linest=3, color=200
 oplot, teff_bin_vec_ctr, asp_ll_99, linest=3, color=200


;;; New plot
set_plot, 'ps'
device, filename = 'fast_rot_frac.eps', /encapsulated

device, xs = 13, ys= 8, /inches
loadct, 13


err_upper = sg_ul_68-sg_frac
err_lower = sg_frac-sg_ll_68

;teff_rev = reverse(teff_bin_vec_ctr)
teff_rev = teff_bin_vec_ctr

plot, teff_rev, sg_frac, ps=6, yr=[0,1], xr = [4000,2500], xtit='T_eff', ytit='Fraction of Fast Rotators (Vsini > 5 km/s) per T_eff bin', charsize=1.5, charthick=2
oplot, [3200,3200], [0,1], linest=2
oploterror, teff_rev, sg_frac, err_upper, /hibar, ps=6
oploterror, teff_rev, sg_frac, err_lower, /lobar, ps=6

; plot, teff_rev, reverse(sg_frac), ps=6, yr=[0,1], xr = [4000,2500]
; oploterror, teff_rev, reverse(sg_frac), reverse(err_upper), /hibar
; oploterror, teff_rev, reverse(sg_frac), reverse(err_lower), /lobar




sptype = ['M0','M1','M2','M3','M4','M5','M6','M7']
sptype_temp = [3850,3700,3550,3400,3200,3050,2800,2650]
axis, xaxis=1, xticks=8, xtickv=sptype_temp, xtickn=sptype, charthick=2, charsize=1.5


device, /close_file

; oploterror, teff_bin_vec_ctr, sg_frac, sg_ul_99-sg_frac, /hibar, ps=6
; oploterror, teff_bin_vec_ctr, sg_frac, sg_frac - sg_ll_99, /lobar, ps=6

; oploterror, teff_bin_vec_ctr, asp_frac, asp_ul_99-asp_frac, /hibar, errcolor='red', ps=1
; oploterror, teff_bin_vec_ctr, asp_frac, asp_frac-asp_ll_99, /lobar, errcolor='red', ps=1



stop






















;;; Doing it the dumb way:
;;; Run the fit 100 times on my vsini, then literally copy the code
;;; and do it for ASP

;;; Also looping over threshold 5-10
param_arr_sg = dblarr(3,100,6)
param_arr_asp = dblarr(3,100,6)

for thresh_i = 5,10 do begin

    thresh = double(thresh_i)

    for bin = 0, 0 do begin
        
        ; Retrieve bin info
        bin_info = data[bin]
        
        vsini_base = bin_info.vsini_sg
        
        nvel = n_elements(vsini_base) 
        
        ; Make vectors of vsini_detected and vsini_nondetected
        det_idx = where(vsini_base gt thresh, ndet, complement=nondet_idx, $
                        ncomplement=nnondet)
        
        vsini_det = vsini_base[det_idx]
        if nnondet gt 0 then vsini_nondet = vsini_base[nondet_idx]
        
        ; Loop 100 times
        for iter = 0, 99 do begin
            refit1:
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
                goto, refit1
                ;param_arr[*,bin,iter] = replicate(!values.d_nan,3)
                ;status_arr[bin,iter] = -1
                
            endif
            print, result
        endfor
        print, "Teff bin: " + strtrim(bin,2) + " finished."
        
    endfor
    
    ; save my result before the next loop
    param_arr_temp = reform(param_arr, 3, 100)
    param_arr_sg[*,*,thresh_i-5] = param_arr_temp
    
    
;;; Copied loop with vsini_base changed to aspcap vsini
    for bin = 0, 0 do begin
        
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
            refit2:
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
                goto, refit2
                ;param_arr[*,bin,iter] = replicate(!values.d_nan,3)
                ;status_arr[bin,iter] = -1
                
            endif
            print, result
        endfor
        print, "Teff bin: " + strtrim(bin,2) + " finished."
        
    endfor
    
    ; save this result
    param_arr_temp = reform(param_arr, 3, 100)
    param_arr_asp[*,*,thresh_i-5] = param_arr_temp
    
endfor
    

param_avg_sg = mean(param_arr_sg, dimension=2)
param_avg_asp = mean(param_arr_asp, dimension=2)

param_sig_sg = stddev(param_arr_sg, dimension=2)
param_sig_asp = stddev(param_arr_asp, dimension=2)

help, param_avg_sg
print, param_avg_sg
print, ''
help, param_avg_asp
print, param_avg_asp


threshvec = lindgen(6)+5
!p.multi = [0,1,3]
v_sg = param_avg_sg[0,*]
v_asp = param_avg_asp[0,*]
k_sg = param_avg_sg[1,*]
k_asp = param_avg_asp[1,*]
s_sg = param_avg_sg[2,*]
s_asp = param_avg_asp[2,*]

vs_sg = param_sig_sg[0,*]
vs_asp = param_sig_asp[0,*]
ks_sg = param_sig_sg[1,*]
ks_asp = param_sig_asp[1,*]
ss_sg = param_sig_sg[2,*]
ss_asp = param_sig_asp[2,*]

plot, threshvec, v_sg, ps=6, xr = [4,11], /xs, yr=[0,5], title = "Best-fit parameters for VFIT (square) and ASPCAP (triangle) data vs. detection threshold", ytit='Vsini_0', charsize=2.0
oplot, threshvec, v_asp, ps=5
oploterror, threshvec, v_sg, vs_sg, /nohat
oploterror, threshvec, v_asp, vs_asp, /nohat

plot, threshvec, k_sg, ps=6, xr = [4,11], /xs, yr=[0.5,1.5], ytit='Kappa', charsize=2
oplot, threshvec, k_asp, ps=5
oploterror, threshvec, k_sg, ks_sg, /nohat
oploterror, threshvec, k_asp, ks_asp, /nohat

plot, threshvec, s_sg, ps=6, xr = [4,11], /xs, yr=[0,2], ytit='Sigma', xtit='Detection Threshold (km/s)', charsize=2
oplot, threshvec, s_asp, ps=5
oploterror, threshvec, s_sg, ss_sg, /nohat
oploterror, threshvec, s_asp, ss_asp, /nohat

stop


    ;save, /variables, filename='vsini_kan_vfit.sav'

plot_the_data:

;restore, 'vsini_kan_vfit.sav'

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
;snr_fit_xvec = maken(0d0,1000d0,4000)

teff_fit_coeff = poly_fit(teff_bin_vec_ctr, vsini_mean, 1, $
                          measure_errors=vsini_stddev, yfit=teff_fit)
teff_fit_yvec = poly(teff_fit_xvec, teff_fit_coeff)
; co1_str = string(round(teff_fit_coeff[1]*1d3)/1d3, format='(D6.3)')
; co0_str = string(round(teff_fit_coeff[0]*1d1)/1d1, format='(D4.1)')

; snr_fit_coeff = poly_fit(snr_bin_vec_ctr, vsini_mean, 1, $
;                          measure_errors=vsini_stddev, yfit=snr_fit)
;snr_fit_coeff = robust_poly_fit(snr_bin_vec_ctr, vsini_mean, 1, $
;                         snr_fit, /double)
;snr_fit_yvec = poly(snr_fit_xvec, snr_fit_coeff)
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

set_plot, 'ps'
device, filename = 'vsini_kaniadakis_vfit.eps', /encapsulated
device, /color, bits=8
device, xs = 13, ys= 8, /inches
loadct, 13


;  plot, snr_bin_vec_ctr, vsini_mean, ps=6, xtit = "snr", ytit = "Vsini_0", $
;    tit=eq_string, charsize=2.5
;  oplot, snr_fit_xvec, snr_fit_yvec
;  oploterr, snr_bin_vec_ctr, vsini_mean, vsini_stddev

 plot, teff_bin_vec_ctr, vsini_mean, ps=6, xtit = "SNR", ytit = "Peak Vsini", $
   tit = "Peak Vsini vs. Teff (VFIT)", charsize=2.5
 ;oplot, snr_fit_xvec, snr_fit_yvec
 oploterr, teff_bin_vec_ctr, vsini_mean, vsini_stddev

device, /close_file

stop

outstr = {teff_bin_ctr:teff_bin_vec_ctr, $
          vsini_mean:vsini_mean, $
          vsini_stddev:vsini_stddev}

;mwrfits, outstr, 'vsini_teff_vfit.fits', /create

stop


end


