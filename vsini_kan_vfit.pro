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
    wait, 0.0001
endif

; Take log likelihood
logl_vec = alog(bigl_vec)
; Sum the vector of log likelihoods
logl = total(logl_vec)

if vsini_0 lt 0 then v_penalty = abs(logl)*0.5d0*exp(-1d0*vsini_0) $
else v_penalty=0
if kappa lt 0 then klow_penalty = abs(logl)*0.5d0*exp(-1d0*kappa) $
else klow_penalty=0
if kappa gt 1 then khigh_penalty = abs(logl)*0.5d0*exp(kappa-1d0) $
else khigh_penalty=0
if sigma lt 0 then s_penalty = abs(logl)*0.5d0*exp(-1d0*sigma) $
else s_penalty=0

return, -1d0*logl+v_penalty+klow_penalty+khigh_penalty+s_penalty

end







pro vsini_kan_vfit, vis=vis

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

; !p.multi=0

; ;window, 0, xs=1550, ys=1000
; set_plot, 'ps'
; device, filename = 'ASPCAP_PHOENIX_vsini.eps', /encapsulated

; device, xs = 9, ys=8, /inches
; device, /color, decomposed = 0


; plot, vsini, vsini_sg, ps=6, symsize=0.5, xr=[0,65], yr=[0,65], /xs, /ys, tit="Comparison of ASPCAP Vsini and PHOENIX-based Vsini", $
;   xtit = "Vsini_ASPCAP (km/s)", ytit="Vsini_PHOENIX (km/s)"
; oplot, 5+lindgen(70), 5+lindgen(70), linest=3
; oplot, lindgen(70), replicate(5d0, 70), linest=2
; oplot, [5,5], [0,70], linest=2

; device, /close_file

; stop

; set_plot, 'x'

; !p.multi=[0,2,3]
; plot, teff, vsini, xr = [2500,4000], ps=6, charsize=2, yr=[0,70], tit="ASPCAP and SG vsini values"
; oplot, teff, vsini_sg, ps=2, color=200
; plot, teff, vsini, xr = [2500,4000], ps=6, charsize=2, yr=[0,20], tit="ASPCAP and SG vsini values (zoom)"
; oplot, teff, vsini_sg, ps=2, color=200



; ;plot, [teff[nondet_asp]], [vsini[nondet_asp]], xr = [2500,4000], ps=6, charsize=2
; ;oplot, [teff[nondet_asp]], [vsini_sg[nondet_asp]], ps=6, color=200

; plot, teff[nondet_sg], vsini[nondet_sg], xr = [2500,4000], ps=6, charsize=2, yr=[0,70], tit="ASPCAP:detection, SG:non-detection"
; oplot, teff[nondet_sg], vsini_sg[nondet_sg], ps=2, color=200

; plot, teff[det_both], vsini[det_both], xr = [2500,4000], ps=6, charsize=2, yr=[0,70], tit="Detections in both analyses"
; oplot, teff[det_both], vsini_sg[det_both], ps=2, color=200

; plot, teff[nondet_both], vsini[nondet_both], xr = [2500,4000], ps=6, charsize=2, yr=[0,70], tit="Non-detections in both analyses"
; oplot, teff[nondet_both], vsini_sg[nondet_both], ps=2, color=200

; res = vsini[det_both]-vsini_sg[det_both]
; plot, teff[det_both], res, ps=6, title="Residuals between detections in both analyses. RMSE="+strtrim(stddev(res),2), charsize = 2, yr=[-35,35]
; ;plot, vsini[det_both]-vsini_sg[det_both], ps=6, title=strtrim(stddev(vsini[det_both]-vsini_sg[det_both]),2), charsize=2

; !p.multi=0

; stop




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

for bin = 0, nbins-1 do begin
;for bin = 12, 12 do begin
if bin eq 10 then continue
if bin eq 11 then continue
if bin eq 13 then continue
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
        
        ;boot_okay = 1
        
        ;while boot_okay do begin
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
            
            refit:
            vsini_guess = randomu(seed,1)*median(vsini_boot)
            kappa_guess = randomu(seed,1)*0.5d0
            sigma_guess = randomu(seed,1)*4d0
            
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

save, /variables, filename='vsini_kan_vfit2.sav'

plot_the_data:

restore, 'vsini_kan_vfit2.sav'

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
device, filename = 'vsini_kaniadakis_vfit2.eps', /encapsulated
device, /color, bits=8
device, xs = 13, ys= 8, /inches
loadct, 13


;  plot, snr_bin_vec_ctr, vsini_mean, ps=6, xtit = "snr", ytit = "Vsini_0", $
;    tit=eq_string, charsize=2.5
;  oplot, snr_fit_xvec, snr_fit_yvec
;  oploterr, snr_bin_vec_ctr, vsini_mean, vsini_stddev

 plot, teff_bin_vec_ctr, vsini_mean, ps=6, xtit = "T_eff (K)", ytit = "Peak Vsini (km/s)", $
   tit = "Kaniadakis Distribtion Peak Vsini vs. T_eff", charsize=2, charthick=2
 ;oplot, snr_fit_xvec, snr_fit_yvec
 oploterr, teff_bin_vec_ctr, vsini_mean, vsini_stddev

device, /close_file

stop

outstr = {teff_bin_ctr:teff_bin_vec_ctr, $
          vsini_mean:vsini_mean, $
          vsini_stddev:vsini_stddev}

mwrfits, outstr, 'vsini_teff_vfit2.fits', /create

stop


end
