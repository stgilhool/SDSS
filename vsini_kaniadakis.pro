


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

function kaniadakis_fit, params, $
                        vsini_vec=vsini_vec, $
                        vsini_xvec=vsini_xvec, $
                        vsini_hist=vsini_hist

; Get params
vsini_0 = params[0]
kappa = params[1]
sigma = params[2]

; Make distribution
nbins = n_elements(vsini_hist) 
vsini_max = max(vsini_xvec)
y = dindgen(vsini_max)
yp = y-vsini_0

sqrt_term = sqrt(1d0 + (kappa^2)*(yp^4)/(sigma^4))
other_term = kappa * (-1d0 * (yp^2)/(sigma^2))
combined_term = (sqrt_term + other_term)^(1/kappa)

kan_dist = y * combined_term

; test other dist
;q_dist = y * (1d0 + (1d0 - kappa)*(-1d0*(y^2)/(sigma^2)))^(1d0/(1d0-kappa))

;kan_dist = q_dist

; Integrate from 0-5 in the kan dist
non_det_kan = total(kan_dist[0:4])

kan_dist_copy = kan_dist

kan_dist = kan_dist[5:-1]
kan_dist = [non_det_kan, kan_dist]

kan_dist = kan_dist[0:nbins]

; Normalize
;kan_dist_normed = kan_dist/total(kan_dist,/double)
nvsini = total(vsini_hist)

;kan_dist_final = kan_dist_normed * nvsini
kan_dist_final = kan_dist * nvsini

; test
kan_dist_final = kan_dist_final[1:-1]
vsini_hist_final = vsini_hist[1:-1]
vsini_xvec = vsini_xvec[1:-1]

;stop
; Calculate chi2
residual = kan_dist_final-vsini_hist_final

dev = residual^2

; Plot
plot, vsini_xvec, vsini_hist, ps=10
oplot, vsini_xvec, kan_dist_final, ps=10, color=200
;wait,0.05


return, dev

end











pro vsini_kaniadakis

; Read in the ASPCAP data
aspcap_file = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'
aspcap_info = mrdfits(aspcap_file,1)

; Make list of Teff b/w 2600 and 4000, and vsini gt 0
teff_all = aspcap_info.teff
teff_err_all = aspcap_info.teff_err
vsini_all = aspcap_info.vsini

sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0, nsel)
; take just our selection
sel_info = aspcap_info[sel_idx]
teff = sel_info.teff
teff_err = sel_info.teff_err
vsini = sel_info.vsini


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


for bin = 0, nbins-1 do begin

    ; Retrieve bin info
    bin_info = data[bin]
    
    vsini_base = bin_info.vsini
    
    nvel = n_elements(vsini_base) 

    ; Loop 100 times
    for iter = 0, 99 do begin

        ; Add error (of 2 km/s)
        std_error = 2d0
        vsini_err = randomn(seed, nvel)*std_error
        
        ; Make sure detections stay above 5km/s
        det_idx = where(vsini_base gt 5, ndet, complement=nondet_idx, $
                        ncomplement=nnondet)
        
        vsini_boot = vsini_base
        
        ;vsini_boot[det_idx] = (vsini_base[det_idx] + vsini_err[det_idx])
        ;> 5d0
        vsini_boot = (vsini_base + vsini_err)
        
        
        ; Make custom histogram
        vsini_xvec = lindgen(max(vsini_base)) + 5
        vsini_xvec = [0,vsini_xvec]
        vsini_locate_vec = value_locate(vsini_xvec, vsini_boot)
        vsini_hist = histogram(vsini_locate_vec, binsize = 1, min=0)
        
        
        ; Make functargs
        functargs = {vsini_vec:vsini_boot, $
                     vsini_xvec:vsini_xvec, $
                     vsini_hist:vsini_hist}
        
        ; Make parinfo
        parinfo_samp = {parname:'', $
                        value:0d0, $
                        limited:[0,0], $
                        limits:[0d0,0d0], $
                        fixed:0}
        parinfo = replicate(parinfo_samp, 3)
        
        parinfo[0].parname = 'vsini'
        parinfo[0].value = median(vsini_boot)
        parinfo[0].limited = [1,0]
        
        
        parinfo[1].parname = 'kappa'
        parinfo[1].value = 0.5d0
        parinfo[1].limited = [1,0]
        parinfo[1].limits = [0d0, 0d0]
        
        parinfo[2].parname = 'sigma'
        parinfo[2].value = 3.5d0
        parinfo[2].limited = [1,0]
        
        
        ; Fit
        result = mpfit('kaniadakis_fit', parinfo=parinfo, functargs=functargs, $
                       status=status, bestnorm=chi2)

        param_arr[*,bin,iter] = result
        status_arr[bin,iter] = status
        chi2_arr[bin,iter] = chi2

    endfor
    print, "Teff bin: " + strtrim(bin,2) + " finished."

endfor


; Plot results
; Get mean and stddev of vsini_0
vsini_arr = param_arr[0,*,*]
vsini_arr = reform(vsini_arr)
vsini_mean = mean(vsini_arr, dimension=2, /nan)
vsini_stddev = stddev(vsini_arr, dimension=2, /nan)

plot, teff_bin_vec_ctr, vsini_mean, ps=6, xtit = "Teff", ytit = "Vsini"
oploterr, teff_bin_vec_ctr, vsini_mean, vsini_stddev


stop
end
