; Shows the effect of changing detection threshold for stars with Teff
; 3500-4000, for my vsini and ASPCAP vsini

pro vsini_threshtest_vfit_vs_asp, vis=vis

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

!p.multi=0


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


; plot, teff[nondet_sg], vsini[nondet_sg], xr = [2500,4000], ps=6, charsize=2, yr=[0,70], tit="ASPCAP:detection, SG:non-detection"
; oplot, teff[nondet_sg], vsini_sg[nondet_sg], ps=2, color=200

; plot, teff[det_both], vsini[det_both], xr = [2500,4000], ps=6, charsize=2, yr=[0,70], tit="Detections in both analyses"
; oplot, teff[det_both], vsini_sg[det_both], ps=2, color=200

; plot, teff[nondet_both], vsini[nondet_both], xr = [2500,4000], ps=6, charsize=2, yr=[0,70], tit="Non-detections in both analyses"
; oplot, teff[nondet_both], vsini_sg[nondet_both], ps=2, color=200

; res = vsini[det_both]-vsini_sg[det_both]
; plot, teff[det_both], res, ps=6, title="Residuals between detections in both analyses. RMSE="+strtrim(stddev(res),2), charsize = 2, yr=[-35,35]
;plot, vsini[det_both]-vsini_sg[det_both], ps=6, title=strtrim(stddev(vsini[det_both]-vsini_sg[det_both]),2), charsize=2

!p.multi=0






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bin by Teff
;binsize = 100
binsize = 500
;min_temp = 2600
min_temp = 3500

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


;;; Make histograms of each data set
sg_vel = data[0].vsini_sg
asp_vel = data[0].vsini

for thresh = 1, 10 do begin
    ; group detections and non-detections
    sgdet_idx = where(sg_vel gt thresh, ndetsg, complement=nondetsg_idx, $
                      ncomplement=nnondetsg)
    aspdet_idx = where(asp_vel gt thresh, ndetasp, complement=nondetasp_idx, $
                       ncomplement=nnondetasp)

    ; make preliminary histograms
    bs = 1
    sg_hist = histogram(sg_vel, binsize = bs, min=minvsini, max=maxvsini)
    asp_hist = histogram(asp_vel, binsize = bs, min=minvsini, max=maxvsini)
    
    ; replace the nondet bins with the average number of non-dets per unit vel
    ; (assuming threshold is 5)
    binmin = fix(minvsini)
    binmax = fix(thresh - minvsini - 1)
    sg_hist[binmin:binmax] = nnondetsg/float(binmax - binmin + 1)
    asp_hist[binmin:binmax] = nnondetasp/float(binmax - binmin + 1)

    distmax = (max(sg_hist) > max(asp_hist)) + 1

    nb = n_elements(asp_hist)
    xvec = lindgen(nb) + minvsini
    plot, xvec, asp_hist, ps=10, title = "Comparison with threshold: "+strtrim(thresh,2)+" SG: "+strtrim(nnondetsg,2)+" | ASP: "+strtrim(nnondetasp,2)+" | DIFF: "+strtrim(abs(nnondetsg-nnondetasp),2), yr=[0,distmax]
    oplot, xvec, sg_hist, ps=10, color=200
    dxstop
endfor

stop
end
