; Plot number unreported teff vs snr

pro teff_unreported

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

; Take just GKd and Md models
mod_class = asp_info.aspcap_class
gk_idx = where(strmatch(mod_class,'GKd*') eq 1)
md_idx = where(strmatch(mod_class,'Md*') eq 1)
gkmd_idx = where(strmatch(mod_class,'GKd*') or strmatch(mod_class,'Md*'))


; Fit a line to each and to combined
teff_gk = teff[gk_idx]
vsini_gk = vsini[gk_idx]
snr_gk = snr[gk_idx]

teff_md = teff[md_idx]
vsini_md = vsini[md_idx]
snr_md = snr[md_idx]

teff_gkmd = teff[gkmd_idx]
vsini_gkmd = vsini[gkmd_idx]
snr_gkmd = snr[gkmd_idx]

;;; Get number unreported vsini per snr bin
; Set up bin vectors
snr_max = max(snr_all)
snr_min = min(snr_all)
snr_min_100 = floor(snr_min/100)*100
snr_max_100 = ceil(snr_max/100)*100
;snr_diff = snr_max-snr_min
;snr_diff_100 = snr_diff/100
snr_diff = snr_max_100-snr_min_100
snr_diff_100 = snr_diff/100 + 1

snr_vec = lindgen(snr_diff_100)*100+snr_min_100

snr_idx = value_locate(snr_vec, snr_all)

; Loop through each snr bin and find number of reported and unreported
; vsini's
ndet_vec = []
nnan_vec = []
for bin = 0, snr_diff_100-2 do begin
    ;print, bin
    ;print, snr_vec[bin], snr_vec[bin+1]
    vnan_bin_idx = where(snr_idx eq bin)
    teff_samp = teff_all[vnan_bin_idx]
    det_idx = where(teff_samp gt 0, ndet, complement=nodet_idx, ncomplement=nnodet)
    ndet_vec = [ndet_vec, ndet]
    nnan_vec = [nnan_vec, nnodet]
    ;print, snr_all[vnan_bin_idx]
    ;print, ndet
    ;print, nnodet
    ;stop
endfor

; Plot
snr_bin_ctr = snr_vec[0:-2] + (snr_vec[1:-1]-snr_vec[0:-2])/2
nspec_vec = ndet_vec + nnan_vec
nanfrac = double(nnan_vec)/double(nspec_vec)

!p.multi = [0,1,2]

window, 0, xs = 1400, ys=900

plot, snr_bin_ctr, nanfrac, ps=10, xtit = "SNR", ytit = "Fraction teff that are NaN", charsize = 2

plot, snr_bin_ctr, nspec_vec, ps=10, xtit = "SNR", ytit = "Number of teff measurements", charsize = 2, /ylog
oplot, snr_bin_ctr, nnan_vec, ps =10, color=200

!p.multi = 0
stop
    




end
