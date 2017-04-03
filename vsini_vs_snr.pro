; Plot vsini vs snr, also plot teff vs snr

pro vsini_vs_snr

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


window, 0, xs = 1200, ys = 900
; Plot vsini vs snr
plot, snr_gkmd, vsini_gkmd, ps = 6, xtit = 'SNR', ytit = 'Vsini', charsize = 2, yr=[0,20], /ys
oplot, snr_gk, vsini_gk, ps = 6, color=200
oplot, snr_md, vsini_md, ps = 6, color=99999

stop
; Plot teff vs snr
plot, snr_gkmd, teff_gkmd, ps=6, yr=[2500, 4100], /ys, xtit = 'SNR', ytit = 'Teff', charsize = 2
oplot, snr_gk, teff_gk, ps=6, color=200
oplot, snr_md, teff_md, ps=6, color=99999


stop


; gkfit_co = robust_poly_fit(teff_gk, vsini_gk, 1, gkfit)
; mdfit_co = robust_poly_fit(teff_md, vsini_md, 1, mdfit)
; gkmdfit_co = robust_poly_fit(teff_gkmd, vsini_gkmd, 1, gkmdfit)

; plot, teff, vsini, ps=6
; oplot, teff_gk, vsini_gk, ps=6, color=200
; oplot, teff_gk, gkfit, color=200
; oplot, teff_md, vsini_md, ps=6, color=99999
; oplot, teff_md, mdfit, color=99999

; plot, teff, vsini, ps=6, yr= [0,30]

; lit_idx = [180,286,520,584,622,740,840,953,1038,1076,281,1132,143,904,958,485,622,1142,339,944,424,869,513,523,572,573,592,826,884,915,949]

; lit_vsini = [12.7,8.8,6.7,26.5,5.6,13.5,14.0,7.0,16.1,13.2,4.0,4.0,4.0,2.5,4.0,2.0,2.5,22.0,4.5,2.5,0.9,6.8,3.0,3.0,3.0,3.7,1.0,6.4,5.3,2.63,2.3]

; teff_lit = teff_all[lit_idx]

; oplot, teff_lit, lit_vsini, ps=6, color=200

stop

end
