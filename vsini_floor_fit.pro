; Fit a line to the vsini vs Teff trend for each model

pro vsini_floor_fit

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

; Take just GKd and Md models
mod_class = asp_info.aspcap_class
gk_idx = where(strmatch(mod_class,'GKd*') eq 1)
md_idx = where(strmatch(mod_class,'Md*') eq 1)
gkmd_idx = where(strmatch(mod_class,'GKd*') or strmatch(mod_class,'Md*'))


; Fit a line to each and to combined
teff_gk = teff[gk_idx]
vsini_gk = vsini[gk_idx]
teff_md = teff[md_idx]
vsini_md = vsini[md_idx]
teff_gkmd = teff[gkmd_idx]
vsini_gkmd = vsini[gkmd_idx]

gkfit_co = robust_poly_fit(teff_gk, vsini_gk, 1, gkfit)
mdfit_co = robust_poly_fit(teff_md, vsini_md, 1, mdfit)
gkmdfit_co = robust_poly_fit(teff_gkmd, vsini_gkmd, 1, gkmdfit)

plot, teff, vsini, ps=6
oplot, teff_gk, vsini_gk, ps=6, color=200
oplot, teff_gk, gkfit, color=200
oplot, teff_md, vsini_md, ps=6, color=99999
oplot, teff_md, mdfit, color=99999

;window, 0, xs=1200, ys=900

set_plot, 'ps'
device, filename = 'vsini_teff_lit4.eps', /encapsulated
device, /color, bits=8
device, xs = 13, ys= 8, /inches
loadct, 13


plotsym, 0, /fill

plot, teff, vsini, ps=6, symsize = 0.25, yr= [0,30], xtit = 'Teff', ytit = 'Vsini', charsize=2.0, title = "ASPCAP Vsini vs. Teff with literature comparison"

;lit_idx =[180,286,520,584,622,740,840,953,1038,1076,281,1132,143,904,958,485,622,1142,339,944,424,869,513,523,572,573,592,826,884,915,949]
lit_idx = [180,286,840,281,1132,143,904,958,485,622,1142,339,944,424,869,513,523,572,573,592,826,884,915,949,664,666,667,668,674,680,686, $
           694,695,698,700,702,707,710,711,712,715,719,729,732,737,738,742,746,747,750,756,759,869,944]

;lit_vsini =
;[12.7,8.8,6.7,26.5,5.6,13.5,14.0,7.0,16.1,13.2,4.0,4.0,4.0,2.5,4.0,2.0,2.5,22.0,4.5,2.5,0.9,6.8,3.0,3.0,3.0,3.7,1.0,6.4,5.3,2.63,2.3]
lit_vsini = [12.7,8.8,14.0,4.0,4.0,4.0,2.5,4.0,2.0,2.5,22.0,4.5,2.5,0.9,6.8,3.0,3.0,3.0,3.7,1.0,6.4,5.3,2.63,2.3, $
             4.5,2.8,3.0,1.8,1.1,9.0,1.4,6.8,20.3,2.8,20.9,0.8,17.7,4.7,26.1,4.3,5.0,1.3,16.5,6.4,6.8,1.0,3.5,3.2,2.3,1.9,3.5,2.7,6.8,2.5]


teff_lit = teff_all[lit_idx]

oplot, teff_lit, lit_vsini, ps=8, color=255

; These guys are from Deshpande 2013
;oplot, teff_lit[[2,3,4,5,7,8,9]], lit_vsini[[2,3,4,5,7,8,9]], ps=5, color=99999, symsize=2
device, /close_file


;;;;;;;;;;;
set_plot, 'ps'
;device, filename = 'vsini_teff_snr.eps', /encapsulated


device, /color, decomposed = 0    
;device, xs = 13, ys= 8, /inches
loadct, 3
loadct, 3, rgb_table=rgb_tab


gamma_ct, 0.4, /current

rgb_tab = [[R_curr], [G_curr], [B_curr]]

minsnr = 0
maxsnr = 1799
binsize = 100
snr_hist = histogram(snr, binsize=binsize, min=minsnr, max=maxsnr,reverse_indices=ri)
nbins = n_elements(snr_hist) 
null = where(snr_hist gt 0, nbins_filled)
snr_node_vec = maken(minsnr,maxsnr+1,nbins+1)
snr_ctr_vec = snr_node_vec[0:-2]+(binsize/2)

color_vec_orig = lindgen(nbins_filled+2)*(256/(nbins_filled+2))

color_vec = byte(color_vec_orig[1:-2])
;color_vec = reverse(color_vec)

;plot, teff, vsini, yr= [0,30], xtit = 'Teff', ytit = 'Vsini',title =
;"ASPCAP Vsini vs. Teff with SNR", /nodata
p1 = plot(teff, vsini, yr= [0,30], xtit = 'Teff', ytit = 'Vsini',title = "ASPCAP Vsini vs. Teff color-coded by SNR", /nodata, /buffer)

color_idx = 0
for bin=0, nbins-1 do begin
    ; Recover the data from the bin
    if ri[bin] ne ri[bin+1] then begin
        bin_idx = ri[ri[bin]:ri[bin+1]-1]
        teff_bin = teff[bin_idx]
        vsini_bin = vsini[bin_idx]
        color_i = color_vec[color_idx]
        color_i_vec = replicate(color_i, n_elements(teff_bin)) 
        
        ;0.13437500      0.13281250      0.92187500      0.88164062

        ;oplot, teff_bin, vsini_bin, ps=8, color=color_vec[color_idx]
        p1 = plot(teff_bin, vsini_bin, rgb_table = rgb_tab, symbol='square', sym_size = 0.5, linestyle=6, $
                  vert_colors=color_i_vec,/overplot, position=[0.1, 0.16, 0.89, 0.91])
        
        color_idx++
    endif
endfor
    

cb = colorbar(target=p1, orientation=0, position=[0.25, 0.05, 0.76, 0.08], $
              tickvalues=color_vec_orig, tickname=['0','100','200','300','400','500','600','700','800','900','1000','1100','1200','1400','1500','1800','1900'], $
              text_orientation=90)

;snr_mag_scale = (max(snr)-min(snr))/7.
;snr_mag = (snr-min(snr))/snr_mag_scale
;snr_mag = bytscl(snr, top=99)
;test = scatterplot(teff, vsini, yr= [0,30], xtit = 'Teff', ytit = 'Vsini',title = "ASPCAP Vsini vs. Teff with SNR", magnitude=snr_mag, rgb_table=3)



; These guys are from Deshpande 2013
;oplot, teff_lit[[2,3,4,5,7,8,9]], lit_vsini[[2,3,4,5,7,8,9]], ps=5, color=99999, symsize=2
;device, /close

p1.save, "vsini_teff_snr.eps"
stop


stop

end
