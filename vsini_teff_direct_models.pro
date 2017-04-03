; Fit a line to the vsini vs Teff trend for each model

pro vsini_teff_direct_models

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
;Load houdebine data as well
houdfile = '/home/stgilhool/APOGEE/vsini_literature/houdebine/houdebine_all.fits'
houdstr = mrdfits(houdfile,1)
htemp = houdstr.teff
hvsini = houdstr.vsini
hupper = houdstr.upper_limit

hupper_idx = where(hupper eq 1, nul, complement=hdet_idx, ncomplement=nhdet)

htempdet = htemp[hdet_idx]
hvsinidet = hvsini[hdet_idx]
htempul = htemp[hupper_idx]
hvsiniul = hvsini[hupper_idx]

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

;color_vec_orig = lindgen(nbins_filled+2)*(256/(nbins_filled+2))
nbins=8
color_vec_orig = lindgen(10)*(200/(10))

color_vec = byte(color_vec_orig[1:-2])
;color_vec = reverse(color_vec)

;plot, teff, vsini, yr= [0,30], xtit = 'Teff', ytit = 'Vsini',title =
;"ASPCAP Vsini vs. Teff with SNR", /nodata
p1 = plot(teff, vsini, xr=[2600,4000], yr= [0,30], xtit = 'Teff', ytit = 'Vsini',title = "ASPCAP Vsini vs. Teff color-coded by ASPCAP_CLASS", /nodata, /buffer)

model_type = ['GKd_a','GKd_b','GKd_c','GKd_d','Md_a','Md_b','Md_c','Md_d']

color_idx = 0
for bin=0, nbins-1 do begin
    ; Recover the data from the bin
    bin_idx = where(strtrim(model,2) eq model_type[bin])
    teff_bin = teff[bin_idx]
    vsini_bin = vsini[bin_idx]
    color_i = color_vec[color_idx]
    color_i_vec = replicate(color_i, n_elements(teff_bin)) 
        
        ;0.13437500      0.13281250      0.92187500      0.88164062

        ;oplot, teff_bin, vsini_bin, ps=8, color=color_vec[color_idx]
    p1 = plot(teff_bin, vsini_bin, rgb_table = rgb_tab, symbol='square', sym_size = 0.5, linestyle=6, $
              vert_colors=color_i_vec,/overplot, name=model_type[bin]) 
    ;position=[0.1, 0.16, 0.89, 0.91], 
    color_idx++

endfor
    

;  cb = colorbar(target=p1, orientation=0, position=[0.25, 0.05, 0.76, 0.08], $
;                tickvalues=color_vec_orig, tickname=['test','GKd_a','GKd_b','GKd_c','GKd_d','Md_a','Md_b','Md_c','Md_d','test'], minor=0)

leg = legend(font_size=8, linestyle=0, position=[0.9,0.9], transparency=50, shadow=0)
for i = 0, nbins-1 do begin
    print, leg[i].text_color
    newcolor = reform(rgb_tab[color_vec[i],*])
    leg[i].text_color = newcolor
    print, leg[i].text_color
    print, ''
endfor




;t1 = text(0.05, 0.06, 'Model_type')

; loadct, 1, rgb_table=rgb_tab

; p1 = plot(htempul, hvsiniul, rgb_table = rgb_tab, symbol='triangle', sym_size = 0.5, linestyle=6, $
;           vert_colors=color_i_vec,/overplot, position=[0.1, 0.16, 0.89, 0.91]) 
; p1 = plot(htempdet, hvsinidet, rgb_table = rgb_tab, symbol='square', sym_size = 0.5, linestyle=6, $
;           vert_colors=color_i_vec,/overplot, position=[0.1, 0.16, 0.89, 0.91]) 

 loadct, 0, rgb_table=rgb_tab
 p1 = plot([2000,5000], [5,5], rgb_table = rgb_tab, linestyle=2, $
           /overplot, position=[0.1, 0.16, 0.89, 0.91]) 
 p1 = plot([3500,3500], [0,100], rgb_table = rgb_tab, linestyle=2, $
           /overplot, position=[0.1, 0.16, 0.89, 0.91]) 

;snr_mag_scale = (max(snr)-min(snr))/7.
;snr_mag = (snr-min(snr))/snr_mag_scale
;snr_mag = bytscl(snr, top=99)
;test = scatterplot(teff, vsini, yr= [0,30], xtit = 'Teff', ytit = 'Vsini',title = "ASPCAP Vsini vs. Teff with SNR", magnitude=snr_mag, rgb_table=3)



; These guys are from Deshpande 2013
;oplot, teff_lit[[2,3,4,5,7,8,9]], lit_vsini[[2,3,4,5,7,8,9]], ps=5, color=99999, symsize=2
;device, /close

p1.save, "vsini_teff_model1.eps"
stop


stop

end
