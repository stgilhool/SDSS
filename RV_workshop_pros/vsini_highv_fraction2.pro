pro vsini_highv_fraction2, fast_rot_thresh, vis=vis

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
sgfile = '/home/stgilhool/APOGEE/master_info6.fits'
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
;load rainbow color table with discrete color tags	      
loadct,39,/silent
setcolors,/system_variables,decomposed=0,/silent  ;requires setcolors.pro

;plotting parameters to set axes and fonts
!p.background=0    ;default background color
!p.color=255
!p.charsize=1.7		;text default size
!p.charthick=6		;text default thickness
!x.thick=5		;thicken x-axis 
!y.thick=5		;thicken y-axis
!p.font=1		;set default font
!p.thick=5		;set default plotting line thickness
!x.margin=[7,4]

;set psym=8 to circle by default
circind=findgen(40) * (!pi*2/39.)
defsysv,'!circ',transpose( [[cos(circind)],[sin(circind)]])  ;user symbol vertices
usersym,!circ(0,*),!circ(1,*),/fill  ;circle plot


PSOPEN,'fast_rot_frac_new.eps',encapsulated=1,/INCHES,XSIZE=10,YSIZE=6,/COLOR

sg_frac = reform(sg_frac)


err_upper = sg_ul_68-sg_frac
err_lower = sg_frac-sg_ll_68

;teff_rev = reverse(teff_bin_vec_ctr)
teff_rev = teff_bin_vec_ctr

plot, teff_rev, sg_frac, ps=8, yr=[0,1], xr = [4000,2500], $
  xtit='Effective Temperature (K)', ytit='Fraction of Rotators', $
  charsize=2, charthick=2
oplot, [3200,3200], [0,1], linest=2
oploterror, teff_rev, sg_frac, err_upper, /hibar, ps=8
oploterror, teff_rev, sg_frac, err_lower, /lobar, ps=8

; plot, teff_rev, reverse(sg_frac), ps=6, yr=[0,1], xr = [4000,2500]
; oploterror, teff_rev, reverse(sg_frac), reverse(err_upper), /hibar
; oploterror, teff_rev, reverse(sg_frac), reverse(err_lower), /lobar




sptype = ['M0','M1','M2','M3','M4','M5','M6','M7']
sptype_temp = [3850,3700,3550,3400,3200,3050,2800,2650]
axis, xaxis=1, xticks=8, xtickv=sptype_temp, xtickn=sptype, charthick=6, charsize=2


psclose

; oploterror, teff_bin_vec_ctr, sg_frac, sg_ul_99-sg_frac, /hibar, ps=6
; oploterror, teff_bin_vec_ctr, sg_frac, sg_frac - sg_ll_99, /lobar, ps=6

; oploterror, teff_bin_vec_ctr, asp_frac, asp_ul_99-asp_frac, /hibar, errcolor='red', ps=1
; oploterror, teff_bin_vec_ctr, asp_frac, asp_frac-asp_ll_99, /lobar, errcolor='red', ps=1



stop




end
