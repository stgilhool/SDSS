pro mearth_sini_test

; Readin results from mearth_period_comp
mearth_file = '/home/stgilhool/APOGEE/vsini_literature/mearth_period_comp.fits'
m_info = mrdfits(mearth_file,1)

; Estimate radii for the stars
	; Read in mean stellar parameters from table
;mean_param_file =
;'/home/stgilhool/APOGEE/misc/EEM_dwarf_UBVIJHK_colors_Teff.txt'
;mean_param_file = '/home/stgilhool/APOGEE/misc/EEM_just_table.txt'
mean_param_file = '/home/stgilhool/APOGEE/misc/EEM_just_what_i_want.txt'

; readcol, mean_param_file, spt, teff, logt, BC_V, M_V, logL, B_V, Bt_Vt, U_B, $
;   V_Rc, V_Ic, V_Ks, J_H, H_K, Ks_W1, M_sun, log_age, b_y, spt1, M_J, M_Ks, Mbol, $
;   i_z, z_Y, W1_W2, $
;   format = 'A,L,X,X,X,F,X,X,X,X,X,X,X,X,X,X,X,X,A,X,X,X,X,X,X', skipline=1
;   ;format = 'A,L,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,A,F,F,F,F,F,F', skipline=21
readcol, mean_param_file, spt, teff, logl, format='A,F,F', skipline=1

	; Interpolate to get Luminosities for the temperatures
luminosity = 10d0^logl

myteff = m_info.teff

est_luminosity = interpol(luminosity, teff, myteff)

	; Calculate radii
teff_solar = 5771.8
radius = ((teff_solar/myteff)^2) * sqrt(est_luminosity)

; Use the radii to get v_eq for the newton stars
day_sec = 24d0*3600d0
radius_solar = 695660

v_eq = 2d0*!pi*radius_solar*radius/(m_info.period*day_sec)

vrot = m_info.vrot

vsini = m_info.vsini

type = m_info.type
; Get sini for all the stars
sini = vsini/v_eq

;pass_idx = where(sini gt 0 and sini le 1 and vsini gt 0, npass)
npass = n_elements(type) 
pass_idx = lindgen(npass)

;gd_idx = where(vsini gt 0)
apg_idx = m_info.apg_idx[pass_idx]
nwt_idx = m_info.nwt_idx[pass_idx]
objid = m_info.objid[pass_idx]
period = m_info.period[pass_idx]
teff = m_info.teff[pass_idx]


 vsini=vsini[pass_idx]
 v_eq = v_eq[pass_idx]
 vrot = vrot[pass_idx]
 type = type[pass_idx]


; Get alternative vrot and vsini for some of these guys
alt_vrot = vrot
alt_vsini = vsini

alt1 = where(objid eq '2M09301445+2630250', nullcheck)
if nullcheck eq 0 then message, 'typo alert' else begin
    alt_vrot[alt1] = 22.3d0
    alt_vsini[alt1] = 6.7d0
endelse

alt2 = where(objid eq '2M11474074+0015201', nullcheck)
if nullcheck eq 0 then message, 'typo alert' else begin
    alt_vrot[alt2] = 16.1d0
    alt_vsini[alt2] = 5.6d0
endelse

alt3 = where(objid eq '2M18452147+0711584', nullcheck)
if nullcheck eq 0 then message, 'typo alert' else begin
    alt_vrot[alt3] = 17.2d0
    alt_vsini[alt3] = 16.1d0
endelse

alt4 = where(objid eq '2M19173151+2833147', nullcheck)
if nullcheck eq 0 then message, 'typo alert' else begin
    alt_vrot[alt4] = 16.2d0
    alt_vsini[alt4] = 13.2d0
endelse

alt5 = where(objid eq '2M11005043+1204108', nullcheck)
if nullcheck eq 0 then message, 'typo alert' else begin
    alt_vsini[alt5] = 26.5d0
endelse

alt6 = where(objid eq '2M12265737+2700536', nullcheck)
if nullcheck eq 0 then message, 'typo alert' else begin
    alt_vsini[alt6] = 13.5d0
endelse

alt7 = where(objid eq '2M16370146+3535456', nullcheck)
if nullcheck eq 0 then message, 'typo alert' else begin
    alt_vsini[alt7] = 7.0d0
endelse

; Consider erros and calculate a sini with errors

; Read in the ASPCAP data
aspcap_file = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'
aspcap_info = mrdfits(aspcap_file,1)

aspcapflag = aspcap_info.aspcapflag
aspcapflag = aspcapflag[apg_idx]

;cut the aspcap VSINI_WARN and VSINI_BAD guys
vsini_warn_bit = 14
vsini_bad_bit = 30
star_bad_bit = 23
flag_bits = [vsini_warn_bit, vsini_bad_bit]
;flag_bits = [vsini_warn_bit, vsini_bad_bit, star_bad_bit]
flag_dec = long(total(2L^flag_bits))

;flag_idx = where((aspcapflag and flag_dec) ne 0, nflg, complement=unflag_idx, ncomplement=nunflg)


;det_idx = where(type eq 'A' or type eq 'B', ndet)
;print, ndet
det_idx = where(type eq 'A' or type eq 'B' and (aspcapflag and flag_dec) eq 0, ndet)
;print, ndet
;stop

vr = vrot[det_idx]
e_vr = 0.14d0*vr
avr = alt_vrot[det_idx]
e_avr = 0.14d0 * avr

vs = vsini[det_idx]
avs = alt_vsini[det_idx]
e_vs = replicate(3d0, ndet)

si = vs/avr

e_si = abs(si)*sqrt((e_vs/vs)^2 + (e_avr/avr)^2)

; Take just the detections
drot_idx = where(vs gt 8d0, ndrot, complement=ndrot_idx, ncomplement=nndrot)

vrdet = vr[drot_idx]
e_vrdet = e_vr[drot_idx]
avrdet = avr[drot_idx]
e_avrdet = e_avr[drot_idx]
vsdet = vs[drot_idx]
avsdet = avs[drot_idx]
e_vsdet = e_vs[drot_idx]

sidet = si[drot_idx]
e_sidet = e_si[drot_idx]

; Treat non-detections
vrndet = vr[ndrot_idx]
e_vrndet = e_vr[ndrot_idx]
avrndet = avr[ndrot_idx]
e_avrndet = e_avr[ndrot_idx]
vsndet = vs[ndrot_idx]
vsndet_ul = vsndet
vsndet_ll = replicate(0d0, nndrot)

avsndet = avs[ndrot_idx]
e_vsndet = e_vs[ndrot_idx]

sindet = si[ndrot_idx]
e_sindet = e_si[ndrot_idx]
sindet_ul = sindet
sindet_ll = replicate(0d0, nndrot)

; Plot our vsini against vrot (taking care for upperlimits)
window, 0, xs = 1400, ys=1000
plot, vr, vs, xtit = 'Newton V_eq', ytit = 'Vsini', xr=[0,75], yr=[0,75], /xs, /ys, charsize=2, /nodata
oplot, dindgen(100), dindgen(100), linest=2
;oplot, avr, vs, color=200, ps=6
oplot, vrdet, vsdet, ps=6
oploterror, vrdet, vsdet, e_vrdet, e_vsdet
oplot, vrndet, vsndet, ps=6
oploterror, vrndet, vsndet, e_vrndet, vsndet_ul-vsndet, /hibar
oploterror, vrndet, vsndet, e_vrndet, vsndet_ll-vsndet, /lobar

stop

; Plot our vsini against vrot
plot, vr, vs, ps=6, xtit = 'Newton V_eq', ytit = 'Vsini', xr=[0,100], yr=[0,100], /xs, /ys, charsize=2
oplot, dindgen(100), dindgen(100), linest=2
oplot, avr, vs, color=200, ps=6
oplot, vr, vs, ps=6

stop
; Plot sini 
vsort = sort(vr)
;plot, vsort, si[vsort], ps = 6, yr=[0,2], xr = [-1, 29], /xs
plot, lindgen(n_elements(si)), si, ps = 6, xr = [-1, 29], /xs
;vdetsort = sort(vrdet)
;oploterror, vdetsort, sidet[vdetsort], e_sidet[vdetsort]
oploterror, drot_idx, sidet, e_sidet
oplot, [-1, 30], [1d0, 1d0], linest=2
;vndetsort = sort[vrndet]
;oploterror, vndetsort, sindet[vndetsort],
;sindet_ll[vndetsort]-sindet[vndetsort], /lobar
oploterror, ndrot_idx, sindet, sindet_ll-sindet, /lobar

stop

; Plot sini 
vsort = sort(vr)
;plot, vsort, si[vsort], ps = 6, yr=[0,2], xr = [-1, 29], /xs
plot, vr, si, ps = 6, xr = [-1, 75], /xs, yr=[0,2]
;vdetsort = sort(vrdet)
;oploterror, vdetsort, sidet[vdetsort], e_sidet[vdetsort]
oploterror, vrdet, sidet, e_sidet
oplot, [-1, 75], [1d0, 1d0], linest=2
;vndetsort = sort[vrndet]
;oploterror, vndetsort, sindet[vndetsort],
;sindet_ll[vndetsort]-sindet[vndetsort], /lobar
oploterror, vrndet, sindet, sindet_ll-sindet, /lobar

stop

; Make and output the structure

tagnames = tag_names(m_info)

tagnames = [tagnames, 'v_eq','sini']

minfo_out = create_struct(tagnames, apg_idx, nwt_idx, objid, type, period, vrot, vsini, teff, v_eq, si)

mwrfits, minfo_out, '/home/stgilhool/APOGEE/misc/sini_data2.fits', /create



print, npass

nwt_det_idx = where(type eq 'A' or type eq 'B', ndet)




stop


; Test if the distribution matches random i

end
