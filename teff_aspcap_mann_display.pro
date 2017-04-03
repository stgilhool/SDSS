pro teff_aspcap_mann_display

; Show a comparison between our Teffs calculated from the Mann 2015?
; paper empirical relations with the new ASPCAP Teffs

; These Mann Teffs come from polynomials:
; Relation 1: 
; Teff/3500 = a + bx + cx^2 + dx^3 + ex^4 + f[Fe/H]
; Relation 2: 
; Teff/3500 = a + bx + cx^2 + dx^3 + ex^4 + fz + gz^2
; where x = (r-J), and z = (J-H)

; Read in Mann results
mann = mrdfits('/home/stgilhool/APOGEE/teff_results/teff_mann_empirical.fits',1)              

;help, mann
;** Structure <1c020a48>, 4 tags, length=7408, data length=7404, refs=1:
;   REL1_IDX        LONG      Array[57]
;   TEFF_REL1       DOUBLE    Array[57]
;   REL2_IDX        LONG      Array[560]
;   TEFF_REL2       DOUBLE    Array[560]

; Read in ASPCAP values
aspcap = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits',1)

; Make teff vectors for the two relations
teff_aspcap = aspcap.teff
; get rid of non-results in the aspcap vector
nan_idx = where(teff_aspcap lt 0, nnan)
teff_aspcap[nan_idx] = !values.d_nan
; make the two vectors
teff_aspcap1 = teff_aspcap[mann.rel1_idx]
teff_aspcap2 = teff_aspcap[mann.rel2_idx]
teff_rel1 = mann.teff_rel1
teff_rel2 = mann.teff_rel2

; Get rid of bad Mann Teff values
nan1_idx = where(teff_rel1 gt 6000, nnan1)
nan2_idx = where(teff_rel2 gt 6000, nnan2)
; replace bad ones (2 in rel1 and 11 in rel2)
if nnan1 gt 0 then teff_rel1[nan1_idx] = !values.d_nan
if nnan2 gt 0 then teff_rel2[nan2_idx] = !values.d_nan

; Check goodness of fit
res1 = teff_aspcap1 - teff_rel1
rms1 = stddev(res1, /nan)
res2 = teff_aspcap2 - teff_rel2
rms2 = stddev(res2, /nan)

; Plot the fit
; Vector for 1:1 comparison line
xxx = dindgen(6000)

; Plot relation 1
window, 0, xs = 1500, ys = 1000
plot, teff_aspcap1, teff_rel1, ps =6, /yno, title = $
  "Mann (r-J), [Fe/H] empirical relation | RMSE = "+strtrim(rms1,2), $
  xtitle = 'Teff_ASPCAP', ytitle = 'Teff_Mann', charsize = 2.0
oplot, xxx, xxx

dxstop

; Plot relation 2
window, 1, xs = 1500, ys = 1000
plot, teff_aspcap2, teff_rel2, ps =6, /yno, title = $
  "Mann (r-J), (J-H) empirical relation | RMSE = "+strtrim(rms2,2),$
  xtitle = 'Teff_ASPCAP', ytitle = 'Teff_Mann', charsize = 2.0
oplot, xxx, xxx

dxstop

; EXTRA: removing the offset in relation 2

high_idx = where(teff_rel2 gt teff_aspcap2, nhigh)
if nhigh gt 0 then begin
    offset_vec = teff_rel2[high_idx]-teff_aspcap2[high_idx]
    offset_mean = mean(offset_vec)
endif else message, "There definitely should be some high_idx"

teff_rel2_calib = teff_rel2 - offset_mean
res2calib = teff_rel2_calib - teff_aspcap2
rms2calib = stddev(res2calib, /nan)

; Re-Plot relation 2
window, 2, xs = 1500, ys = 1000
plot, teff_aspcap2, teff_rel2_calib, ps =6, /yno, title = $
  "Mann (r-J), (J-H) empirical relation (offset removed)| RMSE = "+strtrim(rms2calib,2),$
  xtitle = 'Teff_ASPCAP', ytitle = 'Teff_Mann', charsize = 2.0
oplot, xxx, xxx

stop

end
