; Put together Teff estimates from Mann and ASPCAP

pro teff_for_goldsample

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

dxstop


;;; Readin gold sample
gs = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/gold_sample_nohot.fits',0)
gs_idx = where(gs eq 1, ngs)

struct = {APOGEE_ID:aspcap[0].apogee_id, $
          TEFF_ASPCAP:aspcap[0].teff, $
          TEFF_MANN:teff_rel2_calib[0]}

outstruct = replicate(struct, ngs)

foreach gsindex, gs_idx, idx do begin

    ;Add aspcap Teff
    tasp = aspcap[gsindex].teff
    if tasp le 2000 or tasp gt 4200 then tasp = !values.d_nan

    ;Add Mann Teff
    mannidx = where(mann.rel2_idx eq gsindex, nmann)
    if nmann eq 1 then tmann = teff_rel2_calib[mannidx] $
      else tmann = !values.d_nan

    str = struct
    str.apogee_id = aspcap[gsindex].apogee_id
    str.teff_aspcap = tasp
    str.teff_mann = tmann

    outstruct[idx] = str

endforeach

mwrfits, outstruct, '/home/stgilhool/APOGEE/APOGEE_data/teff_for_goldsample.fits', /create

stop
    


end
