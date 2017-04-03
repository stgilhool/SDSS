function hist_conf, vector, binsize=binsize, binmin=binmin
; Make histogram to determine mode, and confidence limits

hist = histogram(vector, binsize=binsize, min=binmin, /nan)

xvec = dindgen(n_elements(hist))*binsize + binmin + (binsize/2d0)

; Get confidence interval
pdf = hist/total(hist,/double)
cdf = total(pdf, /cum)

conf = interpol(xvec, cdf, [0.25,0.75])

; Get mode (this is not robust)
smoothpdf = gauss_smooth(pdf, 3, /edge_wrap)

mode_idx = where(smoothpdf eq max(smoothpdf))
mode = mean(xvec[mode_idx])

; output
params = [mode, conf]

return, params

end

pro veq_vs_teff_step_tgrid_display

restore, 'veq_realization_data.sav'
;info = mrdfits('veq_vs_teff_stepfit_results.fits',1)
info = mrdfits('veq_vs_teff_stepfit_tgrid_results2.fits',1)

vbinsize = 0.05d0
tbinsize = 25d0

vbinmin = 0d0
tbinmin = 2500d0


;v1young = hist_conf(info.yresults[0,*], binsize=vbinsize, binmin=vbinmin)
;v2young = hist_conf(info.yresults[1,*], binsize=vbinsize, binmin=vbinmin)
;v1old = hist_conf(info.oresults[0,*], binsize=vbinsize, binmin=vbinmin)
;v2old = hist_conf(info.oresults[1,*], binsize=vbinsize, binmin=vbinmin)
v1young = hist_conf(info.best_veq_young_cool, binsize=vbinsize, binmin=vbinmin)
v2young = hist_conf(info.best_veq_young_hot, binsize=vbinsize, binmin=vbinmin)
v1old = hist_conf(info.best_veq_old_cool, binsize=vbinsize, binmin=vbinmin)
v2old = hist_conf(info.best_veq_old_hot, binsize=vbinsize, binmin=vbinmin)

;tcrit_young = hist_conf(info.yresults[2,*], binsize=tbinsize, binmin=tbinmin)
;tcrit_old = hist_conf(info.oresults[2,*], binsize=tbinsize, binmin=tbinmin)
tcrit_young = hist_conf(info.best_tcrit_young, binsize=tbinsize, binmin=tbinmin)
tcrit_old = hist_conf(info.best_tcrit_old, binsize=tbinsize, binmin=tbinmin)


; Also, prepare the binned averages for overplotting
young_binned_avg = mean(info.veq_bin_avg_young, dim=2, /nan)
old_binned_avg = mean(info.veq_bin_avg_old, dim=2, /nan)
nbins_teff = n_elements(young_binned_avg) 
;young_binned_err = sqrt(total(info.veq_bin_stderr_young^2d0,2,/nan))
;old_binned_err = sqrt(total(info.veq_bin_stderr_old^2d0,2,/nan))

; Loop through and get the confidence interval for each point
young_binned_errbars = dblarr(2, nbins_teff)
old_binned_errbars = dblarr(2, nbins_teff)
young_binned_conf = dblarr(3, nbins_teff)
old_binned_conf = dblarr(3, nbins_teff)


for binnum = 0, nbins_teff-1 do begin

    young_avg = young_binned_avg[binnum]
    old_avg = old_binned_avg[binnum]
    
    young_binned_err = hist_conf(info.veq_bin_avg_young[binnum], $
                                 binsize=vbinsize, binmin=vbinmin)
    old_binned_err = hist_conf(info.veq_bin_avg_old[binnum], $
                               binsize=vbinsize, binmin=vbinmin)

    ; Save the errorbars
    young_binned_errbars[*,binnum] = abs(young_binned_err[1:2]-young_avg)
    old_binned_errbars[*,binnum] = abs(old_binned_err[1:2]-old_avg)
    
    ; Save the confidence interval info for comparison
    young_binned_conf[*,binnum] = young_binned_err
    old_binned_conf[*,binnum] = old_binned_err

endfor

teff_bin_vec = lindgen(nbins_teff)*100 + 2650
;;; This was just for visualization

; hist_ty = histogram(info.yresults[2,*], binsize=tbinsize, min=tbinmin)

; xvec_ty = dindgen(n_elements(hist_ty))*tbinsize + tbinmin + (tbinsize/2d0)
; xplotvec_ty = dindgen(n_elements(hist_ty))*tbinsize + tbinmin

; ; Get confidence interval
; pdf_ty = hist_ty/total(hist_ty,/double)


; hist_to = histogram(info.oresults[2,*], binsize=tbinsize, min=tbinmin)

; xvec_to = dindgen(n_elements(hist_to))*tbinsize + tbinmin + (tbinsize/2d0)
; xplotvec_to = dindgen(n_elements(hist_to))*tbinsize + tbinmin

; ; Get confidence interval
; pdf_to = hist_to/total(hist_to,/double)

; plot, xplotvec_ty, pdf_ty, ps=10, yr=[0, max(pdf_ty) > max(pdf_to)]
; oplot, xplotvec_to, pdf_to, ps=10, color=200


;;; Start plotting
; Filled circle definition
usym = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
USERSYM, COS(usym), SIN(usym), /FILL


qage = age_values[0,*]

yng_ndet_idx = where(qage le 2 and vsini le 5, nyndet)
yng_det_idx = where(qage le 2 and vsini gt 5, nydet)
old_ndet_idx = where(qage gt 2 and vsini le 5, nondet)
old_det_idx = where(qage gt 2 and vsini gt 5, nodet)




;  set_plot, 'ps'
;  device, filename = "veq_vs_teff_agemulti_stepfit_tgrid.eps"
;  device, /color, bits=8
;  device, xs=11,ys=8, /inches
;  loadct, 13


; !p.multi=[0,1,2]

; ;Young Detections
; plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Teff", $
;   ytit="V_eq", title = "Age < 2 Gyrs", charsize=1.5, $
;   xr=[2500,4000], yr=[3,200], /ylog, /ys, /nodata
; oplot, teff[yng_det_idx], veq_values[0,yng_det_idx], ps=8
; oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], replicate(70d0,nydet), veq_values[1,yng_det_idx], /lobar,/nohat, ps=3
; oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], replicate(70d0,nydet), veq_values[2,yng_det_idx], /hibar,/nohat, ps=3
; ;Young non-detections
; plotsym, 5, /fill ;Upside-down triangle
; oplot, teff[yng_ndet_idx], veq_values[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
; ; Add the fit
; ;lower v1
; oplot, [min(teff),tcrit_young[1]], replicate(v1young[1],2), linest=3, thick=2, color=400
; ;upper v1
; oplot, [min(teff),tcrit_young[2]], replicate(v1young[2],2), linest=3, thick=2, color=400
; ;mode v1
; oplot, [min(teff),tcrit_young[0]], replicate(v1young[0],2), linest=2, thick=2, color=400
; ;lower v2
; oplot, [tcrit_young[1],max(teff)], replicate(v2young[1],2), linest=3, thick=2, color=400
; ;upper v2
; oplot, [tcrit_young[2],max(teff)], replicate(v2young[2],2), linest=3, thick=2, color=400
; ;mode v2
; oplot, [tcrit_young[0],max(teff)], replicate(v2young[0],2), linest=2, thick=2, color=400
; ;lower tcrit
; oplot, replicate(tcrit_young[1],2), [v1young[1],v2young[1]], linest=3, thick=2, color=400
; ;upper tcrit
; oplot, replicate(tcrit_young[2],2), [v1young[2],v2young[2]], linest=3, thick=2, color=400
; ;mode tcrit
; oplot, replicate(tcrit_young[0],2), [v1young[0],v2young[0]], linest=2, thick=2, color=400
; ; Add the binned errors
; plotsym, 8, /fill
; oplot, teff_bin_vec, young_binned_avg, ps=8, color=200
; oploterror, teff_bin_vec, young_binned_avg, young_binned_err, ps=3, color=200
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles


; ;Old Detections
; plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Teff", $
;   ytit="V_eq", title="Age > 2 Gyrs", charsize = 1.5, $
;   xr=[2500,4000], yr=[3,200], /ylog, /ys, /nodata
; oplot, teff[old_det_idx], veq_values[0,old_det_idx], ps=8
; oploterror, teff[old_det_idx], veq_values[0,old_det_idx], replicate(70d0,nodet), veq_values[1,old_det_idx], /lobar,/nohat, ps=3
; oploterror, teff[old_det_idx], veq_values[0,old_det_idx], replicate(70d0,nodet), veq_values[2,old_det_idx], /hibar,/nohat, ps=3
; ;Old non-detections
; plotsym, 5, /fill ;Upside-down triangle
; oplot, teff[old_ndet_idx], veq_values[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
; ; Add the fit
; ;lower v1
; oplot, [min(teff),tcrit_old[1]], replicate(v1old[1],2), linest=3, thick=2, color=400
; ;upper v1
; oplot, [min(teff),tcrit_old[2]], replicate(v1old[2],2), linest=3, thick=2, color=400
; ;mode v1
; oplot, [min(teff),tcrit_old[0]], replicate(v1old[0],2), linest=2, thick=2, color=400
; ;lower v2
; oplot, [tcrit_old[1],max(teff)], replicate(v2old[1],2), linest=3, thick=2, color=400
; ;upper v2
; oplot, [tcrit_old[2],max(teff)], replicate(v2old[2],2), linest=3, thick=2, color=400
; ;mode v2
; oplot, [tcrit_old[0],max(teff)], replicate(v2old[0],2), linest=2, thick=2, color=400
; ;lower tcrit
; oplot, replicate(tcrit_old[1],2), [v1old[1],v2old[1]], linest=3, thick=2, color=400
; ;upper tcrit
; oplot, replicate(tcrit_old[2],2), [v1old[2],v2old[2]], linest=3, thick=2, color=400
; ;mode tcrit
; oplot, replicate(tcrit_old[0],2), [v1old[0],v2old[0]], linest=2, thick=2, color=400
; ; Add the binned errors
; plotsym, 8, /fill
; oplot, teff_bin_vec, old_binned_avg, ps=8, color=200
; oploterror, teff_bin_vec, old_binned_avg, old_binned_err, ps=3, color=200
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles


; !p.multi = 0

; device, /close_file


;  set_plot, 'ps'
;  device, filename = "veq_vs_feh.eps"
;  device, /color, bits=8
;  device, xs=11,ys=8, /inches
;  loadct, 13


; !p.multi=[0,1,2]

; ;Young Detections
; plot, feh[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="[Fe/H]", $
;   ytit="V_eq", title = "Age < 2 Gyrs", charsize=1.5, $
;   xr=[-1.5,1], yr=[3,200], /ylog, /ys, /xs, /nodata
; oplot, feh[yng_det_idx], veq_values[0,yng_det_idx], ps=8
; oploterror, feh[yng_det_idx], veq_values[0,yng_det_idx], replicate(0.5d0,nydet), veq_values[1,yng_det_idx], /lobar,/nohat, ps=3
; oploterror, feh[yng_det_idx], veq_values[0,yng_det_idx], replicate(0.5d0,nydet), veq_values[2,yng_det_idx], /hibar,/nohat, ps=3
; ;Young non-detections
; plotsym, 5, /fill ;Upside-down triangle
; oplot, feh[yng_ndet_idx], veq_values[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

; ;Old Detections
; plot, feh[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="[Fe/H]", $
;   ytit="V_eq", title="Age > 2 Gyrs", charsize = 1.5, $
;   xr=[-1.5,1], yr=[3,200], /ylog, /ys, /xs, /nodata
; oplot, feh[old_det_idx], veq_values[0,old_det_idx], ps=8
; oploterror, feh[old_det_idx], veq_values[0,old_det_idx], replicate(0.5d0,nodet), veq_values[1,old_det_idx], /lobar,/nohat, ps=3
; oploterror, feh[old_det_idx], veq_values[0,old_det_idx], replicate(0.5d0,nodet), veq_values[2,old_det_idx], /hibar,/nohat, ps=3
; ;Old non-detections
; plotsym, 5, /fill ;Upside-down triangle
; oplot, feh[old_ndet_idx], veq_values[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

; !p.multi=0

; device, /close_file


; set_plot, 'ps'
; device, filename = "veq_vs_teff_agemulti_clear2.eps"
; device, /color, bits=8
; device, xs=12,ys=10, /inches
; loadct, 13


; !p.multi=[0,1,2]

; errstring = textoidl('\sigma_{Teff}')

; ;Young Detections
; plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Effective Temperature (K)", $
;   ytit="Equatorial Velocity (km/s)", title = "Age < 2 Gyrs", charsize=1.5, $
;   xr=[2600,4000], yr=[4,100], /ys, /nodata, /ylog, /xs
; oplot, teff[yng_det_idx], veq_values[0,yng_det_idx], ps=8
; oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], veq_values[1,yng_det_idx], /lobar,/nohat, ps=3
; oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], veq_values[2,yng_det_idx], /hibar,/nohat, ps=3
; oploterror, [2700], [6], replicate(68d0,2), replicate(0,2), /nohat, ps=8
; xyouts, [2690], [6.5], errstring
; ;Young non-detections
; plotsym, 5, /fill ;Upside-down triangle
; oplot, teff[yng_ndet_idx], veq_values[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
; ; Add the fit
; ;lower v1
; oplot, [2500,tcrit_young[1]], replicate(v1young[1],2), linest=2, thick=2, color=400
; ;upper v1
; oplot, [2500,tcrit_young[2]], replicate(v1young[2],2), linest=2, thick=2, color=400
; ;mode v1
; oplot, [2500,tcrit_young[0]], replicate(v1young[0],2), thick=2, color=400
; ;lower v2
; oplot, [tcrit_young[1],max(teff)], replicate(v2young[1],2), linest=2, thick=2, color=400
; ;upper v2
; oplot, [tcrit_young[2],max(teff)], replicate(v2young[2],2), linest=2, thick=2, color=400
; ;mode v2
; oplot, [tcrit_young[0],max(teff)], replicate(v2young[0],2), thick=2, color=400
; ;lower tcrit
; oplot, replicate(tcrit_young[1],2), [v1young[1],v2young[1]], linest=2, thick=2, color=400
; ;upper tcrit
; oplot, replicate(tcrit_young[2],2), [v1young[2],v2young[2]], linest=2, thick=2, color=400
; ;mode tcrit
; oplot, replicate(tcrit_young[0],2), [v1young[0],v2young[0]], thick=2, color=400



; ;Old Detections
; plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Effective Temperature (K)", $
;   ytit="Equatorial Velocity (km/s)", title="Age > 2 Gyrs", charsize = 1.5, $
;   xr=[2600,4000], yr=[4,100], /ys, /nodata, /ylog, /xs
; oplot, teff[old_det_idx], veq_values[0,old_det_idx], ps=8
; oploterror, teff[old_det_idx], veq_values[0,old_det_idx], veq_values[1,old_det_idx], /lobar,/nohat, ps=3
; oploterror, teff[old_det_idx], veq_values[0,old_det_idx], veq_values[2,old_det_idx], /hibar,/nohat, ps=3
; oploterror, [2700], [6], replicate(68d0,2), replicate(0,2), /nohat, ps=8
; xyouts, [2690], [6.5], errstring
; ;Old non-detections
; plotsym, 5, /fill ;Upside-down triangle
; oplot, teff[old_ndet_idx], veq_values[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
; usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
; ; Add the fit
; ;lower v1
; oplot, [2500,tcrit_old[1]], replicate(v1old[1],2), linest=2, thick=2, color=400
; ;upper v1
; oplot, [2500,tcrit_old[2]], replicate(v1old[2],2), linest=2, thick=2, color=400
; ;mode v1
; oplot, [2500,tcrit_old[0]], replicate(v1old[0],2), thick=2, color=400
; ;lower v2
; oplot, [tcrit_old[1],max(teff)], replicate(v2old[1],2), linest=2, thick=2, color=400
; ;upper v2
; oplot, [tcrit_old[2],max(teff)], replicate(v2old[2],2), linest=2, thick=2, color=400
; ;mode v2
; oplot, [tcrit_old[0],max(teff)], replicate(v2old[0],2), thick=2, color=400
; ;lower tcrit
; oplot, replicate(tcrit_old[1],2), [v1old[1],v2old[1]], linest=2, thick=2, color=400
; ;upper tcrit
; oplot, replicate(tcrit_old[2],2), [v1old[2],v2old[2]], linest=2, thick=2, color=400
; ;mode tcrit
; oplot, replicate(tcrit_old[0],2), [v1old[0],v2old[0]], thick=2, color=400

; !p.multi = 0

; device, /close_file



set_plot, 'ps'
device, filename = "veq_vs_teff_agemulti_tgrid_clear.eps"
device, /color, bits=8
device, xs=12,ys=10, /inches
loadct, 13


!p.multi=[0,1,2]

errstring = textoidl('\sigma_{Teff}')

;Young Detections
plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Effective Temperature (K)", $
  ytit="Equatorial Velocity (km/s)", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2600,4000], yr=[4,100], /ys, /nodata, /ylog, /xs
oplot, teff[yng_det_idx], veq_values[0,yng_det_idx], ps=8
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], veq_values[1,yng_det_idx], /lobar,/nohat, ps=3
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], veq_values[2,yng_det_idx], /hibar,/nohat, ps=3
oploterror, [2700], [6], replicate(68d0,2), replicate(0,2), /nohat, ps=8
xyouts, [2690], [6.5], errstring
;Young non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[yng_ndet_idx], veq_values[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
; Add the fit
;lower v1
oplot, [2500,tcrit_young[1]], replicate(v1young[1],2), linest=2, thick=2, color=400
;upper v1
oplot, [2500,tcrit_young[2]], replicate(v1young[2],2), linest=2, thick=2, color=400
;mode v1
oplot, [2500,tcrit_young[0]], replicate(v1young[0],2), thick=2, color=400
;lower v2
oplot, [tcrit_young[1],max(teff)], replicate(v2young[1],2), linest=2, thick=2, color=400
;upper v2
oplot, [tcrit_young[2],max(teff)], replicate(v2young[2],2), linest=2, thick=2, color=400
;mode v2
oplot, [tcrit_young[0],max(teff)], replicate(v2young[0],2), thick=2, color=400
;lower tcrit
oplot, replicate(tcrit_young[1],2), [v1young[1],v2young[1]], linest=2, thick=2, color=400
;upper tcrit
oplot, replicate(tcrit_young[2],2), [v1young[2],v2young[2]], linest=2, thick=2, color=400
;mode tcrit
oplot, replicate(tcrit_young[0],2), [v1young[0],v2young[0]], thick=2, color=400
; Add the binned errors
plotsym, 8, /fill
oplot, teff_bin_vec, young_binned_avg, ps=8, color=100
oplot, teff_bin_vec, young_binned_conf[0,*], ps=8, color=200
oploterror, teff_bin_vec, young_binned_avg, young_binned_errbars[0,*], /lobar,$
  ps=3, color=100
oploterror, teff_bin_vec, young_binned_avg, young_binned_errbars[1,*], /hibar,$
  ps=3, color=100
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles



;Old Detections
plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Effective Temperature (K)", $
  ytit="Equatorial Velocity (km/s)", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2600,4000], yr=[4,100], /ys, /nodata, /ylog, /xs
oplot, teff[old_det_idx], veq_values[0,old_det_idx], ps=8
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], veq_values[1,old_det_idx], /lobar,/nohat, ps=3
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], veq_values[2,old_det_idx], /hibar,/nohat, ps=3
oploterror, [2700], [6], replicate(68d0,2), replicate(0,2), /nohat, ps=8
xyouts, [2690], [6.5], errstring
;Old non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[old_ndet_idx], veq_values[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
; Add the fit
;lower v1
oplot, [2500,tcrit_old[1]], replicate(v1old[1],2), linest=2, thick=2, color=400
;upper v1
oplot, [2500,tcrit_old[2]], replicate(v1old[2],2), linest=2, thick=2, color=400
;mode v1
oplot, [2500,tcrit_old[0]], replicate(v1old[0],2), thick=2, color=400
;lower v2
oplot, [tcrit_old[1],max(teff)], replicate(v2old[1],2), linest=2, thick=2, color=400
;upper v2
oplot, [tcrit_old[2],max(teff)], replicate(v2old[2],2), linest=2, thick=2, color=400
;mode v2
oplot, [tcrit_old[0],max(teff)], replicate(v2old[0],2), thick=2, color=400
;lower tcrit
oplot, replicate(tcrit_old[1],2), [v1old[1],v2old[1]], linest=2, thick=2, color=400
;upper tcrit
oplot, replicate(tcrit_old[2],2), [v1old[2],v2old[2]], linest=2, thick=2, color=400
;mode tcrit
oplot, replicate(tcrit_old[0],2), [v1old[0],v2old[0]], thick=2, color=400
; Add the binned errors
plotsym, 8, /fill
oplot, teff_bin_vec, old_binned_avg, ps=8, color=100
oplot, teff_bin_vec, old_binned_conf[0,*], ps=8, color=200
oploterror, teff_bin_vec, old_binned_avg, old_binned_errbars[0,*], /lobar,$
  ps=3, color=100
oploterror, teff_bin_vec, old_binned_avg, old_binned_errbars[1,*], /hibar,$
  ps=3, color=100
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles



!p.multi = 0

device, /close_file

stop

end
