function hist_conf, vector, binsize=binsize, binmin=binmin, conf_bounds=conf_bounds
; Make histogram to determine mode, and confidence limits

if n_elements(conf_bounds) eq 0 then conf_bounds = [0.25,0.75]

hist = histogram(vector, binsize=binsize, min=binmin, /nan)

xvec = dindgen(n_elements(hist))*binsize + binmin + (binsize/2d0)

; Get confidence interval
pdf = hist/total(hist,/double)
cdf = total(pdf, /cum)

conf = interpol(xvec, cdf, conf_bounds)

; Get mode (this is not robust)
smoothpdf = gauss_smooth(pdf, 3, /edge_wrap)

mode_idx = where(smoothpdf eq max(smoothpdf))
mode = mean(xvec[mode_idx])

; output
params = [mode, conf]

return, params

end

pro veq_vs_teff_smooth_display

restore, 'veq_realization_data.sav'
;info = mrdfits('veq_vs_teff_stepfit_results.fits',1)
info = mrdfits('veq_vs_teff_smoothfit_results.fits',1)

vbinsize = 0.05d0
tbinsize = 25d0
cbinsize = 0.005

vbinmin = 0d0
tbinmin = 2500d0
cbinmin = -10d0

; v1young = hist_conf(info.yresults[0,*], binsize=vbinsize, binmin=vbinmin)
; v2young = hist_conf(info.yresults[1,*], binsize=vbinsize, binmin=vbinmin)
; v1old = hist_conf(info.oresults[0,*], binsize=vbinsize, binmin=vbinmin)
; v2old = hist_conf(info.oresults[1,*], binsize=vbinsize, binmin=vbinmin)

; tcrit_young = hist_conf(info.yresults[2,*], binsize=tbinsize, binmin=tbinmin)
; tcrit_old = hist_conf(info.oresults[2,*], binsize=tbinsize, binmin=tbinmin)

a_young = hist_conf(info.yresults[0,*], binsize=vbinsize, binmin=vbinmin)
b_young = hist_conf(info.yresults[1,*], binsize=tbinsize, binmin=tbinmin)
c_young = hist_conf(info.yresults[2,*], binsize=cbinsize, binmin=cbinmin)
d_young = hist_conf(info.yresults[3,*], binsize=vbinsize, binmin=vbinmin)


a_old = hist_conf(info.oresults[0,*], binsize=vbinsize, binmin=vbinmin)
b_old = hist_conf(info.oresults[1,*], binsize=tbinsize, binmin=tbinmin)
c_old = hist_conf(info.oresults[2,*], binsize=cbinsize, binmin=cbinmin)
d_old = hist_conf(info.oresults[3,*], binsize=vbinsize, binmin=vbinmin)

ntstep = 1400
teff_vec = dindgen(ntstep) + 2600
teff_arr = rebin(reform(teff_vec, 1, ntstep), 3, ntstep)

ayoung = rebin(a_young, 3, ntstep)
byoung = rebin(b_young, 3, ntstep)
cyoung = rebin(c_young, 3, ntstep)
dyoung = rebin(d_young, 3, ntstep)

aold = rebin(a_old, 3, ntstep)
bold = rebin(b_old, 3, ntstep)
cold = rebin(c_old, 3, ntstep)
dold = rebin(d_old, 3, ntstep)

gomp_young = ayoung*exp(-1d0*exp(-1d0*cyoung*(teff_arr-byoung)))+dyoung
gomp_old = aold*exp(-1d0*exp(-1d0*cold*(teff_arr-bold)))+dold


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






info2 = mrdfits('veq_vs_teff_stepfit_tgrid_results2.fits',1)

vbinsize = 0.05d0
tbinsize = 25d0

vbinmin = 0d0
tbinmin = 2500d0

v1young = hist_conf(info2.best_veq_young_cool, binsize=vbinsize, binmin=vbinmin)
v2young = hist_conf(info2.best_veq_young_hot, binsize=vbinsize, binmin=vbinmin)
v1old = hist_conf(info2.best_veq_old_cool, binsize=vbinsize, binmin=vbinmin)
v2old = hist_conf(info2.best_veq_old_hot, binsize=vbinsize, binmin=vbinmin)

tcrit_young = hist_conf(info2.best_tcrit_young, binsize=tbinsize, binmin=tbinmin)
tcrit_old = hist_conf(info2.best_tcrit_old, binsize=tbinsize, binmin=tbinmin)


; Also, prepare the binned averages for overplotting
young_binned_avg = mean(info2.veq_bin_avg_young, dim=2, /nan)
old_binned_avg = mean(info2.veq_bin_avg_old, dim=2, /nan)
nbins_teff = n_elements(young_binned_avg) 

; Loop through and get the confidence interval for each point
young_binned_errbars = dblarr(2, nbins_teff)
old_binned_errbars = dblarr(2, nbins_teff)
young_binned_conf = dblarr(3, nbins_teff)
old_binned_conf = dblarr(3, nbins_teff)


for binnum = 0, nbins_teff-1 do begin

    young_avg = young_binned_avg[binnum]
    old_avg = old_binned_avg[binnum]
    
    young_binned_err = hist_conf(info2.veq_bin_avg_young[binnum], $
                                 binsize=vbinsize, binmin=vbinmin)
    old_binned_err = hist_conf(info2.veq_bin_avg_old[binnum], $
                               binsize=vbinsize, binmin=vbinmin)

    ; Save the errorbars
    young_binned_errbars[*,binnum] = abs(young_binned_err[1:2]-young_avg)
    old_binned_errbars[*,binnum] = abs(old_binned_err[1:2]-old_avg)
    
    ; Save the confidence interval info2 for comparison
    young_binned_conf[*,binnum] = young_binned_err
    old_binned_conf[*,binnum] = old_binned_err

endfor

teff_bin_vec = lindgen(nbins_teff)*100 + 2650








set_plot, 'ps'
device, filename = "veq_vs_teff_agemulti_smooth_clear.eps"
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
;mode
oplot, teff_vec, reform(gomp_young[0,*]), thick=2, color=400
;lower bound
oplot, teff_vec, reform(gomp_young[1,*]), linest=2, thick=2, color=400
;upper bound
oplot, teff_vec, reform(gomp_young[2,*]), linest=2, thick=2, color=400
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
;mode
oplot, teff_vec, reform(gomp_old[0,*]), thick=2, color=400
;lower bound
oplot, teff_vec, reform(gomp_old[1,*]), linest=2, thick=2, color=400
;upper bound
oplot, teff_vec, reform(gomp_old[2,*]), linest=2, thick=2, color=400
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
