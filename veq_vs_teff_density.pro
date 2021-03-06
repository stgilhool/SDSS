function hist_conf_ul, vector, binsize=binsize, binmin=binmin, conf_limit=conf_limit
; Make histogram to determine mode, and confidence limits for upper
; limit cases

if n_elements(conf_limit) eq 0 then conf_limit = 0.9

hist = histogram(vector, binsize=binsize, min=binmin, /nan)

xvec = dindgen(n_elements(hist))*binsize + binmin + (binsize/2d0)

; Get confidence interval
pdf = hist/total(hist,/double)
cdf = total(pdf, /cum)

conf = interpol(xvec, cdf, conf_limit)

conf = [0d0, conf]

; Get mode (this is not robust)
smoothpdf = gauss_smooth(pdf, 3, /edge_wrap)

mode_idx = where(smoothpdf eq max(smoothpdf))
mode = mean(xvec[mode_idx])

; output
params = [mode, conf]

return, params

end


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

pro veq_vs_teff_density

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
;USERSYM, COS(usym), SIN(usym), /FILL, psize=0.1
plotsym, 0, 0.1, /fill

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

teff_yreal_all = teff_realizations[*,young_star_idx]
veq_yreal_all = veq_realizations[*,young_star_idx]
nyoung_tot = nyoung*nsamples
teff_yreal_vec = reform(teff_yreal_all, nyoung_tot)
veq_yreal_vec = reform(veq_yreal_all, nyoung_tot)

;teff_bvec = lindgen(nbins_teff+1)*100+2600
tbinsize=50

ty_hist = histogram(teff_yreal_all, min=tbinmin, binsize=tbinsize, max=3999.999, reverse_indices=ri_ty)


teff_oreal_all = teff_realizations[*,old_star_idx]
veq_oreal_all = veq_realizations[*,old_star_idx]
nold_tot = nold*nsamples
teff_oreal_vec = reform(teff_oreal_all, nold_tot)
veq_oreal_vec = reform(veq_oreal_all, nold_tot)

to_hist = histogram(teff_oreal_all, min=tbinmin, binsize=tbinsize, max=3999.999, reverse_indices=ri_to)




ntbins=n_elements(ty_hist) 

tbin_vec = long(lindgen(ntbins+1)*tbinsize+tbinmin)
nvel_elem = n_elements(veq_vec) 
young_psf_arr = dblarr(nvel_elem, ntbins)
old_psf_arr = dblarr(nvel_elem, ntbins)


veq_vec_2 = dindgen(1001)


; !p.multi = [0,5,6]

;  set_plot, 'ps'
;  device, filename = "veq_vs_teff_pdf_estimate.eps"
;  device, /color, bits=8
;  device, xs=12,ys=10, /inches
;  loadct,13

set_plot, 'x'
!p.multi=[0,1,2]
window, 0, xs=1600, ys=900

yconf_bounds = dblarr(3, ntbins)
oconf_bounds = dblarr(3, ntbins)
yconf_bounds2 = dblarr(3, ntbins)
oconf_bounds2 = dblarr(3, ntbins)
yconf_bounds3 = dblarr(3, ntbins)
oconf_bounds3 = dblarr(3, ntbins)
yconf_boundsp = dblarr(2, ntbins)
oconf_boundsp = dblarr(2, ntbins)

for tbin=0, ntbins-1 do begin

    ; Get the indices for the temperature bin
    ybin_idx = ri_ty[ri_ty[tbin]:ri_ty[tbin+1]-1]
    obin_idx = ri_to[ri_to[tbin]:ri_to[tbin+1]-1]

    ; Get the velocities
    yveq = veq_yreal_vec[ybin_idx]
    oveq = veq_oreal_vec[obin_idx]

    ; Histogram those bitches up
    yveq_hist = histogram(yveq, min=vbinmin, binsize=1)
    oveq_hist = histogram(oveq, min=vbinmin, binsize=1)

    xvec_y = dindgen(n_elements(yveq_hist))*1d0 +vbinmin
    xvec_o = dindgen(n_elements(oveq_hist))*1d0 +vbinmin

    ;young_psf_arr[*,tbin] = yveq_hist/total(yveq_hist, /double)
    ;old_psf_arr[*,tbin] = oveq_hist/total(oveq_hist, /double)

    ypdf = yveq_hist/total(yveq_hist, /double)
    opdf = oveq_hist/total(oveq_hist, /double)

    yconf = hist_conf(yveq, binsize=0.01, binmin=0, conf_bounds=[0.16,0.84])
    oconf = hist_conf(oveq, binsize=0.01, binmin=0, conf_bounds=[0.16,0.84])

    yconfp = confidence_interval(yveq, min=0d0, conf_interval=0.68)
    oconfp = confidence_interval(oveq, min=0d0, conf_interval=0.68)

    yconf2 = hist_conf(yveq, binsize=0.01, binmin=0, conf_bounds=[0.025,0.975])
    oconf2 = hist_conf(oveq, binsize=0.01, binmin=0, conf_bounds=[0.025,0.975])

    if tbin ge 3 then oconf2 = hist_conf_ul(oveq, binsize=0.01, binmin=0, conf_limit=0.95)
    if tbin ge 10 then yconf2 = hist_conf_ul(yveq, binsize=0.01, binmin=0, conf_limit=0.95)

    yconf3 = hist_conf(yveq, binsize=0.01, binmin=0, conf_bounds=[0.005,0.995])
    oconf3 = hist_conf(oveq, binsize=0.01, binmin=0, conf_bounds=[0.005,0.995])



    yconf_bounds[*,tbin] = yconf
    oconf_bounds[*,tbin] = oconf

    yconf_bounds2[*,tbin] = yconf2
    oconf_bounds2[*,tbin] = oconf2
    
    yconf_bounds3[*,tbin] = yconf3
    oconf_bounds3[*,tbin] = oconf3

    yconf_boundsp[*,tbin] = yconfp
    oconf_boundsp[*,tbin] = oconfp

   titstr = "Teff: "+strtrim(tbin_vec[tbin],2)+" - "+strtrim(tbin_vec[tbin+1],2)+" K"

    ;plot, veq_vec_2, ypdf, ps=10, xr=[0,75], yr =
    ;[0,max(ypdf)>max(opdf)], tit=titstr
   ;plot, veq_vec_2, ypdf, ps=10, xr=[0,100], yr = [0,0.14],
   ;tit=titstr
   plot, xvec_y, ypdf, ps=10, xr=[0,100], yr = [0,0.14], tit=titstr
    oplot, replicate(yconfp[0],2), [0,1], linest=2
    oplot, replicate(yconfp[1],2), [0,1], linest=2
    oplot, replicate(yconf[1],2), [0,1], linest=2, color=200
    oplot, replicate(yconf[2],2), [0,1], linest=2, color=200
    
    ;plot, veq_vec_2, opdf, ps=10, tit = strtrim(oconf2,2)+" |
    ;"+strtrim(oconf,2)
    plot, xvec_o, opdf, ps=10, xr = [0,100], yr=[0,0.14]
    oplot, replicate(oconfp[0],2), [0,1], linest=2
    oplot, replicate(oconfp[1],2), [0,1], linest=2
    oplot, replicate(oconf[1],2), [0,1], linest=2, color=200
    oplot, replicate(oconf[2],2), [0,1], linest=2, color=200
    
    print, titstr
    print, "Young: "+strtrim(yconfp,2)+" | "+strtrim(yconf,2)
    print, "Old:   "+strtrim(oconfp,2)+" | "+strtrim(oconf,2)
    print, ''

    

    dxstop

endfor
stop   
;device, /close_file

!p.multi=[0,1,2]
 set_plot, 'ps'
 device, filename = "veq_vs_teff_agemulti_density_conf_2sig.eps"
 device, /color, bits=8
 device, xs=12,ys=10, /inches

;Young population
; plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Effective Temperature (K)", $
;   ytit="Equatorial Velocity (km/s)", title = "Age < 2 Gyrs", charsize=1.5, $
;   xr=[2600,4000], yr=[1,300], /ys, /nodata, /ylog, /xs




density, teff_realizations[*,young_star_idx], veq_realizations[*,young_star_idx], $
  /ylog, yr=[1,300], xr=[2600,4000], /ys, /xs, /dlog, drange=[5,1000], ct=-1, $
  ytit="Equatorial Velocity (km/s)", $
  title = "Age < 2 Gyrs", charsize=1.5
;overplot the confidence intervals
loadct, 13
for tbin = 0, ntbins-1 do begin
    
    xtemp = tbin_vec[tbin] + (tbinsize/2)
    if tbin le 9 then begin
        oplot, replicate(xtemp,2), yconf_bounds2[1:2,tbin], color=400, thick=3
        oplot, [-5,5]+xtemp, replicate(yconf_bounds2[1,tbin],2), color=400, thick=3
        oplot, [-5,5]+xtemp, replicate(yconf_bounds2[2,tbin],2), color=400, thick=3
        plotsym, 8, /fill
        oplot, [xtemp], [yconf_bounds2[0,tbin]], color=400, ps=8, thick=3
    endif else begin
        plotsym, 5, /fill
        oplot, [-5,5]+xtemp, replicate(yconf_bounds2[2,tbin],2), color=400, thick=3
        oplot, replicate(xtemp,2), [1d0,yconf_bounds2[2,tbin]], color=400, thick=3
        oplot, [xtemp], [1.1d0], color=400, ps=8, thick=3
    endelse
endfor


density, teff_realizations[*,old_star_idx], veq_realizations[*,old_star_idx], /ylog, $
  yr=[1,300], xr=[2600,4000], /ys, /xs, /dlog, drange=[5,1000], ct=-1, $
  xtit="Effective Temperature (K)", $
  title = "Age > 2 Gyrs", charsize=1.5
;overplot the confidence intervals
loadct, 13
for tbin = 0, ntbins-1 do begin
    
    xtemp = tbin_vec[tbin] + (tbinsize/2)
    if tbin le 2 then begin
        oplot, replicate(xtemp,2), oconf_bounds2[1:2,tbin], color=400, thick=3
        oplot, [-5,5]+xtemp, replicate(oconf_bounds2[1,tbin],2), color=400, thick=3
        oplot, [-5,5]+xtemp, replicate(oconf_bounds2[2,tbin],2), color=400, thick=3
        plotsym, 8, /fill
        oplot, [xtemp], [oconf_bounds2[0,tbin]], color=400, ps=8, thick=3
    endif else begin
        plotsym, 5, /fill
        ; draw hat
        oplot, [-5,5]+xtemp, replicate(oconf_bounds2[2,tbin],2), color=400, thick=3
        ; draw vertical line
        oplot, replicate(xtemp,2), [1d0,oconf_bounds2[2,tbin]], color=400, thick=3
        ; draw downwards arrow
        oplot, [xtemp], [1.1d0], color=400, ps=8, thick=3
    endelse
endfor

device, /close_file

stop

;;;;;;;;;;;;;;;;;;;;;;;
set_plot, 'ps'
device, filename = "veq_vs_teff_agemulti_density.eps"
device, /color, bits=8
device, xs=12,ys=10, /inches
loadct, 13


!p.multi=[0,1,2]

errstring = textoidl('\sigma_{Teff}')

;Young Detections
plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Effective Temperature (K)", $
  ytit="Equatorial Velocity (km/s)", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2600,4000], yr=[1,300], /ys, /nodata, /ylog, /xs
for rnum = 0, nsamples/2-1 do begin
    ;oplot, teff_realizations[rnum,yng_det_idx],
    ;veq_realizations[rnum,yng_det_idx], ps=8
    oplot, teff_realizations[rnum,young_star_idx], veq_realizations[rnum,young_star_idx], ps=3
    ;Young non-detections
    plotsym, 5, 0.1, /fill ;Upside-down triangle
    ;oplot, teff_realizations[rnum,yng_ndet_idx],
    ;veq_realizations[rnum,yng_ndet_idx], ps=8
    ;oplot, teff_realizations[rnum,yng_ndet_idx], veq_realizations[rnum,yng_ndet_idx], ps=3
endfor
;usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
plotsym, 0, /fill
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
;usersym, cos(usym), sin(usym), /fill, psize=0.1 ;Return to filled in circles
plotsym, 0, 0.1, /fill

;Old Detections
plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Effective Temperature (K)", $
  ytit="Equatorial Velocity (km/s)", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2600,4000], yr=[1,300], /ys, /nodata, /ylog, /xs
for rnum = 0, nsamples/2-1 do begin
    ;oplot, teff_realizations[rnum,old_det_idx],
    ;veq_realizations[rnum,old_det_idx], ps=8
    oplot, teff_realizations[rnum,old_star_idx], veq_realizations[rnum,old_star_idx], ps=3
    ;Old non-detections
    plotsym, 5, 0.1, /fill ;Upside-down triangle
    ;oplot, teff_realizations[rnum,old_ndet_idx], veq_realizations[rnum,old_ndet_idx], ps=8
endfor

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
;usersym, cos(usym), sin(usym), /fill ;Return to filled in circles
plotsym, 0, /fill

!p.multi = 0

device, /close_file

stop

end
