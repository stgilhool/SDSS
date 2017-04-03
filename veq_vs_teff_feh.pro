; Find a relation between vsini, teff and age

function sini_probability, sample
common pdf_info, age_pdf, age_vec
; Give the probability for each element in sample

sini_prob = sample/sqrt(1d0-sample^2)

return, sini_prob

end

function age_probability, sample
common pdf_info, age_pdf, age_vec
; Give the probability for each element in sample

prob_vec = interpol(age_pdf, age_vec, sample)

return, prob_vec

end

pro veq_vs_teff_feh

common pdf_info, age_pdf, age_vec

; Read in metallicity
minfo = mrdfits('/home/stgilhool/APOGEE/master_info6.fits',1)
feh = minfo.feh_rgrss

; Read in aspcap data
asp_str = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits',1)
teff_asp = asp_str.teff


; Read in Vsini, Teff, and U,V,W (although we're just using W)
info_str = mrdfits('spacevelocity_wflags.fits',1)
uarr = info_str.u
uarr = -1d0*uarr ;Flip the sign due to silly convention
varr = info_str.v
warr = info_str.w
vsini = info_str.vsini
teff = info_str.teff ;these are teff_vfit
p_membership = info_str.thickthin_prob
absk = info_str.absk_model


; Read in gold-sample flags
flags = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/gold_sample.fits',0)
; Save flag indices
g_idx = where(flags eq 0, ngood, complement=b_idx, ncomplement=nbad)
napo = n_elements(flags) 
nstars = n_elements(g_idx) 

; Read in the age PDFS
stellar_age_estimator, age_str

age_vec = age_str.age
age_pdf_arr = age_str.prob

;;; Use Mann 2012 to get M and R
rcoeff = [16.77d0, -54.321d0, 57.6627d0, -19.6994d0]
rfehco = 0.4565d0
mcoeff = [0.5858d0, 0.3872d0, -0.1217d0, 0.0106d0, -2.7262d-4]

radius = poly(teff/3500d0, rcoeff) * (1d0 + rfehco*feh)
mass = poly(absk, mcoeff)

msun = 1.989d30 ;kg 
rsun = 6.957d5 ;km


moment_i = (1d0/5d0)*mass*radius^2 ; approximately 1/2 of the moment for sphere, see moment_of_inertia.pro and cullen's email
;moment_permass = moment_i/mass
moment_permass = (1d0/5d0) * radius^2
ang_factor = (1d0/5d0) * radius

; Use ASPCAP teff only, except for where they are undefined
;teff_asp = teff_asp[g_idx]
newg_idx = where(teff_asp gt 0 and teff_asp le 4100 and flags eq 0, nstars)
; modify vectors
uarr = uarr[newg_idx]
varr = uarr[newg_idx]
warr = uarr[newg_idx]
vsini = vsini[newg_idx]
teff = teff_asp[newg_idx] ; using aspcap teff
age_pdf_arr = age_pdf_arr[*,newg_idx]
p_membership = p_membership[newg_idx]

absk = absk[newg_idx]
feh = feh[newg_idx]

mass = mass[newg_idx]
radius = radius[newg_idx]

moment_i = moment_i[newg_idx]
moment_permass = moment_permass[newg_idx]
ang_factor = ang_factor[newg_idx]

print, nstars
dxstop



;;;;;; Generate the sini pdf
inclination = dblarr(1001)/1000d0*!dpi/2d0
sini_dist = tan(inclination)


; Generate N samples for each star, based on their age, temperature,
; and vsini pdfs
nsamples = 1d5

age_realizations = dblarr(nsamples,nstars)
teff_realizations = dblarr(nsamples, nstars)
veq_realizations = dblarr(nsamples, nstars)

veq_values = dblarr(3, nstars)
age_values = dblarr(3, nstars)
;sini_realizations = dblarr(nsamples, nstars)
;vsini_realizations = dblarr(nsamples, nstars)

window, 0, xs=1200,ys=1000
!p.multi=[0,1,2]

for starnum = 0, nstars-1 do begin

    ; PDF
    age_pdf = age_pdf_arr[*,starnum]
    age_mode_idx = where(age_pdf eq max(age_pdf), nagemode)
    if nagemode then age_mode = age_vec[age_mode_idx] $
      else message, "Problem with age_mode_idx"
    ; save age
    age_values[0,starnum] = age_mode
    
    ; Get confidence interval for age
     age_cdf = total(age_pdf, /cumulative)
     ; Interpolate to find the age values where the CDF = 5 and 95%
     age_confidence = interpol(age_vec, age_cdf, [0.05,0.95])
     if age_confidence[0] gt age_mode then age_confidence[0]=age_mode
     if age_confidence[1] lt age_mode then age_confidence[1]=age_mode
     
     ;save age conf
     age_values[1,starnum] = abs(age_confidence[0]-age_mode)
     age_values[2,starnum] = abs(age_confidence[1]-age_mode)

     ;test
     ; plot, age_vec, age_pdf
;      oplot, replicate(age_confidence[0],2), [0,1], linest=2
;      oplot, replicate(age_confidence[1],2), [0,1], linest=2
;      oplot, replicate(age_mode,2), [0,1], linest=3
;      plot, age_vec, age_cdf
;      oplot, replicate(age_confidence[0],2), [0,1], linest=2
;      oplot, replicate(age_confidence[1],2), [0,1], linest=2
;      oplot, replicate(age_mode,2), [0,1], linest=3
;      dxstop
    
; Create a sample for this star
    ;sample_pop_age = create_sample_distribution('age_probability',minmax(age_vec),nsamples)
    ;sample_pop_teff =
    ;create_sample_distribution('teff_probability',minmax(teff),nsamples)
    ; Get a sample of realistic sini values 
    sample_pop_sini = create_sample_distribution('sini_probability',[0d0,1d0],nsamples)
    ;add 2km/s uncertainty on measured val
    if vsini[starnum] le 5 then begin
        ; non-detection
        sample_pop_vsini = randomu(seed, nsamples) * 5d0
    endif else if vsini[starnum] gt 5 then begin
        ; detection
        vsini_uncertainty = randomn(seed,nsamples)* 2d0
        sample_pop_vsini = (vsini[starnum] + vsini_uncertainty)
        ; if values lt 5 arrise, redo them
        low_idx = where(sample_pop_vsini le 5, nlow)
        while nlow gt 0 do begin
            new_uncertainty = randomn(seed,nlow) * 2d0
            sample_pop_vsini[low_idx] = vsini[starnum] + new_uncertainty
            low_idx = where(sample_pop_vsini le 5, nlow)
        endwhile
    endif


    ;temporary
    ;vsini_realizations[*,starnum] = sample_pop_vsini
    ;sini_realizations[*,starnum] = sample_pop_sini

    sample_pop_veq = sample_pop_vsini/sample_pop_sini

    ; Determine the 95% confidence interval for the V_eq
     vbinsize = 0.1d0
    
     veq_hist = histogram(sample_pop_veq, binsize=vbinsize, min=0)
     
     veq_plotvec = dindgen(n_elements(veq_hist))*vbinsize
     veq_vec = veq_plotvec + vbinsize/2
     ;Make the CDF and PDF
     veq_pdf = veq_hist/total(veq_hist,/double)
     veq_cdf = total(veq_pdf, /cumulative)
     ; Interpolate to find the Veq values where the CDF = 5 and 95%
     veq_confidence = interpol(veq_vec, veq_cdf, [0.05,0.95])
    

     

    ; Determine the mode
     ;smooth the pdf
     smoothpdf = gauss_smooth(veq_pdf,3,/edge_wrap)
     ; find the mode by finding the index of the max of smoothpdf
     mode1_idx = where(smoothpdf eq max(smoothpdf))
     mode1 = mean(veq_vec[mode1_idx])
     
     
     veq_error = abs(veq_confidence-mode1)
     
     ; find the mode by d(pdf)/dveq=0
     ; dpdf = (smoothpdf[1:*]-smoothpdf[0:-2])
;      ; Nevermind, mode1 works fine

;      ;Plot the pdf
;      plot, veq_plotvec, veq_pdf, ps=10, xr=[0,100]
;      oplot, replicate(veq_confidence[0],2), [0,1], linesty=2
;      oplot, replicate(veq_confidence[1],2), [0,1], linesty=2
;      oplot, replicate(mode1,2), [0,1], linesty=3
;      oplot, replicate(mode2,2), [0,1], linesty=3, color=200
;      ;Plot the cdf
;      plot, veq_plotvec, veq_cdf, ps=10, xr=[0,100]
;      oplot, replicate(veq_confidence[0],2), [0,1], linesty=2
;      oplot, replicate(veq_confidence[1],2), [0,1], linesty=2
;      oplot, replicate(mode1,2), [0,1], linesty=3
;      oplot, replicate(mode2,2), [0,1], linesty=3, color=200
;      ;Plot the smoothed pdf
;      plot, veq_plotvec[0:-2], dpdf, ps=10, xr=[0,100]
;      oplot, veq_plotvec[0:-2], replicate(0d0,n_elements(dpdf)), linest=4 
;      oplot, replicate(veq_confidence[0],2), [0,1], linesty=2
;      oplot, replicate(veq_confidence[1],2), [0,1], linesty=2
;      oplot, replicate(mode1,2), [0,1], linesty=3
     


;      stop

    teff_uncertainty = randomn(seed,nsamples)*70d0
    sample_pop_teff = teff[starnum] + teff_uncertainty



    ; Save the sample
    ;age_realizations[*,starnum] = sample_pop_age
    teff_realizations[*,starnum] = sample_pop_teff
    veq_realizations[*,starnum] = sample_pop_veq
    ;veq_values[*,starnum] = [mode1,veq_confidence]
    veq_values[*,starnum] = [mode1,veq_error]
    print, "Starnum: "+strtrim(starnum,2)+" finished"

endfor


;;; Start plotting
; Filled circle definition
usym = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
USERSYM, COS(usym), SIN(usym), /FILL

!p.multi=0

young_star_idx = where(age_values[0,*] le 2, nyoung, $
                       complement=old_star_idx, ncomplement=nold)

set_plot, 'ps'
device, filename = "veq_vs_teff_agecolor_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


plot, teff, veq_values[0,*], ps=8, xtit="Teff", ytit="V_eq", yr=[1,300], xr=[2500,4000], /ylog, /ys
oploterror, teff, veq_values[0,*], replicate(70d0,nstars), veq_values[1,*], /lobar, /nohat, ps=3
oploterror, teff, veq_values[0,*], replicate(70d0,nstars), veq_values[2,*], /hibar, /nohat, ps=3
oplot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, color=200
oploterror, teff[young_star_idx], veq_values[0,young_star_idx], replicate(70d0,nyoung), veq_values[1,young_star_idx], /lobar, errcolor='red', /nohat, ps=3
oploterror, teff[young_star_idx], veq_values[0,young_star_idx], replicate(70d0,nyoung), veq_values[2,young_star_idx], /hibar, errcolor='red', /nohat, ps=3

device, /close_file

stop

set_plot, 'ps'
device, filename = "veq_vs_teff_agemulti_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi=[0,1,2]



;Young
plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2500,4000], yr=[1,300], /ylog, /ys
oploterror, teff[young_star_idx], veq_values[0,young_star_idx], replicate(70d0,nyoung), veq_values[1,young_star_idx], /lobar,/nohat, ps=3
oploterror, teff[young_star_idx], veq_values[0,young_star_idx], replicate(70d0,nyoung), veq_values[2,young_star_idx], /hibar,/nohat, ps=3
;Old
plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2500,4000], yr=[1,300], /ylog, /ys
oploterror, teff[old_star_idx], veq_values[0,old_star_idx], replicate(70d0,nold), $
  veq_values[1,old_star_idx], /lobar, ps=3,/nohat
oploterror, teff[old_star_idx], veq_values[0,old_star_idx], replicate(70d0,nold), $
  veq_values[2,old_star_idx], /hibar, ps=3,/nohat

!p.multi = 0

device, /close_file

; plot, teff, veq_values[0,*], ps=6, xtit="Teff", ytit="V_eq", yr=[1,300], xr=[2500,4000]
; oploterror, teff, veq_values[0,*], replicate(70d0,nstars), veq_values[1,*], /lobar
; oploterror, teff, veq_values[0,*], replicate(70d0,nstars), veq_values[2,*], /hibar

stop

set_plot, 'ps'
device, filename = "ang_mom_vs_teff_agecolor_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


;ang_mom = veq_values * moment_permass/radius
ang_mom = veq_values * rebin(reform(ang_factor,1,nstars), 3, nstars)

plot, teff, ang_mom[0,*], ps=8, xtit="Teff", ytit="Ang. Momentum per unit mass", $
  xr=[2500,4000], yr=[0.1,10], /ylog, /ys
oploterror, teff, ang_mom[0,*], replicate(70d0,nstars), ang_mom[1,*], /lobar, /nohat, ps=3
oploterror, teff, ang_mom[0,*], replicate(70d0,nstars), ang_mom[2,*], /hibar, /nohat, ps=3
oplot, teff[young_star_idx], ang_mom[0,young_star_idx], ps=8, color=200
oploterror, teff[young_star_idx], ang_mom[0,young_star_idx], replicate(70d0,nyoung), ang_mom[1,young_star_idx], /lobar, errcolor='red', /nohat, ps=3
oploterror, teff[young_star_idx], ang_mom[0,young_star_idx], replicate(70d0,nyoung), ang_mom[2,young_star_idx], /hibar, errcolor='red', /nohat, ps=3

device, /close_file

stop


set_plot, 'ps'
device, filename = "ang_mom_vs_teff_agemulti_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi=[0,1,2]
;Young
plot, teff[young_star_idx], ang_mom[0,young_star_idx], ps=8, xtit="Teff", $
  ytit="L per unit mass", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2500,4000], yr=[0.1,10], /ylog, /ys
oploterror, teff[young_star_idx], ang_mom[0,young_star_idx], replicate(70d0,nyoung), ang_mom[1,young_star_idx], /lobar, /nohat, ps=3
oploterror, teff[young_star_idx], ang_mom[0,young_star_idx], replicate(70d0,nyoung), ang_mom[2,young_star_idx], /hibar, /nohat, ps=3
;Old
plot, teff[old_star_idx], ang_mom[0,old_star_idx], ps=8, xtit="Teff", $
  ytit="L per unit mass", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2500,4000], yr=[0.1,10], /ylog, /ys
oploterror, teff[old_star_idx], ang_mom[0,old_star_idx], replicate(70d0,nold), $
  ang_mom[1,old_star_idx], /lobar, /nohat, ps=3
oploterror, teff[old_star_idx], ang_mom[0,old_star_idx], replicate(70d0,nold), $
  ang_mom[2,old_star_idx], /hibar, /nohat, ps=3

!p.multi = 0

device, /close_file

stop

tbinsz = 200
tbinmin = 2600
ntbins = 7
thist = histogram(teff, binsize=tbinsz, min=tbinmin, nbins=ntbins, reverse_indices=ri)
tbin_ends = lindgen(n_elements(thist)+1)*tbinsz + tbinmin
tbin_ctrs = lindgen(n_elements(thist))*tbinsz + tbinmin + (tbinsz/2)


set_plot, 'ps'
device, filename = "veq_vs_age_teff_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi = [0,2,4]

for bin = 0, ntbins-1 do begin
    ;Check if there is data in the bin
    bintitle = "Teff: "+strjoin(strtrim(tbin_ends[bin:bin+1],2),'-')
    if ri[bin] ne ri[bin+1] then begin
        data_idx = ri[ri[bin]:ri[bin+1]-1]
        bin_vvec = veq_values[*,data_idx]
        bin_avec = age_values[*,data_idx]

        binv = bin_vvec[0,*]
        binv_errlo = bin_vvec[1,*]
        binv_errhi = bin_vvec[2,*]
        
        bina = bin_avec[0,*]
        bina_errlo = bin_avec[1,*]
        bina_errhi = bin_avec[2,*]
        
        plot, bina, binv, ps=8, xtit="Age", $
          ytit="V_eq", title = bintitle, charsize=1.5, $
          xr=[0,11], yr=[1,300], /ylog, /ys, /xs
        oploterror, bina, binv, bina_errlo, binv_errlo, /lobar, /nohat, ps=3
        oploterror, bina, binv, bina_errhi, binv_errhi, /hibar, /nohat, ps=3
    endif else begin
        plot, age_values[0,*], veq_values[0,*], xtit="Age", $
          ytit="V_eq", title = bintitle, charsize=1.5, $
          xr=[0,11], yr=[1,300], /nodata, /ylog, /ys, /xs
    endelse
endfor

device, /close_file

stop

set_plot, 'ps'
device, filename = "ang_mom_vs_age_teff_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi = [0,2,4]

for bin = 0, ntbins-1 do begin
    ;Check if there is data in the bin
    bintitle = "Teff: "+strjoin(strtrim(tbin_ends[bin:bin+1],2),'-')
    if ri[bin] ne ri[bin+1] then begin
        data_idx = ri[ri[bin]:ri[bin+1]-1]
        bin_mvec = ang_mom[*,data_idx]
        bin_avec = age_values[*,data_idx]
        
        plot, bin_avec[0,*], bin_mvec[0,*], ps=8, xtit="Age", $
          ytit="L per M", title = bintitle, charsize=1.5, $
          xr=[0,10], yr=[0.1,10], /ylog, /ys
        oploterror, bin_avec[0,*], bin_mvec[0,*], bin_avec[1,*], bin_mvec[1,*], /lobar, /nohat, ps=3
        oploterror, bin_avec[0,*], bin_mvec[0,*], bin_avec[2,*], bin_mvec[2,*], /hibar, /nohat, ps=3
    endif else begin
        plot, age_values[0,*], ang_mom[0,*], xtit="Age", $
          ytit="V_eq", title = bintitle, charsize=1.5, $
          xr=[0,10], yr=[0.1,10], /nodata, /ylog, /ys
    endelse
endfor



device, /close_file

stop


;;; For the _limits_ plots, we replace the non-detections with upside
;;; triangles to denote upper limits

;ndet_idx = where(vsini le 5, nndet, complement=det_idx, ncomplement=ndet)

qage = age_values[0,*]

yng_ndet_idx = where(qage le 2 and vsini le 5, nyndet)
yng_det_idx = where(qage le 2 and vsini gt 5, nydet)
old_ndet_idx = where(qage gt 2 and vsini le 5, nondet)
old_det_idx = where(qage gt 2 and vsini gt 5, nodet)

set_plot, 'ps'
device, filename = "veq_vs_teff_agemulti_limits_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi=[0,1,2]



;Young Detections
plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2500,4000], yr=[1,300], /ylog, /ys, /nodata
oplot, teff[yng_det_idx], veq_values[0,yng_det_idx], ps=8
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], replicate(70d0,nydet), veq_values[1,yng_det_idx], /lobar,/nohat, ps=3
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], replicate(70d0,nydet), veq_values[2,yng_det_idx], /hibar,/nohat, ps=3
;Young non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[yng_ndet_idx], veq_values[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

;Old Detections
plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2500,4000], yr=[1,300], /ylog, /ys, /nodata
oplot, teff[old_det_idx], veq_values[0,old_det_idx], ps=8
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], replicate(70d0,nodet), veq_values[1,old_det_idx], /lobar,/nohat, ps=3
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], replicate(70d0,nodet), veq_values[2,old_det_idx], /hibar,/nohat, ps=3
;Old non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[old_ndet_idx], veq_values[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

!p.multi = 0

device, /close_file



set_plot, 'ps'
device, filename = "ang_mom_vs_teff_agemulti_limits_ylog.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi=[0,1,2]


;Young Detectiona
plot, teff[young_star_idx], ang_mom[0,young_star_idx], ps=8, xtit="Teff", $
  ytit="L per unit Mass", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2500,4000], yr=[0.1,10], /ylog, /ys, /nodata
oplot, teff[yng_det_idx], ang_mom[0,yng_det_idx], ps=8
oploterror, teff[yng_det_idx], ang_mom[0,yng_det_idx], replicate(70d0,nydet), ang_mom[1,yng_det_idx], /lobar,/nohat, ps=3
oploterror, teff[yng_det_idx], ang_mom[0,yng_det_idx], replicate(70d0,nydet), ang_mom[2,yng_det_idx], /hibar,/nohat, ps=3
;Young non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[yng_ndet_idx], ang_mom[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

;Old Detections
plot, teff[old_star_idx], ang_mom[0,old_star_idx], ps=8, xtit="Teff", $
  ytit="L per unit Mass", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2500,4000], yr=[0.1,10], /ylog, /ys, /nodata
oplot, teff[old_det_idx], ang_mom[0,old_det_idx], ps=8
oploterror, teff[old_det_idx], ang_mom[0,old_det_idx], replicate(70d0,nodet), ang_mom[1,old_det_idx], /lobar,/nohat, ps=3
oploterror, teff[old_det_idx], ang_mom[0,old_det_idx], replicate(70d0,nodet), ang_mom[2,old_det_idx], /hibar,/nohat, ps=3
;Old non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[old_ndet_idx], ang_mom[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

!p.multi = 0

device, /close_file



set_plot, 'ps'
device, filename = "veq_vs_teff_agemulti_limits.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi=[0,1,2]



;Young Detections
plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2500,4000], yr=[1,300], /ys, /nodata
oplot, teff[yng_det_idx], veq_values[0,yng_det_idx], ps=8
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], replicate(70d0,nydet), veq_values[1,yng_det_idx], /lobar,/nohat, ps=3
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], replicate(70d0,nydet), veq_values[2,yng_det_idx], /hibar,/nohat, ps=3
;Young non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[yng_ndet_idx], veq_values[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

;Old Detections
plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2500,4000], yr=[1,300], /ys, /nodata
oplot, teff[old_det_idx], veq_values[0,old_det_idx], ps=8
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], replicate(70d0,nodet), veq_values[1,old_det_idx], /lobar,/nohat, ps=3
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], replicate(70d0,nodet), veq_values[2,old_det_idx], /hibar,/nohat, ps=3
;Old non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[old_ndet_idx], veq_values[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

!p.multi = 0

device, /close_file



set_plot, 'ps'
device, filename = "veq_vs_teff_fehmulti_limits.eps"
device, /color, bits=8
device, xs=11,ys=8, /inches
loadct, 13


!p.multi=[0,1,2]



lfeh_star_idx = where(feh le -0.086, nlfeh, $
                       complement=hfeh_star_idx, ncomplement=nhfeh)


lfeh_ndet_idx = where(feh le -0.086 and vsini le 5, nlndet)
lfeh_det_idx = where(feh le -0.086 and vsini gt 5, nldet)
hfeh_ndet_idx = where(feh gt -0.086 and vsini le 5, nhndet)
hfeh_det_idx = where(feh gt -0.086 and vsini gt 5, nhdet)


;Low Fe/H Detections
plot, teff[lfeh_star_idx], veq_values[0,lfeh_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title = "[Fe/H] < -0.08 (median [Fe/H])", charsize=1.5, $
  xr=[2500,4000], yr=[3,300], /ys, /nodata, /ylog
oplot, teff[lfeh_det_idx], veq_values[0,lfeh_det_idx], ps=8
oploterror, teff[lfeh_det_idx], veq_values[0,lfeh_det_idx], replicate(70d0,nldet), veq_values[1,lfeh_det_idx], /lobar,/nohat, ps=3
oploterror, teff[lfeh_det_idx], veq_values[0,lfeh_det_idx], replicate(70d0,nldet), veq_values[2,lfeh_det_idx], /hibar,/nohat, ps=3
;Young non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[lfeh_ndet_idx], veq_values[2,lfeh_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

;High feh Detections
plot, teff[hfeh_star_idx], veq_values[0,hfeh_star_idx], ps=8, xtit="Teff", $
  ytit="V_eq", title="[Fe/H] > -0.08", charsize = 1.5, $
  xr=[2500,4000], yr=[3,300], /ys, /nodata, /ylog
oplot, teff[hfeh_det_idx], veq_values[0,hfeh_det_idx], ps=8
oploterror, teff[hfeh_det_idx], veq_values[0,hfeh_det_idx], replicate(70d0,nhdet), veq_values[1,hfeh_det_idx], /lobar,/nohat, ps=3
oploterror, teff[hfeh_det_idx], veq_values[0,hfeh_det_idx], replicate(70d0,nhdet), veq_values[2,hfeh_det_idx], /hibar,/nohat, ps=3
;High feh non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[hfeh_ndet_idx], veq_values[2,hfeh_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

!p.multi = 0

device, /close_file



set_plot, 'ps'
device, filename = "veq_vs_teff_agemulti_clear.eps"
device, /color, bits=8
device, xs=12,ys=10, /inches
loadct, 13


!p.multi=[0,1,2]

errstring = textoidl('\sigma_{Teff}')

;Young Detections
plot, teff[young_star_idx], veq_values[0,young_star_idx], ps=8, xtit="Effective Temperature (K)", $
  ytit="Equatorial Velocity (km/s)", title = "Age < 2 Gyrs", charsize=1.5, $
  xr=[2600,4000], yr=[5,300], /ys, /nodata, /ylog, /xs
oplot, teff[yng_det_idx], veq_values[0,yng_det_idx], ps=8
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], veq_values[1,yng_det_idx], /lobar,/nohat, ps=3
oploterror, teff[yng_det_idx], veq_values[0,yng_det_idx], veq_values[2,yng_det_idx], /hibar,/nohat, ps=3
oploterror, [2700], [6], replicate(68d0,2), replicate(0,2), /nohat, ps=8
xyouts, [2690], [7], errstring
;Young non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[yng_ndet_idx], veq_values[2,yng_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

;Old Detections
plot, teff[old_star_idx], veq_values[0,old_star_idx], ps=8, xtit="Effective Temperature (K)", $
  ytit="Equatorial Velocity (km/s)", title="Age > 2 Gyrs", charsize = 1.5, $
  xr=[2600,4000], yr=[5,300], /ys, /nodata, /ylog, /xs
oplot, teff[old_det_idx], veq_values[0,old_det_idx], ps=8
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], veq_values[1,old_det_idx], /lobar,/nohat, ps=3
oploterror, teff[old_det_idx], veq_values[0,old_det_idx], veq_values[2,old_det_idx], /hibar,/nohat, ps=3
oploterror, [2700], [6], replicate(68d0,2), replicate(0,2), /nohat, ps=8
xyouts, [2690], [7], errstring
;Old non-detections
plotsym, 5, /fill ;Upside-down triangle
oplot, teff[old_ndet_idx], veq_values[2,old_ndet_idx], ps=8 ;Plot the upper limit of the confidence interval
usersym, cos(usym), sin(usym), /fill ;Return to filled in circles

!p.multi = 0

device, /close_file



print, "DONE!!!"

stop


; help, veq_realizations
; help, sini_realizations
; help, vsini_realizations

; window,0,xs=1600,ys=900
; !p.multi = [0,1,3]

; for starnum = 0, nstars-1 do begin

; ;     v_vec = dblarr(600)/2.
; ;     s_vec = dblarr(101)/100.

;     veq_i = veq_realizations[*,starnum]
;     vsini_i = vsini_realizations[*,starnum]
;     sini_i = sini_realizations[*,starnum]
    
;     vbinsize = 0.5d0
;     sbinsize = 0.01d0
    
;     veq_hist = histogram(veq_i, binsize=vbinsize, min=0)
;     vsini_hist = histogram(vsini_i, binsize=vbinsize, min=0)

;     sini_hist = histogram(sini_i, binsize=sbinsize, min=0)

;     ve_vec = dindgen(n_elements(veq_hist))*vbinsize
;     vs_vec = dindgen(n_elements(vsini_hist))*vbinsize
;     si_vec = dindgen(n_elements(sini_hist))*sbinsize

;     plot, si_vec, sini_hist, ps=10, ytit = "N", xtit="sini"
;     plot, vs_vec, vsini_hist, ps=10, ytit="N", xtit="Vsini"
;     plot, ve_vec, veq_hist, ps=10, ytit="N", xtit="Veq"
    
;     dxstop
; endfor
; stop

;;;
; Now we will have N realizations of the data with age, teff and vsini
; for each star

; For each realization, calculate the variance/covariance matrix

covar = dblarr(3,3,nsamples)
corr  = dblarr(3,nsamples)


for trial = 0, nsamples-1 do begin

    
    trial_age = reform(age_realizations[trial,*])
    mean_age = mean(trial_age)
    trial_age = trial_age-mean_age

    trial_teff = reform(teff_realizations[trial,*])
    mean_teff = mean(trial_teff)
    trial_teff = (trial_teff-mean_teff)

    trial_veq = reform(veq_realizations[trial,*])
    mean_veq = mean(trial_veq)
    trial_veq = trial_veq-mean_veq

    matrix = [[trial_age],[trial_teff],[trial_veq]]

    covar_i = matrix_multiply(matrix,matrix,/atranspose)

    ; Divide by N
    covar_i = covar_i/(nstars-1)

    ; Get the diagonal entries
    covar_diag = diag_matrix(covar_i)
    
    ; Get the correlations
    corr_01 = covar_i[0,1]/sqrt(covar_diag[0]*covar_diag[1]) ;age/temp
    corr_02 = covar_i[0,2]/sqrt(covar_diag[0]*covar_diag[2]) ;age/veq
    corr_12 = covar_i[1,2]/sqrt(covar_diag[1]*covar_diag[2]) ;temp/veq

    corr[*,trial] = [corr_01,corr_02,corr_12]


 ;   covar_i = covar_i/total(sqrt(covar_diag))

    help, covar_i
    ;dxstop

    covar[*,*,trial] = covar_i

    ;dxstop

endfor

help, corr
help, covar

set_plot, 'ps'
device, filename = "vsini_teff_age_correlation_both.eps"
device, /color, bits=8
device, xs=6,ys=4, /inches
loadct, 13

;window, 0, xs=1600, ys=900
!p.multi = [0,1,3]
;plot, corr[0,*], ps=6, title="Age-Teff correlation", charsize=1.5
;plot, corr[1,*], ps=6, title="Age-Veq correlation", charsize=1.5
;plot, corr[2,*], ps=6, title="Teff-Veq correlation", charsize=1.5

binsize = 0.005d0
at_hist = histogram(corr[0,*], binsize=binsize, omin=atom)
av_hist = histogram(corr[1,*], binsize=binsize, omin=avom)
tv_hist = histogram(corr[2,*], binsize=binsize, omin=tvom)

at_vec = dindgen(n_elements(at_hist))*binsize + atom
av_vec = dindgen(n_elements(av_hist))*binsize + avom
tv_vec = dindgen(n_elements(tv_hist))*binsize + tvom

plot, at_vec, at_hist, ps=10, title="Distribution of Age-Teff correlation", charsize=1.5, xr=[-0.2,0.2],/xs
plot, av_vec, av_hist, ps=10, title="Distribution of Age-Veq correlation", charsize=1.5, xr=[-0.2,0.2],/xs
plot, tv_vec, tv_hist, ps=10, title="Distribution of Teff-Veq correlation", charsize=1.5, xr=[-0.2,0.2],/xs

!p.multi=0

device, /close_file

stop
end
