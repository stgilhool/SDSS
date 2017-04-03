; Find a relation between vsini, teff and age

function fit_stepfunction, p, veq_data=veq, teff_data=teff, tcrit=crit_t

high_v = p[0]
low_v = p[1]
;crit_t = p[2]


; Sort the data into low and high temp
low_temp_idx = where(teff le crit_t, nlow, $
                     complement=high_temp_idx, ncomplement=nhigh)
; Deviates
dev = dblarr(n_elements(veq)) 
dev[low_temp_idx] = veq[low_temp_idx]-high_v
dev[high_temp_idx] = veq[high_temp_idx]-low_v

; Plot
;  tvec = findgen(1501) + 2600
;  lt_plot_idx = where(tvec le crit_t, nlt, $
;                      ncomplement=nht, complement=ht_plot_idx)
;  plot, teff, veq, /nodata, yr=[0,300]
;  oplot, tvec[lt_plot_idx], replicate(high_v,nlt)
;  oplot, tvec[ht_plot_idx], replicate(low_v,nht)
;  oplot, teff, veq, ps=6
;  wait, 0.01

return, dev

end

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

pro veq_vs_teff_stepfit_tgrid, load_data=load_data

if n_elements(load_data) eq 0 then load_data = 0
if load_data then goto, loaddata

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

save, /variables, filename = "veq_realization_data.sav"
loaddata:
restore, "veq_realization_data.sav"



young_star_idx = where(age_values[0,*] le 2, nyoung, $
                       complement=old_star_idx, ncomplement=nold)

qage = age_values[0,*]

yng_ndet_idx = where(qage le 2 and vsini le 5, nyndet)
yng_det_idx = where(qage le 2 and vsini gt 5, nydet)
old_ndet_idx = where(qage gt 2 and vsini le 5, nondet)
old_det_idx = where(qage gt 2 and vsini gt 5, nodet)


;ang_mom = veq_values * moment_permass/radius
ang_mom = veq_values * rebin(reform(ang_factor,1,nstars), 3, nstars)

sample_parinfo = {value:0d0, $
                  limited:[1,1], $
                  limits:[0,0], $
                  step:0}

parinfo = replicate(sample_parinfo, 2)

; pars are v1, v2
guess = [15d0, 5d0]
parinfo.value = guess

parinfo[0].limits = [0d0,200d0]
parinfo[1].limits = [0d0,200d0]

;young_results = dblarr(2,nsamples)
;old_results = dblarr(2,nsamples)

; Set up vector for stepping through tcrits
ntrials_tcrit = 25
teff_trial = lindgen(ntrials_tcrit)*25 + 3000
; Set up vector for binning by teff and averaging
nbins_teff = 14
teff_bin_vec = lindgen(nbins_teff+1)*100 + 2600

; Set up structure to store results of fits
result_str_sample = {veq_young_cool:dblarr(ntrials_tcrit), $
                     veq_young_hot:dblarr(ntrials_tcrit), $
                     chi2_young:dblarr(ntrials_tcrit), $
                     status_young:dblarr(ntrials_tcrit), $
                     veq_old_cool:dblarr(ntrials_tcrit), $
                     veq_old_hot:dblarr(ntrials_tcrit), $
                     chi2_old:dblarr(ntrials_tcrit), $
                     status_old:dblarr(ntrials_tcrit), $
                     best_young_idx:0, $
                     best_tcrit_young:0L, $
                     best_veq_young_cool:0d0, $
                     best_veq_young_hot:0d0, $
                     best_old_idx:0, $
                     best_tcrit_old:0L, $
                     best_veq_old_cool:0d0, $
                     best_veq_old_hot:0d0, $
                     veq_bin_avg_young:dblarr(nbins_teff), $
                     veq_bin_stddev_young:dblarr(nbins_teff), $
                     veq_bin_stderr_young:dblarr(nbins_teff), $
                     veq_bin_avg_old:dblarr(nbins_teff), $
                     veq_bin_stddev_old:dblarr(nbins_teff), $
                     veq_bin_stderr_old:dblarr(nbins_teff)}

result_str = replicate(result_str_sample, nsamples)

; SET UP STRUCTURE TO STORE RESULTS OF AVERAGE VEQ BY TEFF BIN


tcrit = 0 ;placeholder


; Loop through each iteration and fit a step function
for fitnum = 0, nsamples-1 do begin

    ysamp = veq_realizations[fitnum, young_star_idx]
    osamp = veq_realizations[fitnum, old_star_idx]
    ysamp = reform(ysamp)
    osamp = reform(osamp)

    yteff = teff_realizations[fitnum, young_star_idx]
    oteff = teff_realizations[fitnum, old_star_idx]
    yteff = reform(yteff)
    oteff = reform(oteff)


    yfunctargs = {veq_data:ysamp, $
                 teff_data:yteff, $
                 tcrit:tcrit}
    
    ofunctargs = {veq_data:osamp, $
                  teff_data:oteff, $
                 tcrit:tcrit}
    
             
    ; Loop through the tcrits and fit step functions
    foreach tcrit, teff_trial, teff_trial_idx do begin
        
        ; Assign the critical temperature
        yfunctargs.tcrit = tcrit
        ofunctargs.tcrit = tcrit
        
        ; Fit the velocities for the young pop
        yresult = mpfit('fit_stepfunction', parinfo=parinfo, $
                        functargs=yfunctargs, bestnorm=ychi, status=ystat, /quiet)

        oresult = mpfit('fit_stepfunction', parinfo=parinfo, $
                        functargs=ofunctargs, bestnorm=ochi, status=ostat, /quiet)

        ; Check and store results
        if n_elements(yresult) ne 1 then begin
            result_str[fitnum].veq_young_cool[teff_trial_idx] = yresult[0]
            result_str[fitnum].veq_young_hot[teff_trial_idx] = yresult[1]
            result_str[fitnum].chi2_young[teff_trial_idx] = ychi
            result_str[fitnum].status_young[teff_trial_idx] = ystat
        endif else begin
            result_str[fitnum].veq_young_cool[teff_trial_idx] = !values.d_nan
            result_str[fitnum].veq_young_hot[teff_trial_idx] = !values.d_nan
            result_str[fitnum].chi2_young[teff_trial_idx] = !values.d_nan
            result_str[fitnum].status_young[teff_trial_idx] = !values.d_nan
        endelse

        if n_elements(oresult) ne 1 then begin
            result_str[fitnum].veq_old_cool[teff_trial_idx] = oresult[0]
            result_str[fitnum].veq_old_hot[teff_trial_idx] = oresult[1]
            result_str[fitnum].chi2_old[teff_trial_idx] = ochi
            result_str[fitnum].status_old[teff_trial_idx] = ostat
        endif else begin
            result_str[fitnum].veq_old_cool[teff_trial_idx] = !values.d_nan
            result_str[fitnum].veq_old_hot[teff_trial_idx] = !values.d_nan
            result_str[fitnum].chi2_old[teff_trial_idx] = !values.d_nan
            result_str[fitnum].status_old[teff_trial_idx] = !values.d_nan
        endelse
    endforeach
    
        
    ; Determine best tcrits
    chi2_young = result_str[fitnum].chi2_young
    chi2_old = result_str[fitnum].chi2_old
    
    min_chi2_young = min(chi2_young, /nan)
    min_chi2_old = min(chi2_old, /nan)
    
    ;get the best of the young population
    best_young_idx = where(chi2_young eq min_chi2_young, youngchk)
    if youngchk gt 0 then begin
        best_idx = best_young_idx[0]
        result_str[fitnum].best_young_idx = best_idx
        result_str[fitnum].best_tcrit_young = teff_trial[best_idx]
        result_str[fitnum].best_veq_young_cool = $
          result_str[fitnum].veq_young_cool[best_idx]
        result_str[fitnum].best_veq_young_hot = $
          result_str[fitnum].veq_young_hot[best_idx]
    endif else begin
        result_str[fitnum].best_young_idx = !values.d_nan
        result_str[fitnum].best_tcrit_young = !values.d_nan
        result_str[fitnum].best_veq_young_cool = !values.d_nan
        result_str[fitnum].best_veq_young_hot = !values.d_nan
    endelse

    ;get the best of the old population
    best_old_idx = where(chi2_old eq min_chi2_old, oldchk)
    if oldchk gt 0 then begin
        best_idx = best_old_idx[0]
        result_str[fitnum].best_old_idx = best_idx
        result_str[fitnum].best_tcrit_old = teff_trial[best_idx]
        result_str[fitnum].best_veq_old_cool = $
          result_str[fitnum].veq_old_cool[best_idx]
        result_str[fitnum].best_veq_old_hot = $
          result_str[fitnum].veq_old_hot[best_idx]
    endif else begin
        result_str[fitnum].best_old_idx = !values.d_nan
        result_str[fitnum].best_tcrit_old = !values.d_nan
        result_str[fitnum].best_veq_old_cool = !values.d_nan
        result_str[fitnum].best_veq_old_hot = !values.d_nan
    endelse


    ;;; Also, bin by Teff and calculate mean Veq's
    young_bin_idx = value_locate(teff_bin_vec, yteff)
    old_bin_idx = value_locate(teff_bin_vec, oteff)

    for tbin = 0, nbins_teff-1 do begin
        
        yidx = where(young_bin_idx eq tbin, nyidx)
        oidx = where(old_bin_idx eq tbin, noidx)

        if nyidx ne 0 then begin
            result_str[fitnum].veq_bin_avg_young[tbin] = mean(ysamp[yidx])
            result_str[fitnum].veq_bin_stddev_young[tbin] = $
              total((mean(ysamp[yidx])-ysamp[yidx])^2)/(double(nyidx)-1d0)
            result_str[fitnum].veq_bin_stderr_young[tbin] = $
              result_str[fitnum].veq_bin_stddev_young[tbin]/sqrt(double(nyidx))
        endif else begin
            result_str[fitnum].veq_bin_avg_young[tbin] = !values.d_nan
            result_str[fitnum].veq_bin_stddev_young[tbin] = !values.d_nan
            result_str[fitnum].veq_bin_stderr_young[tbin] = !values.d_nan
        endelse 

        if noidx ne 0 then begin
            result_str[fitnum].veq_bin_avg_old[tbin] = mean(osamp[oidx])
            result_str[fitnum].veq_bin_stddev_old[tbin] = $
              total((mean(ysamp[oidx])-ysamp[oidx])^2)/(double(noidx)-1d0)
            result_str[fitnum].veq_bin_stderr_old[tbin] = $
              result_str[fitnum].veq_bin_stddev_old[tbin]/sqrt(double(noidx))
        endif else begin
            result_str[fitnum].veq_bin_avg_old[tbin] = !values.d_nan
            result_str[fitnum].veq_bin_stddev_old[tbin] = !values.d_nan
            result_str[fitnum].veq_bin_stderr_old[tbin] = !values.d_nan
        endelse 
        
    endfor

    print, fitnum
    ;dxstop

endfor




mwrfits, result_str, 'veq_vs_teff_stepfit_tgrid_results2.fits', /create
print, "Results output"


stop

;;; Start plotting
; Filled circle definition
usym = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points, 
; and set the filled flag:
USERSYM, COS(usym), SIN(usym), /FILL


v1ysort = sort(young_results[0,*])
v2ysort = sort(young_results[1,*])
tcysort = sort(young_results[2,*])
v1osort = sort(old_results[0,*])
v2osort = sort(old_results[1,*])
tcosort = sort(old_results[2,*])

v1y_lower_idx = v1ysort[1d3]
v1y_upper_idx = v1ysort[nsamples-1d3-1]


vbinsize = 0.05d0
tbinsize = 25d0

vbinmin = 0d0
tbinmin = 2500d0

test = histogram(young_results[0,*], binsize=vbinsize, min=vbinmin)
vhist_plotvec = dindgen(n_elements(test))*vbinsize + vbinmin
vhist_xvec = dindgen(n_elements(test))*vbinsize + vbinmin + (vbinsize/2)

plot, vhist_plotvec, test, ps=10

stop

end
