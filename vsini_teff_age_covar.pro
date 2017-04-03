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

pro vsini_teff_age_covar

common pdf_info, age_pdf, age_vec

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
teff = info_str.teff
p_membership = info_str.thickthin_prob


; Read in gold-sample flags
flags = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/gold_sample.fits',0)
; Save flag indices
g_idx = where(flags eq 0, ngood, complement=b_idx, ncomplement=nbad)
napo = n_elements(flags) 
nstars = n_elements(g_idx) 

; Generate the age PDFS
stellar_age_estimator, age_str

age_vec = age_str.age
age_pdf_arr = age_str.prob

help, age_pdf_arr

; Use ASPCAP teff only, except for where they are undefined
;teff_asp = teff_asp[g_idx]
newg_idx = where(teff_asp gt 0 and teff_asp le 4100 and flags eq 0, nstars)
; modify vectors
uarr = uarr[newg_idx]
varr = uarr[newg_idx]
warr = uarr[newg_idx]
vsini = vsini[newg_idx]
teff = teff_asp[newg_idx]
age_pdf_arr = age_pdf_arr[*,newg_idx]

help, age_pdf_arr

p_membership = p_membership[newg_idx]


; Generate the sini pdf
inclination = dblarr(1001)/1000d0*!dpi/2d0
sini_dist = tan(inclination)

print, nstars
dxstop


; Generate N samples for each star, based on their age, temperature,
; and vsini pdfs
nsamples = 5d4

age_realizations = dblarr(nsamples,nstars)
teff_realizations = dblarr(nsamples, nstars)
veq_realizations = dblarr(nsamples, nstars)

;sini_realizations = dblarr(nsamples, nstars)
;vsini_realizations = dblarr(nsamples, nstars)

for starnum = 0, nstars-1 do begin

    ; PDF
    age_pdf = age_pdf_arr[*,starnum]

    ; Create a sample for this star
    sample_pop_age = create_sample_distribution('age_probability',minmax(age_vec),nsamples)
    ;sample_pop_teff =
    ;create_sample_distribution('teff_probability',minmax(teff),nsamples)
    ; Get a sample of realistic sini values 
    sample_pop_sini = create_sample_distribution('sini_probability',[0d0,1d0],nsamples)
    ;add 2km/s uncertainty on measured val
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

    teff_uncertainty = randomn(seed,nsamples)*70d0
    sample_pop_teff = teff[starnum] + teff_uncertainty



    ; Save the sample
    age_realizations[*,starnum] = sample_pop_age
    teff_realizations[*,starnum] = sample_pop_teff
    veq_realizations[*,starnum] = sample_pop_veq

    print, "Starnum: "+strtrim(starnum,2)+" finished"

endfor

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
device, filename = "vsini_teff_age_correlation_aspcap.eps"
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
