; Find a relation between vsini, teff and age

function age_probability, sample
common pdf_info, age_pdf, age_vec
; Give the probability for each element in sample

prob_vec = interpol(age_pdf, age_vec, sample)

return, prob_vec

end

pro vsini_teff_age

common pdf_info, age_pdf, age_vec

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

; Generate N samples for each star, based on their PDFs
nsamples = 1d5

age_realizations = dblarr(nsamples,nstars)

for starnum = 0, nstars-1 do begin

    ; PDF
    age_pdf = age_pdf_arr[*,starnum]

    ; Create a sample for this star
    sample_pop = create_sample_distribution('age_probability',minmax(age_vec),nsamples)

    ; Save the sample
    age_realizations[*,starnum] = sample_pop

endfor

; Now we will have N realizations of the data with age, teff and vsini
; for each star

; For each realization, 3d plot Vsini
; Xaxis is Teff, Yaxis is Age, Z is Vsini

for trial = 0, nsamples-1 do begin

    trial_age = reform(age_realizations[trial,*])

    p1 = plot3d(teff, trial_age, vsini, 'o', /sym_filled)

    dxstop

endfor

end
