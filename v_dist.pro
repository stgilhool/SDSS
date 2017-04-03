; Procedure to reconstruct underlying v-distribution

;;; FUNCTION: mcvdist_optimize
;;; mpfit calls this function, which does a monte carlo simulation
;;; given parameters that determine the model v distribution

;;; FUNC: model_v_generate
;;; generate a model v-distribution from input parameters

;;; FUNC: model_sini_generate
;;; generate a model sini population based on known sini distribution 

;;; FUNC: model_vsini_generate
;;; use model v and model sini distributions to make a model sini
;;; population

;;; Determine the likelihood of 

function v_samp_prob, v_samples, v_pdf

; The v_pdf is in bins of 1 km/s for something like 0-100 km/s
; Make a histogram of the input v sample
v_hist = histogram(v_samples, binsize = 1, min=0, omin = om, reverse_indices = ri)

nh = n_elements(v_hist) 
nsamp = n_elements(v_samples) 


v_prob = dblarr(nsamp)
; Loop through the ri vector, and recover which bin each v_sample went
; in, and give each v_sample the corresponding prob
for vbin = 0, nh-1 do begin

    bin_prob = v_pdf[vbin]

    if ri[vbin+1] ne ri[vbin] then begin
        vsamp_idx = ri[ri[vbin]:ri[vbin+1]-1]
        v_prob[vsamp_idx] = bin_prob
    endif
endfor

return, v_prob

end
    



function create_sample_distribution, pdf_function, pdf_xminmax, n_samples
; Create a sample distribution for a general pdf function
; INPUTS:
; pdf_function - string containing the name of a function
;                which will calculate the probability of
;                a random value
; pdf_xminmax  - 2-element vector which gives the range
;                of random values to generate
; n_samples    - Desired number of samples


; Initialize return variable
sample_pop = []

repeat begin

    ; create a bunch of random numbers between xmin and xmax
    pdf_xrange = pdf_xminmax[1] - pdf_xminmax[0]
    samp_gen = randomu(seed, 10e5) * pdf_xrange + pdf_xminmax[0]
    ; probability of each number is its sin
    prob_x = call_function(pdf_function, samp_gen)
    
    ; create more random numbers which will serve to check if
    ; we accept or reject the num_gen numbers
    chk_gen = randomu(seed, 10e5)
    ; accept or reject
    accept_idx = where(prob_x ge chk_gen, acpcnt)
    
    if acpcnt gt 0 then begin
        keepers = samp_gen[accept_idx]
        ; add the keepers to the sini_vec
        sample_pop = [sample_pop, keepers]
    endif
    
endrep until n_elements(sample_pop) ge n_samples

; Keep only the first n_samples
sample_pop = sample_pop[0:n_samples - 1]

return, sample_pop

end

function v_pdf_generate, params

if n_elements(params) eq 1350 then begin
    model_v_prelim = params
    ;Normalize
    model_v_dist = model_v_prelim/total(model_v_prelim)
endif else message, "not equipped to deal with alternate v-parameterizations"

return, model_v_dist
end

function v_pop_generate, params, n_samples

if n_params() eq 0 then begin
    print, "Calling sequence:"
    print, "result = v_pop_generate(params,n_samples)"
    message, "Syntax Error"
endif

; First, generate the v_pdf
v_pdf = v_pdf_generate(params)

; Now, draw samples from that distribution
v_pdf_func = 'v_samp_prob'
v_minmax = [0d0,100d0]

;v_pop = draw_samples(v_pdf_func, v_minmax, n_samples)
v_pop = create_sample_distribution(v_pdf_func, v_minmax, n_samples)

return, v_pop

end



function sini_pop_generate, vec_length

if n_params() eq 0 then begin
    print, "Calling sequence:"
    print, "result = sini_pop_generate(vec_length)"
    message, "Syntax Error"
endif

; Make tan(i) distribution of sini
sini_vec = []

repeat begin

    ; create a bunch of random numbers between 0 and pi/2
    i_gen = randomu(seed, 10e5) * (!pi/2.0)
    ; probability of each number is its sin
    prob_sini = tan(i_gen)
    
    ; create more random numbers which will serve to check if
    ; we accept or reject the num_gen numbers
    chk_gen = randomu(seed, 10e5)
    ; accept or reject
    accept_idx = where(prob_sini ge chk_gen, acpcnt)
    
    if acpcnt gt 0 then begin
        keepers = sin(i_gen[accept_idx])
        ; add the keepers to the sini_vec
        sini_vec = [sini_vec, keepers]
    endif
    
endrep until n_elements(sini_vec) ge vec_length

; Keep only the first vec_length entries
sini_vec = sini_vec[0:vec_length - 1]

return, sini_vec

end




function mcvdist_optimize, p, KEY=key

n_obs = 1350L
; Make trial underlying v-distribution
; don't need anymore, but can use as a check
v_dist_prob = v_pdf_generate(p)

; Make data histogram for comparison
binsize = 4
inmin = 0
inmax = 100
vsini_vfit_hist = histogram(data, $
                            binsize = binsize, min=inmin, max=inmax, $
                            reverse_indices=vfitri)

; Do MC simulation of many model populations
n_mc = 10e5
vsini_result = dblarr(n_elements(vsini_vfit_hist), n_mc)
for trial = 0, n_mc-1 do begin
    
    ; Generate a population of vsini's where p(sini) = tan(i)
    sini_pop = sini_pop_generate(n_obs)
    
    ; Generate a population of v's according to trial v_pdf
    v_pop = v_pop_generate(p,n_obs)
    
    ; Combine
    vsini_pop = v_pop * sini_pop
    
    ; Bin it up
    vsini_model_hist = histogram(vsini_pop, $
                                 binsize = binsize, min=inmin, max=inmax, $
                                 reverse_indices=vmodri)

    ; Store
    vsini_result[*,trial] = vsini_model_hist
endfor

; Compare result

end






; isrt = sort(i_vec)
; i_vec_srt = i_vec[isrt]

; i_hist = histogram(i_vec_srt, binsize = 0.01, omin=om)

; i_hist = i_hist/total(i_hist)

; i_vals = (dindgen(n_elements(i_hist))*0.01) + om

; plot, i_vals, i_hist, ps=10
; oplot, i_vals, sin(i_vals)/total(sin(i_vals)), color=200
; stop
; ; Check that sini is distributed as tani
; sini_vec = sin(i_vec_srt)
; sini_hist = histogram(sini_vec, binsize = 0.01, omin=om)
; sini_hist = sini_hist/total(sini_hist)
; sini_vals = asin((dindgen(n_elements(sini_hist))*0.01) + om)
; plot, sini_vals, sini_hist, ps=10
; tani = tan(sini_vals)
; oplot, sini_vals, tani/total(tani), color=200
; stop
; ; Another check
; ; cosi is uniformly distributed
; cosi_vec = randomu(seed, 10e6)
; i_acos_vec = acos(cosi_vec)
; i_acos_hist = histogram(i_acos_vec, binsize = 0.01, omin=om)

; i_acos_hist = i_acos_hist/total(i_acos_hist)

; i_acos_vals = (dindgen(n_elements(i_acos_hist))*0.01) + om

; plot, i_vals, i_hist, ps = 10, xr = [0, max(i_vals) > max(i_acos_vals)]
; oplot, i_acos_vals, i_acos_hist, color=200, ps = 10

; stop

