;;; Return the RMS error based on the simulation from RMS_ERROR_SIMULATION

function rms_error_optimize, p

common amoeba_info, rms_data, npts_data, sig_vec_rough, rms_vec_rough

sigma = p[0]

; Do default 100 draws of npts_data at this sigma
ndraws = 100
rms_vec = []

for draw = 0, ndraws-1 do begin

    sample = randomn(seed, npts_data, /double) * sigma

    rms_sample_run = sqrt(mean(total(sample^2)))

    rms_vec = [rms_vec, rms_sample_run]

endfor

rms = mean(rms_vec)
;stderr = stddev(rms_vec)

residual = (rms-rms_data)^2

;plot, sig_vec_rough, abs(rms_vec_rough-rms_data)
;oplot, [sigma], [rms], ps=6
;wait, 0.1

return, residual

end


function rms_error0, rms_data, npts_data

common amoeba_info, rms_data1, npts_data1, sig_vec_rough, rms_vec_rough
rms_data1 = rms_data
npts_data1 = npts_data


rms = 0
sigma = 0.5
iter = 0

; Result vectors
sig_vec = []
rms_vec = []

while rms lt rms_data*2 do begin

    rms_sample = []
    ; Get 3 rms values
    for run = 0, 20 do begin
        sample = randomn(seed, npts_data, /double) * sigma

        rms_sample_run = sqrt(mean(total(sample^2)))

        rms_sample = [rms_sample, rms_sample_run]
    endfor

    rms = mean(rms_sample)

    rms_vec = [rms_vec, rms]

    sig_vec = [sig_vec, sigma]

    sigma=sigma+0.5
    iter++

endwhile

sig_vec_rough = sig_vec
rms_vec_rough = rms_vec

diff_vec = abs(rms_vec-rms_data)





sigma_min_rough_idx = where(diff_vec eq min(diff_vec))

sigma_guess = sig_vec[sigma_min_rough_idx]
sigma_scale = 0.5*sigma_guess
ftol = 0.9


print, "Sigma guess = " + strtrim(sigma_guess,2)
amoeba_full = tic()

result = amoeba3(ftol, function_name='rms_error_optimize', p0=sigma_guess, scale=sigma_scale, nmax=200)

toc, amoeba_full

; Plot result
;plot, sig_vec, diff_vec
;oplot, [result,result], [0,rms_data*2], linest=2


return, result

end
