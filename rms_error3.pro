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


pro rms_error3

rms = 0
sigma = 1d0
iter = 0
npts = 1

; Result vectors
npts_vec = []
rms_vec = []
sig_vec = []

while npts lt 1000 do begin

    rms_sample = []
    ; Get 1000 rms values
    for run = 0, 999 do begin
        sample = randomn(seed, npts, /double) * sigma

        rms_sample_run = sqrt(mean(total(sample^2)))

        rms_sample = [rms_sample, rms_sample_run]
    endfor

    rms = mean(rms_sample)
    sig = stddev(rms_sample)

    rms_vec = [rms_vec, rms]
    npts_vec = [npts_vec, npts]
    sig_vec = [sig_vec, sig]
    ;sigma=sigma+0.5
    
    iter++
    npts++

endwhile



plot, npts_vec, sig_vec

stop



end
