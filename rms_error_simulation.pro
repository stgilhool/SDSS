;;; Find the error for a given RMS by simulating Gaussians of
;;; different sigma and calculating their RMS's
;;; WHOOPS, there's an N-dependence, duh.  Okay, making it a function.

pro rms_error_simulation

rms = 0
sigma = 1
npoints = 1000
iter = 0

; Result vectors
sig_vec = []
rms_vec = []

while rms lt 1000 do begin

    rms_sample = []
    ; Get 3 rms values
    for run = 0, 2 do begin
        sample = randomn(seed, npoints, /double) * sigma

        rms_sample_run = sqrt(mean(total(sample^2)))

        rms_sample = [rms_sample, rms_sample_run]
    endfor

    rms = mean(rms_sample)

    rms_vec = [rms_vec, rms]

    sig_vec = [sig_vec, sigma]

    print, 'Iteration Number: ' + strtrim(iter,2) + ' completed.'
    print, 'Sigma    : ' + strtrim(sigma,2)
    print, 'RMS (all): ' + strtrim(rms_sample,2)
    print, 'RMS      : ' + strtrim(rms,2)
    print, ''

    sigma++
    iter++

endwhile
    
outstructure = {Sigma:sig_vec, RMS:rms_vec}
mwrfits, outstructure, 'rms_error_simulation_data.fits', /create
    

stop

end
