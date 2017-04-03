;;; Return the RMS error based on a simulation
;;; First, we find a rough value of a Gaussian sigma which yields the
;;; RMS of the data for a given number of data points, N.  Then, we
;;; optimize sigma using AMOEBA to get the exact value of sigma which
;;; corresponds to the data.  Finally, we take lots of draws from that
;;; Gaussian, and calculate the RMS for each draw.  The standard
;;; deviation of that set of RMS's is our simulated error in the RMS


; Return the vector of RMS's for npoints, sigma, and ntrials
function rms_error_simulate, sigma, npts, ntrials=ntrials

if n_elements(ntrials) eq 0 then ntrials = 1000

; Make an array of random numbers with sigma of sigma
sample = randomn(seed, ntrials, npts, /double) * sigma

; Get the RMS of each column (each trial)
rms_sample = sqrt(mean(sample^2, dimension=2))

;Return the vector of RMS's
return, rms_sample

; rms_sample = []
; ; Get 1000 rms values
; for run = 0, 999 do begin
;     sample = randomn(seed, npts, /double) * sigma
    
;     rms_sample_run = sqrt(mean(total(sample^2)))
    
;     rms_sample = [rms_sample, rms_sample_run]
; endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Use AMOEBA to optimize the Gaussian sigma for a given RMS and number
; of data points
function rms_error_optimize, p

common amoeba_info, rms_data, npts_data, sig_vec_rough, rms_vec_rough, vis

sigma = p[0]

; Do default 100 draws of npts_data at this sigma
ndraws = 100
; rms_vec = []

; for draw = 0, ndraws-1 do begin

;     sample = randomn(seed, npts_data, /double) * sigma

;     rms_sample_run = sqrt(mean(total(sample^2)))

;     rms_vec = [rms_vec, rms_sample_run]

; endfor

rms_vec = rms_error_simulate(sigma, npts_data, ntrials=ndraws)

rms = mean(rms_vec)
;stderr = stddev(rms_vec)

residual = (rms-rms_data)^2

; Overplot the fitting on the simulated results from before
if vis then begin
    plot, sig_vec_rough, abs(rms_vec_rough-rms_data)
    oplot, [sigma], [abs(rms-rms_data)], ps=6
    wait, 0.1
endif

return, residual

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;; MAIN ;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function rms_error, rms_data, npts_data, vis=vis

if n_elements(vis) eq 0 then vis = 0

common amoeba_info, rms_data1, npts_data1, sig_vec_rough, rms_vec_rough, vis1
vis1 = vis
rms_data1 = rms_data
npts_data1 = npts_data


; Initialize values for the while loop
rms = 0
sigma = 1
iter = 0

; Result vectors
sigma_vec = []
rms_vec = []


while rms lt rms_data*2 do begin

    rms_sample_vec = rms_error_simulate(sigma, npts_data1, ntrials=25)

    rms = mean(rms_sample_vec)

    rms_vec = [rms_vec, rms]

    sigma_vec = [sigma_vec, sigma]

    sigma=sigma+1
    iter++

endwhile

; Save vectors into common block for plotting later
sig_vec_rough = sigma_vec
rms_vec_rough = rms_vec

; Make a vector whose minimum corresponds to the best sigma for given
; RMS, npts
diff_vec = abs(rms_vec-rms_data)
; Get the index where sigma is optimized
sigma_min_rough_idx = where(diff_vec eq min(diff_vec))

; Use that index to find the best sigma, roughly, and initialize
; amoeba with guess, scale and ftol
sigma_guess = sigma_vec[sigma_min_rough_idx]
sigma_scale = 0.5*sigma_guess
ftol = 0.9

; Give some output to user
print, "Sigma guess = " + strtrim(sigma_guess,2)
amoeba_full = tic()

; Perform AMOEBA optimization of sigma
result = amoeba3(ftol, function_name='rms_error_optimize', p0=sigma_guess, scale=sigma_scale, nmax=200)
; Check the result and save it
if result[0] eq -1 then begin
    final_sigma = sigma_guess
    final_sigma = rms_data
    ;message, "AMOEBA failed to converge"
endif else begin
    final_sigma = result[0]
endelse

final_sigma = rms_data

toc, amoeba_full

;;; Finally, use this optimized sigma to draw 1000 samples and compute
;;; their RMS's

final_rms_vec = rms_error_simulate(final_sigma, npts_data, ntrials=1000)
rms_error = stddev(final_rms_vec)



; Plot result
if vis then begin
    plot, sigma_vec, diff_vec, title = 'RMS = '+strtrim(rms_data,2)+' | Npts = '+strtrim(npts_data,2)+ $
      ' | RMS error = '+strtrim(rms_error,2), charsize = 2, xtit = 'Gaussian Sigma', $
      ytit = 'Difference between RMS (data) and simulated RMS'
    oplot, [final_sigma,final_sigma], [0,rms_data*2], linest=2
endif    

return, rms_error

end
