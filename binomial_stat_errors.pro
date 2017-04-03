;+
; NAME: binomial_stat_errors
;
; PURPOSE:
;    To calculate upper and lower limits for data subject 
;    to binomial statistics, given some confidence level
;    Based on equations 15 and 16 in Gehrels, 1986
;
; CATEGORY:
;    Statistics and uncertainty
;
; CALLING SEQUENCE:
;    limits = binomial_stat_errors(n, k, CL)
;
; INPUTS:
;    n - total number of events
;    k - number of 'successes'
;    CL - desired confidence level (percent)
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;    limits - 2-element array with [lower_limit, upper_limit]
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;    fit_info - contains input variables passed to fitting functions
;
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;-


; Calculate the confidence level, given an upper limit
; This is used to fit the prescribed confidence level and
; solve for the upper limit by optimization
function cl_upper, p_upper

common fit_info, n_in, k_in, cl_in

sum_vec = []

for x = 0, k_in do begin

    nchoosek_i = binomial_coefficient(n_in, x)
    term1 = p_upper^x
    term2 = (1d0 - p_upper)^(n_in-x)
    sum_i = nchoosek_i * term1 * term2
    ; save the summation term
    sum_vec = [sum_vec, sum_i]

endfor

; Sum up the terms
sum_tot = total(sum_vec, /double)
    
; Compute the difference between 1 - confidence level and the
; summation
diff = abs(sum_tot - (1d0-cl_in))

diff1 = diff + 1d0


; Add penalty if p_upper gt 1.0 or lt 0.0
if p_upper gt 1 then begin

    degree_bad = abs(p_upper - 1d0)
    penalty = 10d0*exp(10d0*degree_bad)

endif else if p_upper lt 0 then begin

    degree_bad = abs(p_upper)
    penalty = 10d0*exp(10d0*degree_bad)

endif else penalty = 0d0

chi = (diff1+penalty)/0.1d0

return, chi

end


; Calculate the confidence level, given an lower limit
; This is used to fit the prescribed confidence level and
; solve for the lower limit by optimization
function cl_lower, p_lower

common fit_info, n_in, k_in, cl_in

sum_vec = []

for x = 0, k_in-1 do begin

    nchoosek_i = binomial_coefficient(n_in, x)
    term1 = p_lower^x
    term2 = (1d0 - p_lower)^(n_in-x)
    sum_i = nchoosek_i * term1 * term2
    ; save the summation term
    sum_vec = [sum_vec, sum_i]

endfor

; Sum up the terms
sum_tot = total(sum_vec, /double)
    
;help, sum_vec
;print, sum_vec
;help, sum_tot
;print, sum_tot
;dxstop

; Compute the difference between the confidence level and the
; summation
diff = abs(sum_tot - cl_in)

diff1 = diff + 1d0



; Add penalty if p_lower gt 1.0 or lt 0.0
if p_lower gt 1 then begin

    degree_bad = abs(p_lower - 1d0)
    penalty = 10d0*exp(10d0*degree_bad)

endif else if p_lower lt 0 then begin

    degree_bad = abs(p_lower)
    penalty = 10d0*exp(10d0*degree_bad)

endif else penalty = 0d0

chi = (diff1+penalty)/0.1d0

return, chi

end


function binomial_stat_errors, n, k, conf_level

common fit_info, n_in, k_in, cl_in

; Check inputs (later)


; Prepare inputs for amoeba
n_in = n
k_in = k
cl_in = (conf_level/100d0)

ftol = 1d-8

initial_frac = double(k)/double(n)
if initial_frac eq 0 then begin
    guess_low = 0d0
    guess_high = 0.5d0
endif else begin

    guess_low = 0.5d0 * initial_frac
    guess_high = (1.5d0 * initial_frac) < 1d0
endelse

;scale_low = 0.5d0
;scale_high = 0.5d0
scale_low = 1d0
scale_high = 1d0


;;; Compute lower limit
; Check if k=0
if k eq 0 then begin
    
    lower_limit = 0d0

endif else begin
    iter = 0
    repeat begin

        low_guess = 0d0 > randomu(seed,1)*guess_low*3 < initial_frac

        ;help, low_guess
        ;print, low_guess[0]
        ;help, initial_frac


        result = amoeba3(ftol, function_name='cl_lower', p0=[low_guess], $
                         scale=[scale_low], nmax=10)

        ; check result
        ;if result[0] ne -1 then lower_limit = result[0] else message,
        ;"bad fit low"
        if result[0] ne -1 then begin
            ; Check that result works
            chi = cl_lower(result[0])
            ;print, chi
            if abs(chi-10d0) lt 1d-5 then terminate = 1 else terminate = 0
        endif else begin
            ;message, "bad fit low"
            print, "Bad fit low, trying to calculate upper limit for inverse scenario"

            nfail = n_in - k_in
            guess = [double(nfail/n_in)]
            scale = [0.5d0]
            
            k_in_copy = k_in
            k_in = nfail
            result = amoeba3(ftol, function_name='cl_upper', p0=guess, $
                         scale=scale, nmax=10)

            if result[0] ne -1 then begin
                ; Check result
                chi = cl_upper(result[0])
                
                if abs(chi-10d0) lt 1d-5 then begin
                    ; Take 1 - result to get the lower limit
                    result[0] = 1d0-result[0]
                    terminate = 1 
                endif else begin
                    ; Replace k_in with the old value before reiterating
                    k_in = k_in_copy
                    
                    print, "Bad fit low, reiterating... Iteration: "+strtrim(iter,2)
                    iter++
                    terminate = 0
                endelse
            endif else begin
                ; Replace k_in with the old value before reiterating
                k_in = k_in_copy
                ; Prepare to reiterate
                print, "Bad fit low, reiterating... Iteration: "+strtrim(iter,2)
                iter++
                terminate = 0
            endelse
                
        endelse
        if terminate eq 0 then begin
            print, "USING DUMB STATISTICS"
            nfail = n_in-k_in
            lower_lim = (1d0/n_in)*(double(k_in)-1.96d0*sqrt(nfail*k_in/double(n_in)))
            result[0] = lower_lim
            terminate = 1
            print, "result = " +strtrim(result,2)
        endif

        ;print, iter
        ;iter++
    endrep until terminate eq 1
    
    lower_limit = result[0]

endelse

; Compute upper limit
; check if k=n
if k eq n then begin
    
    upper_limit = 1d0

endif else begin

    iter = 0
    repeat begin
        
        high_guess = initial_frac > randomu(seed,1)*guess_high*3 < 1d0
        result = amoeba3(ftol, function_name='cl_upper', p0=[high_guess], $
                         scale=[scale_high], nmax=10)
        ; check result
        ;if result[0] ne -1 then upper_limit = result[0] else message,
        ;"bad fit high"
        if result[0] ne -1 then begin
            ; Check that result works
            chi = cl_upper(result[0])
            if abs(chi-10d0) lt 1d-5 then terminate = 1 else terminate = 0
            ;print, chi
        endif else begin
            ;message, "bad fit high"
            print, "Bad fit high, reiterating... Iteration: "+strtrim(iter,2)
            iter++
            terminate=0
        endelse
        
        if terminate eq 0 then begin
            print, "USING DUMB STATISTICS - UPPER"
            nfail = n_in-k_in
            upper_lim = (1d0/n_in)*(double(k_in)+1.96d0*sqrt(nfail*k_in/double(n_in)))
            result[0] = upper_lim
            terminate = 1
            print, "result = " +strtrim(result,2)
        endif

        ;print, iter
        ;iter++
    endrep until terminate eq 1
    
    upper_limit = result[0]
    


; result = amoeba3(ftol, function_name='cl_upper', p0=guess_high, scale=scale_high)
    
;     ; check result
;     if result[0] ne -1 then upper_limit = result[0] else message, "bad fit high"
    
endelse

; Return result of [lower_limit, upper_limit]
binomial_limits = [lower_limit, upper_limit]

return, binomial_limits

end
