; Procedure to calculate confidence limits using binomial statistics
; Using Gehrels, 1986

; Calculate the confidence level, given an upper limit
; This is used to fit the prescribed confidence level and
; solve for the upper limit by optimization
function cl_upper, p_upper, n_in=n_in, k_in=k_in, cl_in=cl_in

sum_vec = []

for x = 0, n_in do begin

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

chi = diff1/0.1d0

return, chi

end


; Calculate the confidence level, given an lower limit
; This is used to fit the prescribed confidence level and
; solve for the lower limit by optimization
function cl_lower, p_lower, n_in=n_in, k_in=k_in, cl_in=cl_in

sum_vec = []

for x = 0, n_in-1 do begin

    nchoosek_i = binomial_coefficient(n_in, x)
    term1 = p_lower^x
    term2 = (1d0 - p_lower)^(n_in-x)
    sum_i = nchoosek_i * term1 * term2
    ; save the summation term
    sum_vec = [sum_vec, sum_i]

endfor

; Sum up the terms
sum_tot = total(sum_vec, /double)
    
; Compute the difference between the confidence level and the
; summation
diff = abs(sum_tot - cl_in)

diff1 = diff + 1d0

chi = diff1/0.1d0

return, chi

end


function binomial_stat_errors_backup, n, k, conf_level

; Check inputs (later)


; Prepare inputs for mpfit
functargs = {n_in:n, $  
             k_in:k, $
             cl_in:conf_level}

;;; Compute lower limit
; Check if k=0
if k eq 0 then begin
    
    lower_limit = 0d0

endif else begin

    parinfo_low = {value:0.5d0, $
                   fixed:0, $
                   step:0.001d0, $
                   limited:[1,1], $
                   limits:[0d0,1d0]}

    
    result = mpfit('cl_lower', functargs=functargs, parinfo=parinfo_low, status=status)

    ; check result
    if status gt 0 then lower_limit = result[0]

endelse

; Compute upper limit
; check if k=n
if k eq n then begin
    
    upper_limit = 1d0

endif else begin

    parinfo_high = {value:0.5d0, $
                    fixed:0, $
                    step:0.001d0, $
                    limited:[1,1], $
                    limits:[0d0,1d0]}
    
    result = mpfit('cl_upper', functargs=functargs, parinfo=parinfo_high, status=status_high)

    ; check result
    if status_high gt 0 then upper_limit = result[0]

endelse

; Return result of [lower_limit, upper_limit]
binomial_limits = [lower_limit, upper_limit]

return, binomial_limits

end
