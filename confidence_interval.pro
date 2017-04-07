;+
; NAME: CONFIDENCE_INTERVAL
;
; PURPOSE: Estimate confidence intervals for
;          arbitrary distributions
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-

function ci_optimize, lower_limit, data_vec=data_vec, $
                      cdf=cdf, ci=ci, xplotvec=xplotvec, $
                      plotpdf=plotpdf, visualize=visualize

bounds = lower_limit+[0d0,ci]

data_bounds = interpol(data_vec, cdf, bounds)

data_range = abs(data_bounds[1]-data_bounds[0])

if visualize then begin
    plot, data_vec, cdf, title = strtrim(data_range,2)
    oplot, replicate(data_bounds[0],2), [0,1], linest=2
    oplot, replicate(data_bounds[1],2), [0,1], linest=2
    oplot, minmax(data_vec), replicate(bounds[0],2), linest=2
    oplot, minmax(data_vec), replicate(bounds[1],2), linest=2
    
    plot, xplotvec, plotpdf, ps=10
    oplot, replicate(data_bounds[0],2), [0,1], linest=2
    oplot, replicate(data_bounds[1],2), [0,1], linest=2
    wait, 0.1
endif

return, data_range

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function confidence_interval, vector, min=binmin, $
                              conf_interval=conf_interval,$
                              visualize=visualize

if n_elements(visualize) eq 0 then visualize = 0
if n_elements(conf_interval) eq 0 then conf_interval=0.9

; Estimate binsize
vsort = vector[sort(vector)]

dvector = vsort[1:*]-vsort[0:-2]

minstep = min(dvector)
maxstep = max(dvector)
medstep = median(dvector)
meanstep = mean(dvector)

smallstep = 0.1*minstep
bigstep = 7*meanstep


; Construct the CDF from the data vector
n_data = n_elements(vector) 
cdf_data = (dindgen(n_data)+1)/n_data
cdf = [0d0, cdf_data]

vsort_prime = interpol(vsort, cdf_data, cdf)
if n_elements(binmin) eq 0 then begin
    binmin = vsort_prime[0]
endif 
    
; Estimate PDF for plotting
plothist = histogram(vector, binsize=bigstep, min=binmin)
plotpdf = plothist/total(plothist, /double)
xplotvec = dindgen(n_elements(plothist))*bigstep + binmin


; Search for the smallest interval in vector which corresponds to
; conf_interval
lims = [0d0,1d0-conf_interval]
guess = mean(lims)
parinfo = {value:guess, limited:[1,1], limits:lims, step:medstep}

functargs = {data_vec:vsort_prime, $
             cdf:cdf, $
             ci:conf_interval, $
             xplotvec:xplotvec, $
             plotpdf:plotpdf, $
             visualize:visualize}

!p.multi=[0,1,2]
lower_bound = mpfit('ci_optimize', parinfo=parinfo, $
                    functargs=functargs, status=status, /quiet)
!p.multi=0

if status gt 0 then begin
    final_bounds = lower_bound + [0d0,conf_interval]
    conf_range = interpol(vsort_prime, cdf, final_bounds)
endif else message, "Optimization of confidence iterval failed"

return, conf_range

end
