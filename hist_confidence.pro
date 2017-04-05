;+
; NAME: HIST_CONFIDENCE
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

pro hist_confidence, vector, min=binmin, max=binmax, binsize=binsize

if n_elements(binmin) eq 0 then binmin = min(vector, /nan)
if n_elements(binmax) eq 0 then binmax = max(vector, /nan)

; Estimate binsize
vsort = vector[sort(vector)]

dvector = vsort[1:*]-vsort[0:-2]

minstep = min(dvector)

maxstep = max(dvector)

medstep = median(dvector)

meanstep = mean(dvector)

stepsize = 7*meanstep
testhist = histogram(vector, binsize=stepsize, min=binmin)
pdf = testhist/total(testhist, /double)
cdf = total(pdf, /cum)

xvec = dindgen(n_elements(testhist))*stepsize+binmin
xxvec = dindgen(n_elements(testhist)*7)*stepsize/7d0+binmin

!p.multi = [0,1,2]
plot, xvec, pdf, ps=10
plot, xvec, cdf, ps=10
oplot, xxvec, interpol(cdf, xvec, xxvec), color=200

stop



end
