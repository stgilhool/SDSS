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

pro hist_confidence, vector, min=binmin, max=binmax, conf_interval=conf_interval

if n_elements(binmin) eq 0 then binmin = min(vector, /nan)
if n_elements(binmax) eq 0 then binmax = max(vector, /nan)

; Estimate binsize
vsort = vector[sort(vector)]

dvector = vsort[1:*]-vsort[0:-2]

minstep = min(dvector)
maxstep = max(dvector)
medstep = median(dvector)
meanstep = mean(dvector)

smallstep = 0.1*minstep
bigstep = 7*meanstep

; Estimate PDF and CDF
hist = histogram(vector, binsize=smallstep, min=binmin)

pdf = testhist/total(testhist, /double)
cdf = total(pdf, /cum)

; xvec with points in the middle of bins
xvec = dindgen(n_elements(hist))*smallstep + binmin + (smallstep/2d0)

; make sure the cdf and xvec are bounded
cdf = [0d0, cdf]
xvec = [binmin, xvec]

plothist = histogram(vector, binsize=bigstep, min=binmin)
xplotvec = dindgen(n_elements(plothist))*bigstep + binmin


!p.multi = [0,1,2]
plot, xplotvec, plothist, ps=10
plot, xvec, cdf
oplot, xxvec, interpol(cdf, xvec, xxvec), color=200

stop



end
