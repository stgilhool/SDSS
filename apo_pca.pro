;+
; NAME:
;
; PURPOSE:
;   To perform a priciple component analysis on APOGEE
;   M-dwarf spectra
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
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
;       Thu Feb 4 05:46:57 2016,
;       <stgilhool@iroquois.physics.upenn.edu>
;
;		
;
;-

pro apo_pca

; Readin spectra
nfil = 1000
apostr = readin_apo(/flux, nfiles=nfil)


; Continuum normalize

; find lowest common blue pixel
; find highest common blue pixel
minpix = lonarr(nfil)
maxpix = lonarr(nfil)
for i = 0, nfil - 1 do begin
    minpix[i] = fxpar(*(apostr[i].header), "BOVERMIN")
    maxpix[i] = fxpar(*(apostr[i].header), "BOVERMAX")
endfor

;fp = max(minpix)
lp = min(maxpix)
;fp = lp - 801L
fp = 2600L
np = lp-fp+1
;

; Truncate spectra
spectra = apostr.spectra[fp:lp, *]
norm = spectra ; initialize to same dims only
help, apostr.spectra
help, spectra
help, np

xvc = dindgen(np)

for i = 0, nfil-1 do begin
    norm[0,i] = continuum_fit(xvc, spectra[*,i])
    spectra[0,i] = spectra[*,i]/norm[*,i]
endfor

; Convert to optical depth space
medspec = median(spectra, dimension=2)
;tau = spectra
tau = -1d0*alog(spectra)

; remove mean
meantau = mean(tau, dimension=2, /nan)
meantau_arr = rebin(meantau, np, nfil)

;meantau = mean(tau, dimension=1, /nan)
;meantau_arr = rebin(reform(meantau, 1, nfil), np, nfil)
tau = tau-meantau_arr

taunan = where(finite(tau) ne 1, nancnt)
if nancnt gt 0 then begin
    tau[taunan] = 0
endif
;stop

;tau=transpose(tau)
;meantau_arr = transpose(meantau_arr)


; Perform PCA
result = pcomp(tau, coefficients=coeffs, /double, eigenvalues=evals, variances=variances)

;evecs = coeffs/rebin(sqrt(evals), np, np)
evecs = coeffs/rebin(sqrt(evals), np, np)

;print, evecs[*,0]
testdata = evecs ## transpose(tau)
testdata = transpose(testdata)

r = testdata + meantau_arr

testresult = result/rebin(sqrt(transpose(evals)), np, nfil)
tr = testresult+meantau_arr

;evecs=evecs + rebin(meantau, np, np)


!p.multi=[0,1,4]
window, 0, xs=1500, ys=800
for i=0,nfil/3 - 1 do begin
plot, medspec, yr =[0.5,1.1], /ys, title = 'Median'
;plot, exp(-1d0*r[*,3*i]),yr =[0.5,1.1], /ys, title = strtrim(3*i,2)
;plot, exp(-1d0*r[*,3*i+1]),yr =[0.5,1.1], /ys, title = strtrim(3*i+1,2)
;plot, exp(-1d0*r[*,3*i+2]),yr =[0.5,1.1], /ys, title = strtrim(3*i+2,2)
;plot, evecs[*,3*i],yr =[0.5,1.1], /ys, title = strtrim(3*i,2)
;plot, evecs[*,3*i+1],yr =[0.5,1.1], /ys, title = strtrim(3*i+1,2)
;plot, evecs[*,3*i+2],yr =[0.5,1.1], /ys, title = strtrim(3*i+2,2
plot, exp(-1d0*evecs[*,3*i]),yr =[0.5,1.1], /ys, title = strtrim(3*i,2)
plot, exp(-1d0*evecs[*,3*i+1]),yr =[0.5,1.1], /ys, title = strtrim(3*i+1,2)
plot, exp(-1d0*evecs[*,3*i+2]),yr =[0.5,1.1], /ys, title = strtrim(3*i+2,2)
dxstop
endfor
!p.multi=0
stop

end
