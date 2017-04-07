function hist_conf, vector, binsize=binsize, binmin=binmin, conf_bounds=conf_bounds
; Make histogram to determine mode, and confidence limits

if n_elements(conf_bounds) eq 0 then conf_bounds = [0.25,0.75]

hist = histogram(vector, binsize=binsize, min=binmin, /nan)

xvec = dindgen(n_elements(hist))*binsize + binmin + (binsize/2d0)

; Get confidence interval
pdf = hist/total(hist,/double)
cdf = total(pdf, /cum)

conf = interpol(xvec, cdf, conf_bounds)

; Get mode (this is not robust)
smoothpdf = gauss_smooth(pdf, 3, /edge_wrap)

mode_idx = where(smoothpdf eq max(smoothpdf))
mode = mean(xvec[mode_idx])

; output
params = [mode, conf]

return, params

end
