function param_fit, p

common amoeba_info, data, nnondet, thresh, xvec, full_data, vis, quiet

mu = p[0]
sigma = p[1]




ndata = n_elements(full_data) 

 ; make a gaussian pdf
delta_x = xvec[1]-xvec[0]
 gaussvec = gaussian(xvec, [1d0,mu,sigma])
 gaussnorm = total(gaussvec)*delta_x
 gaussvec = gaussvec/gaussnorm

; Make CDF
cdf = total(gaussvec, /cumulative)*delta_x

; renormalization const is 1/(cdf_xmax-cdf_xmin)
;cdfmax = double(interpol(cdf, xvec, max(xvec)))
;cdfmin = double(interpol(cdf, xvec, thresh))
;renorm = (cdfmax-cdfmin)^(-1d0)


; help, cdfmax
; print, cdfmax
; help, cdfmin
; print, cdfmin
; help, renorm
; print, renorm

 ; integrate the pdf up to 3
 ;nondet_idx = where(xvec lt thresh)
 ;nondet_gauss = gaussvec[nondet_idx]
 ;nondet_prob = total(nondet_gauss)*delta_x


;;; Calculate 1-S(x) = F(thresh) for the censored guys
censored_prob = double(interpol(cdf, xvec, thresh))


;bigl_vec =
;renorm*(1d0/(sqrt(2d0*!pi)*sigma))*exp(-1d0*((data-mu)^2d0)/(2d0*sigma^2d0))

; Likelihood for the detections
bigl_vec = (1d0/(sqrt(2d0*!pi)*sigma))*exp(-1d0*((data-mu)^2d0)/(2d0*sigma^2d0))

; Add the likelihood for nnondet non-detections
if nnondet gt 0 then begin
    bigl_nondet = replicate(censored_prob, nnondet)
    bigl_vec = [bigl_nondet, bigl_vec]
endif

 ;bigl_vec = [nondet_prob, bigl_vec]

 logl_vec = alog(bigl_vec)

 logl = total(logl_vec)
;likelihood = exp(logl)

;logl = 1d2 - total(((data-mu)^2)/(2d0*(sigma^2)))
;logl = 1d2 - total(((data_norm-mu)^2)/(2d0*(sigma^2)))
;  realmu = mean(data)
; realsig = stddev(data)

; print, mu, realmu
; print, sigma, realsig
; print, logl
; print, likelihood
; print, ''
; wait, 0.05


if vis then begin
    data_hist = histogram(full_data, min=min(xvec), max=max(xvec), binsize=delta_x)
    plot, xvec, data_hist, ps=10, xr=[-5,25], yr=[0,10]
    oplot, xvec, gaussvec*ndata, color=200
    oplot, [thresh,thresh], [!y.crange[0],!y.crange[1]], linest=2, color=99999
    wait, 0.0001
endif

return, -1d0*logl

end


pro censored_gaussian_test, nsamp=nsamp, mu=mu, sigma=sigma, thresh=thresh, vis=vis, quiet=quiet, result=outstr

if n_elements(nsamp) eq 0 then nsamp = 1d2
if n_elements(mu) eq 0 then mu = 10d0
if n_elements(sigma) eq 0 then sigma = 5d0
if n_elements(thresh) eq 0 then thresh = 3d0
if n_elements(vis) eq 0 then vis=0
if n_elements(quiet) eq 0 then quiet=1


common amoeba_info, data, nnondet, thresh1, xvec, full_data, vis1, quiet1
thresh1 = thresh
vis1 = vis
quiet1 = quiet



; Make Gaussian distribution with sigma=5 and mu=10
gdist = randomn(seed, 1d5)*sigma + mu
;gdist = randomn(seed, 1d5)

; Make an xvector for the gaussian
 npts = 1001
 delta_x = 0.1d0

 xvec = (dindgen(npts)-(npts/2)) * delta_x

mu_vec = dblarr(50)
sig_vec = dblarr(50)

for loop = 0, 49 do begin
    ; Draw samples
    samp_idx = long(randomu(seed, nsamp-1)*1d5 - 1)
        
    full_data = gdist[samp_idx]
    
    ; sample_dist = histogram(full_data, binsize=0.25, omin=om)

    ;     sample_x = findgen(n_elements(sample_dist))*0.25 + om
    
    ;     plot, sample_x, sample_dist, ps=10
    ;     wait, 0.5
    
    ;data = gdist[0:1000]
    ;data = gdist[0:50]
    ;full_data= gdist[0:1d4]
    
    ; test = call_function('param_fit', [1d0, 1d0])
    ; test = call_function('param_fit', [-1d0, 1d0])
    ; test = call_function('param_fit', [0d0, 1.5d0])
    ; test = call_function('param_fit', [0d0, 0.5d0])
    ; test = call_function('param_fit', [0d0, 1d0])
    ;  stop
    
    ;binsize = 0.5d0
    ;data = histogram(gdist, binsize = binsize, omin=om)
    ;data_norm = data/(total(data)*binsize)
    ;data_x = findgen(n_elements(data))*binsize + om
    ;plot, xvec, h, ps=10
    ;stop
    
    ;data = data/total(data)
    ; fit for sigma and mu
    ; parinfo_samp = {parname:'', $
    ;                 value:0d0, $
    ;                 limited:[0,0], $
    ;                 limits:[0d0,0d0], $
    ;                 fixed:0}
    
    ; parinfo = replicate(parinfo_samp, 2)
    
    ; parinfo[0].parname = 'mu'
    ; parinfo[0].value = 10d0
    
    ; parinfo[1].parname = 'sigma'
    ; parinfo[1].value = 5d0
    
    ; ; Functargs
    ; functargs = {data:data}
    
    ;result = mpfit('param_fit', parinfo = parinfo, functargs=functargs, $
    ;               status=status, bestnorm=neg_likelihood)
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;    Treat data points < 3 as non-detections
    
    
    det_idx = where(full_data ge thresh, ndet, complement=nondet_idx, ncomplement=nnondet)
    
    data = full_data[det_idx]
    
    ;data=full_data
    
    
    
    
    mu_randomfactor = double(randomu(seed,1))*0.1d0*mu
    sigma_randomfactor = double(randomu(seed,1))*0.1d0*sigma


    guess = [mu+mu_randomfactor,sigma+sigma_randomfactor]
    scale = [mu_randomfactor*3d0, sigma_randomfactor*3d0]
    
    result = amoeba3(1d-10, function_name = 'param_fit', $
                     p0 = guess, scale = scale)
    
    if quiet eq 0 then print, result
    
    
    if n_elements(result) gt 1 then begin
        mu_vec[loop]=result[0]
        sig_vec[loop]=result[1]
    endif else begin
        mu_vec[loop]=!values.d_nan
        sig_vec[loop]=!values.d_nan
    endelse
endfor

mu_out = mean(mu_vec, /nan)
emu_out = stddev(mu_vec, /nan)

sig_out = mean(sig_vec, /nan)
esig_out = stddev(sig_vec, /nan)

outstr = {mu:mu_out, $
          mu_err:emu_out, $
          sigma:sig_out, $
          sigma_err:esig_out}

print, "Final result is: "
print, "Mean mu: " + strtrim(mean(mu_vec, /nan), 2)
print, "Mean sigma: " + strtrim(mean(sig_vec, /nan), 2)

;    stop
    
end
