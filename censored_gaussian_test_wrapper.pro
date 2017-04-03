pro censored_gaussian_test_wrapper, _EXTRA=ex

;initialize
; mu_vec = []
; mu_err_vec = []
; sigma_vec = []
; sigma_err_vec = []

result_vec = []
thresh_vec = []

; THIS NEEDS TO CHANGE IF THE DEFAULTS OF THE OTHER PRO CHANGE
if tag_exist(ex,'mu') then mu_sim = ex.mu else mu_sim=10d0
if tag_exist(ex,'sigma') then sigma_sim = ex.sigma else sigma_sim=5d0
; if finite(ex.vis) then sigma_sim = ex.sigma else sigma_sim=5d0


for i = 0, 20 do begin

    thresh = i
    ;censored_gaussian_test, thresh=thresh, mu=mu_sim,
    ;sigma=sigma_sim, result=result

    ;censored_gaussian_test, thresh=thresh, result=result, _EXTRA=ex
    truncated_gaussian_test, thresh=thresh, result=result, _EXTRA=ex

    result_vec = [result_vec, result]
    thresh_vec = [thresh_vec, thresh]
;     mu_vec = [mu_vec, mu]
;     mu_err_vec = [mu_err_vec, mu_err]
;     sigma_vec = [sigma_vec, sigma]
;     sigma_err_vec = [sigma_err_vec, sigma_err]

endfor

mu_vec = result_vec.mu
mu_err_vec = result_vec.mu_err
sigma_vec = result_vec.sigma
sigma_err_vec = result_vec.sigma_err

window,0,xs=1200,ys=900

!p.multi = [0,1,2]
mu_yr0 = mean(mu_vec)-(max(mu_err_vec)*1.5d0)
mu_yr1 = mean(mu_vec)+(max(mu_err_vec)*1.5d0)

sigma_yr0 = mean(sigma_vec)-(max(sigma_err_vec)*1.5d0)
sigma_yr1 = mean(sigma_vec)+(max(sigma_err_vec)*1.5d0)

nthresh = n_elements(thresh_vec)
maxthresh = max(thresh_vec) 

plot, thresh_vec, mu_vec, ps=6, xr = [-1,maxthresh+1], $
  yr =[mu_yr0,mu_yr1], $
  title = 'mu parameter vs. upper limit threshold', $
  xtit='threshold', ytit='mu', charsize=1.5
oplot, thresh_vec, replicate(mu_sim, n_elements(thresh_vec)), linest=2
oploterr, thresh_vec, mu_vec, mu_err_vec

plot, thresh_vec, sigma_vec, ps=6, xr = [-1,maxthresh+1], $
  yr=[sigma_yr0,sigma_yr1], $
  title = 'sigma parameter vs. upper limit threshold', $
  xtit='threshold', ytit='sigma', charsize=1.5
oplot, thresh_vec, replicate(sigma_sim, n_elements(thresh_vec)), linest=2
oploterr, thresh_vec, sigma_vec, sigma_err_vec

stop

end
