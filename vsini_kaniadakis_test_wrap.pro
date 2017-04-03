pro vsini_kaniadakis_test_wrap, _EXTRA=ex

std_err_vec = [1d-4, 1d-3, 1d-2, 1d-1, 1d0, 2d0, 3d0]

; these need to change if the defaults change in the underlying pro
if tag_exist(ex,'vsini_0') then vsini_0 = ex.vsini_0 else vsini_0 = 7d0
if tag_exist(ex,'kappa') then kappa = ex.kappa else kappa = 0.5d0
if tag_exist(ex,'sigma') then sigma = ex.sigma else sigma = 3d0

; loop
result_mean_arr = []
result_stddev_arr = []
for iter = 0, n_elements(std_err_vec)-1  do begin
    
    std_err_i = std_err_vec[iter]

    vsini_kaniadakis_test, std_err = std_err_i, result_mean=result_mean, result_stddev=result_stddev, _EXTRA=ex

    result_mean_arr = [[result_mean_arr],[result_mean]]
    result_stddev_arr = [[result_stddev_arr],[result_stddev]]


endfor

; Plot results
window, 0, xs = 1200, ys = 900
!p.multi = [0,1,3]

vsini_mean = result_mean_arr[0,*]
kappa_mean = result_mean_arr[1,*]
sigma_mean = result_mean_arr[2,*]

vsini_error = result_stddev_arr[0,*]
kappa_error = result_stddev_arr[1,*]
sigma_error = result_stddev_arr[2,*]

vsini_yr0 = min(vsini_mean)-(max(vsini_error)*1.5d0)
vsini_yr1 = max(vsini_mean)+(max(vsini_error)*1.5d0)

kappa_yr0 = min(kappa_mean)-(max(kappa_error)*1.5d0)
kappa_yr1 = max(kappa_mean)+(max(kappa_error)*1.5d0)

sigma_yr0 = min(sigma_mean)-(max(sigma_error)*1.5d0)
sigma_yr1 = max(sigma_mean)+(max(sigma_error)*1.5d0)


plot, std_err_vec, vsini_mean, yr=[vsini_yr0,vsini_yr1], ps=6, tit='Vsini result', ytit='vsini_ctr', charsize = 2.0, /xlog
oploterr, std_err_vec, vsini_mean, vsini_error
oplot, std_err_vec, replicate(vsini_0, n_elements(std_err_vec)), linest=2

plot, std_err_vec, kappa_mean, yr=[kappa_yr0,kappa_yr1], ps=6, tit='Kappa result', ytit='kappa', charsize = 2.0, /xlog
oploterr, std_err_vec, kappa_mean, kappa_error
oplot, std_err_vec, replicate(kappa, n_elements(std_err_vec)), linest=2

plot, std_err_vec, sigma_mean, yr=[sigma_yr0,sigma_yr1], ps=6, tit='Sigma result', ytit='sigma', xtit='Assumed Error in Vsini measurements', charsize = 2.0, /xlog
oploterr, std_err_vec, sigma_mean, sigma_error
oplot, std_err_vec, replicate(sigma, n_elements(std_err_vec)), linest=2

stop

end
