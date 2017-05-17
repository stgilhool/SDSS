; Testing cross correlation techniques to determine Vsini
; FT(obs) = FT(template*broadening kernel) = FFT(template) x FFT(broadening kernel)
function make_broadening_kernel, vsini, OVERSAMP=oversamp, _EXTRA=ex
;;; ROTATION
if n_elements(oversamp) eq 0 then oversamp=77
if n_elements(ex) eq 0 then epsilon=0.6 

bigc = 3d5 ;km/s
dlogwl = 6d-6
; change in velocity equal to 1 oversampled logwl pixel
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc

;help, deltav

; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid, epsilon=epsilon)
;help, velgrid
; Convert velocity-space x-axis to logwl-space
velgrid_kernel_x = alog10(1d0 + (velgrid/bigc))
;help, velgrid_kernel_x
; Make a similar vector which matches the logwl-grid
nlsfrot = n_elements(velgrid)+6 ;pad to make sure we get the whole kernel
lsf_rot_kernel_x = (dindgen(nlsfrot)-(nlsfrot/2)) * dlogwl/oversamp
lsf_x = (dindgen(nlsfrot)-(nlsfrot/2))/oversamp
;help, lsf_rot_kernel_x
; Spline lsf_rot onto the logwl_grid vector
lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x, /spline)

; Make sure interpolation doesn't introduce negatives
rot_neg_idx = where(lsf_rot lt 0d0, ncnt)
if ncnt gt 0 then lsf_rot[rot_neg_idx]=0d0
; Normalize

lsf_rot = lsf_rot/total(lsf_rot, /double)

; Recover velgrid
nvelgrid = n_elements(velgrid) 
lsf_velgrid_x = (dindgen(nvelgrid)-(nvelgrid/2))/oversamp
velgrid_out = interpol(velgrid, lsf_velgrid_x, lsf_x)

return, {lsf_rot:lsf_rot, kernel_logwl:lsf_rot_kernel_x, kernel_x:lsf_x, velgrid:velgrid_out}

end



pro vsini_correlate_test_fft

bigc = 3d5
d_logwl = 6d-6

; Load the APOGEE spectra
files = [289,1142]
apg_info = readin_apo2(files=files, hdu=8)

xmin = 1300L
xmax = 2300L

aidx = 0
wl_grid_full = apg_info[aidx].wl_grid
wl_grid = apg_info[aidx].wl_grid[xmin:xmax]

;Oversample
oversamp = 11L
npix = n_elements(wl_grid) 
npix_over = npix*oversamp
x = dindgen(npix)
xx = (dindgen(npix_over)-(oversamp/2))/oversamp
wl_grid_over = interpol(wl_grid, x, xx)

xnorm = x-npix/2
xxnorm = xx-npix_over/(2*oversamp)

; apo_spec = apg_info[0].spectra[xmin:xmax]
; apo_spec_over = interpol(apo_spec, wl_grid, wl_grid_over)
; aso_norm = sqrt(total(apo_spec_over^2))

; Make APOGEE LSF
lsfco = *apg_info[aidx].output
npix_lsf = 51L * oversamp
lsfx = (dindgen(npix_lsf)-(npix_lsf/2))/oversamp
lsf0 = (xmax-xmin)/2 + xmin
lsf = lsf_gh(lsfx+lsf0,lsf0,lsfco)
lsf = lsf/total(lsf,/double)


; Lag vectors
lag_over = lindgen(npix_over) - npix_over/2

vel_x_over = bigc*xxnorm*d_logwl

; Load the PHOENIX template
pho_str = readin_phoenix([4200],['+0.5'],['4.5'],wl_grid)
pho = pho_str.spec
pho_over_str = readin_phoenix([4200],['+0.5'],['4.5'],wl_grid_over)
pho_over = pho_over_str.spec

;pho_over_convol = convol(pho_over, lsf, /edge_wrap)
pho_over_convol = pho_over


; change APG spec to template with lsf and rot_lsf
g1_kern = make_broadening_kernel(60d0, oversamp=oversamp)
g1 = g1_kern.lsf_rot
g1 = g1/total(g1,/double)

g1_mean0 = (g1-mean(g1))
g1_normed = g1_mean0/(sqrt(total(g1_mean0^2,/double)))

g1_velgrid = g1_kern.velgrid

; put G1 on the same x scale
g1_compare = interpol(g1_normed, g1_velgrid, vel_x_over)

pho_mean0 = pho_over-mean(pho_over)
pho_normed = pho_mean0/sqrt(total(pho_mean0^2, /double))

;apo_spec_over = convol(pho_over_convol, g1, /edge_wrap)
;apo_spec_over = convol(pho_over, g1, /edge_wrap)
apo_spec_over = convol(pho_mean0, g1_mean0, /edge_wrap)
apo_spec_over_2 = c_correlate(g1_compare, pho_over, lag_over)
aso_norm = sqrt(total(apo_spec_over^2))

plot, c_correlate(g1_compare, pho_over, lag_over)
oplot,c_correlate(pho_over,g1_compare, lag_over), color=200
dxstops
plot, vel_x_over, apo_spec_over
oplot, vel_x_over, apo_spec_over_2, color=200
dxstop


plot, wl_grid_over, apo_spec_over
oplot, wl_grid_over, pho_over, color=200
dxstop




; Take FT of apo_spec_over, pho_over and g1



ft_obs = fft(apo_spec_over, /double, /center)
ft_pho = fft(pho_normed, /double, /center)
ft_g1 = fft(g1_compare, /double, /center)

plot, ft_obs
help, ft_obs
dxstop
plot, ft_pho
help, ft_pho
dxstop
plot, ft_g1
help, ft_g1
dxstop
plot, ft_obs
oplot, ft_pho*ft_g1, color=200
dxstop


dxstop

stop



end
