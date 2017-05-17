; Testing the preliminary steps in the Vsini correlation approach
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro vsini_correlate_test

;load rainbow color table with discrete color tags	      
loadct,39,/silent
setcolors,/system_variables,/silent,decomposed=0  ;requires setcolors.pro

;plotting parameters to set axes and fonts
!p.background=0    ;default background color
!p.color=255
!p.charsize=1.7		;text default size
!p.charthick=6		;text default thickness
!x.thick=5		;thicken x-axis 
!y.thick=5		;thicken y-axis
!p.font=1		;set default font
!p.thick=5		;set default plotting line thickness
!x.margin=[7,4]

;set psym=8 to circle by default
circind=findgen(40) * (!pi*2/39.)
defsysv,'!circ',transpose( [[cos(circind)],[sin(circind)]])  ;user symbol vertices
usersym,!circ(0,*),!circ(1,*),/fill  ;circle plot




bigc = 3d5
d_logwl = 6d-6

; Load the APOGEE spectra
files = [326, 793, 397,289,1142]
apg_info = readin_apo2(files=files, hdu=8)

xmin = 350L
xmax = 1300L

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

apo_spec = apg_info[aidx].spectra[xmin:xmax]
apo_spec_over = interpol(apo_spec, wl_grid, wl_grid_over)
aso_norm = sqrt(total(apo_spec_over^2))

; Make APOGEE LSF
lsfco = *apg_info[aidx].output
npix_lsf = 51L * oversamp
lsfx = (dindgen(npix_lsf)-(npix_lsf/2))/oversamp
lsf0 = (xmax-xmin)/2 + xmin
lsf = lsf_gh(lsfx+lsf0,lsf0,lsfco)
lsf = lsf/total(lsf,/double)


apo_spec_over = convol(apo_spec_over, lsf, /edge_wrap)

; Load the PHOENIX template
pho_str = readin_phoenix([2800],['-0.0'],['4.5'],wl_grid)
pho = pho_str.spec
pho_over_str = readin_phoenix([2800],['-0.0'],['4.5'],wl_grid_over)
pho_over = pho_over_str.spec

pho_over_convol = convol(pho_over, lsf, /edge_wrap)
;pho_over_convol = pho_over
;pho_over = pho_over_convol




;pho_over_convol = pho_over_convol+0.05


poc_norm = sqrt(total(pho_over_convol^2))

; change APG spec to template with lsf and rot_lsf
vsini_known_vec = [60d0, 30d0, 25d0, 20d0, 15d0, 10d0, 8d0, 5d0]

vsini_known_vec = [32.25d0]

g_results = dblarr(n_elements(xxnorm),n_elements(vsini_known_vec))



for i = 0, n_elements(vsini_known_vec) -1 do begin
vsini_known = vsini_known_vec[i]

vsini_estimate = vsini_known
g1_kern = make_broadening_kernel(vsini_estimate, oversamp=oversamp)
g1 = g1_kern.lsf_rot
g1 = g1/total(g1,/double)

g1_velgrid = g1_kern.velgrid


;apo_spec_over = convol(pho_over_convol, g1, /edge_wrap)
;signal_to_noise = 50d0
;noise = randomn(seed, n_elements(apo_spec_over))/signal_to_noise
;apo_spec_over = apo_spec_over + noise
;aso_norm = sqrt(total(apo_spec_over^2))



; Plot Fig1.1
psopen, 'vcorr_test_spectra.eps', /encapsulated, xs=10, ys=10, /inches, /color

plot, wl_grid_over, apo_spec_over+0.3, yr=[0.5,1.5]
oplot, wl_grid_over, pho_over_convol

psclose

;plot, wl_grid_over, apo_spec_over
;oplot, wl_grid_over, pho_over, color=200
;dxstop
;plot, wl_grid_over, -1d0*alog(apo_spec_over), tit="Optical Depth"
;oplot, wl_grid_over, -1d0*alog(pho_over_convol), color=200
;dxstop
plot, wl_grid_over, apo_spec_over, thick=4, tit='APG spec and Template spec, both with LSF'
oplot, wl_grid_over, pho_over_convol, thick=3, color=!red
dxstop



; Make CCF_TT
lag = xnorm
lag_over = lindgen(npix_over) - npix_over/2

vel_x = bigc*xnorm*d_logwl
;vel_x_over = bigc*xxnorm*d_logwl
vel_x_over = (10d0^(xxnorm*d_logwl)-1d0)*bigc

;ccf_tt = a_correlate(pho, lag, /double)
ccf_tt = c_correlate(pho, pho, lag, /double)
ccf_tt_over = a_correlate(pho_over, lag_over, /double)
;ccf_tt_over = a_correlate(pho_over_convol, lag_over, /double)
;ccf_tt_over = c_correlate(pho_over_convol,pho_over_convol, lag_over, /double)
;ccf_tt_over = convol(pho_over_convol,pho_over_convol,poc_norm^2, /edge_wrap)

;ccf_tt_over = ccf_tt_over/(poc_norm^2)

; shift up
;if min(ccf_tt_over) lt 0 then ccf_tt_over = ccf_tt_over + abs(min(ccf_tt_over))


; get derivative in order to locate the central peak
delta_vel_over = vel_x_over[1] - vel_x_over[0]
vel_over_plus = vel_x_over + (0.5d0*delta_vel_over)
vel_over_minus = vel_x_over - (0.5d0*delta_vel_over)

ccf_tt_plus = interpol(ccf_tt_over, vel_x_over, vel_over_plus)
ccf_tt_minus = interpol(ccf_tt_over, vel_x_over, vel_over_minus)

ccf_tt_deriv = (ccf_tt_plus - ccf_tt_minus)/delta_vel_over

zero_idx = where(abs(ccf_tt_deriv) le 1d-3*max(ccf_tt_over), nzero)

vel_zeroes = vel_x_over[zero_idx]
; plot the ccf with lines at the minima and maxima from zeroes of the
; derivative
;plot, vel_x_over, ccf_tt_over, xr=[-50,50], ps=6, symsize=0.2
plot, vel_x_over, ccf_tt_over, thick=3, /yno, tit='CCF_TT with min/maxima'
dxstop
for z = 0, n_elements(zero_idx)-1 do begin
    oplot, replicate(vel_zeroes[z],2), [-100,100], thick=0.5, linesty=2
endfor
dxstop
; get the velocity of the first minima
first_vel_zero_idx = (where(vel_zeroes gt 5))[0]
first_vel_zero = vel_zeroes[first_vel_zero_idx]

; other approach
first_vel_zero_idx_plus = zero_idx[first_vel_zero_idx]
ccf_x0 = npix_over/2
ccf_peak_width = 2L*(first_vel_zero_idx_plus-ccf_x0)+1L

first_vel_zero_idx_minus = first_vel_zero_idx_plus-ccf_peak_width
first_vel_zero_minus = vel_x_over[first_vel_zero_idx_minus]
first_vel_zero_plus = vel_x_over[first_vel_zero_idx_plus]


; plot to test
plot, vel_x_over, ccf_tt_over, xr=[-60,60], ps=6, symsize=0.2, tit="Zoomed CCF_TT with boundaries of central peak marked"
oplot, replicate(first_vel_zero,2), [-100,100], linest=2, thick=2
oplot, replicate(first_vel_zero_minus,2), [-100,100], linest=2, color=200, thick=2
oplot, replicate(first_vel_zero_plus,2), [-100,100], linest=2, color=200, thick=2


dxstop

; Remove first peak
ccf_tt_peak = replicate(0d0, npix_over)
ccf_tt_peak[first_vel_zero_idx_minus:first_vel_zero_idx_plus] = $
  ccf_tt_over[first_vel_zero_idx_minus:first_vel_zero_idx_plus] - $
  mean(ccf_tt_over[[first_vel_zero_idx_minus,first_vel_zero_idx_plus]])

;peak_idx = [first_vel_zero_idx_minus+35:first_vel_zero_idx_plus-35]
peak_idx = where(abs(vel_x_over) le 5)
ccf_tt_peak_fit = mpfitpeak(vel_x_over[peak_idx], ccf_tt_over[peak_idx], ccftt_coeff, nterms=3, estimates=[max(ccf_tt_over),0.0,5])

ccf_tt_peak_y = gaussian(vel_x_over, ccftt_coeff)

plot, vel_x_over, ccf_tt_over, xr=[-100,100], tit ="Central CCF_TT peak with gaussian fit"
oplot, vel_x_over[peak_idx], ccf_tt_peak_fit, color=200, thick=6
oplot, vel_x_over, ccf_tt_peak_y, color=99999, thick=3
dxstop

;ccf_tt2 = ccf_tt_over - ccf_tt_peak
ccf_tt2 = ccf_tt_over - ccf_tt_peak_y


plot, vel_x_over, ccf_tt_over, xr = [-100, 100], tit = "ccf_tt with ccf_tt2"
oplot, vel_x_over, ccf_tt2, color=200

dxstop


; Make CCF_DT

ccf_dt_over = c_correlate(apo_spec_over, pho_over, lag_over, /double)
;ccf_dt_over = c_correlate(apo_spec_over, pho_over_convol, lag_over,/double)
;ccf_dt_over = convol(apo_spec_over, pho_over_convol, /edge_wrap)
;ccf_dt_over = ccf_dt_over/(aso_norm*poc_norm)
plot, vel_x_over, ccf_tt_over, xr=[-100,100], tit='CCF_TT with CCF_DT'
oplot, vel_x_over, ccf_dt_over, color=200

; shift up
;if min(ccf_dt_over) lt 0 then ccf_dt_over = ccf_dt_over + abs(min(ccf_dt_over))

dxstop
; Fit the central peak of CCF_DT with a gaussian to get G_1
;g1 = mpfitpeak(vel_x_over, ccf_dt_over, fitco, nterms=4,
;estimates=[0.5,0.0,5.0,0.0])
g1_fit = mpfitpeak(vel_x_over, ccf_dt_over, fitco, nterms=3, estimates=[0.5,0.0,5.0])
g1_peak_fit = mpfitpeak(vel_x_over[first_vel_zero_idx_minus+35:first_vel_zero_idx_plus-35], $
                            ccf_dt_over[first_vel_zero_idx_minus+35:first_vel_zero_idx_plus-35], $
                            ccfdt_coeff, nterms=3, estimates=[max(ccf_dt_over),0.0,60])

g1_y = gaussian(vel_x_over, ccfdt_coeff)

;g1_kern = make_broadening_kernel(fitco[2], oversamp=oversamp)

;g1_kern = make_broadening_kernel(60d0, oversamp=oversamp)
;g1 = g1_kern.lsf_rot
;g1_velgrid = g1_kern.velgrid
;g1 = g1_y
; normalize g1
;g12 = g1/total(g1,/double)

plot, vel_x_over, ccf_dt_over, xr=[-100,100], tit='CCF_DT'
;oplot, g1_velgrid, g12*max(ccf_dt_over)/max(g12), color=200
;oplot, vel_x_over, g1, color=99999
dxstop
; Calculate G from G=CCF_DT-CCF2_TT*G1

plot, vel_x_over, ccf_dt_over, tit='comparison of ccf_dt and ccf_tt2'
oplot, vel_x_over, ccf_tt2, color=200
dxstop

sidelobes = convol(ccf_tt2, g1, /edge_mirror)

;scale sidelobes
sidelobe_idx = where(abs(vel_x_over) ge 100)

; fit sidelobes and ccf_dt sidelobes
ccf_tt_side = sidelobes[sidelobe_idx]
ccf_dt_side = ccf_dt_over[sidelobe_idx]

sidelobe_scaled = (sidelobes-mean(ccf_tt_side)+mean(ccf_dt_side))*stddev(ccf_dt_side)/stddev(ccf_tt_side)

plot, ccf_dt_over, tit="CCF_DT with scaled CCF_TT2*G1"
oplot, sidelobe_scaled, thick=1, color=!red
dxstop

g = ccf_dt_over-sidelobe_scaled

plot, vel_x_over, ccf_dt_over, tit='comparison of ccf_dt and ccf_tt2*g1'
oplot, vel_x_over, sidelobe_scaled, color=200
dxstop


; plot, vel_x_over, g, xr=[-75,75]

; dxstop

; Mask over the one weird pixel in the middle
; bigneg_idx = where(g le -0.5, nbigneg)
; if nbigneg gt 0 then begin
;     g[bigneg_idx] = !values.d_nan
;     g = interpol(g, vel_x_over, vel_x_over, /nan)
; endif

; oplot, vel_x_over, g, color=200

; dxstop
; Test G against known rotation kernel
;vsini_known = 60d0

g_known = make_broadening_kernel(vsini_known, oversamp=99)

g_test = g_known.lsf_rot
g_test_vel = g_known.velgrid

g_results[*,i] = g


endfor

psopen, 'vsini_correlate_mearth_test.eps', /encapsulated, xs=10, ys=10, /inches, /color

    g_known = make_broadening_kernel(vsini_known, oversamp=99)

    g_test = g_known.lsf_rot
    g_test_vel = g_known.velgrid
    plot, vel_x_over, g_results, thick=6, $
      xr = [-1.5*vsini_known, 1.5*vsini_known], $
      tit='Vsini = '+strtrim(fix(vsini_known),2), $
      xtit = 'Delta V (km/s)', $
      yr = [-0.1, max(g_results)*1.1]
    oplot, g_test_vel, g_test*max(g_results)/max(g_test), color=!red, thick=4

    psclose
stop

psopen, 'vsini_correlate_noiseless_halflsf_test.eps', /encapsulated, xs=10, ys=5, /inches, /color
!p.multi = [0,4,2]
covec = [!magenta, !blue, !cyan, !yellow, !forest, !green, !orange, !red]
for i = 0, n_elements(vsini_known_vec) -1 do begin
    
    vsini_known = vsini_known_vec[i]
    g_known = make_broadening_kernel(vsini_known, oversamp=99)

    g_test = g_known.lsf_rot
    g_test_vel = g_known.velgrid
    plot, vel_x_over, g_results[*,i], thick=5, $
      xr = [-1.5*vsini_known, 1.5*vsini_known], $
      tit='Vsini = '+strtrim(fix(vsini_known),2), $
      xtit = 'Delta V (km/s)', $
      yr = [-0.1, max(g_results[*,i])*1.1]
    oplot, g_test_vel, g_test*max(g_results[*,i])/max(g_test), color=covec[i], thick=3
endfor

psclose

stop

plot, ccf_dt_over, tit='test of CCF_DT = CCF_TT1*G + CCF_TT2*G'
oplot, convol(ccf_tt_peak_y,g1,/edge_wrap)+convol(ccf_tt2,g1,/edge_wrap), color=200

dxstop

plot, ccf_tt_over, tit='test of CCF_TT = CCF_TT1 + CCF_TT2'
oplot, ccf_tt_peak_y + ccf_tt2, color=200

dxstop

plot, convol(ccf_tt_over, g1, /edge_wrap), tit='test of CCF_TT*G = CCF_TT1*G + CCF_TT2*G'
oplot, convol(ccf_tt_peak_y,g1,/edge_wrap)+convol(ccf_tt2,g1,/edge_wrap), color=200

dxstop

psopen, 'vcorr_test_fig1.eps', /encapsulated, xs=10, ys=10, /inches, /color

plot, vel_x_over, ccf_tt_over+2
oplot, vel_x_over, sidelobes+2, color=!red

oplot, vel_x_over, sidelobes+1.5

oplot, vel_x_over, ccf_dt_over+1

oplot, vel_x_over, g+0.5

psclose

psopen, 'vcorr_test_fig2.eps', /encapsulated, xs=10, ys=10, /inches, /color

plot, vel_x_over, g, xr=[-100,100]
oplot, g_test_vel, g_test*max(g)/max(g_test), color=!red, ps=8

psclose


stop

end
