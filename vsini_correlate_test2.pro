; Testing the preliminary steps in the Vsini correlation approach, and
; adding the FT of G to calculate vsini (rather than just compare with
; the theoretical kernel
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

pro vsini_correlate_test2, display=display

if n_elements(display) eq 0 then display = 0

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
logwl0 = 4.179d0

xmin = 350L
xmax = 1300L

; Make wl grid
xpixels = lindgen(xmax-xmin+1)+xmin
logwl_grid = poly(xpixels,[logwl0,d_logwl])
wl_grid = 10d0^logwl_grid

;Oversample
oversamp = 11L
npix = n_elements(wl_grid) 
npix_over = npix*oversamp
x = dindgen(npix)
xx = (dindgen(npix_over)-(oversamp/2))/oversamp
wl_grid_over = interpol(wl_grid, x, xx)

xnorm = x-npix/2
xxnorm = xx-npix_over/(2*oversamp)

; Make APOGEE LSF out of a Gaussian with width 1.34 logwl pixels
lsf_sigma = 1.34d0
npix_lsf = 35L * oversamp
lsfx = (dindgen(npix_lsf)-(npix_lsf/2))/oversamp
lsf = gaussian(lsfx, [1d0/(oversamp*lsf_sigma*sqrt(2d0*!dpi)),0d0,lsf_sigma], /double)


; Load the PHOENIX template
pho_over_str = readin_phoenix([3000],['-0.0'],['4.5'],wl_grid_over)
pho_over = pho_over_str.spec

pho_over_convol = convol(pho_over, lsf, /edge_wrap)

; Make fake APOGEE spectrum from PHOENIX template
;apo_spec_over_unbroadened = pho_over_convol
apo_spec_over_unbroadened = pho_over


; Step through Vsini values and do the correlations
vsini_known_vec = [60d0, 30d0, 25d0, 20d0, 15d0, 10d0, 8d0, 5d0]

nvsini = n_elements(vsini_known_vec)
g_results = dblarr(n_elements(xxnorm),nvsini)



for i = 0, nvsini -1 do begin

    vsini_known = vsini_known_vec[i]

    vsini_estimate = vsini_known
    
    ; Make broadening kernel for the data (this should be estimated
    ; from the data in the final analysis, but we will use prior
    ; knowledge here for testing)
    g1_kern = make_broadening_kernel(vsini_estimate, oversamp=oversamp)
    ; normalized broadening kernel
    g1 = g1_kern.lsf_rot
    g1 = g1/total(g1,/double)
    ; abscissa velocity values
    g1_velgrid = g1_kern.velgrid

    ; Fake broadening the data
    apo_spec_over = convol(apo_spec_over_unbroadened, g1, /edge_wrap)
    ;signal_to_noise = 50d0
    ;noise = randomn(seed, n_elements(apo_spec_over))/signal_to_noise
    ;apo_spec_over = apo_spec_over + noise
    

    if display then begin
        plot, wl_grid_over, apo_spec_over, thick=4, $
          tit='APG spec and Template spec, both with LSF'
        oplot, wl_grid_over, pho_over_convol, thick=3, color=!red
        dxstop
    endif


    ; Make CCF_TT
    lag_over = lindgen(npix_over) - npix_over/2
    vel_x_over = (10d0^(xxnorm*d_logwl)-1d0)*bigc

    ccf_tt_over = a_correlate(pho_over, lag_over, /double)


    ; Remove first peak
    peak_idx = where(abs(vel_x_over) le 5)
    ccf_tt_peak_fit = mpfitpeak(vel_x_over[peak_idx], $
                                ccf_tt_over[peak_idx], ccftt_coeff, $
                                nterms=3, estimates=[max(ccf_tt_over),0.0,5])

    ccf_tt_peak_y = gaussian(vel_x_over, ccftt_coeff)

    if display then begin
        plot, vel_x_over, ccf_tt_over, xr=[-100,100], $
          tit ="Central CCF_TT peak with gaussian fit"
        oplot, vel_x_over[peak_idx], ccf_tt_peak_fit, color=200, thick=6
        oplot, vel_x_over, ccf_tt_peak_y, color=99999, thick=3
        dxstop, 'vel_x_over', 'ccf_tt_peak_y'
    endif
    
    ; CCF_TT2
    ccf_tt2 = ccf_tt_over - ccf_tt_peak_y

    if display then begin
        plot, vel_x_over, ccf_tt_over, xr = [-100, 100], tit = "ccf_tt with ccf_tt2"
        oplot, vel_x_over, ccf_tt2, color=200
        dxstop
    endif

    ;;; Make CCF_DT

    ccf_dt_over = c_correlate(apo_spec_over, pho_over, lag_over, /double)

    if display then begin
        plot, vel_x_over, ccf_tt_over, xr=[-100,100], tit='CCF_TT with CCF_DT'
        oplot, vel_x_over, ccf_dt_over, color=200
        dxstop
    
    

    plot, vel_x_over, ccf_dt_over, tit='comparison of ccf_dt and ccf_tt2'
    oplot, vel_x_over, ccf_tt2, color=200
    dxstop
    endif
    
    ; Calculate G from G=CCF_DT-CCF2_TT*G1
    sidelobes = convol(ccf_tt2, g1, /edge_mirror)

    ;scale sidelobes
    sidelobe_idx = where(abs(vel_x_over) ge 100)

    ; fit sidelobes and ccf_dt sidelobes
    ccf_tt_side = sidelobes[sidelobe_idx]
    ccf_dt_side = ccf_dt_over[sidelobe_idx]

    sidelobe_scaled = (sidelobes-mean(ccf_tt_side)+mean(ccf_dt_side)) * $
      stddev(ccf_dt_side)/stddev(ccf_tt_side)

    if display then begin
        plot, vel_x_over, ccf_dt_over, tit="CCF_DT with scaled CCF_TT2*G1"
        oplot, vel_x_over, sidelobe_scaled, thick=1, color=!red
        dxstop
    endif
    g = ccf_dt_over-sidelobe_scaled

    ; Test G against known rotation kernel

    g_known = make_broadening_kernel(vsini_known, oversamp=99)

    g_test = g_known.lsf_rot
    g_test_vel = g_known.velgrid

    g_results[*,i] = g


endfor

; Try the FT for the vsini=60 case
ti = 0
vsini_known = vsini_known_vec[ti]

; velx_fft_idx = where(abs(vel_x_over) le 1201,nfft)
; help, nfft
; if nfft mod 2 eq 0 then message, "Even number of elements in rotation kernels!"
; velx_fft = vel_x_over[velx_fft_idx]
velx_fft_idx = lindgen(n_elements(vel_x_over)) 
velx_fft = vel_x_over

; retrieve the derived kernel
g_emp_long = g_results[*,ti]
g_emp = g_emp_long[velx_fft_idx]

; Make derived kernel fall to zero exponentially above 100 km/s
decay_vec = replicate(1d0,n_elements(vel_x_over))
velx_zero_idx = where(abs(vel_x_over) ge 100,nfft,complement=mid_idx)
velx_decay_vec = abs(vel_x_over)-100
velx_decay_vec[mid_idx] = 0d0
velx_decay_vec = velx_decay_vec/50d0
decay_vec = exp(-1d0*velx_decay_vec)

plot, vel_x_over, decay_vec
dxstop
plot, vel_x_over, g_emp
oplot, vel_x_over, g_emp*decay_vec, thick=2, color=200
dxstop




g_emp[velx_zero_idx] = 0d0

; make theoretical kernel
g_now = make_broadening_kernel(vsini_known, oversamp=oversamp)
g_comp_norm = g_now.lsf_rot
g_comp_vel = g_now.velgrid
; put on same abscissa vector
g_comp_norm_long = interpol(g_comp_norm, g_comp_vel, velx_fft)
;scale g_comp to g_emp size
g_comp = g_comp_norm_long*max(g_emp)/max(g_comp_norm_long)

plot, vel_x_over, g_emp_long, xr=[-150,150]
oplot, velx_fft, g_comp, color=!red, thick = 2
oplot, velx_fft, g_emp, color=!blue, thick=4
dxstop, 'vel_x_over', 'g_emp_long', 'velx_fft', 'g_comp', 'g_emp'

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

;;; Take the FFTs
g_emp_fft = fft(g_emp, /double);, /center)
g_comp_fft = fft(g_comp, /double);, /center)

nstep = n_elements(velx_fft) 
bigc = 3d5 ;km/s
dlogwl = 6d-6
; change in velocity equal to 1 oversampled logwl pixel
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc

sigvec_x = dindgen((nstep - 1)/2) + 1
is_N_even = (nstep MOD 2) EQ 0
if (is_N_even) then $
  sigvec = double([0.0, sigvec_x, nstep/2, -nstep/2 + sigvec_x]/(nstep*deltav)) $
else $
  sigvec = double([0.0, sigvec_x, -(nstep/2 + 1) + sigvec_x]/(nstep*deltav))

vsini_sigvec = k1/(sigvec)



plot, sigvec*1d3, alog10(g_comp_fft), xr=[0,40], ps=-6, thick=1
plot, sigvec*1d3, alog10(g_emp_fft), xr=[0,40], ps=-6, thick=1
dxstop
plot, sigvec*1d3, vsini_sigvec, xr=[0,40], ps=-6, thick=1
;oplot, sigvec*1d3, alog10(g_emp_fft), color=!red
;axis, 
dxstop, 'g_emp_fft', 'g_comp_fft'

; Take amplitude of FT

g_emp_fft_amp = alog10(sqrt(real_part(g_emp_fft)^2 + imaginary(g_emp_fft)^2))
g_comp_fft_amp = alog10(sqrt(real_part(g_comp_fft)^2 + imaginary(g_comp_fft)^2))

; smooth
g_smooth = gauss_smooth(g_emp_fft_amp, 2)
g_emp_fft_amp = g_smooth


; Find the location of the first peak using the third derivative
d_gfft_minus = shift(g_emp_fft_amp, -1) - g_emp_fft_amp
d_gfft_plus = g_emp_fft_amp - shift(g_emp_fft_amp, 1)
d_gfft = (d_gfft_minus + d_gfft_plus)/2d0

; 2nd derivative
dd_gfft = shift(d_gfft, -1) - d_gfft

; other 2nd derivative
dd_gfft_minus = shift(d_gfft, -1) - d_gfft
dd_gfft_plus = d_gfft - shift(d_gfft, 1)
dd_gfft_other = (dd_gfft_minus + dd_gfft_plus)/2d0


; 3rd derivative
ddd_gfft = dd_gfft - shift(dd_gfft, 1)

; Find zeroes

plot, sigvec*1d3, g_smooth, xr=[9,13], ps=3, symsize=1
dxstop
plot, sigvec*1d3, d_gfft, xr=[9,13], ps=3, symsize=1, yr=[-0.1,0.1], tit = "First"
dxstop
plot, sigvec*1d3, dd_gfft_other, xr=[9,13], ps=3, symsize=1, yr=[-0.1,0.1], tit="Second"
dxstop
plot, sigvec*1d3, ddd_gfft, xr=[9,13], ps=3, symsize=1, yr=[-0.1,0.1], tit="Third"

dxstop


zidx = []
sig = []

plot, sigvec*1d3, alog10(g_emp_fft), xr=[0,40]
for i = 10, 500 do begin

    ; First deriv should be neg, x, pos
    fder = d_gfft[[i-1,i+1]]
    if fder[0] lt 0 and fder[1] gt 0 then begin
        ; Second should be pos, pos, pos
        fder2 = dd_gfft_other[i-1:i+1]
        if total(fder2 gt 0) eq 3 then begin
            ; Third should be pos, x, neg
            fder3 = ddd_gfft[[i-1,i+1]]
            if fder3[0] gt 0 and fder3[1] lt 0 then begin
                sigval = interpol(sigvec[i-1:i+1],ddd_gfft[i-1:i+1],0d0)
                sig = [sig, sigval]
                zidx = [zidx, i]
                oplot, replicate(sigval*1d3,2), [-100,100], linest=2, thick=1
                dxstop
            endif
        endif
    endif

endfor




stop

end
