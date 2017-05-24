; Using code from vsini_correlate_test and vsini_correlate_test2, we
; run the cross-correlation procedure from Diaz et al, and test the
; recoverability of Vsini for PHOENIX templates as a function of temperature
function vct_find_zeroes_2, xvec, yvec


show=0

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

; Slide a window along and find local minima

ysmooth = gauss_smooth(yvec, 1)

window_width = 21L
window_idx_x = lindgen(window_width)-(window_width/2)

npix_search = 1000L

minidx_vec = dblarr(npix_search)
max_vec = dblarr(npix_search)

for i = 0, npix_search-1 do begin

    window_idx = window_idx_x + i

    ; Get the minima
    minval = min(ysmooth[window_idx], minidx)
    min_idx = window_idx[minidx]
    ; Take the max on either side of the window and average
    maxval_minus = max(ysmooth[window_idx[0:(window_width/2)]])
    maxval_plus = max(ysmooth[window_idx[(window_width/2)+1:*]])
    maxval_avg = mean([maxval_minus,maxval_plus])
    
    minidx_vec[i] = min_idx
    max_vec[i] = maxval_avg

endfor

hist = histogram(minidx_vec, min=0, reverse_indices=ri)
minima_idx = where(hist gt (window_width/2)+1, nminima)

if nminima gt 0 then begin
   
    minima = ysmooth[minima_idx]
    maxima = max_vec[minima_idx]
    
    ; Calculate the "depth"
    peak_depth = (maxima-minima)/abs(maxima)

    ; Take only peaks that are 10%
;     plot, xvec*1d3, ysmooth+10, xr=[0,100]
;     oplot, xvec[0:n_elements(hist)-1]*1d3, hist, ps=10, $
;       thick=2, co=!red
;     for w = 0, nminima-1 do begin
;         oplot, replicate(xvec[minima_idx[w]],2)*1d3, [-100,100], thick=1, co=!blue
;         print, peak_depth[w]
;         dxstop
;     endfor
; dxstop
    real_minima_idx = where(peak_depth gt 0.1, nreal)
    if nreal gt 0 then begin
        first_peak_idx = minima_idx[real_minima_idx[0]]
        sigpeak = xvec[first_peak_idx]
        vsini = k1/sigpeak
        ;print, vsini
        ;dxstop
    endif
    
endif



;plot, xvec[0:1000]*1d3, ysmooth[0:1000] + 10, xr=[0,100]

;oplot, xvec[om:om+n_elements(hist)-1]*1d3, histogram(minidx_vec), ps=10
;dxstop
return, sigpeak

end


function vct_find_zeroes, xvec, yvec

show=0

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

; Find the location of the first minimum using 3 derivatives
ysmooth = gauss_smooth(yvec, 0)

; 1st deriv
d_minus = shift(ysmooth, -1) - ysmooth
d_plus = ysmooth - shift(ysmooth, 1)
d1 = (d_minus + d_plus)/2d0

; 2nd derivative
dd = shift(d1, -1) - d1

; other 2nd derivative
dd_minus = shift(d1, -1) - d1
dd_plus = d1 - shift(d1, 1)
dd_other = (dd_minus + dd_plus)/2d0


; 3rd derivative
ddd = dd - shift(dd, 1)

; Find zeroes
zidx = []
sig = []

if show then begin
    vs_vec = k1/xvec

    vs_tickvec = [60,30,25,20,15,10,8,5]
    vs_tickval = (k1/vs_tickvec)*1e3

    plot, xvec*1d3, ysmooth, xr=[0,150], ps=-8, xs=8
    axis, xaxis=1, xticks=n_elements(vs_tickvec), $
      xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
endif

for i = 10, 1000 do begin

    ; First deriv should be neg, x, pos
    fder = d1[[i-2,i-1,i+1,i+2]]
    if fder[0] lt 0 and fder[1] lt 0 and fder[2] gt 0 and fder[3] gt 0 then begin
        ; Second should be pos, pos, pos
        fder2 = dd_other[i-1:i+1]
        if total(fder2 gt 0) eq 3 then begin
            ; Third should be pos, x, neg
            fder3 = ddd[[i-2,i-1,i+1,i+2]]
            if fder3[0] gt 0 and fder3[1] gt 0 $
              and fder3[2] lt 0 and fder3[3] lt 0 then begin
                sigval = interpol(xvec[i-2:i+2],ddd[i-2:i+2],0d0)
                sig = [sig, sigval]
                zidx = [zidx, i]
                if show then begin
                    vsini = k1/sigval
                    oplot, replicate(sigval*1d3,2), [-10,10], linest=2, thick=1
                endif
            endif
        endif
    endif

endfor

dxstop

nsig = n_elements(sig) 

if nsig eq 0 then return, -1 else return, sig[0:9<(nsig-1)]

end

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


pro vsini_correlate_templatetest

display = 0
ps=1

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

lag_over = lindgen(npix_over) - npix_over/2
vel_x_over = (10d0^(xxnorm*d_logwl)-1d0)*bigc
deltav = ((10d0^(d_logwl/oversamp))-1d0)*bigc

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

; Make APOGEE LSF out of a Gaussian with width 1.34 logwl pixels
lsf_sigma = 1.34d0
npix_lsf = 35L * oversamp
lsfx = (dindgen(npix_lsf)-(npix_lsf/2))/oversamp
lsf = gaussian(lsfx, [1d0/(oversamp*lsf_sigma*sqrt(2d0*!dpi)),0d0,lsf_sigma], /double)

; for later use in making result go to zero quickly        
velx_zero_idx = where(abs(vel_x_over) ge 100,nfft,complement=mid_idx)
velx_decay_vec = abs(vel_x_over)-100
velx_decay_vec[mid_idx] = 0d0
velx_decay_vec = velx_decay_vec/50d0
decay_vec = exp(-1d0*velx_decay_vec)

; for later use.  the abscissa vector for the fft
sigvec_x = dindgen((npix_over - 1)/2) + 1
is_N_even = (npix_over MOD 2) EQ 0
if (is_N_even) then $
  sigvec = [0d0, sigvec_x, npix_over/2, -npix_over/2 + sigvec_x]/(npix_over*deltav) $
else $
  sigvec = [0d0, sigvec_x, -(npix_over/2 + 1) + sigvec_x]/(npix_over*deltav)

vsini_sigvec = k1/(sigvec)


; Step through Temps and Vsini values and do the correlations
;vsini_known_vec = [25d0, 20d0, 15d0, 10d0, 8d0, 5d0]
vsini_known_vec = [15d0, 14d0, 13d0, 12d0, 11d0, 10d0, 9d0, 8d0, 7d0, 6d0, 5d0]

nvsini = n_elements(vsini_known_vec)

nteff = 15
teff_vec = lindgen(nteff)*100 + 2600 ;teff=2600-4000

; Initialize results arrays
g_results = dblarr(nteff, nvsini, npix_over)
gfft_results = g_results
sig_results = dblarr(nteff, nvsini)


if ps then begin
    psopen, 'recoverability_test4.eps', /encapsulated, xs=10, ys=7, /inches, /color
    !p.multi=[0,5,3]
endif




foreach teff, teff_vec, tidx do begin

    print, teff
    ; Load the PHOENIX template
    pho_over_str = readin_phoenix([teff],['-0.0'],['4.5'],wl_grid_over)
    pho_over = pho_over_str.spec
    pho_over_convol = convol(pho_over, lsf, /edge_wrap)
    if ps then begin
        ;plot, vsini_known_vec, vsini_known_vec, xr=[0,26],
        ;yr=[0,26], /ys, /xs, $
        plot, vsini_known_vec, vsini_known_vec, xr=[5,16], yr=[5,16], /ys, /xs, $
          /nodata, $
          xtit = "True Vsini", ytit="Recovered Vsini", $
          tit = "Vsini Recovery for Teff = "+strtrim(teff,2)
    endif


    ; Step through the vsini's
    for i = 0, nvsini - 1 do begin

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
        apo_spec = convol(pho_over, lsf, /edge_wrap)
        apo_spec_over = convol(apo_spec, g1, /edge_wrap)
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
        ccf_tt_over = a_correlate(pho_over, lag_over, /double)
        ;ccf_tt_over = a_correlate(pho_over_convol, lag_over, /double)

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
        
        ccf_dt_over = c_correlate(apo_spec_over, pho_over, lag_over,/double)
        ;ccf_dt_over = c_correlate(apo_spec_over, pho_over_convol, lag_over, /double)
        
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

        g_interp = interpol(g_test, g_test_vel, vel_x_over)

        
        ; Clean G by making it fall to 0 exponentially past 100 km/s
        g_out = g * decay_vec
        
        g_results[tidx, i, *] = g_out

        ;;; TAKE FFT
        g_fft = fft(g_out, /double)
        
        g_fft_amp = sqrt(real_part(g_fft)^2 + imaginary(g_fft)^2)
        
        glog_fft = alog10(g_fft_amp)

        gfft_results[tidx, i, *] = glog_fft


        ; Compare with known kernel
        gcomp = g_interp * max(g_out)/max(g_interp)
        gcomp_fft = fft(gcomp, /double)
        gcomp_fft_amp = sqrt(real_part(gcomp_fft)^2 + imaginary(gcomp_fft)^2)
        gcomplog_fft = alog10(gcomp_fft_amp)

        ; put LSF on velgrid
        lsf_over = interpol(lsf, lsfx, xxnorm)
        ;lsf_comp = lsf_over * max(g_out)/max(lsf_over)
        lsf_comp = lsf_over
        lsf_fft = fft(lsf_comp, /double)
        lsf_fft_amp = sqrt(real_part(lsf_fft)^2 + imaginary(lsf_fft)^2)
        
        lsflog_fft = alog10(lsf_fft_amp)

        complog_fft = alog10(lsf_fft_amp*gcomp_fft_amp)

        

        ; plot, sigvec*1d3, glog_fft, xr=[0,100], yr=[-10,-2]
;         oplot, sigvec*1d3, gcomplog_fft, color=!red
;         oplot, sigvec*1d3, lsflog_fft, color=!orange
;         oplot, sigvec*1d3, gcomplog_fft + lsflog_fft, color=!blue
;         dxstop

        ;;; Find the location (and Vsini) of first zero
        print, vsini_known
        sig_0 = vct_find_zeroes_2(sigvec, glog_fft)
        
        sig_results[tidx,i] = sig_0[0]

        vsini_recovered = k1/sig_0

        

        ;;; PLOT
        ;oplot, [vsini_known], [vsini_recovered], ps=8, color=!red
        if ps then begin
            oplot, replicate(vsini_known,n_elements(sig_0)) , [k1/sig_0], ps=8, color=!red
        endif
    endfor
    
    if ps then oplot, lindgen(100), lindgen(100), thick=1, linest=2
    
endforeach

if ps then psclose        

stop

end