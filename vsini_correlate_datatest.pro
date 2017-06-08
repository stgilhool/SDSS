; Using code from vsini_correlate_test and vsini_correlate_test2, we
; run the cross-correlation procedure from Diaz et al, and test the
; recoverability of Vsini for PHOENIX templates as a function of temperature
function vcorr_model_template, params, wl_grid=wl_grid, temp_spec=temp_spec

wl_shift = params[0]
c0 = params[1]
c1 = params[2]
if n_elements(params) eq 4 then scale = params[3]

; Shift wavelength of template
wl_grid_shifted = wl_grid * (1d0 + (wl_shift/(!const.c/1d3)))

spec_shifted = interpol(temp_spec, wl_grid_shifted, wl_grid)

; Adjust the continuum
npix = n_elements(spec_shifted) 
npix_half =  npix/2
continuum_orig = replicate(1d0, npix) 
xnorm = (dindgen(n_elements(spec_shifted))-npix_half)/npix_half

cont_corr = poly(xnorm, [c0,c1])
spec_shifted_2 = spec_shifted * cont_corr

; If necessary, scale the template
if n_elements(params) eq 4 then begin
    spec_scaled = spec_shifted_2^scale 
endif else begin
    spec_scaled = spec_shifted_2
endelse
; Downsample

spec_model = spec_scaled

return, spec_model

end


function fit_template_mpfit, params, apo_spec=apo_spec, apo_err=apo_err, temp_spec=temp_spec, mask=mask, wl_grid=wl_grid

model_spec = vcorr_model_template(params, temp_spec=temp_spec, wl_grid=wl_grid)

; Deviates
dev = (apo_spec-model_spec)/apo_err

; Mask bad pixels
mask_pix = where(mask eq 1, nmask)
if nmask gt 0 then dev[mask_pix] = 0d0

; Mask the ends, as well
dev[0:9] = 0d0
dev[-10:-1] = 0d0

comp_spec = model_spec
if nmask gt 0 then comp_spec[mask_pix]=!values.d_nan

plot, wl_grid, apo_spec, /xs, yr=[0.5,1.1], /ys, thick=1
oplot, wl_grid, comp_spec, color=!red, ps=-8, thick=1, symsize=0.2
wait, 0.1


; Return the deviates
return, dev

end

function fit_template_amoeba, params
common amoeba_info


end

function vct_find_zeroes_2, xvec, yvec


show=1

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

; Slide a window along and find local minima

ysmooth = gauss_smooth(yvec, 1.5)
;ysmooth = smooth(yvec, 3)

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
minima_idx_all = where(hist gt (window_width/2)+1, nminima_all)
minima_idx_idx = where(minima_idx_all gt 10, nminima)
minima_idx = minima_idx_all[minima_idx_idx]

dxstop, 'minima_idx_all', 'minima_idx'


if nminima gt 0 then begin
   
    minima = ysmooth[minima_idx]
    maxima = max_vec[minima_idx]
    
    ; Calculate the "depth"
    peak_depth = (maxima-minima)/abs(maxima)

    ;Take only peaks that are 10%
    if show then begin
        vs_tickvec = [60,30,25,20,15,10,8,5]	
        vs_tickval = (k1/vs_tickvec)*1e3

        plot, xvec*1d3, ysmooth+10, xr=[0,150], ps=-8, xs=8
        
        axis, xaxis=1, xticks=n_elements(vs_tickvec), $
          xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
        
        oplot, xvec[0:n_elements(hist)-1]*1d3, hist, ps=10, $
          thick=2, co=!red
        for w = 0, (5<(n_elements(minima_idx)-1)) do begin
            oplot, replicate(xvec[minima_idx[w]],2)*1d3, [-100,100], thick=1, co=!blue
            print, peak_depth[w]
            dxstop
        endfor
        
    endif
    real_minima_idx = where(peak_depth gt 0.06, nreal)
    if nreal gt 0 then begin
        first_peak_idx = minima_idx[real_minima_idx[0]]
        sigpeak = xvec[first_peak_idx]
        vsini = k1/sigpeak
        ;print, vsini
        ;dxstop

        if show then begin
            oplot, replicate(sigpeak,2)*1d3, [-100,100], thick=3, co=!orange
            xyouts, sigpeak*1d3 + 1, 8, strtrim(vsini,2)
            dxstop
        endif
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


pro vsini_correlate_datatest

display = 1
ps=0

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

files = [289,1142]
teff_vec = [3100,3200]
vsini_known_vec = [25d0,25d0]
apg_info = readin_apo2(files=files, hdu=8)
;lsf_info = *apg_info[files].output
;dxstop, 'lsf_info'
nfiles = n_elements(files) 

xmin=6428
xmax=8156

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



; Initialize results arrays

if ps then begin
    psopen, 'vsc_datatest.eps', /encapsulated, xs=10, ys=7, /inches, /color
    !p.multi=[0,5,3]
endif



for specnum = 0, nfiles-1 do begin

    ; Make APOGEE LSF
    ;lsf_sigma = 1.34d0
    npix_lsf = 35L * oversamp
    lsfx = (dindgen(npix_lsf)-(npix_lsf/2))/oversamp
    lsfco = *apg_info[specnum].output
    lsf0 = (xmax-xmin)/2 + xmin
    lsf = lsf_gh(lsfx+lsf0, lsf0, lsfco)
    lsf = lsf/total(lsf,/double)
    
    teff = teff_vec[specnum]

    ; Load the PHOENIX template
    pho_spec_str = readin_phoenix([teff],['-0.0'],['4.5'],wl_grid)
    pho_spec = pho_spec_str.spec
    pho_over_str = readin_phoenix([teff],['-0.0'],['4.5'],wl_grid_over)
    pho_over = pho_over_str.spec
    pho_over_convol = convol(pho_over, lsf, /edge_wrap)

    pho_spec_convol = downsample_tophat(pho_over_convol, oversamp)

    if ps then begin
        ;plot, vsini_known_vec, vsini_known_vec, xr=[0,26],
        ;yr=[0,26], /ys, /xs, $
        plot, vsini_known_vec, vsini_known_vec, xr=[4,16], yr=[4,16], /ys, /xs, $
          /nodata, $
          xtit = "True Vsini", ytit="Recovered Vsini", $
          tit = "Vsini Recovery for Teff = "+strtrim(teff,2)
    endif


    ; use vsini
    vsini_known = vsini_known_vec[specnum]
        
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
    
    ; Data
    apo_spec = apg_info[specnum].spectra[xmin:xmax]
    apo_err = apg_info[specnum].errors[xmin:xmax]
    apo_spec_over_raw = interpol(apo_spec, x, xx)
    
    apo_mask_str = readin_apo2(files=[files[specnum]], hdu=3)
    apo_mask_arr = *apo_mask_str.output
    apo_mask = apo_mask_arr[xmin:xmax,0]
    ; Change mask to just the bits we want
    flag_bit = [0,1,2,3,4,5,6,12,13,14]
    flag_big = indgen(17)
    flag_decimal = long(total(2L^flag_bit))

    mask = flag_decimal and apo_mask
    mask_idx = where(mask gt 0, nmask)
    mask_cp = mask
    if nmask gt 0 then mask[mask_idx] = 1
    mask = byte(mask)
    
    apo_spec_cp = apo_spec
    apo_spec_nan = apo_spec
    apo_spec_nan[mask_idx] = !values.d_nan
    apo_err_cp = apo_err
    apo_err[mask_idx] = 1d8
    ;Interpolate
    apo_spec = interpol(apo_spec_nan, x, x, /nan)

    ; NOW, oversample
    apo_spec_over = interpol(apo_spec, x, xx)
    
    plot, wl_grid_over, apo_spec_over_raw, /xs, yr=[0.6, 1.1], /ys, thick=7
    oplot, wl_grid, apo_spec_nan, ps=8, color=!red
    oplot, wl_grid_over, apo_spec_over, color=!blue, thick=1
    dxstop
   
    if display then begin
        plot, wl_grid_over, apo_spec_over, thick=2, $
          tit='APG spec and Template spec, both with LSF'
        oplot, wl_grid_over, pho_over_convol, thick=2, color=!red
        dxstop
    endif
    

    ;;; FIT the data to the model (scale, wl_shift and normalize)
    ;;; before cross_correlating
    ;; Inputs: apo_spec, mask, apo pho_over, params(s, wl, n0, n1)
    functargs = {apo_spec:apo_spec_cp, $
                 apo_err:apo_err, $
                 temp_spec:pho_spec_convol, $
                 wl_grid:wl_grid, $
                 mask:mask}

    parinfo_str = {value:0d0, $
                   limited:[0,0], $
                   limits:[0d0,0d0], $
                   step:0d0}

    nfree = 4
    parinfo = replicate(parinfo_str, nfree)
    parinfo[0].limited=[1,1]
    parinfo[0].limits=[-5d0,5d0]
    parinfo[0].step = 0.05d0
    parinfo[1].limited=[1.1]
    parinfo[1].limits=[0.9,1.1]
    parinfo[1].value = 1d0
    parinfo[2].step = 0.1d0
    parinfo[3].value=1d0
    parinfo[3].step = 0.05d0

                   
    r = mpfit('fit_template_mpfit', functargs=functargs, parinfo=parinfo, $
              status=status, bestnorm=chi2)
    
    spec_model = vcorr_model_template(r, wl_grid=wl_grid, temp_spec=pho_spec)

    ; Mask areas of mismatch
    res = apo_spec - spec_model
    res[mask_idx] = !values.d_nan
    median_res = median(res)
    bad_res_idx = where(res gt 0.9*max(res,/nan), nbad)
    if nbad gt 0 then begin
        spec_model[bad_res_idx] = apo_spec[bad_res_idx]
        ; don't know why this next line is here...
        apo_spec_over = interpol(apo_spec, x, xx)
    endif

        spec_model_over = interpol(spec_model, x, xx)

    if display then begin
        plot, apo_spec_over, tit="masked spectra"
        oplot, spec_model_over, color=!red, thick=2
        dxstop
    endif
        
    ;dxstop, 'r', 'status', 'chi2'
    ;;; Make CCF_DT
    
    ;ccf_dt_over_0 = c_correlate(apo_spec_over, pho_over,
    ;lag_over,/double)
    ccf_dt_over_0 = c_correlate(apo_spec_over, spec_model_over, lag_over,/double)
    ;ccf_dt_over = c_correlate(apo_spec_over, pho_over_convol, lag_over, /double)
    
            
    
    ;;; Recenter
    lagmax = where(ccf_dt_over_0 eq max(ccf_dt_over_0))
    pho_over_shift = shift(pho_over, npix_over/2 - lagmax)

    ccf_dt_over = c_correlate(apo_spec_over, pho_over_shift, lag_over,/double)
    ccf_dt_over2 = ccf_dt_over_0
    if display then begin
        plot, vel_x_over, ccf_dt_over_0, xr=[-100,100], $
          tit='CCF_DT and recentered CCF_DT'
        oplot, vel_x_over, ccf_dt_over, color=200
        oplot, shift(vel_x_over, npix_over/2 + lagmax), ccf_dt_over2, color=!red
        dxstop
    endif

    ccf_dt_over = ccf_dt_over2

    ; Make CCF_TT
    ccf_tt_over = a_correlate(pho_over, lag_over, /double)
    ;ccf_tt_over = a_correlate(pho_over_convol, lag_over, /double)
    ; Add a convolution with LSF
    ccf_tt_over_cp = ccf_tt_over
    ccf_tt_over = convol(ccf_tt_over_cp, lsf, /edge_wrap)
    
    if display then begin
        plot, vel_x_over, ccf_tt_over_cp, xr=[-50,50], tit="CCF_TT and CCF_TT*LSF"
        oplot, vel_x_over, ccf_tt_over, color=!red, thick=3
        nolsfmax = max(ccf_tt_over_cp, nolsfidx)
        lsfmax = max(ccf_tt_over, lsfidx)
        dxstop, 'nolsfidx', 'lsfidx'
    endif 
        
    
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
    
    if display then begin
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
    g_theory = convol(g_interp, lsf, /edge_wrap)
    
    ; Clean G by making it fall to 0 exponentially past 100 km/s
    g_out = g * decay_vec
    
    g_out = g_out - min(g_out)
    
    g_out = g_out * decay_vec
    
    ;g_results[tidx, i, *] = g_out
    
    
    g_theory = g_theory*max(g_out)/max(g_theory)
    
    if display then begin
        plot, vel_x_over, g_out, tit="Derived broadening kernel with theory", $
          xr = [-500,500]
        oplot, vel_x_over, g_theory, color=!red, thick=3
        dxstop
    endif
    
    
        ;;; TAKE FFT
    g_fft = fft(g_out, /double)
    
    g_fft_amp = sqrt(real_part(g_fft)^2 + imaginary(g_fft)^2)
    
    glog_fft = alog10(g_fft_amp)
    
    ;gfft_results[tidx, i, *] = glog_fft
    
    
    ; Compare with known kernel
    gcomp = g_theory
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
    
    
    
    plot, sigvec*1d3, gauss_smooth(glog_fft,1), xr=[-5,25], yr=[-15,-2], xs=8
    oplot, sigvec*1d3, gauss_smooth(gcomplog_fft,1), color=!red
    vs_tickvec = [60,30,25,20,15,10,8,5]	
    vs_tickval = (k1/vs_tickvec)*1e3
    axis, xaxis=1, xticks=n_elements(vs_tickvec), $
      xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
    ;         oplot, sigvec*1d3, lsflog_fft, color=!orange
    ;         oplot, sigvec*1d3, gcomplog_fft + lsflog_fft, color=!blue
    
    
        ;;; Find the location (and Vsini) of first zero
    print, vsini_known
    ;sig_0 = vct_find_zeroes_2(sigvec, glog_fft)
    
    sig_0 = vct_find_zeroes_2(sigvec, glog_fft)
    
    ;         sig_results[tidx,i] = sig_0[0]
    
    vsini_recovered = k1/sig_0
    print, vsini_recovered
    
    

        ;;; PLOT
    ;oplot, [vsini_known], [vsini_recovered], ps=8, color=!red
    if ps then begin
        oplot, replicate(vsini_known,n_elements(sig_0)) , [k1/sig_0], ps=8, color=!red
    endif
endfor

if ps then oplot, lindgen(100), lindgen(100), thick=1, linest=2



if ps then psclose        

stop

end
