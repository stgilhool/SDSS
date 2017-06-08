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
;wait, 0.05


; Return the deviates
return, dev

end


function vct_find_zeroes_2, xvec, yvec


show=0

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

vsini_sigvec = k1/xvec
; define number of pixels to search
null = where(vsini_sigvec gt 3, npix_search)
fft_domain = where(vsini_sigvec le 100 and vsini_sigvec gt 3, npix_domain)
first_pix = min(fft_domain)
last_pix = max(fft_domain)

; Slide a window along and find local minima
ysmooth = gauss_smooth(yvec, 1.5)

window_width = 21L
window_idx_x = lindgen(window_width)-(window_width/2)

minidx_vec = dblarr(npix_search)
max_vec = dblarr(npix_search)

for i = first_pix, last_pix do begin

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

;dxstop, 'minima_idx_all', 'minima_idx'


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
        sigpeak_vec = xvec[minima_idx[real_minima_idx]]
        first_peak_idx = minima_idx[real_minima_idx[0]]
        sigpeak = xvec[first_peak_idx]
        vsini_vec = k1/sigpeak_vec
        vsini = k1/sigpeak
        ;print, vsini
        ;dxstop

        if show then begin
            oplot, replicate(sigpeak,2)*1d3, [-100,100], thick=3, co=!orange
            xyouts, sigpeak*1d3 + 1, 8, strtrim(vsini,2)
            dxstop
        endif
    endif else begin
        sigpeak_vec = [-99999]
        sigpeak = -99999
        vsini = 0d0
    endelse
    
endif

sigpeak_str = {minima_idx:minima_idx,$
               minima_depth:peak_depth ,$
               minima_sigval:xvec[minima_idx], $
               minima_vsini:k1/xvec[minima_idx], $
               sigvec:xvec[fft_domain],$
               fftvec:yvec[fft_domain]}

return, sigpeak_str

end


function ccf_dt_kernel_estimate, params, ccf_dt=ccf_dt, vel_x=vel_x, oversamp=oversamp

vsini = params[0]

; Make broadening kernel
kernel_str = make_broadening_kernel(vsini, oversamp=oversamp)
kernel = kernel_str.lsf_rot
velgrid = kernel_str.velgrid

; Scale broadening kernel
g1 = interpol(kernel, velgrid, vel_x)
g1_scaled = g1*max(ccf_dt)/max(g1)

; Return deviates
dev = ccf_dt - g1_scaled

return, dev

end


function vsini_correlate, lsf=lsf, $
                          display=display, $
                          pho_spec=pho_spec, $
                          pho_over=pho_over, $
                          oversamp=oversamp, $
                          apo_spec=apo_spec, $
                          apo_err=apo_err, $
                          mask=mask, $
                          pix_x=pix_x, $
                          xx=xx, $
                          vel_x_over=vel_x_over, $
                          wl_grid=wl_grid, $
                          lag_over=lag_over, $
                          decay_vec=decay_vec, $
                          sigvec=sigvec, $
                          vsini_vfit=vsini_vfit

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



; Vars
npix = n_elements(apo_spec) 
npix_over = npix*oversamp
epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3


; Make the necessary PHOENIX spectra
pho_over_convol = convol(pho_over, lsf, /edge_wrap)

pho_spec_convol = downsample_tophat(pho_over_convol, oversamp)

; Data
apo_spec_over_raw = interpol(apo_spec, pix_x, xx)

; Change mask to just the bits we want
mask_idx = where(mask eq 1, nmask)
apo_spec_cp = apo_spec
apo_spec_nan = apo_spec
apo_spec_nan[mask_idx] = !values.d_nan
apo_err_cp = apo_err
apo_err[mask_idx] = 1d8
;Interpolate
apo_spec = interpol(apo_spec_nan, pix_x, pix_x, /nan)

; NOW, oversample
apo_spec_over = interpol(apo_spec, pix_x, xx)
wl_grid_over = interpol(wl_grid, pix_x, xx)

; Show the raw and masked APOGEE spectra
if display then begin
    plot, wl_grid_over, apo_spec_over_raw, /xs, yr=[0.6, 1.1], /ys, thick=7, $
      tit = "Raw and Masked APOGEE spectra"
    oplot, wl_grid, apo_spec_nan, ps=8, color=!red
    oplot, wl_grid_over, apo_spec_over, color=!blue, thick=1
    dxstop
endif

; Show the APOGEE spectra with the LSF-Broadened PHOENIX spectrum
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
;parinfo[3].limited=[1,1]
;parinfo[3].limits=[2d0, 4d0]
parinfo[3].step = 0.05d0


r = mpfit('fit_template_mpfit', functargs=functargs, parinfo=parinfo, $
          status=status, bestnorm=chi2, maxiter=50)

if n_elements(r) eq 1 then begin
    mpfiterror=1
    spec_model = pho_spec
endif else begin
    mpfiterror=0
    spec_model = vcorr_model_template(r, wl_grid=wl_grid, temp_spec=pho_spec)
endelse

; Mask areas of mismatch
res = apo_spec - spec_model
res[mask_idx] = !values.d_nan
median_res = median(res)
bad_res_idx = where(res gt 0.95*max(res,/nan), nbad)
nbad = 0
if nbad gt 0 then begin
    spec_model[bad_res_idx] = mean(spec_model)
    apo_spec[bad_res_idx] = mean(apo_spec, /nan)
    
    apo_spec_over = interpol(apo_spec, pix_x, xx)
endif

spec_model_over = interpol(spec_model, pix_x, xx)

if display then begin
    plot, apo_spec_over, tit="APOGEE spec with fitted (and poss. masked) model"
    oplot, spec_model_over, color=!red, thick=2
    dxstop
endif

;;;;;;;;;;;;;;;
;;; Make CCF_DT
;;;;;;;;;;;;;;;

ccf_dt_over_0 = c_correlate(apo_spec_over, spec_model_over, lag_over,/double)

;;; Recenter
lagmax = where(ccf_dt_over_0 eq max(ccf_dt_over_0))
vel_x_shifted = shift(vel_x_over, npix_over/2 - lagmax)
ccf_dt_over = interpol(ccf_dt_over_0, vel_x_over, vel_x_shifted)

if display then begin
    plot, vel_x_over, ccf_dt_over_0, xr=[-100,100], $
      tit='CCF_DT and recentered CCF_DT'
    oplot, vel_x_over, ccf_dt_over, color=!red
    dxstop
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ESTIMATE G FROM CCF_DT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

peak_idx = where(abs(vel_x_over) le 8)

functargs_ccf = {ccf_dt:ccf_dt_over[peak_idx], $
                 vel_x:vel_x_over[peak_idx]}

parinfo_ccf = {value:15d0, $
               limited:[1,1], $
               limits:[0.1d0,100d0], $
               step:1d0}

; g1_vel = mpfit('ccf_dt_kernel_estimate', parinfo=parinfo_ccf, $
;                         functargs=functargs_ccf, status=stat_ccf, $
;                         bestnorm=chi2_ccf)
; ; check that the fit worked
; if stat_ccf gt 0 then g1_vsini = g1_vel[0] $
;   else g1_vsini = vsini_vfit
g1_vsini = vsini_vfit

; make kernel from estimated vsini
g1_str = make_broadening_kernel(g1_vsini, oversamp=oversamp)

g1 = g1_str.lsf_rot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Make CCF_TT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

lsf_option=0

case lsf_option of 
    
    1: begin
    ; Hit CCF_TT with LSF after autocorrelation
        ccf_tt_over = a_correlate(pho_over, lag_over, /double)
        ccf_tt_over_cp = ccf_tt_over
        ccf_tt_over = convol(ccf_tt_over_cp, lsf, /edge_wrap)
        plottitle = "CCF_TT and CCF_TT*LSF"
        end

    2: begin
    ; Autocorrelate LSF-broadened template
        ccf_tt_over = a_correlate(pho_over_convol, lag_over, /double)
        ccf_tt_over_cp = ccf_tt_over
        plottitle = "CCF_TT with LSF-convolved template"
        end

    else: begin
    ; No LSF
        ccf_tt_over = a_correlate(pho_over, lag_over, /double)
        ccf_tt_over_cp = ccf_tt_over
        plottitle = "CCF_TT with no LSF"
        end

endcase

;;; Recenter
lagmax_tt = where(ccf_tt_over eq max(ccf_tt_over))
vel_x_shifted_tt = shift(vel_x_over, npix_over/2 + lagmax_tt)
ccf_tt_over_recenter = interpol(ccf_tt_over, vel_x_over, vel_x_shifted_tt)
ccf_tt_over = ccf_tt_over_recenter

if display then begin
    plot, vel_x_over, ccf_tt_over_cp, xr=[-50,50], tit=plottitle
    oplot, vel_x_over, ccf_tt_over, color=!red, thick=3
    nolsfmax = max(ccf_tt_over_cp, nolsfidx)
    lsfmax = max(ccf_tt_over, lsfidx)
    dxstop, 'nolsfidx', 'lsfidx'
endif 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MAKE CCF_TT2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Scale and remove sidelobes 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Clean G by making it fall to 0 exponentially past 100 km/s
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; g_out = g * decay_vec

; g_out = g_out - min(g_out)

; g_out = g_out * decay_vec
g_out = g

if display then begin
    plot, vel_x_over, g_out, tit="Derived broadening kernel with estimate", $
      xr = [-500,500]
    oplot, vel_x_over, g1*max(g_out)/max(g1), color=!red, thick=3
    dxstop
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; TAKE FFT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

g_fft = fft(g_out, /double)

g_fft_amp = sqrt(real_part(g_fft)^2 + imaginary(g_fft)^2)

glog_fft = alog10(g_fft_amp)

if display then begin
    plot, sigvec*1d3, gauss_smooth(glog_fft,1), xr=[-5,25], yr=[-15,-2], xs=8
    vs_tickvec = [60,30,25,20,15,10,8,5]	
    vs_tickval = (k1/vs_tickvec)*1e3
    axis, xaxis=1, xticks=n_elements(vs_tickvec), $
      xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Find the location (and Vsini) of first zero
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

sigpeak_str_0 = vct_find_zeroes_2(sigvec, glog_fft)

sigpeak_str = create_struct(sigpeak_str_0, 'vel_x_over', vel_x_over, 'kernel', g_out)

; plot, sigvec*1d3, glog_fft, xr=[0,20], yr=[-15,-2], xs=8, /nodata
; oplot, sigpeak_str.sigvec*1d3, sigpeak_str.fftvec
; axis, xaxis=1, xticks=n_elements(vs_tickvec), $
;   xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
; for i = 0, n_elements(sigpeak_str.minima_idx)-1 do begin
;     oplot, replicate(sigpeak_str.minima_sigval[i],2), [-100,100], linest=2, thick=1
; endfor




;if sigpeak_vec[0] gt 0 then vsini_vec = k1/sigpeak_vec $
;  else vsini_vec = 0

;return, vsini_vec

return, sigpeak_str

end
