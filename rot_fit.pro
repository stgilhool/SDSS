;+
; NAME:
; PURPOSE:
;    To measure vsini in APOGEE spectra
;    by looking at blended lines
;
; CATEGORY:
; CALLING SEQUENCE:
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;
;       Thu Jan 7 13:41:25 2016,
;       <stgilhool@iroquois.physics.upenn.edu>
;
;		
;
;-


function model_lines, x_min, x_max, amp, ctr, sig, vsini 

common globals, bigc, oversamp, dlogwl
common apo_data, lsf_apo_coeff, wl_apo, logwl_apo


lsf_coeff = lsf_apo_coeff
logwl = logwl_apo
;function model_lines, x_pixel, amp, sig, ctr, vsini

; Make pixel position vectors
nx = x_max-x_min+1
nxx = nx * oversamp

x = dindgen(nx) + x_min
xx= (dindgen(nxx) - (oversamp/2))/oversamp + x_min

nlines = n_elements(ctr)
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc

; Oversample the logwl grid
wl_log_over = interpol(logwl[x[0]:x[0]+nx-1], x, xx)

;;;;;;;;;;;;;;;;;;;
;;; Make model spectrum from Gaussians
;;;;;;;;;;;;;;;;;;;

gauss_pars = [transpose(amp), transpose(ctr), transpose(sig)]
tau_arr = dblarr(nxx, nlines)

; make the gaussians
for i = 0, nlines-1 do begin
    tau_arr[*,i] = gaussian(xx, gauss_pars[*,i])
endfor

opt_depth = total(tau_arr, 2)
spec = exp(-1d0*opt_depth)

;;;;;;;;;;;;;;;;
;;; Convolve with rotation kernel
;;;;;;;;;;;;;;;;

; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid)
nlsfrot = n_elements(velgrid)+10 ;pad to make sure we get the whole kernel
; Convert velocity-space x-axis to logwl-space
velgrid_kernel_x = alog10(1d0 + (velgrid/bigc))
; Make a similar vector which matches the logwl-grid
lsf_rot_kernel_x = (dindgen(nlsfrot)-(nlsfrot/2)) * dlogwl/oversamp
; Spline lsf_rot onto the logwl_grid vector
;lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x,
;/spline)
lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x)

; Make sure interpolation doesn't introduce negatives
rot_neg_idx = where(lsf_rot lt 0d0, ncnt)
if ncnt gt 0 then lsf_rot[rot_neg_idx]=0d0
; Normalize
lsf_rot = lsf_rot/total(lsf_rot, /double)

; Convolve with the spectrum
rot_spec = convol(spec, lsf_rot, /edge_wrap)

;;;;;;;;;;;;;;;;
;;; Convolve with APOGEE lsf
;;;;;;;;;;;;;;;;

;npix_lsf = 25L*oversamp
npix_lsf = nxx
;pad0 = 10L  ;must be even
pad0 = 0
pad = pad0*oversamp

lsf_y = dblarr(nxx, nxx-pad)

convol_spec = rebin(rot_spec, nxx, nxx-pad)

;lsf_xx = dindgen(npix_lsf) + x[0]
lsf_x0 = xx[pad/2:nxx-1-(pad/2)]

for i = 0, nxx-1-pad do begin
    ; make lsf centered at pixel xx[i] and save it
    lsf_y[0,i] = lsf_gh(xx, lsf_x0[i], lsf_coeff)/oversamp
    
endfor

; convolve   
apo_spec = total(lsf_y*convol_spec, 1, /double)
apo_spec[0:4*oversamp-1] = 1d0
apo_spec[-4*oversamp:-1] = 1d0
;apo_spec = convol_spec[*,0]
; Downsample
spec_down = downsample_tophat(apo_spec, oversamp)

;xfinal_test = x[pad0/2:nx-1-(pad0/2)]

return, spec_down

end



function fit_function, p, CTR=ctr, X_MIN=x_min, X_MAX=x_max, LINES_PER_GROUP=lines_per_group, AMPIDX0=ampidx0, SIGIDX0=sigidx0, OBS_SPEC=obs_spec, OBS_ERR=obs_err, VIS=vis

common globals, bigc, oversamp, dlogwl

ngroups = n_elements(lines_per_group)

vsini = p[-1]

for grp = 0, ngroups-1 do begin
    
    amp = p[ampidx0[grp]:ampidx0[grp]+lines_per_group[grp]-1]
    sig = p[sigidx0[grp]:sigidx0[grp]+lines_per_group[grp]-1]

    if grp eq 0 then begin
        model_spec = model_lines(x_min[grp], x_max[grp], amp, ctr, sig, vsini)
    endif else begin
        new_spec = model_lines(x_min[grp], x_max[grp], amp, ctr, sig, vsini)
        model_spec = [model_spec, new_spec]
    endelse
endfor

res = model_spec-obs_spec

err = obs_err

dev = res/err

dev = dev[5:-6]

xvec = lindgen(x_max[0]-x_min[0]+1) + x_min[0]

if vis eq 1 then begin
    plot,xvec, obs_spec, /xs, yr = [0.3,1.2]
    oplot, xvec, model_spec, color=200
    plot, xvec, res, ps=3
endif

return, dev

end




pro rot_fit, fnum

common globals, bigc, oversamp, dlogwl
bigc = 3d5 ;km/s
oversamp = 7L
dlogwl = 6d-6

vis = 1

common apo_data, lsf_apo_coeff, wl_apo, logwl_apo

; Read in a structure containing APOGEE lsf coefficients
; and other data
apo_str = readin_apo(hdu=8, /irtf, nfiles=20)
lsf_apo_coeff_ptr = apo_str[1].output
lsf_apo_coeff = *lsf_apo_coeff_ptr

wl_apo = apo_str[1].wl_grid
logwl_apo = apo_str[1].logwl_grid
obs_spectrum = apo_str[fnum].spectra
obs_error = apo_str[fnum].errors


; Change later to be read in
lines_per_group = [3]
nlinestot = total(lines_per_group)
ctr = [2908, 2915.5, 2927]
x_min = [2898L]
x_max = [2937L]
amp_guess = replicate(1d0, nlinestot)
sig_guess = replicate(0.5d0, nlinestot)
vsini_guess = 10d0
;goto, skip_over
ctr = [5782.5, 5793, 5798]
x_min = [5773L]
x_max = [5809L]

amp_guess = replicate(1d0, nlinestot)
sig_guess = replicate(0.5d0, nlinestot)
vsini_guess = 10d0
;skip_over

;ctr = [2908, 2915.5, 2927]
;x_min = [2898L]
;x_max = [2937L]

;ctr = [2908, 2915.5, 2927]
;x_min = [2898L]
;x_max = [2937L]





ngroups = n_elements(lines_per_group)


ampidx1 = total(lines_per_group, /cumulative)-1
if ngroups gt 1 then begin
    ampidx0 = [0, (ampidx1[0:-2]+1)]
endif else ampidx0 = [0]
sigidx0 = ampidx0 + nlinestot
sigidx1 = ampidx1 + nlinestot

obs_spec = obs_spectrum[x_min:x_max]
obs_err = obs_error[x_min:x_max]


; FIT
functargs = {ctr:ctr, $
              x_min:x_min, $
              x_max:x_max, $
              lines_per_group:lines_per_group, $
              ampidx0:ampidx0, $
              sigidx0:sigidx0, $
              obs_spec:obs_spec, $
              obs_err:obs_err, $
             vis:vis $
             }

parinfo0 = {value:0d0, $
            limited:[1,1], $
            limits:[0d0,0d0], $
            mpside:2, $
            fixed:0 $
           }

parinfo = replicate(parinfo0, 2*nlinestot+1)

;amp_guess = replicate(2.4d0, nlines)
;sig_guess = replicate(0.35d0, nlines)

parinfo[0:nlinestot-1].value = amp_guess
parinfo[nlinestot:2*nlinestot-1].value = sig_guess
parinfo[2*nlinestot].value = vsini_guess

;parinfo[0:2*nlinestot-1].fixed=1

parinfo[0:nlinestot-1].limits[0] = 0d0
parinfo[0:nlinestot-1].limits[1] = 10d0
parinfo[nlinestot:2*nlinestot-1].limits[0] = 1d-1
parinfo[nlinestot:2*nlinestot-1].limits[1] = 1.5d0
parinfo[2*nlinestot].limits = [0.5d0, 1d3]

window, 0, xsize = 1300, ysize = 900
!p.multi=[0,1,2]
r = mpfit('fit_function', parinfo=parinfo, functargs=functargs, bestnorm=chi2, status=status, nfev=nfev, dof=dof)

print, "Amp (fit): " + strtrim(r[0:nlinestot-1],2)
print, "Sig (fit): " + strtrim(r[nlinestot:2*nlinestot-1], 2)
print, "vsini (fit): " + strtrim(r[2*nlinestot],2)
print, ''
;print, "Amp (sim): " + strtrim(amp,2)
;print, "Sig (sim): " + strtrim(sig, 2)
;print, "vsini (sim): " + strtrim(vsini, 2)
;print, ''
print, "Chi2: " + strtrim(chi2,2)
print, "Chi2dof: " + strtrim(chi2/dof, 2)
print, "Status: " + strtrim(status, 2)

stop
 


;;; Make fake data

; Get wl soln from SDSS
; Get spectrum for comparison
; Get LSF for modeling
; Make Gaussian absorption features parameterized by sigma and amp
; Put them about 10-15 pixels apart
; Convolve with given vsini rotation kernel
; Convolve with APOGEE LSF
; Add noise (1%)

;;; Make modeling function
; Take line centroids as inputs
; Make model Gaussians with trial sigmas and amps
; Convolve with trial vsini rotation kernel
; Convolve with APOGEE LSF
; Return CHI2 

end
