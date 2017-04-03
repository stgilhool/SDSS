;+
; NAME:
; PURPOSE:
;    To make fake data and test if measuring vsini
;    by looking at blended lines is feasible
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

function model_data, p, CTR=ctr, LSF_COEFF=lsf_coeff, LOGWL=logwl, OBS_SPEC=obs_spec, OBS_ERR=obs_err

; Constants
bigc = 3d5 ;km/s
oversamp = 7L
nx = 1000L
nxx = nx * oversamp
nlines = n_elements(ctr)
dlogwl = 6d-6
deltav = 0.1d0

; Inputs
amp = p[0:nlines-1]
sig = p[nlines:2*nlines-1]
vsini = double(p[2*nlines])

; Make pixel position vectors
x = dindgen(nx)
xx= (dindgen(nxx) - (oversamp/2))/oversamp

; Oversample the logwl grid
wl_log_over = interpol(logwl[0:nx-1], x, xx)

;;; Make model spectrum from Gaussians
gauss_pars = [transpose(amp), transpose(ctr), transpose(sig)]

tau_arr = dblarr(nxx, nlines)

; make the gaussians
for i = 0, nlines-1 do begin
    tau_arr[*,i] = gaussian(xx, gauss_pars[*,i])
endfor

opt_depth = total(tau_arr, 2)
spec = exp(-1d0*opt_depth)

;plot, spec
;spec = 1d0 - opt_depth ;not sure which is correct 

;;; Convolve with rotation kernel

; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid)

nlsfrot = n_elements(velgrid)
; Convert velocity-space x-axis to logwl-space
velgrid_kernel_x = alog10(1d0 + (velgrid/bigc))
; Make a similar vector which matches the logwl-grid
lsf_rot_kernel_x = (dindgen(nlsfrot)-(nlsfrot/2)) * dlogwl/oversamp
; Spline lsf_rot onto the logwl_grid vector
lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x)


; Normalize
lsf_neg_idx = where(lsf_rot lt 0, negcnt)
if negcnt gt 0 then begin
    lsf_rot[lsf_neg_idx] = 0d0
endif
lsf_rot = lsf_rot/total(lsf_rot, /double)

; Convolve with the spectrum
rot_spec = convol(spec, lsf_rot)

;;; Convolve with APOGEE lsf
npix_lsf = 25L*oversamp

apo_spec = rot_spec

for i = 0, nxx-1-npix_lsf do begin

    lsf_xx = xx[i:i+npix_lsf-1]
    lsf_x0 = xx[i+(npix_lsf/2)]

    lsf_y = lsf_gh(lsf_xx, lsf_x0, lsf_coeff)/oversamp
   
    convolspec = total(lsf_y*rot_spec[i:i+npix_lsf-1], /double)
    apo_spec[i+(npix_lsf/2)] = convolspec
   
endfor

; Downsample
spec_down = downsample_tophat(apo_spec, oversamp)


res = spec_down-obs_spec

err = obs_err

dev = res/err

return, dev

end


function rt_make_data, amp, ctr, sig, vsini, LSF=lsf, NOISE=noise, WL=wl
; Create fake data for this test given a vector of absorption line
; centroids, intrinsic widths and amplitudes, a stellar vsini, an
; LSF and sigma of gaussian noise to be added
common apo_data, lsf_apo_coeff, wl_apo, logwl_apo
; ADD check inputs

; Constants
bigc = 3d5 ;km/s
oversamp = 7L
nx = 1000L
nxx = nx * oversamp
nlines = n_elements(ctr)

; Convert log wl scale to linear
wl_log = logwl_apo
wl_lin = wl_apo

; Make pixel position vectors
x = dindgen(nx)
xx= (dindgen(nxx) - (oversamp/2))/oversamp

; Oversample the logwl grid
wl_log_over = interpol(wl_log[0:nx-1], x, xx)
dlogwl = 6d-6


;;; Make model spectrum from Gaussians
gauss_pars = [transpose(amp), transpose(ctr), transpose(sig)]

tau_arr = dblarr(nxx, nlines)

; make the gaussians
for i = 0, nlines-1 do begin
    tau_arr[*,i] = gaussian(xx, gauss_pars[*,i])
endfor

opt_depth = total(tau_arr, 2)
spec = exp(-1d0*opt_depth)
;spec = 1d0 - opt_depth ;not sure which is correct 

;;; Convolve with rotation kernel
; FIX I think I can set deltav such that no interpolation is necessary later
deltav = 0.1d0
;vsini = 30
; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid)
nlsfrot = n_elements(velgrid)
; Convert velocity-space x-axis to logwl-space
velgrid_kernel_x = alog10(1d0 + (velgrid/bigc))
; Make a similar vector which matches the logwl-grid
lsf_rot_kernel_x = (dindgen(nlsfrot)-(nlsfrot/2)) * dlogwl/oversamp
; Spline lsf_rot onto the logwl_grid vector
lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x)

; Normalize
lsf_neg_idx = where(lsf_rot lt 0, negcnt)
if negcnt gt 0 then begin
    lsf_rot[lsf_neg_idx] = 0d0
endif
lsf_rot = lsf_rot/total(lsf_rot, /double)

; Convolve with the spectrum
rot_spec = convol(spec, lsf_rot)

;;; Convolve with APOGEE lsf
npix_lsf = 25L*oversamp

apo_spec = rot_spec

for i = 0, nxx-1-npix_lsf do begin

    lsf_xx = xx[i:i+npix_lsf-1]
    lsf_x0 = xx[i+(npix_lsf/2)]

    lsf_y = lsf_gh(lsf_xx, lsf_x0, lsf_apo_coeff)/oversamp
   
    convolspec = total(lsf_y*rot_spec[i:i+npix_lsf-1], /double)
    apo_spec[i+(npix_lsf/2)] = convolspec
   
endfor

; Downsample
spec_down = downsample_tophat(apo_spec, oversamp)

if n_elements(noise) ne 0 then begin
    if noise gt 0 then begin
        final_spec = spec_down * (1d0+(randomn(seed, nx)*noise*0.01d0))
    endif else final_spec = spec_down
endif else final_spec = spec_down


return, final_spec

end

pro rot_test, vsini

common apo_data, lsf_apo_coeff, wl_apo, logwl_apo

restore, 'correlation_v4.sav'

; Read in a structure containing APOGEE lsf coefficients
; and other data
lsf_apo_str = readin_apo(hdu=8, /irtf, nfiles=2)
lsf_apo_coeff_ptr = lsf_apo_str[1].output
lsf_apo_coeff = *lsf_apo_coeff_ptr
wl_apo = lsf_apo_str[1].wl_grid
logwl_apo = lsf_apo_str[1].logwl_grid


amp = [2.74d0, 2.25d0, 2.5d0]
ctr = [485,500,525]
sig = [0.4d0,0.3d0,0.33d0]
nlines = n_elements(ctr)
;vsini=5d0


noise = 1d-2
;noise = 0
fake_data = rt_make_data(amp, ctr, sig, vsini, noise=noise)

obs_err = noise*0.01d0 > 0.01d0

; FIT
functargs = {ctr:ctr, $
             lsf_coeff:lsf_apo_coeff, $
             logwl:logwl_apo, $
             obs_spec:fake_data, $
             obs_err:obs_err $
             }
parinfo0 = {value:0d0, $
            limited:[1,1], $
            limits:[0d0,0d0], $
            mpside:2, $
            fixed:0 $
           }

parinfo = replicate(parinfo0, 2*nlines+1)

;amp_guess = replicate(2.4d0, nlines)
;sig_guess = replicate(0.35d0, nlines)
amp_guess = amp
sig_guess = sig
vsini_guess = 15d0

parinfo[0:nlines-1].value = amp_guess
parinfo[nlines:2*nlines-1].value = sig_guess
parinfo[2*nlines].value = vsini_guess

parinfo[0:2*nlines-1].fixed=1

parinfo[0].limits = [0d0, 10d0]
parinfo[1].limits = [0d0, 10d0]
parinfo[2].limits = [0d0, 10d0]
parinfo[3].limits = [0d0, 1.5d0]
parinfo[4].limits = [0d0, 1.5d0]
parinfo[5].limits = [0d0, 1.5d0]
parinfo[2*nlines].limits = [0.5d0, 1d3]

!p.multi=[0,1,2]
r = mpfit('model_data', parinfo=parinfo, functargs=functargs, bestnorm=chi2, status=status, nfev=nfev, dof=dof)

print, "Amp (fit): " + strtrim(r[0:nlines-1],2)
print, "Sig (fit): " + strtrim(r[nlines:2*nlines-1], 2)
print, "vsini (fit): " + strtrim(r[2*nlines],2)
print, ''
print, "Amp (sim): " + strtrim(amp,2)
print, "Sig (sim): " + strtrim(sig, 2)
print, "vsini (sim): " + strtrim(vsini, 2)
print, ''
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
