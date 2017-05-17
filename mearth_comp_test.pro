; Procedure to test the Vsini measurements against MEarth rotation
; periods
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

pro mearth_comp_test, specnum

input_file = '/home/stgilhool/APOGEE/vsini_results/mearth_comparison.fits'
instr = mrdfits(input_file, 1)

; Get the indices for the overlaps
apg_idx = instr.apg_idx
n_overlaps = n_elements(apg_idx) 

; Get apg info
apg_info = readin_apo2(files=apg_idx, hdu=8)
; Get master info
master_info = mrdfits('/home/stgilhool/APOGEE/master_info6.fits',1)
m_info = master_info[apg_idx]

; Save LSF coefficients
lsf_co = dblarr(26,n_overlaps)
for i = 0, n_overlaps-1 do lsf_co[*,i] = *apg_info[i].output

; Take a chunk
xmin = 1000
xmax = 3000

; Load a PHOENIX template
wl_grid_chunk = apg_info[0].wl_grid[xmin:xmax]
; oversamp the wlgrid
oversamp = 21L
npix = n_elements(wl_grid_chunk) 
npix_over = npix*oversamp
x = dindgen(npix)
xx = (dindgen(npix_over)-(oversamp/2))/oversamp
wl_grid_over = interpol(wl_grid_chunk, x, xx)
pho = readin_phoenix([2800],['-0.0'],['4.5'],wl_grid_over)

; Take a chunk, broaden template, and cross-correlate

apgspec = apg_info[specnum].spectra[xmin:xmax]
phospec = pho.spec

; oversample apg spec
apgspec_over = interpol(apgspec, wl_grid_chunk, wl_grid_over)

; make lsf

npix_lsf = 31*oversamp
lsfx = (dindgen(npix_lsf)-(npix_lsf/2))/double(oversamp)
lsf0 = (xmax-xmin)/2 + xmin

lsf = lsf_gh(lsfx+lsf0,lsf0,lsf_co[*,specnum])
lsf = lsf/total(lsf, /double)

; Broaden the PHOENIX spectrum
; Create rotation kernel
vsini_vec = dindgen(19)+1d0
nvsini = n_elements(vsini_vec) 
npix_lag = 31L*oversamp
lag = lindgen(npix_lag) - (npix_lag/2)

ccorr_result = dblarr(npix_lag, nvsini)

!p.multi=[0,1,2]

foreach vsini, vsini_vec, vsini_idx do begin

    ; Doppler broaden
    rot_kern = make_broadening_kernel(vsini, oversamp=oversamp)
    lsf_rot = rot_kern.lsf_rot
        
    rot_spec = convol(phospec, lsf_rot, /edge_wrap)
    
    ; LSF broaden
    phospec_convol = convol(rot_spec, lsf, /edge_wrap)
    
    phospec_final = phospec_convol+0.05

    ; Cross-correlate

    ccorr_result[*,vsini_idx] = c_correlate(apgspec_over, phospec_final, lag)
    
    plot, wl_grid_over, apgspec_over, yr = [0.8, 1.1], /ys, xr=[1.54d4,1.55d4]
    oplot, wl_grid_over, phospec_final, color=200
    
    plot, lag, ccorr_result[*,vsini_idx]

    print, "VSINI     PEAK"
    print, vsini, max(ccorr_result[*,vsini_idx])
    
endforeach



stop



end
