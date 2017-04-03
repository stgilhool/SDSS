; Brute force attempt to measure vsini


function model_observation, par, KEY=key


common constants, oversamp, bigc, dlogwl, deltav
common function_args, xxnorm, lsf, apo_logwl_over, pho_spec_over, apo_obs, apo_err, vsini
; need nxx or normx0

; Take input pars
spec_scale = par[0]
norm0 = par[1]
norm1 = par[2]


; Adjust normalization
normshift = norm1 * (xxnorm) + norm0
pho_spec_over_norm = pho_spec_over/normshift

; scale model spectrum line strength
pho_spec_scaled = pho_spec_over_norm^spec_scale

;;;;;;;;;;;;;;;;
;;; Convolve with rotation kernel
;;;;;;;;;;;;;;;;

; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid)
; Convert velocity-space x-axis to logwl-space
velgrid_kernel_x = alog10(1d0 + (velgrid/bigc))

; Make a similar vector which matches the logwl-grid
nlsfrot = n_elements(velgrid)+6 ;pad to make sure we get the whole kernel
lsf_rot_kernel_x = (dindgen(nlsfrot)-(nlsfrot/2)) * dlogwl/oversamp
; Spline lsf_rot onto the logwl_grid vector
lsf_rot = interpol(lsf_rot_vel, velgrid_kernel_x, lsf_rot_kernel_x)

; Make sure interpolation doesn't introduce negatives
rot_neg_idx = where(lsf_rot lt 0d0, ncnt)
if ncnt gt 0 then lsf_rot[rot_neg_idx]=0d0
; Normalize

lsf_rot = lsf_rot/total(lsf_rot, /double)

; Convolve with the spectrum
rot_spec = convol(pho_spec_scaled, lsf_rot, /edge_wrap)

;;;;;;;;;;;;;;;;
;;; Convolve with APOGEE lsf
;;;;;;;;;;;;;;;;
model_spec = convol(rot_spec, lsf, /edge_wrap)

; Downsample
final_spec = downsample_tophat(model_spec, oversamp)

; Calc chi2
res = final_spec-apo_obs

dev = (final_spec-apo_obs)/apo_err

dev2 = dev^2

chi2 = total(dev2, /nan, /double)

;penalize chi2 for non-physical models
nan_check = where(finite(dev2) eq 0, nancnt, ncomplement=gdcnt)
if nancnt gt 0 then begin
    avg_chi = chi2/gdcnt
    chi2 = chi2+(3*avg_chi*nancnt)
endif

plot, apo_obs, /xs, yr=[0.6,1.2]
oplot, final_spec, color=200
plot, res, ps =3, yr=[-0.3, 0.3], /xs
;print, p, vsini, chi2

;return, chi2
return, dev

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro vsini_grid

common constants, oversamp, bigc, dlogwl, deltav
common function_args, xxnorm, lsf, apo_logwl_over, pho_spec_over, apo_obs, apo_err, vsini

; Inputs
; x_min and x_max most both be even or both be odd
x_min = 2500
x_max = 3200

; Constants
oversamp=7
bigc = 3d5 ;km/s
dlogwl = 6d-6
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc
nplsf = 141 ; npix lsf kernel
npad = 5 ; npix to mask on edges
df = -8d0 ; factor for converting flux in PHOENIX templates

; Paths
outpath = '/home/stgilhool/APOGEE/vsini_results/mptest/'
ppath = '/home/stgilhool/PHOENIX/logg4.5/'

; Read in apo files for LSF coeffs
napofiles = 1350

lsf_str = readin_apo(nfiles = napofiles, hdu = 8, /flux)

lsf_coeff_arr = dblarr(28, napofiles)

for i = 0, napofiles-1 do begin
    lsf_coeff_arr[0,i] = *(lsf_str[i].output)
endfor

; Read in aspcap files for normalized flux spectra
asp_str = readin_aspcap(nfiles = napofiles, hdu = 3)

; Define indexing vectors, etc.
nx = x_max-x_min+1
nxx = nx * oversamp

x = dindgen(nx) + x_min
xx= (dindgen(nxx) - (oversamp/2))/oversamp + x_min

xxnorm = dindgen(nxx) - (nxx/2)

lsf_0 = xx[nxx/2]
lx = (dindgen(nplsf)-(nplsf/2))/oversamp + lsf_0

nxorig = n_elements(asp_str[0].spectra)

; Get WLs
apo_logwl = asp_str[0].logwl_grid
apo_wl = asp_str[0].wl_grid
; Oversample wavelength grid
apo_logwl_trunc = apo_logwl[x]
apo_logwl_over = interpol(apo_logwl_trunc, x, xx)
    

apo_spec_full = asp_str.spectra
apo_spec = apo_spec_full[x_min:x_max, *]
apo_errors_full = asp_str.errors
apo_errors = apo_errors_full[x_min:x_max, *]

; Metallicity grid
pfeh = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
        '+0.3', '+0.5']

;pfeh = ['-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', '+0.3', '+0.5']
nfeh = n_elements(pfeh)

; Temperature grid
nteff = 20
teff_start = 20
pteff = lindgen(nteff)+teff_start

; Vsini grid
;vsini_vec = dindgen(10)*5 + 4d0
vsini_vec = [4d0, 6d0, 10d0, 16d0, 24d0, 34d0, 46d0, 60d0]
nvel=n_elements(vsini_vec)

; Initialize arrays
pho_spec_arr = ptrarr(nteff, nfeh, /allocate_heap)
pho_wave_arr = ptrarr(nteff, nfeh, /allocate_heap)
chiarr = dindgen(nteff, nfeh, nvel)

;;; Loop to read in and minimally process all PHOENIX files
;;;Loop over teff
foreach teff, pteff, teff_idx do begin
    ;;;Loop over metallicity
    for feh = 0, nfeh-1 do begin
        
        ; Readin PHOENIX
        pfile = ppath+strtrim(teff*100L,2)+pfeh[feh]+'_result.fits'
        pstruct = mrdfits(pfile,1)
        pspec = pstruct.flux
        pwl = pstruct.wave
        npspec = n_elements(pspec)
       
        if npspec eq 1 then begin
            *(pho_spec_arr[teff_idx, feh]) = !values.d_nan
            *(pho_wave_arr[teff_idx, feh]) = !values.d_nan
            continue
        endif else begin
            ; Convert pspec to flux and normalize
            pflux = 10^(pspec+df)
            pho_logwl = alog10(pwl)
        
            pspec_wlrange = where(pho_logwl ge apo_logwl[0] and $
                                  pho_logwl le apo_logwl[-1], npix_pho)
            if npix_pho eq 0 then message, "ERROR: phoenix wl range problem"
            pflux = pflux[pspec_wlrange]
            pho_logwl = pho_logwl[pspec_wlrange]
            
            pflux_norm = continuum_fit(dindgen(npix_pho), pflux)
            pho_spec = pflux/pflux_norm
            *(pho_spec_arr[teff_idx, feh]) = pho_spec
            *(pho_wave_arr[teff_idx, feh]) = pho_logwl
        endelse
    endfor
endforeach

stop
;;; Now do the fitting
;;;Loop over APOGEE files
for afil = 0, napofiles-1 do begin
    
    ; Define outfile
    outname = 'rfile' + strtrim(afil,2) + '.fits'
    outfile = outpath + outname
    outname_temp = 'rfile' + strtrim(afil,2) + '_temp.fits'
    outfile_temp = outpath + outname_temp
    
    ; Check to see if it's being/been processed, and skip if so
    if file_test(outfile) eq 1 or file_test(outfile_temp) eq 1 then continue
    ; Temporarily write the file if we're going to process it
    mwrfits, ['temp'], outfile_temp, /create
    
    ; Start clock
    apo_file_time = tic("Full grid for APOGEE file " + strtrim(afil,2))
    ; Define observation and error vectors
    apo_obs = apo_spec[*,afil]
    apo_err = apo_errors[*,afil]
    apo_err[0:npad-1] = 1d8
    apo_err[-1*(npad):-1] = 1d8
    ; Define LSF
    lsf_coeff = lsf_coeff_arr[*,afil]
    lsf = lsf_gh(lx, lsf_0, lsf_coeff)/oversamp
    
    ;;;Loop over teff
    foreach teff, pteff, teff_idx do begin
        ;;;Loop over metallicity
        for feh = 0, nfeh-1 do begin

            ; Retrieve PHOENIX
            pho_spec  = *(pho_spec_arr[teff_idx, feh])
            pho_logwl = *(pho_wave_arr[teff_idx, feh])
            npspec = n_elements(pho_spec) 
            
            if npspec eq 1 then begin
                chiarr[teff_idx, feh, *] = !values.d_nan
                continue
            endif
            
            pho_spec_over = interpol(pho_spec, pho_logwl, apo_logwl_over)
                        
            ;;;Loop over vsini
            for v = 0, nvel-1 do begin
                
                ; Start clock
                one_run_time = tic()
                
                vsini = vsini_vec[v]

                ; Make initial parameter guesses
                spec_strength = 1d0
                norm0 = 1d0
                norm1 = 0d0
                ; Parameter scales
                sstrength_scale = 0.5d0
                norm0_scale = 0.1d0
                norm1_scale = 0.1d0/nxx
                
                guess = [spec_strength, norm0, norm1]
                scale = [sstrength_scale, norm0_scale, norm1_scale]
                
                ftol = 1d-10
                
                !p.multi = [0,1,2]
                
                ;r = amoeba3(ftol, scale=scale, p0=guess, $
                ;            function_name='model_observation', $
                ;            function_value=fval, nmax=1000, ncalls=ncalls)
                ;chi2 = fval[0]

                r = mpfit('model_observation', guess, maxiter = 300, bestnorm=chi2)
                                
                ; Record chi2
                if n_elements(r) eq 1 then chiarr[teff_idx,feh,v] = !values.d_nan $
                  else chiarr[teff_idx,feh,v] = chi2
                ; Report result
                print, afil, teff*100, pfeh[feh], vsini, $
                  chiarr[teff_idx,feh,v]
                toc, one_run_time

            endfor
        endfor
    endforeach
stop    
    ; Save and write results
    mwrfits, chiarr, outfile, /create, /silent
    file_delete, outfile_temp, /quiet
    print, "Wrote file "+strtrim(afil,2)
    toc, apo_file_time
endfor

print, "Process finished"
stop
end
