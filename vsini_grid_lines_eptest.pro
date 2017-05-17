; Brute force attempt to measure vsini, followed by fitting vsini with
; best results from grid search
; Using full grid
; Using ASPCAP spectra
; Restricting spectrum to select metal lines
; Also varying epsilon (limb-darkening coefficient)


function model_observation, par, KEY=key


common constants, oversamp, bigc, dlogwl, deltav
common function_args, xxnorm, lsf, apo_logwl_over, pho_spec_over, apo_obs, apo_err, vsini, epsilon
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
lsf_rot_vel = lsf_rotate(deltav, vsini, velgrid=velgrid, epsilon=epsilon)
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

; Mask pixels
mask_idx = where(apo_err gt 1d7, nmask, complement=good_idx)
if nmask gt 0 then begin
    res = res[good_idx]
    dev = (final_spec[good_idx]-apo_obs[good_idx])/apo_err[good_idx]
endif

dev2 = dev^2

chi2 = total(dev2, /nan, /double)

;penalize chi2 for non-physical models
nan_check = where(finite(dev2) eq 0, nancnt, ncomplement=gdcnt)
if nancnt gt 0 then begin
    avg_chi = chi2/gdcnt
    chi2 = chi2+(3*avg_chi*nancnt)
endif

;  plot, apo_obs, /xs, yr=[0.6,1.2]
;  oplot, final_spec, color=200
;  plot, res, ps =3, yr=[-0.3, 0.3], /xs
;  wait, 0.1
;print, p, vsini, chi2

;return, chi2
return, dev

end


function model_observation_vfree, par, KEY=key


common constants, oversamp, bigc, dlogwl, deltav
common function_args, xxnorm, lsf, apo_logwl_over, pho_spec_over, apo_obs, apo_err, vsini, epsilon
; need nxx or normx0

; Take input pars
trial_vsini = par[0]
spec_scale = par[1]
norm0 = par[2]
norm1 = par[3]


; Adjust normalization
normshift = norm1 * (xxnorm) + norm0
pho_spec_over_norm = pho_spec_over/normshift

; scale model spectrum line strength
pho_spec_scaled = pho_spec_over_norm^spec_scale

;;;;;;;;;;;;;;;;
;;; Convolve with rotation kernel
;;;;;;;;;;;;;;;;

; Create rotation kernel
lsf_rot_vel = lsf_rotate(deltav, trial_vsini, velgrid=velgrid, epsilon=epsilon)
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

; Mask pixels
mask_idx = where(apo_err gt 1d7, nmask, complement=good_idx)
if nmask gt 0 then begin
    res = res[good_idx]
    dev = (final_spec[good_idx]-apo_obs[good_idx])/apo_err[good_idx]
endif

dev2 = dev^2

chi2 = total(dev2, /nan, /double)

;penalize chi2 for non-physical models
nan_check = where(finite(dev2) eq 0, nancnt, ncomplement=gdcnt)
if nancnt gt 0 then begin
    avg_chi = chi2/gdcnt
    chi2 = chi2+(3*avg_chi*nancnt)
endif

;  plot, apo_obs, /xs, yr=[0.6,1.2]
;  oplot, final_spec, color=200
;  plot, res, ps =3, yr=[-0.3, 0.3], /xs
;  wait, 0.1
;print, p, vsini, chi2

;return, chi2
return, dev

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro vsini_grid_lines_eptest

common constants, oversamp, bigc, dlogwl, deltav
common function_args, xxnorm, lsf, apo_logwl_over, pho_spec_over, apo_obs, apo_err, vsini, epsilon

epsilon = 0.4
; Inputs
; x_min and x_max most both be even or both be odd
x_min = 0L
x_max = 8574L

; Constants
oversamp=7
bigc = 3d5 ;km/s
dlogwl = 6d-6
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc
nplsf = 141 ; npix lsf kernel
npad = 5 ; npix to mask on edges
df = -8d0 ; factor for converting flux in PHOENIX templates

; Paths
;outpath = '/home/stgilhool/APOGEE/vsini_results/mptest/'
outpath = '/home/stgilhool/APOGEE/vsini_results/logg45_justlines_eptest_4/'
ppath = '/home/stgilhool/PHOENIX/logg4.5/'

; Read in apo files for LSF coeffs
napofiles = 1350L

lsf_str = readin_apo(nfiles = napofiles, hdu = 8, /flux)

lsf_coeff_arr = dblarr(28, napofiles)
bovermin_arr = lonarr(napofiles)
bovermax_arr = lonarr(napofiles)

for i = 0, napofiles-1 do begin
    lsf_coeff_arr[0,i] = *(lsf_str[i].output)
    temphead = *(lsf_str[i].header)
    bovermin_temp = fxpar(temphead, 'BOVERMIN')
    bovermax_temp = fxpar(temphead, 'BOVERMAX')
    bovermin_arr[i] = bovermin_temp
    bovermax_arr[i] = bovermax_temp
endfor

; Redefine xmin and xmax, still must both be odd or even
x_min = max(bovermin_arr)
x_max = min(bovermax_arr)
if (x_max-x_min) mod 2 ne 0 then x_max = x_max - 1

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


;;; Lines
nlines=13
linestr0 = {wl:0d0,species:'',chip:''}
linestr = replicate(linestr0, nlines)

linestr[0] = {wl:16723.572d0,species:'Al',chip:'red'}
linestr[1] = {wl:16755.215d0,species:'Al',chip:'red'}
linestr[2] = {wl:16767.951d,species:'Al',chip:'red'}
linestr[3] = {wl:15266.518d,species:'Mn',chip:'blue'}
linestr[4] = {wl:16712.214d,species:'Mn',chip:'red'}
linestr[5] = {wl:16141.281d,species:'Ca',chip:'green'}
linestr[6] = {wl:16099.189d,species:'Si',chip:'green'}
linestr[7] = {wl:16685.453d,species:'Si',chip:'red'}
linestr[8] = {wl:15167.290d,species:'K',chip:'blue'}
linestr[9] = {wl:15172.739d,species:'K',chip:'blue'}
linestr[10] = {wl:15339.246d,species:'Ti',chip:'blue'}
linestr[11] = {wl:15745.048d,species:'Mg',chip:'blue'}
linestr[12] = {wl:15753.317d,species:'Mg',chip:'blue'}


;  window, 0, xs=1800, ys=900
;  plot, 10d0^apo_logwl_trunc, apo_spec[*,0], yr=[0.2,1.5]
;  for l = 0, nlines-1 do oplot, replicate(linestr[l].wl,2), [0,2], linest=2

; Make line mask
linewl_all = linestr.wl
linewl = linewl_all[sort(linewl_all)]
good_idx = []

for l = 0, nlines-1 do begin 
    goodidx = where(abs(10d0^apo_logwl_trunc-linewl[l]) le 5d0, ngood)
    if ngood gt 0 then begin
        good_idx = [good_idx,goodidx]
    endif
endfor
mask = replicate(1, nx)
mask[good_idx] = 0

mask_idx = where(mask eq 1, nmask)

; oplot, 10d0^apo_logwl_trunc[mask_idx], apo_spec[mask_idx,0], color=200


; Metallicity grid
;pfeh = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0',
;'-0.5', '-0.0', $
pfeh = ['-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
        '+0.3', '+0.5']

;pfeh = ['-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', '+0.3', '+0.5']
nfeh = n_elements(pfeh)

; Temperature grid
nteff = 16
teff_start = 24
pteff = lindgen(nteff)+teff_start

; Vsini grid
;vsini_vec = dindgen(10)*5 + 4d0
vsini_vec = [4d0, 6d0, 10d0, 16d0, 24d0, 34d0, 46d0, 60d0, 76d0]
nvel=n_elements(vsini_vec)

; Initialize arrays
pho_spec_arr = ptrarr(nteff, nfeh, /allocate_heap)
pho_wave_arr = ptrarr(nteff, nfeh, /allocate_heap)
chiarr = dblarr(nteff, nfeh, nvel)
chiarr[*,*,*] = !values.d_nan
;best_str_init = {teff:0L, feh:0.0, vsini_best:0d0, chi2_best:0d0}
null_vfree = {teff:0L, feh:'', vsini_best:0d0, chi2_best:0d0, valid_flag:0}
;best_str = replicate(best_str_init, nteff*nfeh)
best_str = replicate(null_vfree, nteff*nfeh)

;;; Loop to read in and minimally process all PHOENIX files
;;;Loop over teff
foreach teff, pteff, teff_idx do begin
    ;;;Loop over metallicity
    for feh = 0, nfeh-1 do begin
        
        ; Readin PHOENIX
        pfile = ppath+strtrim(teff*100L,2)+pfeh[feh]+'_result.fits'
        if file_test(pfile) then begin
            pstruct = mrdfits(pfile,1)
            pspec = pstruct.flux
            pwl = pstruct.wave
            npspec = n_elements(pspec)
        
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
        endif else begin
            ;        if npspec eq 1 then begin
            *(pho_spec_arr[teff_idx, feh]) = !values.d_nan
            *(pho_wave_arr[teff_idx, feh]) = !values.d_nan
            continue
        endelse
    endfor
endforeach

;stop
;;; Now do the fitting
;;;Loop over APOGEE files


;;; Read in ASPCAP stuff and make cuts
aspfile = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'
aspcap_info = mrdfits(aspfile,1)

;;; Read in irtf results
irtf_info = mrdfits('/home/stgilhool/APOGEE/master_info6.fits',1)

; Axe bad spectra
aspcapflag = aspcap_info.aspcapflag
teff_all = aspcap_info.teff
vsini_all = aspcap_info.vsini

;cut the aspcap VSINI_WARN and VSINI_BAD guys
vsini_warn_bit = 14
vsini_bad_bit = 30
star_bad_bit = 23

flag_bits = [vsini_warn_bit, star_bad_bit, vsini_bad_bit]

flag_dec = long(total(2L^flag_bits))

sel_idx = where(teff_all gt 2600 and teff_all lt 4000 and vsini_all gt 0 $
                and (aspcapflag and flag_dec) eq 0, nsel)

asp_info = aspcap_info[sel_idx]
;aspcap_idx = aspcap_idx_all[sel_idx]
;teff_asp = teff_all[sel_idx]
;teff_err_asp = teff_err_all[sel_idx]
vsini_asp = vsini_all[sel_idx]
;snr_asp = snr_all[sel_idx]
model_asp = asp_info.aspcap_class

apoidx_vec = lindgen(napofiles)
; Choose only these files for analysis
;apoidx_vec = [143,180,281,286,339,424,485,513,520,523,572,573,584,592,622,740,826,840,$
;              869,884,904,915,944,949,953,958,1038,1076,1132,1142]


;;test
;pteff = [3600,3700,3800]


foreach afil, apoidx_vec, aidx do begin
    
    ; Define outfile
    outname = 'rfile' + strtrim(afil,2) + '.fits'
    outfile = outpath + outname
    outname_temp = 'rfile' + strtrim(afil,2) + '_temp.fits'
    outfile_temp = outpath + outname_temp

    ; Check if the file has made the cut
    cut_check = where(sel_idx eq afil, ncheck)
    if ncheck eq 0 then begin
        print, "File doesn't make the cut.  Skipping idx: "+strtrim(afil,2)
        continue
    endif
    
    ; Check to see if it's being/been processed, and skip if so
    if file_test(outfile) eq 1 or file_test(outfile_temp) eq 1 then continue
    ; Temporarily write the file if we're going to process it
    mwrfits, ['temp'], outfile_temp, /create
        
    ; Start clock
    apo_file_time = tic("Full grid for APOGEE file " + strtrim(afil,2))
    ; Define observation and error vectors
    apo_obs = apo_spec[*,afil]
    apo_err = apo_errors[*,afil]

    ; Apply linemask
    apo_err[mask_idx] = 1d8
    apo_err[0:npad-1] = 1d8
    apo_err[-1*(npad):-1] = 1d8
    ; Define LSF
    lsf_coeff = lsf_coeff_arr[*,afil]
    lsf = lsf_gh(lx, lsf_0, lsf_coeff)/oversamp


    ;;; Define which teff, and feh's to loop over
    ; teff_asp = aspcap_info[afil].teff
;     tpos = value_locate(pteff*100d0, teff_asp)
;     loop_teff_idx = lindgen(4)+tpos-1
;     tvalid_idx = where(loop_teff_idx ge 0 and loop_teff_idx lt nteff, ntpos)
;     if ntpos eq 0 then message, "Error with determining teff's for looping"
;     teff_loop_idx = loop_teff_idx[tvalid_idx]
    
    

    ; if finite(irtf_info[afil].feh_irtf) then begin
;         pos = value_locate(pfeh, irtf_info[afil].feh_irtf)
;         ; take indices for the 2 above and 2 below
;         loop_feh = lindgen(4)+pos-1
;         ; only keep real (positive) indices
;         valid_idx = where(loop_feh ge 0 and loop_feh lt nfeh, npos)
;         if npos eq 0 then message, "Error with determining feh's for looping"
;         feh_loop = loop_feh[valid_idx]
;     endif else feh_loop = lindgen(n_elements(pfeh))
;     nfeh_loop = n_elements(feh_loop) 
;     feh_0 = feh_loop[0]
;     feh_1 = feh_loop[-1]

    ;;;Loop over teff
    ;foreach teff, pteff, teff_idx do begin
    ;for teff_idx = 16,18 do begin
    for teff_idx = 0, nteff-1 do begin
        teff = pteff[teff_idx]
        ;;;Loop over metallicity
        ;for feh = 0, nfeh-1 do begin
        for feh = 0, nfeh-1 do begin

            ; Retrieve PHOENIX
            pho_spec  = *(pho_spec_arr[teff_idx, feh])
            pho_logwl = *(pho_wave_arr[teff_idx, feh])
            npspec = n_elements(pho_spec) 
            
            if npspec eq 1 then begin
                chiarr[teff_idx, feh, *] = !values.d_nan
                best_str [teff_idx*feh] = null_vfree
                continue
            endif
            
            pho_spec_over = interpol(pho_spec, pho_logwl, apo_logwl_over)
                        
            ;;;Loop over vsini
            nparam = 3
            temporary_results = dblarr(nparam, nvel)
            for v = 0, nvel-1 do begin
            ;for v = 0, 2 do begin
                
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

                r = mpfit('model_observation', guess, maxiter = 300, bestnorm=chi2, /quiet)
                !p.multi = 0                
                ; Record chi2
                if n_elements(r) eq 1 then begin
                    chiarr[teff_idx,feh,v] = !values.d_nan
                    temporary_results[*,v] = !values.d_nan
                endif else begin
                    chiarr[teff_idx,feh,v] = chi2
                    temporary_results[*,v] = r
                endelse
                ; Report result
                print, afil, teff*100, pfeh[feh], vsini, $
                  chiarr[teff_idx,feh,v]
                toc, one_run_time
                
                ;dxstop
            endfor
            
            ; Get the results for the velocities at this Teff and FeH
            chivel = reform(chiarr[teff_idx,feh,*])
            
            best_vsini_guess_idx = where(chivel eq min(chivel, /nan), nmin)
            if nmin ne 0 then begin
                best_vsini_guess = vsini_vec[best_vsini_guess_idx[0]]
                best_param_guess = temporary_results[*,best_vsini_guess_idx[0]]

                vfree_guess = [best_vsini_guess, best_param_guess]

                parinfo_str = {value:0d0, limited:[0,0], limits:[0d0,0d0], relstep:0}
                parinfo = replicate(parinfo_str, n_elements(vfree_guess))

                parinfo.value = vfree_guess
                parinfo[0].limited = [1,1]
                parinfo[0].limits = [0d0, 100d0]
                parinfo[0].relstep = 0.1
                parinfo[1:-1].relstep = 1d-3
                
                vfree_r = mpfit('model_observation_vfree', parinfo=parinfo, maxiter = 750, bestnorm=chi2_vfree, /quiet)
                
                if n_elements(r) ne 1 then begin
                    
                    vfree_str = {teff:teff*100, $
                                 feh:pfeh[feh], $
                                 vsini_best:vfree_r[0], $
                                 chi2_best:chi2_vfree, $
                                 valid_flag:1}
                    best_str[teff_idx*feh] = vfree_str
                endif else begin
                    best_str[teff_idx*feh] = null_vfree
                    print, "fit didn't converge"
                    dxstop
                endelse
            endif
                

        endfor
    endfor

;stop    
    ; Save and write results
    mwrfits, chiarr, outfile, /create, /silent
    ; write best stuff to HDU1
    mwrfits, best_str, outfile, /silent
    file_delete, outfile_temp, /quiet
    print, "Wrote file "+strtrim(afil,2)
    toc, apo_file_time
endforeach

print, "Process finished"
stop
end
