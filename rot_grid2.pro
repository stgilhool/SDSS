; Brute force attempt to measure vsini


function model_observation, x_min, x_max, lsf_apo_coeff, pho_spec, tell_spec, vsini, spec_scale, norm0, norm1, vshift, QUICK=quick

if n_elements(quick) eq 0 then quick=0
common apo_data, apo_logwl
oversamp=7
bigc = 3d5 ;km/s
dlogwl = 6d-6

lsf_coeff = lsf_apo_coeff
logwl = apo_logwl
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
;wl_over = 10d0^wl_log_over

wl_logshift = logwl + alog10(1d0+(vshift/bigc))
;wl_shift = 10d0^wl_logshift
; Shift phoenix spec by small amount

pho_spec_over = interpol(pho_spec, wl_logshift, wl_log_over)
;pho_spec_over = interpol(pho_spec, wl_shift, wl_over)

; Adjust normalization
normx0 = nxx/2.
normshift = norm1 * (xx-normx0) + norm0

pho_spec_over = pho_spec_over*normshift

; scale model spectrum line strength
pho_spec_scaled = pho_spec_over^spec_scale

;window, 1
;plot, pho_spec_scaled, title='pho_spec_scaled', charsize = 2

; and telluric
tell_spec_over = interpol(tell_spec, apo_logwl, wl_log_over)
;if tell_scale lt 0 then tell_scale = 0d0
tell_scale=0d0
tell_spec_scaled = tell_spec_over^tell_scale

;window, 2
;plot, tell_spec_scaled, title='tell_spec_scaled', charsize = 2



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
rot_spec_orig = convol(pho_spec_scaled, lsf_rot, /edge_wrap)

;window, 3
;plot, rot_spec_orig, title='rot_spec_orig', chars = 2


; Add telluric absorption
rot_spec = rot_spec_orig*tell_spec_scaled

;window, 4
;plot, rot_spec, title='full rot_spec', chars = 2

;;;;;;;;;;;;;;;;
;;; Convolve with APOGEE lsf
;;;;;;;;;;;;;;;;

;npix_lsf = 25L*oversamp
npix_lsf = nxx
;pad0 = 10L  ;must be even
pad0 = 0
pad = pad0*oversamp

if quick eq 0 then begin
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
endif else begin

;test
nplsf = 141

lsf_0 = xx[nxx/2]
lx = (dindgen(nplsf)-(nplsf/2))/oversamp + lsf_0
lsf_y = lsf_gh(lx, lsf_0, lsf_coeff)/oversamp
;plot, lx, lsf_y, yr=[0,0.05], chars=2
;oplot, [lsf_0,lsf_0], [0,0.05], linestyle=2

apo_spec = convol(rot_spec, lsf_y, /edge_wrap)


endelse

;apo_spec[0:4*oversamp-1] = 1d0
;apo_spec[-4*oversamp:-1] = 1d0

;window, 5
;plot, apo_spec, title='full model (oversampled)', chars=2
;apo_spec = convol_spec[*,0]
; Downsample
spec_down = downsample_tophat(apo_spec, oversamp)

;xfinal_test = x[pad0/2:nx-1-(pad0/2)]
;window, 6
;plot, spec_down, title='spec_down', chars=2


return, spec_down

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function spec_fit, p, KEY=key

common apo_data, apo_logwl
common funct_args, x_min, x_max, lsf_coeff, pho_spec, tell_spec, vsini, obs, err, vshift

spec_scale = p[0]
;tell_scale = p[1]
norm0 = p[1]
norm1 = p[2]
;vshift = p[3]
;vshift = 1d0

mod_obs = model_observation(x_min, x_max, lsf_coeff, pho_spec, tell_spec, $
                            vsini, spec_scale, norm0, norm1, vshift, /quick)


res = mod_obs-obs

dev = (mod_obs-obs)/err

dev2 = dev^2

chi2 = total(dev2, /nan, /double)

;penalize chi2 for non-physical models
nan_check = where(finite(dev2) eq 0, nancnt, ncomplement=gdcnt)
if nancnt gt 0 then begin
    avg_chi = chi2/gdcnt
    chi2 = chi2+(3*avg_chi*nancnt)
endif

plot, obs, /xs, yr=[0.6,1.2]
oplot, mod_obs, color=200
plot, res, ps =3, yr=[-0.3, 0.3], /xs
;wait, 0.001
print, p, vsini, chi2
return, chi2

end


pro rot_grid2
common apo_data, apo_logwl
common funct_args, x_min, x_max, lsf_coeff, pho_spec, tell_spec, vsini, obs, err, vshift
; Read in apogee files
; Read in telluric from apogee files
; Read in PHOENIX templates
; Dims: FeH, Template depth, Telluric Depth, vsini
; Leave off telluric for now

bigc = 3d5 ;km/s

; Read in apo files
napofiles = 1
apo_str = readin_apo(nfiles = napofiles, hdu = 8, /flux)

lsf = dblarr(28, napofiles)

for i = 0, napofiles-1 do begin
    lsf[0,i] = *(apo_str[i].output)
endfor

apo_str = readin_apo(nfiles = napofiles, hdu = 6, /flux)

asp_str = readin_aspcap(nfiles = napofiles, hdu = 3)


; Truncate all spectra to just blue chip wls
x_min = 2400
x_max = 3200
nx = x_max-x_min+1
nxorig = n_elements(apo_str[0].spectra)

apo_logwl_raw = asp_str[0].logwl_grid
apo_wl_raw = asp_str[0].wl_grid

apo_spec_full = asp_str.spectra
apo_spec = apo_spec_full[x_min:x_max, *]
apo_err_full = asp_str.errors
apo_err = apo_err_full[x_min:x_max, *]

apo_tell = dblarr(nxorig, napofiles)

for i = 0, napofiles-1 do begin
    apo_tell[0,i] = (*(apo_str[i].output))[*,0]
endfor

; "normalize" apo spectra
;apo_norm = apo_str.spectra/apo_spec_full

;apo_tell = apo_tell/apo_norm

; truncate telluric
;apo_tell = apo_tell[x_min:x_max, *]



;;; Read in PHOENIX
;pwlpath ='/home/stgilhool/RV_projects/IRCS_rv/data/supplemental/stellar_template_models/'
ppath = '/home/stgilhool/PHOENIX/'

pfeh = ['+0.5', '+0.3', '-0.0', '-0.5', '-1.0', '-1.5', '-2.0', '-2.5', '-3.0', $
       '-3.5', '-4.0']

;pwlfile = pwlpath + 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'

;pho_wl = mrdfits(pwlfile,0)
;pho_logwl = alog10(pho_wl)


px = dindgen(n_elements(apo_logwl_raw))
npad = 5


nafil = napofiles
pteff = lindgen(20)+20
nfeh = n_elements(pfeh)
vsini_vec = dindgen(40)/2. + 5


chiarr = dindgen(nafil, n_elements(pteff), nfeh, n_elements(vsini_vec))
;;;Loop over APOGEE files
for afil = 0, nafil-1 do begin
    
    asphead = *(asp_str[afil].header)
    ;Shift to heliocentric frame
    vhelio = fxpar(asphead, "VHELIO")
    heliocorr = alog10(1d0+(vhelio/bigc))
    ;apo_logwl = apo_logwl_raw
    apo_logwl = apo_logwl_raw + heliocorr
    apo_wl = 10d0^(apo_logwl)
    ;Shift to air wavelengths
    ;apo_air = apo_wl / (1.0 +  5.792105E-2/(238.0185E0 - (1.E4/apo_wl)^2) + 1.67917E-3/( 57.362E0 - (1.E4/apo_wl)^2))
    
    ;apo_airlogwl = alog10(apo_air)
    ;apo_logwl = apo_airlogwl
    
    ;apo_logwl = apo_logwl[x_min:x_max]

    
    ;;;Loop over teff
    foreach teff, pteff, teff_idx do begin
        ;;;Loop over metallicity
        for feh = 0, nfeh-1 do begin

            pfile = ppath + 'lte0'+strtrim(teff,2)+'-4.5'+pfeh[feh]+'.BT-Settl.7'
            
            if file_test(pfile) eq 0 then begin
                chiarr[afil, teff_idx, feh, *] = !values.d_nan
                continue
            endif
            df = -8d0
            ; Read in PHOENIX model
            readcol, pfile, pwl, c1, c2, c3, c4, format = "D,D,D,D,D"
            pws = sort(pwl)
            pwl = pwl[pws]
            
            pwl_air = pwl
            ;vactoair, pwl, pwl_air
            ;airtovac, pwl, pwl_air
            ;pwl_air=pwl
            ;pwl_air = pwl / (1.0 + 2.735182D-4 + 131.4182 / pwl^2 + 2.76249E8 / pwl^4)

            c1 = c1[pws]
            c2 = c2[pws]
            c3 = c3[pws]
            c4 = c4[pws]
            
            c1 = 10^(c1+df)
            c2 = 10^(c2+df)
            pho_cont = c1/c2

            pho_logwl = alog10(pwl_air)
            pspec = interpol(pho_cont, pho_logwl, apo_logwl)
            ;pspec = interpol(pho_cont, pwl_air, apo_wl)
            

            ;;;Loop over vsini
            for v = 0, n_elements(vsini_vec)-1 do begin
                
                vsini = vsini_vec[v]

                ; Define relevant spectral vectors
                pho_spec = pspec
                tell_spec = apo_tell[*,afil]
                obs = apo_spec[*,afil]
                err = apo_err[*,afil]
                err[0:npad-1] = 1d8
                err[-1*(npad):-1] = 1d8
                lsf_coeff = lsf[*,afil]
                
                ; Make initial parameter guesses
                spec_strength = 10d0
                tell_strength = 1d0
                vshift = 0d0
                norm0 = 1d0
                norm1 = 0d0
                ; Parameter scales
                sstrength_scale = 5d0
                tstrength_scale = 0.9d0
                vshift_scale = 5d0
                norm0_scale = 0.1d0
                norm1_scale = 0.1d0/nx
                
                ;guess = [spec_strength, tell_strength, norm0, norm1]
                guess = [spec_strength, norm0, norm1];, vshift]
                ;scale = [sstrength_scale, tstrength_scale, norm0_scale, norm1_scale]
                scale = [sstrength_scale, norm0_scale, norm1_scale];, vshift_scale]
                
                ftol = 1d-8
                
                !p.multi = [0,1,2]
                term=0
                r = amoeba3(ftol, scale=scale, p0=guess, function_name='spec_fit', $
                            function_value=fval, nmax=200, ncalls=ncalls)
                while term eq 0 do begin
                    dxstop    
                    
                    ;read, term, prompt="enter 0 to continue, 1 to stop"
                    read, vshift, prompt="enter vshift value: "
                    read, spec_strength, prompt="Enter spec_strength: "
                    read, norm0, prompt="Enter norm0: "
                    read, norm1, prompt="Enter norm1: "
                    guess = [spec_strength, norm0, norm1];, vshift]
                    chi2 = spec_fit(guess)	
                    
                endwhile
                
                chi2 = fval[0]
                
                chiarr[afil,teff_idx,feh,v] = chi2
                stop
            endfor
        endfor
    endforeach
endfor

stop
end
