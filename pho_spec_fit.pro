;pro pho_spec_fit, phoenix_spectra, line_info, sdss_spectra

;;; Read in control/options file

; Read in SDSS spectra
; Read in SDSS-related info (ie, LSF, mask)
; For a given SDSS spectrum, find its irtf Fe/H
; Read in PHOENIX spectra

;;; Next, we want to fit a model to the SDSS spectrum

; Take subset of PHOENIX spectra that have similar Fe/H
; Spline PHOENIX spectra onto the SDSS wavelength grid
; Normalize PHOENIX spectra
; Convolve PHOENIX spectra with the SDSS LSF


;;;;;;;;;;;;;;;;;;;;;;;;;
;;;MAKE LSF
;;;;;;;;;;;;;;;;;;;;;;;;;
function make_lsf, lsf_coeffs, x_min, x_max, nplsf, oversamp
; LSF stuff
; Define indexing vectors, etc.
; To make the LSF, we'll call the function read_sdss_lsf_coeffs
; Then we need an x-vector (representing the pixels), which should be
; oversampled
; Then we need a centroid for the lsf.  Simplest way is to take the
; middle pixel and use that lsf for the whole spectrum.  It doesn't
; change that much.  That is what's written below.
; Call lsf_gh(x, x0, coeffs) to create the lsf, and pass that to the
; modeling function.
npix_over = nplsf * oversamp
npix_hwhm = npix_over/2
nx = x_max-x_min+1
nxx = nx * oversamp

x = dindgen(nx) + x_min
xx= (dindgen(nxx) - (oversamp/2))/oversamp + x_min

lsf_0_idx = nxx/2
lsf_0 = xx[nxx/2]
lx = xx[(0)>(lsf_0_idx-npix_hwhm):(n_elements(xx)-1)<(lsf_0_idx+npix_hwhm)]

; Define LSF
lsf = lsf_gh(lx, lsf_0, lsf_coeffs)/oversamp

return, lsf

end
;;;;;;;;;;;;;;;;;;;;;;;;;
;;;END MAKE LSF
;;;;;;;;;;;;;;;;;;;;;;;;;




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Function to fit the model
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; make it easy to call without mpfit
 ; have options to report model, rather than just chi2
; take parameters as input
; adjust PHOENIX model accordingly
; optimize the fit for each PHOENIX model
;;;

;INPUTS:
; single normalized, masked SDSS spectrum
; SDSS error
; single normalized PHOENIX spectrum
 ; if wl_shift is a free par, then we also need PHOENIX_wl
 ; otherwise, we can spline the PHOENIX spec onto the SDSS grid prior
 ; to input
; SDSS wl grid
; SDSS LSF
; Free pars: wl_shift?, norm shift?, optical depth scale
; wl_shift should apply equally to each PHOENIX spec for a given SDSS
; spec
; decided that I will just use a few iterations with small wl_shifts
; ---> not a free parameter

;OUTPUTS:
; deviates vector
; optional: best_fit_model spectrum

function fit_spectrum, free_pars, $
                       pho_wl=pho_wl, $
                       pho_spec=pho_spec, $
                       sdss_wl=sdss_wl, $
                       sdss_spec=sdss_spec, $
                       sdss_err=sdss_err, $
                       lsf=lsf, $
                       oversamp=oversamp, $
                       output_model=output_model, $
                       first_pix=first_pix, $
                       npix_select=npix_select, $
                       vis=vis

; Take input free parameters
npar = n_elements(free_pars) 
tau_scale = free_pars[0]
if npar eq 2 then norm_scale = free_pars[1] $
  else norm_scale = free_pars[1:-1]

; Adjust the optical depth of the PHOENIX spectrum
model_spec_tau = pho_spec^tau_scale

; Convolve with SDSS LSF
model_spec_convol = convol(model_spec_tau, lsf)

; Downsample to SDSS resolution, if necessary
model_spec_downsamp = downsample_tophat(model_spec_convol, oversamp)

; Adjust continuum to account for normalization errors
;FIX first_pix and npix_select
norm_adjust_factor = ircsrv_normcorr(norm_scale, $
                                     first_pix=first_pix, $
                                     npix_select=npix_select)
model_spec_normcorr = model_spec_downsamp/norm_adjust_factor

model_spec = model_spec_normcorr

; Calculate goodness-of-fit
residuals = sdss_spec - model_spec
dev = residuals/sdss_err

; Plot, if necessary
if vis then begin
    yrplotmin = (min(sdss_spec,/nan) < min(model_spec,/nan))
    yrplotmax = (max(sdss_spec,/nan) > max(model_spec,/nan))
    yrplot = [0.7, 1.2]
    plot, sdss_wl, sdss_spec, /xs, yr=yrplot
    oplot, sdss_wl, model_spec, color=200
    ;plot, sdss_wl, residuals, /xs, ps=3
    plot, sdss_wl, dev, /xs, ps=3
endif

;dxstop

; Return result
if output_model then begin
    return, model_spec
endif else begin
    return, dev
endelse

end
;;;;;;;;;;;;;;;;;;;;;;;;
;;;END MODELING FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;
;;;CALL MPFIT
;;;;;;;;;;;;;;;;;;;;;;;;
;;; mpfit each PHOENIX spec to the SDSS spec
; free pars are: optical depth scaling, 
;                normalization fudge factor?
;                wavelength fudge factor?
; output results (parameters, chi2, status, models?)
pro call_mpfit, color

npix = 8575L

oversamp=7


sdss_str = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/apg_data.fits',1)
;make sdss_wl
sdss_wl_coeffs = sdss_str[0].wl_coeffs
sdss_log_wl = poly(dindgen(npix), sdss_wl_coeffs)
sdss_wl = 10d0^sdss_log_wl

; Read in PHOENIX
; Read in phoenix file
;readin if there's no save file.  if there is, just load the save
savfile = 'pho_spec_fit_phoenix_spectra.sav'
if file_test(savfile) eq 0 then begin
    
    ; Make common PHOENIX wl grid
    ; Just oversample the SDSS grid
    npix_oversamp = npix * oversamp
    x = dindgen(npix)
    xx = (dindgen(npix_oversamp)-(oversamp/2))/oversamp
    pho_wl = interpol(sdss_wl, x, xx)

    ; Define the grid point vectors
    teff_vec = lindgen(23)*100 + 2000L
    logg_vec = ['4.5', '5.0']
    feh_vec = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', '+0.3', '+0.5'] 
        
    ; Read
    pho_str = readin_phoenix(teff_vec, feh_vec, logg_vec, pho_wl)
    nspec = n_elements(pho_str)

    ; Save spectra into an array
    pho_arr = pho_str.spec
    ; Save parameters into vectors
    pho_teff = float(pho_str.teff)
    pho_logg = float(pho_str.logg)
    pho_feh  = float(pho_str.feh)
    ;turn -0 to 0
    pho_feh_0_idx = where(pho_feh eq 0, zcnt)
    if zcnt gt 0 then pho_feh[pho_feh_0_idx] = 0.0

    ; Save to reduce read in time
    save, /variables, filename = savfile
endif else restore, savfile
;HARDCODED
;INPUTS
;npar = 30
;mpfit options
maxiter = 200 ;(default = 200)
quiet = 0
myfunct = 'fit_spectrum'

output_model=0

npix_lsf = 49L
vis =0       
npar = 37


; Select chip
case color of
    'blue': begin
        first_pix=0L
        npix_select = 3500L
    end
    'green': begin
        first_pix=3501L
        npix_select = 2700L
    end
    'red': begin
        first_pix = 6201L
        npix_select = 2370L
    end
endcase
;LOAD necessary data
best_sdss_idx = [51,69,197,211,326,356,362,380,$
                 622,740,825,826,982,1060,1061,$
                 1062,1076,1216,1218,1269]
nsdss = n_elements(best_sdss_idx) 

;get flux and then normalize
;mask
pixmask_arr = sdss_str[best_sdss_idx].mask_pixwt
mask_bits = [0,12]
mask_dec = long(total(2L^mask_bits))
mask_arr = pixmask_arr and mask_dec
mask_idx = where(mask_arr ne 0,nmasked)
if nmasked gt 0 then begin
    mask_arr[mask_idx] = 1
endif

sdss_flux = sdss_str[best_sdss_idx].flux_pixwt
sdss_e_flux = sdss_str[best_sdss_idx].e_flux_pixwt
sdss_spec_arr = sdss_flux
sdss_err_arr = sdss_e_flux
; normalize
for snum = 0, n_elements(sdss_flux[0,*])-1 do begin
    spec_str = sdss_str[snum]
    spec = spec_str.flux_pixwt
    err = spec_str.e_flux_pixwt
    ; Remove continuum
    bmin = spec_str.brange[0]
    bmax = spec_str.brange[1]
    gmin = spec_str.grange[0]
    gmax = spec_str.grange[1]
    rmin = spec_str.rrange[0]
    rmax = spec_str.rrange[1]
    
    nblue  = bmax - bmin+1
    ngreen = gmax - gmin+1
    nred   = rmax - rmin+1
    
    bnorm = continuum_fit(dindgen(nblue), spec[bmin:bmax])
    gnorm = continuum_fit(dindgen(ngreen), spec[gmin:gmax])
    rnorm = continuum_fit(dindgen(nred), spec[rmin:rmax])

    ; Normalize and scale errors for each chip
    spec[bmin:bmax] = spec[bmin:bmax]/bnorm
    spec[gmin:gmax] = spec[gmin:gmax]/gnorm
    spec[rmin:rmax] = spec[rmin:rmax]/rnorm
    
    err[bmin:bmax] = err[bmin:bmax]/bnorm
    err[gmin:gmax] = err[gmin:gmax]/gnorm
    err[rmin:rmax] = err[rmin:rmax]/rnorm
    
    ; Zero the regions in between
    spec[0:bmin-1] = 1d0
    spec[bmax+1:gmin-1] = 1d0
    spec[gmax+1:rmin-1] = 1d0
    spec[rmax+1:-1] = 1d0
    
    err[0:bmin-1] = 1d8
    err[bmax+1:gmin-1] = 1d8
    err[gmax+1:rmin-1] = 1d8
    err[rmax+1:-1] = 1d8

    sdss_spec_arr[*,snum] = spec
    sdss_err_arr[*,snum] = err
endfor


for sdss_idx = 0, nsdss-1 do begin

    sdss_overall_idx = best_sdss_idx[sdss_idx]
    sdss_spec = sdss_spec_arr[*,sdss_idx]
    sdss_err = sdss_err_arr[*,sdss_idx]

    ;apply mask
    sdss_mask = mask_arr[*,sdss_idx]
    mask_idx = where(sdss_mask ne 0, nmask)
    if nmask gt 0 then begin
        sdss_err[mask_idx] = 1d8
        sdss_spec[mask_idx] = !values.d_nan
        sdss_spec = interpol(sdss_spec, sdss_wl, sdss_wl, /nan)
    endif





    ; Load LSF
    lsf_file = '/home/stgilhool/APOGEE/APOGEE_data/lsf_coeffs.fits'
    lsf_coeff_arr = mrdfits(lsf_file,0)
    ;make lsf
    lsf_coeffs = lsf_coeff_arr[*,sdss_overall_idx]
    lsf = make_lsf(lsf_coeffs, first_pix, first_pix+npix_select-1, npix_lsf, oversamp)
    
;;; Make parinfo
    parinfo_str = {VALUE:0d0,$
                   FIXED:0, $
                   LIMITED:[0,0], $
                   LIMITS:[0d0,0d0], $
                   PARNAME:'', $
                   STEP:0, $
                   RELSTEP:0, $
                   MPSIDE:0}
    
    parinfo = replicate(parinfo_str, npar)
    ;Vals
    parinfo[0].value = 1d0
    parinfo[1:-1].value = 1d0
    ;Limits
    ;parinfo[0].limited[0] = 1
    ;parinfo[0].limits[0] = 0
    
    parinfo[1:-1].limited[0] = 1
    parinfo[1:-1].limited[1] = 1
    
    parinfo[1:-1].limits[0] = 0.8d0
    parinfo[1:-1].limits[1] = 1.2d0
    
;;; Make Functargs
    ; Adjust for just one piece of spectrum
    pho_wl_idx_0 = first_pix*oversamp
    pho_wl_idx_1 = pho_wl_idx_0 + (npix_select*oversamp)-1
    pho_wl_select = pho_wl[pho_wl_idx_0:pho_wl_idx_1]
    ;pho_spec = pho_spec[pho_wl_idx_0:pho_wl_idx_1]
    
    sdss_wl_idx_0 = first_pix
    sdss_wl_idx_1 = first_pix + npix_select -1
    sdss_wl_select = sdss_wl[sdss_wl_idx_0:sdss_wl_idx_1]
    sdss_spec_select = sdss_spec[sdss_wl_idx_0:sdss_wl_idx_1]
    sdss_err_select = sdss_err[sdss_wl_idx_0:sdss_wl_idx_1]
    
    
    ; Get approprate pho spec and loop
    teff = [2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100]
    
    nteff = n_elements(teff) 
    feh = [-1.0,-0.5,0.0,0.3]
    
    nfeh = n_elements(feh) 
    logg=5.0
    chi2 = []
    result_arr = []
    model_out = []
    status_arr = []
    
    for teff_idx = 0, nteff-1 do begin
        for feh_idx = 0, nfeh-1 do begin
            teff_i = teff[teff_idx]
            feh_i = feh[feh_idx]
            pho_idx = where(pho_teff eq teff_i and pho_logg eq logg and pho_feh eq feh_i, nres)
            if nres eq 1 then pho_spec=pho_arr[*,pho_idx] $
            else begin
                ;message, "No template match"
                print, "WARNING: No template match"
                chi2 = [chi2, !values.d_nan]
                result = replicate(-1,npar)
                model = replicate(1,npix_select)
                status = -1
                result_arr = [[result_arr],[result]]
                model_out = [[model_out],[model]]
                status_arr = [status_arr, status]
                continue
            endelse
            
            pho_spec = pho_spec[pho_wl_idx_0:pho_wl_idx_1]
            
            
            functargs = {pho_wl:pho_wl_select, $
                         pho_spec:pho_spec, $
                         sdss_wl:sdss_wl_select, $
                         sdss_spec:sdss_spec_select, $
                         sdss_err:sdss_err_select, $
                         lsf:lsf, $
                         oversamp:oversamp, $
                         output_model:output_model, $
                         first_pix:first_pix, $
                         npix_select:npix_select, $
                         vis:vis}
            
            
            !p.multi = [0,1,2]
            result = MPFIT(MYFUNCT, $
                           FUNCTARGS=functargs, $
                           PARINFO=parinfo, $
                           MAXITER=maxiter, $
                           QUIET=quiet, $
                           STATUS=status, $
                           BESTNORM=bestnorm, $
                           NFEV=nfev, $
                           NITER=niter)
            !p.multi = 0               
            print, bestnorm
            chi2 = [chi2, double(bestnorm)]
            if n_elements(result) eq 1  then begin
                result = replicate(-1,npar)
                model = replicate(1,npix_select)
            endif else begin
                functargs.output_model = 1
                model = call_function(myfunct, result, _extra=functargs)
                functargs.output_model = 0
            endelse 
            result_arr = [[result_arr],[result]]
            model_out = [[model_out],[model]]
            status_arr = [status_arr, status]
            
        endfor
    endfor
    outstr = {sdss_idx:sdss_overall_idx, $
              model_teff:teff, $
              model_feh:feh, $
              model_logg:logg, $
              model_arr:model_out, $
              sdss_spec:sdss_spec_select, $
              sdss_err:sdss_err_select, $
              sdss_wl:sdss_wl_select, $
              status_arr:status_arr, $
              result_arr:result_arr, $
              chi2:chi2}
    outfile = '/home/stgilhool/APOGEE/teff_results/teff_fit_'+ $
      color+'_'+strtrim(sdss_overall_idx,2)+'.fits'
    mwrfits, outstr, outfile, /create
endfor    
stop
end












;;; Function to interpret results
; This should be callable on its own (without doing the fitting first)
; read in mpfit results
; plot the data and models
; report parameters, chi2 and Teff
;;; Find some way of converging on a final Teff
; easiest: for each logg and Fe/H, do a linear interpolation
;          and pick the one with the lowest chi2



;;; Debugging
; plot the fitting to make sure it's working
; output the data, model, chi2, parameters
; have something to visualize that stuff







pro pho_spec_fit, color

call_mpfit, color

end

