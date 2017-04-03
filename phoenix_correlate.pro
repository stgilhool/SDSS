; Fit the VALD line list to PHOENIX templates
function model_lines, p, nlines = nlines, ngauss=ngauss, pho_wl=pho_wl, pho_spec=pho_spec, vis=vis

common misc, wl0, wl1

; Create continuum
npix = n_elements(pho_spec) 
model = replicate(1d0, npix)

for line = 0, nlines-1 do begin
    ; Get the gaussian parameters
    gparams = p[line*ngauss:(line+1)*ngauss-1]
    ; Create gaussian
    gauss_line = gaussian(pho_wl,gparams, /double)
    ; Subtract from continuum
    model = model - gauss_line
endfor

; Scale the model to fit the PHOENIX flux
continuum_coeffs = p[nlines*ngauss:-1]
model_scale = poly(pho_wl, continuum_coeffs)
model_lines = model*model_scale

; Calculate fit
res = pho_spec - model_lines

err = (1d-2 * exp(1d0-model)) * model_scale

dev = res/err

if vis then begin
    ;plot, pho_wl, pho_spec, /yno, xr=[15154, 16969], /xs
    plot, pho_wl, pho_spec, /yno, xr=[wl0, wl1], /xs
    oplot, pho_wl, model_lines, color=200
    ;plot, pho_wl, dev, ps = 6, symsize = 0.5, xr=[15154, 16969], /xs
    plot, pho_wl, dev, ps = 6, symsize = 0.5, xr=[wl0, wl1], /xs
;    wait, 0.05
endif

return, dev

end

; Read in the PHOENIX files, interpolate, normalize
function readin_phoenix, teff, feh, logg, wl_grid

base_path = '/home/stgilhool/PHOENIX/'
outstruct = []
xvector = dindgen(n_elements(wl_grid))

foreach logg_i, logg, logg_idx do begin
    
    ;logg_path = base_path + 'logg' + strtrim(logg_i,2) + '/'
    logg_path = base_path + 'logg' + logg_i + '/'
    
    foreach teff_i, teff, teff_idx do begin
        
        teff_name = strtrim(teff_i,2)
        
        foreach feh_i, feh, feh_idx do begin
        
            file_name = teff_name + feh_i + '_result.fits'
            file = logg_path + file_name
            fcheck = file_test(file)
            if fcheck then begin
                phospec = mrdfits(file,1)
                wl = phospec.wave
                spec = phospec.flux
                ; ADDED on 6/7 pretty sure it's right
                df = -8d0 ;conversion factor
                realspec = 10d0^(spec+df)
                spec = realspec
                ;end ADDED
                ; Interpolate onto wl_grid
                outspec = interpol(spec, wl, wl_grid)
                ; Normalize
                norm = continuum_fit(xvector, outspec)
                specnorm = outspec/norm

                struct = {spec:specnorm, feh:feh_i, teff:teff_i, logg:logg_i}
                outstruct = [outstruct, struct]
            endif
        endforeach
    endforeach
endforeach
            
return, outstruct

end


pro phoenix_correlate

common misc, wl0, wl1

linelist = '/home/stgilhool/PHOENIX/line_list_phoenix.txt'

phoenix_file = '/home/stgilhool/PHOENIX/logg4.5/3500-0.0_result.fits'



; Options
ngauss = 3
nstellar = 4
npoly = 3
vis = 1


; Read in phoenix file
;readin if there's no save file.  if there is, just load the save
savfile = 'phoenix_spectra.sav'
if file_test(savfile) eq 0 then begin
    ; WL range option
    wl0 = 15154.
    wl1 = 16970.
    ;wl0 = 15250
    ;wl1 = 15450

    ; Make common PHOENIX wl grid
    delta_wl = wl1-wl0
    d_wl = 0.01d0
    npix_pho = delta_wl/d_wl
    pho_wl_grid = dindgen(npix_pho)*d_wl + wl0
    
    ; Define the grid point vectors
    teff_vec = lindgen(23)*100 + 2000L
    logg_vec = ['4.5', '5.0']
    feh_vec = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', '+0.3', '+0.5'] 
        
    ; Read
    pho_str = readin_phoenix(teff_vec, feh_vec, logg_vec, pho_wl_grid)
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

; Read in linelist
readcol, linelist, wl_ctr, species, format = 'D,A', delimiter = "'", nlines=nlines
; truncate
line_idx = where(wl_ctr gt wl0 and wl_ctr lt wl1, nlines)
wl_ctr = wl_ctr[line_idx]
species = species[line_idx]

;;; Calculate EW in the PHOENIX spectra by smoothing with boxcar
bwidth = 25
;pho_depth = d_wl*(1d0-temporary(pho_arr))
pho_depth = d_wl*(1d0-pho_arr)
pho_ew = smooth(pho_depth, bwidth) * bwidth

;;; Calculate Pearson correlations
teff_corr = dblarr(npix_pho)
logg_corr = dblarr(npix_pho)
feh_corr  = dblarr(npix_pho)

for pixel = 0, npix_pho-1 do begin
    ew_i = reform(pho_ew[pixel,*])
    teff_corr[pixel] = find_correlation(ew_i, pho_teff)
    logg_corr[pixel] = find_correlation(ew_i, pho_logg)
    feh_corr[pixel]  = find_correlation(ew_i, pho_feh)
endfor

!p.multi=0
window, 0
plot, teff_corr, ps=3
oplot, logg_corr, ps=3, color=200
oplot, feh_corr, ps=3, color=99999

;;; Analyze
;;; Find where the VALD lines land on the pho_wl_grid
line_idx = []
interact = 1

!p.multi = [0,2,2]
for vline = 0, nlines-1 do begin
    wl_diff = abs(pho_wl_grid - wl_ctr[vline])
    pho_idx = where(wl_diff eq min(wl_diff), nmin)
    ;if nmin eq 1 then line_wl = pho_wl_grid[pho_idx] $
    if nmin eq 1 then line_idx = [line_idx, pho_idx[0]] $
      else if nmine gt 1 then begin
        print, "WARNING:Not a unique match for VALD line wl" $
          + strtrim(wl_ctr[vline],2)
        line_idx = [line_idx, pho_idx[0]]
        dxstop
    endif else message, "No match for VALD line wl = " $
      + strtrim(wl_ctr[vline],2)
    
    ; get line info
    spname = species[vline]
    ctrwl  = wl_ctr[vline]

    ; set up indices to look at features
    wide0 = pho_idx-(5*(bwidth/2))
    wide1 = pho_idx+(5*(bwidth/2))
    nwide = wide1-wide0 + 1
    narr0 = pho_idx-(1*(bwidth/2))
    narr1 = pho_idx+(1*(bwidth/2))
    nnarr = narr1-narr0 + 1

    ;line_idx = [line_idx, line_wl]
    if interact then begin
        ew_line = reform(pho_ew[pho_idx,*])
        window, 1, xs = 1500, ys=900, title = spname + ' at wl = ' + strtrim(ctrwl,2)
        plot, ew_line, pho_teff, ps=3, tit='Teff corr', xtit='ew', ytit='Teff', yr=[1900,4300], /ys
        plot, ew_line, pho_logg, ps=3, tit='logg corr', xtit='ew', ytit='logg', yr=[4,5.5], /ys
        plot, ew_line, pho_feh, ps=3, tit='[Fe/H] corr', xtit='ew', ytit='[Fe/H]'

        teff_corr_plot = teff_corr[pho_idx-bwidth/2:pho_idx+bwidth/2]
        logg_corr_plot = logg_corr[pho_idx-bwidth/2:pho_idx+bwidth/2]
        feh_corr_plot  = feh_corr[pho_idx-bwidth/2:pho_idx+bwidth/2]
        yrmin = min(teff_corr_plot) < min(logg_corr_plot) < min(feh_corr_plot)
        yrmax = max(teff_corr_plot) > max(logg_corr_plot) > max(feh_corr_plot)
        ;plot, teff_corr_plot, tit='Corr Coeffs', yr = [yrmin, yrmax]   
        ;plot, teff_corr_plot, tit='Corr Coeffs', yr = [-1,1]   
        ;oplot, logg_corr[pho_idx-bwidth/2:pho_idx+bwidth/2], linest=2
        ;oplot, feh_corr[pho_idx-bwidth/2:pho_idx+bwidth/2], linest=3
        plot, rebin(pho_wl_grid[wide0:wide1],nwide,nspec), pho_arr[wide0:wide1,*], ps=3, yr = [0.94,1.02], /ys
        oplot, rebin(pho_wl_grid[narr0:narr1],nnarr,nspec), pho_arr[narr0:narr1,*], ps=3, color=200
        
        dxstop
    endif
endfor
!p.multi=0

stop


apg_idx = where(pho_wl ge wl0 and pho_wl le wl1, npix)

; Make line structure
linestr_fid = {LINENUM:0L, $
               WL_CTR:0d0, $
               SPECIES:'', $
               STATUS:0, $
               STELLAR_PAR:replicate(0d0,nstellar), $
               GAUSS_PAR:replicate(0d0,ngauss), $
               PAR_ERR:replicate(0d0,ngauss), $
               FIT_ERR:0d0}
           
linestr = replicate(linestr_fid, nlines)

;Populate
linenum = lindgen(nlines)
stellar_par = [3500, 0.0, 4.5, 0]

linestr.linenum = linenum
linestr.wl_ctr = wl_ctr
linestr.species = species
linestr.stellar_par = stellar_par

; Set up call to mpfit
parinfo_fid = {value:0d0, $
               mpside:0, $
               limited:[0,0], $
               limits:[0d0,0d0], $
               fixed:0, $
               step:0d0, $
               relstep:0d0}

parinfo = replicate(parinfo_fid, ngauss*nlines+npoly)

;;; FILL IN PARINFO
amp_idx = lindgen(nlines)*ngauss
ctr_idx = lindgen(nlines)*ngauss + 1
sig_idx = lindgen(nlines)*ngauss + 2

;parameter value guesses
linestr.gauss_par[0] = replicate(0.01d0, nlines)
linestr.gauss_par[1] = wl_ctr
linestr.gauss_par[2] = replicate(0.1d0, nlines)
gauss_guess = reform(linestr.gauss_par, ngauss*nlines)
poly_guess = [13.62d0, 1d-6, 1d-10]
guess = [gauss_guess, poly_guess]               
parinfo.value = guess

;PARAMETER LIMITS
;amp btw 0 and 1
parinfo[amp_idx].limited = [1,1]
parinfo[amp_idx].limits  = [0d0, 1d0]

;ctr can change by +/-0.3A
parinfo[ctr_idx].limited = [1,1]
ctr_vals = parinfo[ctr_idx].value
ctr_lolim = ctr_vals - 0.3d0
ctr_hilim = ctr_vals + 0.3d0
parinfo[ctr_idx].limits[0] = ctr_lolim
parinfo[ctr_idx].limits[1] = ctr_hilim

;sig is gt 0 and lt 2
parinfo[sig_idx].limited = [1,1]
parinfo[sig_idx].limits = [0d0, 2d0]

;;; functargs
functargs = {nlines:nlines, $
             ngauss:ngauss, $
             pho_wl:pho_wl, $
             pho_spec:pho_spec, $
             vis:vis}

if vis then !p.multi = [0,1,2] else !p.multi = 0



;;;MPFIT
result = mpfit('model_lines', parinfo=parinfo, functargs=functargs, status=status, $
               bestnorm=chi2)

; Save results in linestr
linestr.status = status

param_results = result[0:-1*(npoly+1)] 
param_results = reform(param_results, ngauss, nlines)

linestr.gauss_par = param_results

; Flag lines with bad amps, ctrs, and sigs
bad_amp_idx = where(linestr.gauss_par[0] eq 0, nbdamp)
if nbdamp gt 0 then linestr[bad_amp_idx].status = -1

bad_ctr_idx = where(abs(linestr.gauss_par[1]-wl_ctr) ge 0.20d0, nbdctr)
if nbdctr gt 0 then linestr[bad_ctr_idx].status = -2

bad_sig_idx = where(linestr.gauss_par[2] eq 2.0, nbdsig)
if nbdsig gt 0 then linestr[bad_sig_idx].status = -3

; report good lines
gd_idx = where(linestr.status ge 0, ngd)
if ngd gt 0 then newlinestr = linestr[gd_idx] else newlinestr = linestr[gd_idx]

help, newlinestr, /str

stop

end
