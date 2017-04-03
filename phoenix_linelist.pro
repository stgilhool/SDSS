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






pro phoenix_linelist

common misc, wl0, wl1

linelist = '/home/stgilhool/PHOENIX/line_list_phoenix.txt'
phoenix_file = '/home/stgilhool/PHOENIX/logg4.5/3500-0.0_result.fits'

; Options
ngauss = 3
nstellar = 4
npoly = 3
vis = 1
;wl0 = 15154.
;wl1 = 16970.
wl0 = 15250
wl1 = 15450

; Read in linelist
readcol, linelist, wl_ctr, species, format = 'D,A', delimiter = "'", nlines=nlines
; truncate
line_idx = where(wl_ctr gt wl0 and wl_ctr lt wl1, nlines)
wl_ctr = wl_ctr[line_idx]
species = species[line_idx]

; Read in phoenix file
pho = mrdfits(phoenix_file,1)
pho_wl = pho.wave
apg_idx = where(pho_wl ge wl0 and pho_wl le wl1, npix)
pho_wl = pho_wl[apg_idx]
pho_spec = pho.flux[apg_idx]

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
