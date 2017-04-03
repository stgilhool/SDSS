pro calc_feh, pthresh

inpath = '/home/stgilhool/APOGEE/feh_results/'

; if pthresh (a particular result file) is not specified
; then read them all in and pick the best by chi2dof
if n_params() eq 0 then begin
    pthresh_vec = [0,1,3,4,5,6,7,8,9,10,15,20,25,30,35]
    npthresh = n_elements(pthresh_vec) 
    chi2dof_test = dblarr(npthresh)
    for i = 0, npthresh-1 do begin
        pthresh = pthresh_vec[i]
        infilename = 'feh_fit_'+strtrim(pthresh,2)+'.fits'
        infile = inpath + infilename
        instr = mrdfits(infile,1)
        chi2dof_test[i] = instr.chi2dof
    endfor
    ; just set pthresh to 1 for now FIX
    pthresh = 1
endif

pthresh = 100 ; this corresponds to apo_fehv9 output, which has fixed x0's
                                ;;and doesn't actually care about pthresh
infilename = 'feh_fit_'+strtrim(pthresh,2)+'.fits'
infile = inpath + infilename

instr = mrdfits(infile, 1)

napgspec = 1350
npix_apg = 8575L

nlines = instr.nlines
x0 = instr.x0
fp = instr.fp
lp = instr.lp
npix_win = instr.width
regcoeff = instr.result

; Read in all APOGEE spectra
apg_str = readin_apo(nfiles = napgspec, /flux)
spectra = dblarr(npix_apg, napgspec)
blue_range  = lonarr(3, napgspec)
green_range = lonarr(3, napgspec)
red_range   = lonarr(3, napgspec)

for i = 0, napgspec-1 do begin

    head = *(apg_str[i].header)
    spec = apg_str[i].spectra

    bmin = fxpar(head, 'BOVERMIN')
    bmax = fxpar(head, 'BOVERMAX')
    gmin = fxpar(head, 'GOVERMIN')
    gmax = fxpar(head, 'GOVERMAX')
    rmin = fxpar(head, 'ROVERMIN')
    rmax = fxpar(head, 'ROVERMAX')

    nblue  = bmax - bmin+1
    ngreen = gmax - gmin+1
    nred   = rmax - rmin+1
    
    bnorm = continuum_fit(dindgen(nblue), spec[bmin:bmax])
    gnorm = continuum_fit(dindgen(ngreen), spec[gmin:gmax])
    rnorm = continuum_fit(dindgen(nred), spec[rmin:rmax])

    spec[bmin:bmax] = spec[bmin:bmax]/bnorm
    spec[gmin:gmax] = spec[gmin:gmax]/gnorm
    spec[rmin:rmax] = spec[rmin:rmax]/rnorm

    spectra[*,i] = spec
    blue_range[*,i] = [bmin, bmax, nblue]
    green_range[*,i] = [gmin, gmax, ngreen]
    red_range[*,i] = [rmin, rmax, nred]
endfor

; Get a common chunk for each color
bmin = max(blue_range[0,*], /nan)
bmax = min(blue_range[1,*], /nan)
nblue = bmax - bmin + 1
gmin = max(green_range[0,*], /nan)
gmax = min(green_range[1,*], /nan)
ngreen = gmax - gmin + 1
rmin = max(red_range[0,*], /nan)
rmax = min(red_range[1,*], /nan)
nred = rmax - rmin + 1

; Split spectra up by chip
bspec = spectra[bmin:bmax, *]
gspec = spectra[gmin:gmax, *]
rspec = spectra[rmin:rmax, *]

; Sum up chunks
npix_win = 25L
nwinb= nblue - (npix_win-1)
nwing= ngreen - (npix_win-1)
nwinr= nred - (npix_win-1)

depthb= 1d0-bspec
depthg= 1d0-gspec
depthr= 1d0-rspec

; This gives average value in window centered at each pixel
avg_arrb = smooth(depthb, [npix_win, 1])
avg_arrg = smooth(depthg, [npix_win, 1])
avg_arrr = smooth(depthr, [npix_win, 1])
; This is the coordinates of the cntrs of each window
x_0b  = lindgen(nwinb) + (npix_win/2)
xx0b  = rebin(x_0b, nwinb, napgspec)
yy0b  = rebin(reform(lindgen(napgspec),1,napgspec), nwinb, napgspec)

x_0g  = lindgen(nwing) + (npix_win/2)
xx0g  = rebin(x_0g, nwing, napgspec)
yy0g  = rebin(reform(lindgen(napgspec),1,napgspec), nwing, napgspec)

x_0r  = lindgen(nwinr) + (npix_win/2)
xx0r  = rebin(x_0r, nwinr, napgspec)
yy0r  = rebin(reform(lindgen(napgspec),1,napgspec), nwinr, napgspec)

; Gives sum of pixels in window centered at each pixel
ew_arrb = avg_arrb[xx0b, yy0b] * double(npix_win)
ew_arrg = avg_arrg[xx0g, yy0g] * double(npix_win)
ew_arrr = avg_arrr[xx0r, yy0r] * double(npix_win)

; Locations of feature chunks
corr_x0b = [720L, 2000L]-(npix_win/2)
ew_corr_xb = corr_x0b + (npix_win/2)
ew_corr_yb = ew_arrb[corr_x0b,*]
    
corr_x0g = [300, 400, 510, 1090, 1170, 1290, 1350, 1390, 1580, $
            1610, 1790, 1810, 2080, 2110, 2130, 2200]-(npix_win/2)
ew_corr_xg = corr_x0g + (npix_win/2)
ew_corr_yg = ew_arrg[corr_x0g,*]

corr_x0r = [12, 100, 120, 150, 890, 910, 950, 1050, 1080]-(npix_win/2)
ew_corr_xr = corr_x0r + (npix_win/2)
ew_corr_yr = ew_arrr[corr_x0r,*]

ew_corr_x0 = [ew_corr_xb, ew_corr_xg, ew_corr_xr]
ew_corr_y0 = [ew_corr_yb, ew_corr_yg, ew_corr_yr]
  
n_ewc = n_elements(ew_corr_x0)

ew_corr_fp = ew_corr_x0 - (npix_win/2)
ew_corr_lp = ew_corr_x0 + (npix_win/2)

; Add column for constant term
ew_corr_y = [reform(replicate(1d0,napgspec),1,napgspec),ew_corr_y0]

; Calculate metallicity
feh_result = ew_corr_y ## regcoeff
feh_result = reform(feh_result)

; Check result

feh_irtf = instr.data
feh_model = instr.model

plot, feh_irtf, feh_model, ps = 6, xr=[-2.5,2.5], yr=[-2.5, 2.5]



stop



outfilename = 'feh_results_100.fits'
outfile = inpath + outfilename
mwrfits, feh_result, outfile, /create
stop

end
