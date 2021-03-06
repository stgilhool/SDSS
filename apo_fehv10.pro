;+
; NAME:
;
;
;
; PURPOSE:
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-
pro vline, val,_extra=extra, min=min, max=max
if !y.type eq 1 then yrange = 10^!y.crange else yrange = !y.crange
nv = n_elements(val)
for i = 0, nv-1 do oplot,fltarr(2)+val[i],yrange,_extra=extra
end

function feh_fit, p, ew=ew, feh=feh, err=err, vis=vis

common outputs, model, res

feh_col = reform(feh, 1, n_elements(feh))

model = ew ## p

model = reform(model)

data = feh

res = data - model

dev = res/err

chi2 = sqrt(total(dev^2, /double))

sorti = sort(data)
sortd = data[sorti]
sortm = model[sorti]
sortr = res[sorti]

if vis eq 1 then begin
    ;plot, data, ps = 6, xtitle = "File number", ytitle = "Metalicity", title = strtrim(chi2,2), /xs
    ;oplot, model, ps = 6, color = 200
    ;plot, res, ps = 6, title = "Residuals", /xs
    plot, sortd, ps = 6, xtitle = "File number", ytitle = "Metalicity", title = strtrim(chi2,2), /xs
    oplot, sortm, ps = 6, color = 200
    plot, sortr, ps = 6, title = "Residuals with RMSE: "+strtrim(stddev(res),2), /xs
endif

return, dev

end




pro apo_fehv10, pthresh

if n_params() eq 0 then pthresh = 0.3
;goto, test

; Constants
bigc = 299792.458 ;km/s

; Options
vis = 1

; Paths
root_dir = '/home/stgilhool/APOGEE/'
data_path = root_dir + 'APOGEE_data/'

listfile = data_path + 'mdwarflist.txt'

; Readin file names
readcol, listfile, datafile, format="A", count=nfiles

files = data_path + datafile

files_truncated = strarr(nfiles)
aspfiles = strarr(nfiles)


for i = 0, nfiles-1 do begin

    file_split = strsplit(datafile[i], '/', /extract)
    location = file_split[0]
    files_truncated[i] = file_split[1]

    mass2namesplit = strsplit(files_truncated[i], '-', /extract)
    mass2name = mass2namesplit[-1]
    aspfiles[i] = data_path + location + '/aspcapStar-r5-v603-' + mass2name
    print, files[i]
    print, files_truncated[i]
    print, aspfiles[i]

endfor


; Get ready to read in headers, vhelio, wl coeffs
badspecflag = mrdfits('badspecflag.fits')

; KEEP ALL SPECTRA FOR NOW
badspecflag = replicate(0, n_elements(badspecflag))
gdidx = where(badspecflag eq 0, ngd, complement=bdidx, ncomplement=nbd)


; Readin irtf data
irtf = read_irtf_data()

; Save only the good ones (taken from manual flagging)
irtf = irtf[gdidx]

n_irtf = n_elements(irtf)
irtf_files = strtrim(irtf.APG_file,2)
feh_afk = irtf.AFK
feh_afke= irtf.AFKE

teff_vk = irtf.teff_vk
teff_vkb = irtf.teff_vkb
teff_vks = irtf.teff_vks

v_final = irtf.v_final
v_src = irtf.v_src

kmag = irtf.kmag
logg = irtf.apg_rv_logg

n = n_irtf



header = ptrarr(n, /allocate_heap)
spectra = dblarr(8575L,n)
nx = 8575L
x = dindgen(nx)

vharr = dblarr(n)
wl0arr = dblarr(n)
wl1arr = dblarr(n)

irtf_apo_idx = lonarr(n)
blue_range  = lonarr(3, n)
green_range = lonarr(3, n)
red_range   = lonarr(3, n)
;logwl_grid = dblarr(nx, n)
;log_shift = dblarr(n)

window, 0, xsize = 1400, ysize=1000



for i = 0, n_irtf-1 do begin

    ; Find the matching APOGEE file
    apo_null = where(strmatch(files_truncated, irtf_files[i]) eq 1, matchcnt)
    if matchcnt eq 1 then apo_idx = apo_null else message, "file match messed up"

    irtf_apo_idx[i] = apo_idx

    null = mrdfits(files[apo_idx],0,head)
    ; save header
    *(header[i]) = head

    ;print, head
    ;stop


    ; Read in 1d spectrum array for a file
    specall = mrdfits(files[apo_idx],1)

    ; Save only the first (combined) spectrum
    spec = specall[*,0]
    ; normalize it
    ;snorm = continuum_fit(x, spec)
    ;spec = temporary(spec)/temporary(snorm)

    ; Get vhelio, and wl coeffs for a file
    ; (vhelio is unnecessary, and wl coeffs are same for all)
    vhelio = fxpar(head, 'VHELIO', datatype=0.0d)
    wl0    = fxpar(head, 'CRVAL1', datatype=0.0d)
    wl1    = fxpar(head, 'CDELT1', datatype=0.0d)
  
    ; Remove continuum
    bmin = fxpar(head, 'BOVERMIN')
    bmax = fxpar(head, 'BOVERMAX')
    gmin = fxpar(head, 'GOVERMIN')
    gmax = fxpar(head, 'GOVERMAX')
    rmin = fxpar(head, 'ROVERMIN')
    rmax = fxpar(head, 'ROVERMAX')

    ;bmin = fxpar(head, 'BMIN')
    ;bmax = fxpar(head, 'BMAX')
    ;gmin = fxpar(head, 'GMIN')
    ;gmax = fxpar(head, 'GMAX')
    ;rmin = fxpar(head, 'RMIN')
    ;rmax = fxpar(head, 'RMAX')

    nblue  = bmax - bmin+1
    ngreen = gmax - gmin+1
    nred   = rmax - rmin+1
    
    ;bnorm = continuum_fit(dindgen(nblue), spec[bmin:bmax])
    ;gnorm = continuum_fit(dindgen(ngreen), spec[gmin:gmax])
    ;rnorm = continuum_fit(dindgen(nred), spec[rmin:rmax])

    ;spec[bmin:bmax] = spec[bmin:bmax]/bnorm
    ;spec[gmin:gmax] = spec[gmin:gmax]/gnorm
    ;spec[rmin:rmax] = spec[rmin:rmax]/rnorm

    ; Save quantities
    vharr[i]  = vhelio
    wl0arr[i] = wl0    
    wl1arr[i] = wl1
    spectra[*,i] = spec
    blue_range[*,i] = [bmin, bmax, nblue]
    green_range[*,i] = [gmin, gmax, ngreen]
    red_range[*,i] = [rmin, rmax, nred]

    ;if i eq 0 then plot, spec, /xs else oplot, spec, color=200
    ;wait, 0.1

endfor

; Get a common chunk for each color
; JUST BLUE FOR NOW
bmin = max(blue_range[0,*], /nan)
bmax = min(blue_range[1,*], /nan)
nblue = bmax - bmin + 1
gmin = max(green_range[0,*], /nan)
gmax = min(green_range[1,*], /nan)
ngreen = gmax - gmin + 1
rmin = max(red_range[0,*], /nan)
rmax = min(red_range[1,*], /nan)
nred = rmax - rmin + 1

; Make Blue spectra array
bspec = spectra[bmin:bmax, *]
gspec = spectra[gmin:gmax, *]
rspec = spectra[rmin:rmax, *]

aspec = [bspec, gspec, rspec]

; Sum up chunks

npix_win = 25L
;nwin = nblue - (npix_win-1)
nwinb= nblue - (npix_win-1)
nwing= ngreen - (npix_win-1)
nwinr= nred - (npix_win-1)

;depth = 1d0-bspec
depth = 1d0-aspec
depthb= 1d0-bspec
depthg= 1d0-gspec
depthr= 1d0-rspec

; This gives average value in window centered at each pixel
avg_arr = smooth(depth, [npix_win, 1])

avg_arrb = smooth(depthb, [npix_win, 1])
avg_arrg = smooth(depthg, [npix_win, 1])
avg_arrr = smooth(depthr, [npix_win, 1])
; This is the coordinates of the cntrs of each window
;x_0  = lindgen(nwin) + (npix_win/2)
;xx0  = rebin(x_0, nwin, n)
;yy0  = rebin(reform(lindgen(n),1,n), nwin, n)

x_0b  = lindgen(nwinb) + (npix_win/2)
xx0b  = rebin(x_0b, nwinb, n)
yy0b  = rebin(reform(lindgen(n),1,n), nwinb, n)

x_0g  = lindgen(nwing) + (npix_win/2)
xx0g  = rebin(x_0g, nwing, n)
yy0g  = rebin(reform(lindgen(n),1,n), nwing, n)

x_0r  = lindgen(nwinr) + (npix_win/2)
xx0r  = rebin(x_0r, nwinr, n)
yy0r  = rebin(reform(lindgen(n),1,n), nwinr, n)



;ew_arr  = avg_arr[xx0, yy0] * double(npix_win)
ew_arrb = avg_arrb[xx0b, yy0b] * double(npix_win)
ew_arrg = avg_arrg[xx0g, yy0g] * double(npix_win)
ew_arrr = avg_arrr[xx0r, yy0r] * double(npix_win)

;feh_arr = rebin(reform(feh_afk, 1, n), nwin, n)

    
; Fit the metalicity of each chunk
; Initialize
;pcoeff = dblarr(nwin)

pcoeffb = dblarr(nwinb)
pcoeffg = dblarr(nwing)
pcoeffr = dblarr(nwinr)

; Make median bspec for display
bspec_med = median(bspec, dimension=2)
gspec_med = median(gspec, dimension=2)
rspec_med = median(rspec, dimension=2)




; Calculate Pearson correlation coeff for each window
for i = 0, nwinb-1 do begin
    
    ew = reform(ew_arrb[i,*])
    mean_ew = mean(ew)
    diff_ew = ew - mean_ew

    feh = feh_afk
    mean_feh = mean(feh)
    diff_feh = feh - mean_feh

    r = total(diff_ew*diff_feh, /double)/(sqrt(total(diff_ew^2, /double)) * sqrt(total(diff_feh^2, /double)))
    
    ; Record results
    pcoeffb[i] = r

endfor

for i = 0, nwing-1 do begin
    
    ew = reform(ew_arrg[i,*])
    mean_ew = mean(ew)
    diff_ew = ew - mean_ew

    feh = feh_afk
    mean_feh = mean(feh)
    diff_feh = feh - mean_feh

    r = total(diff_ew*diff_feh, /double)/(sqrt(total(diff_ew^2, /double)) * sqrt(total(diff_feh^2, /double)))
    
    ; Record results
    pcoeffg[i] = r

endfor




for i = 0, nwinr-1 do begin


    
    ew = reform(ew_arrr[i,*])
    mean_ew = mean(ew)
    diff_ew = ew - mean_ew

    feh = feh_afk
    mean_feh = mean(feh)
    diff_feh = feh - mean_feh

    r = total(diff_ew*diff_feh, /double)/(sqrt(total(diff_ew^2, /double)) * sqrt(total(diff_feh^2, /double)))
    
    ; Record results
    pcoeffr[i] = r

endfor

;stop


save, /variables, filename='correlation_v10.sav'
test:
restore, 'correlation_v10.sav'


!p.multi=0
if vis eq 1 then begin
    window, 0, xsize=1300, ysize=1000
    ;!p.multi=[0,1,2]
endif


!p.multi = [0,1,3]
plot, x_0b, pcoeffb, /xs, ps=6, title='Pearson Regression Coefficients (blue)'
oplot, x_0b, replicate(pthresh, n_elements(x_0b))
plot, x_0g, pcoeffg, /xs, ps=6, title='Pearson Regression Coefficients (green)'
oplot, x_0g, replicate(pthresh, n_elements(x_0g))
plot, x_0r, pcoeffr, /xs, ps=6, title='Pearson Regression Coefficients (red)'
oplot, x_0r, replicate(pthresh, n_elements(x_0r))

pcoeffall=dblarr(n_elements(spectra[*,0]))
pcoeffall[bmin+(npix_win/2):bmax-(npix_win/2)] = pcoeffb
pcoeffall[gmin+(npix_win/2):gmax-(npix_win/2)] = pcoeffg
pcoeffall[rmin+(npix_win/2):rmax-(npix_win/2)] = pcoeffr

;stop

!p.multi = 0

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
ew_corr_y = [reform(replicate(1d0,n),1,n),ew_corr_y0]
help, ew_corr_y

; Add column for v-k
;ew_corr_y = [ew_corr_y, reform((v_final-kmag), 1, n)]
;ew_corr_y = [ew_corr_y, reform((teff_vk), 1, n)]

;help, ew_corr_y
;stop
!p.multi=0
if vis eq 1 then begin
    window, 0, xsize=1300, ysize=1000
    !p.multi=[0,1,2]
endif
    
;;; Fit the metalicity of each chunk

common outputs, model, res

functargs = {ew:ew_corr_y, $
             feh:feh_afk, $
             err:feh_afke, $
             vis:vis $
            }
    
parinfo = replicate({value:0.1d0},n_ewc+1)


r = mpfit('feh_fit', parinfo=parinfo, functargs=functargs, status=status, dof=dof, bestnorm=chi2)


print, "result:"
print, r

outlier = where(abs(res) ge 0.2, ocount)
print, ocount
print, outlier
stop

outstr = {result:r, $
          model:model, $
          data:feh_afk, $
          res:res, $
          err:feh_afke, $
          nlines:n_ewc, $
          x0:ew_corr_x0, $
          fp:ew_corr_fp, $
          lp:ew_corr_lp, $
          width:npix_win, $
          status:status, $
          pthresh:pthresh, $
          chi2:chi2, $
          chi2dof:chi2/dof $
          }
         
          
outpath = '/home/stgilhool/APOGEE/feh_results/'
outfile = outpath + 'asp_feh_fit_' +strtrim(fix(pthresh*100),2)+ '.fits'
mwrfits, outstr, outfile, /create


reg = regress(ew_corr_y[1:*,*], feh_afk, measure_errors=feh_afke, sigma=parsigma, const=const, chisq=chisq, /double, correlation=correlation, mcorrelation=mcorrelation, ftest=ftest, status=regstatus, yfit=regfit)

window, 2, xsize=1300, ysize=1000
!p.multi = [0,1,2]

sortfeh = sort(feh_afk)
sfeh = feh_afk[sortfeh]
sregfit = regfit[sortfeh]
sres = sfeh-sregfit

plot, sfeh, ps = 6, xtitle = "File number", ytitle = "Metalicity", title = strtrim(chisq,2), /xs
oplot, sregfit, ps = 6, color = 200
plot, sres, ps = 6, title = "Residuals with RMSE: "+strtrim(stddev(sres),2), /xs


window, 3, xsize=1300, ysize=1000
!p.multi=0
;plot, logg[sortfeh],ps=6, /xs

plot, v_final[sortfeh]-kmag[sortfeh], sres, ps=6, /xs
;plot, teff_vk[sortfeh], ps=6, /xs
;plot, sres, teff_vk[sortfeh], ps = 6, /xs



print, "IDL REGRESS"
print, "Result: ", [const, reform(reg)]
print, "Correlation: "
print, correlation

print, v_src

stop



;if n_elements(r) eq 1 then begin
;    coeff[*,i] = !values.d_nan
;endif else begin
;    coeff[*,i] = r
;endelse

    

; !p.multi=0
; stop

; spec_med = median(spectra, dimension=2)
; sdall = stddev(spectra, dimension=2)
; for i = 0, n_ewc-1 do begin
;     ewx = lindgen(npix_win) + ew_corr_fp[i]
;     spec = spec_med[ew_corr_fp[i]:ew_corr_lp[i]]
;     sd = sdall[ew_corr_fp[i]:ew_corr_lp[i]]
    
;     plot, ewx, spec, /xs, yr=[0.6,1.1], title = "Chunk of spectrum corresponding to x_0="+strtrim(ew_corr_x0[i],2)+" with Rbest="+strtrim(r[i+1],2)+" and P="+strtrim(pcoeffall[ew_corr_x0[i]],2)
;     errplot, ewx, spec-sd, spec+sd

;     dxstop
; endfor



; stop



end
