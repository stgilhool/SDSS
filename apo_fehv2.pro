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

function feh_fit, p, ew=ew, feh=feh, x_0=x_0, npix_win=npix_win, bspec_med=bspec_med, vis=vis

p0 = p[0]
p1 = p[1]

x_fp = x_0 - (npix_win/2)
x = lindgen(npix_win)+x_fp
spec_chunk = bspec_med[x]

model = p0 + p1*ew

data = feh

res = data - model

if vis eq 1 then begin
    plot, x, spec_chunk, title = "Spectrum Window", /xs, yr=[0.7,1.1]

    plot, ew, data, /xs, ps=6, title="Fit for window at pix: "+strtrim(x_0,2)
    oplot, ew, model, color=200
    
endif

return, res

end




pro apo_fehv2

goto, test

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


for i = 0, nfiles-1 do begin

    file_split = strsplit(files[i], '/', /extract)
    files_truncated[i] = file_split[-1]
endfor

; Readin irtf data
irtf = read_irtf_data()

n_irtf = n_elements(irtf)
irtf_files = strtrim(irtf.APG_file,2)
feh_afk = irtf.AFK
    
; Get ready to read in headers, vhelio, wl coeffs
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
    
    bnorm = continuum_fit(dindgen(nblue), spec[bmin:bmax])
    gnorm = continuum_fit(dindgen(ngreen), spec[gmin:gmax])
    rnorm = continuum_fit(dindgen(nred), spec[rmin:rmax])

    spec[bmin:bmax] = spec[bmin:bmax]/bnorm
    spec[gmin:gmax] = spec[gmin:gmax]/gnorm
    spec[rmin:rmax] = spec[rmin:rmax]/rnorm

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

; Make Blue spectra array
bspec = spectra[bmin:bmax, *]

; Sum up chunks

npix_win = 25L
nwin = nblue - (npix_win-1)

depth = 1d0-bspec

; This gives average value in window centered at each pixel
avg_arr = smooth(depth, [npix_win, 1])
; This is the coordinates of the cntrs of each window
x_0  = lindgen(nwin) + (npix_win/2)
xx0  = rebin(x_0, nwin, n)
yy0  = rebin(reform(lindgen(n),1,n), nwin, n)

ew_arr = avg_arr[xx0, yy0] * double(npix_win)

feh_arr = rebin(reform(feh_afk, 1, n), nwin, n)

!p.multi=0
if vis eq 1 then begin
    window, 0, xsize=1300, ysize=1000
    ;!p.multi=[0,1,2]
endif
    
; Fit the metalicity of each chunk
; Initialize
pcoeff = dblarr(nwin)
scoeff = dblarr(nwin)
kcoeff = dblarr(nwin)

; Make median bspec for display
bspec_med = median(bspec, dimension=2)

; Calculate Pearson correlation coeff for each window
; Also calculate Spearman rank, and Kendall rank
n_kpair = (n*(n-1))/2
feh_krank = dblarr(n_kpair)
inc = 0L
for i = 0, n-2 do begin
    for j = i+1, n-1 do begin
        if feh_afk[i] ge feh_afk[j] then feh_krank[inc] = 1 $
          else if feh_afk[i] lt feh_afk[j] then feh_krank[inc] = 0
        inc++
    endfor
endfor
        
    
for i = 0, nwin-1 do begin
    
    ew = reform(ew_arr[i,*])
    mean_ew = mean(ew)
    diff_ew = ew - mean_ew
    ew_rank = sort(ew)



    feh = feh_afk
    mean_feh = mean(feh)
    diff_feh = feh - mean_feh
    feh_rank = sort(feh)

    

    diff_rank = ew_rank - feh_rank

    r = total(diff_ew*diff_feh, /double)/(sqrt(total(diff_ew^2, /double)) * sqrt(total(diff_feh^2, /double)))
    
    s = 1d0 - ((6d0 * total(diff_rank^2))/(n*(n^2-1)))



    einc = 0L
    ew_krank = dblarr(n_kpair)
    for w = 0, n-2 do begin
        for j = w+1, n-1 do begin
            if ew[w] ge ew[j] then ew_krank[einc] = 1 $
            else if ew[w] lt ew[j] then ew_krank[einc] = 0
            einc++
        endfor
    endfor

    concord_check = feh_krank - ew_krank
    concord_idx = where(concord_check eq 0, n_concord, ncomplement=n_discord)
    
    k = (n_concord - n_discord)/double(n_kpair)

    ; Record results
    pcoeff[i] = r
    scoeff[i] = s
    kcoeff[i] = k

endfor

!p.multi = [0,1,3]
plot, x_0, pcoeff, /xs, ps=6, title='Pearson Regression Coefficients'
plot, x_0, scoeff, /xs, ps=6, title='Spearman Coefficients'
plot, x_0, kcoeff, /xs, ps=6, title='Kendall Coefficients'



stop
!p.multi = 0


test:
restore, 'correlation.sav'

; Find the peaks

npix_pwin = 101L
pwin_step = 10L
n_pwin = (nwin-(npix_pwin+1))/pwin_step

;pkx = lonarr(n_pwin)
ppkx_0 = lonarr(n_pwin)
spkx_0 = lonarr(n_pwin)
kpkx_0 = lonarr(n_pwin)

for i = 0, n_pwin-1 do begin
    
    
    
    
    fpix = i*pwin_step
    lpix = fpix + npix_pwin-1
    
    xvec = lindgen(npix_pwin)+x_0[fpix]
    pvec = pcoeff[fpix:lpix]
    ;svec = scoeff[fpix:lpix]
    kvec = kcoeff[fpix:lpix]

    pfit = mpfitpeak(xvec, pvec, ppkcoeff, nterms=3, /positive, status=pstatus)
    kfit = mpfitpeak(xvec, kvec, kpkcoeff, nterms=3, /positive, status=kstatus)
    
    if ppkcoeff[0] ge 0.15 then ppkx_0[i] = fix(ppkcoeff[1]) else ppkx_0[i] = !values.f_nan
    if kpkcoeff[0] ge 0.10 then kpkx_0[i] = fix(kpkcoeff[1]) else kpkx_0[i] = !values.f_nan

    ;ppkx_0[i] = fix(ppkcoeff[1])
    ;kpkx_0[i] = fix(kpkcoeff[1])

    ;if i lt 5 then stop

endfor

bs = 3
phist = histogram(ppkx_0, binsize = bs, min = x_0[0], max=max(x_0), /nan)
khist = histogram(kpkx_0, binsize = bs, min = x_0[0], max=max(x_0), /nan)

phistx = lindgen(n_elements(phist))*bs + x_0[0]
khistx = lindgen(n_elements(khist))*bs + x_0[0]

!p.multi = [0,1,2]
plot, phistx, phist,ps=10, title='Pearson peak pixels'
plot, khistx, khist,ps=10, title='Kendall peak pixels'

ppkpos = where(phist ge 2, ppkcnt)
kpkpos = where(khist ge 2, kpkcnt)

if ppkcnt ne 0 then ppkx = phistx[ppkpos]
if kpkcnt ne 0 then kpkx = khistx[kpkpos]

stop

 
plot, x_0, pcoeff, /xs, ps=6, title='Pearson Regression Coefficients'
for i = 0, ppkcnt-1 do vline, ppkx[i]
plot, x_0, kcoeff, /xs, ps=6, title='Kendall Coefficients'
for i = 0, kpkcnt-1 do vline, kpkx[i]
stop
!p.multi = 0

pkx = [ppkx, kpkx]
pkxu = pkx[uniq(pkx, sort(pkx))]
npk = n_elements(pkxu)

bpkx = dblarr(npix_win, npk)
bpky = dblarr(npix_win, npk)
bpkp = dblarr(npk)
bpkk = dblarr(npk)

window, 1

for i = 0, npk - 1 do begin
    pidx = where(ppkx eq pkxu[i], pidxcnt)
    kidx = where(kpkx eq pkxu[i], kidxcnt)

    bfp = pkxu[i] - (npix_win/2)
    blp = bfp + npix_win - 1
    bpkx[*,i] = lindgen(npix_win) + bfp
    bpky[*,i] = bspec_med[bfp:blp]

    if pidxcnt gt 0 then bpkp[i] = pcoeff[bfp] else bpkp[i] = 0d0
    if kidxcnt gt 0 then bpkk[i] = kcoeff[bfp] else bpkk[i] = 0d0
    
    plot, bpkx[*,i], bpky[*,i], /xs, yr = [0.7,1.1], $
      title = "Chunk of median spectrum P = "+ $
      strtrim(bpkp[i],2)+ " | K = "+strtrim(bpkk[i],2)
    dxstop
endfor



stop

end
