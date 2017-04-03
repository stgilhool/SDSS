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




pro apo_fehv1

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
    !p.multi=[0,1,2]
endif
    
; Fit the metalicity of each chunk
; Initialize
coeff = dblarr(2,nwin)

; Make median bspec for display
bspec_med = median(bspec, dimension=2)


for i = 0, nwin-1 do begin
    
    ew  = reform(ew_arr[i,*])

    functargs = {ew:ew, $
                 feh:feh_afk, $
                 x_0:x_0[i], $
                 npix_win:npix_win, $
                 bspec_med:bspec_med, $
                 vis:vis $
                 }
    
    parinfo = replicate({value:1d0},2)

    
    r = mpfit('feh_fit', parinfo=parinfo, functargs=functargs, status=status)

    if n_elements(r) eq 1 then begin
        coeff[*,i] = !values.d_nan
    endif else begin
        coeff[*,i] = r
    endelse

    
endfor

window, 1, xsize=1300, ysize=1000
plot, x_0, coeff[1,*], /xs, ps=6, title='Linear Coeff'

window, 2
plot, x_0, coeff[0,*], /xs, ps=6, title='Const Coeff'

stop


end
