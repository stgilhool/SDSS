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

pro apogee_start, n=n

if n_elements(n) eq 0 then n=50

bigc = 299792.458 ;km/s

root_dir = '/home/stgilhool/APOGEE/'
data_path = root_dir + 'APOGEE_data/'

listfile = data_path + 'mdwarflist.txt'

; Readin file names
readcol, listfile, datafile, format="A", count=nfiles

files = data_path + datafile

; Read in headers

header = ptrarr(n, /allocate_heap)
spectra = dblarr(8575L,n)
nx = 8575L
x = dindgen(nx)

vharr = dblarr(n)
wl0arr = dblarr(n)
wl1arr = dblarr(n)

;logwl_grid = dblarr(nx, n)
;log_shift = dblarr(n)

for i = 0, n-1 do begin
    ; Read in header for a file
    null = mrdfits(files[i],0,head)
    ; save header
    *(header[i]) = head
    
    ; Read in 1d spectrum array for a file
    specall = mrdfits(files[i],1)

    ; Save only the first (combined) spectrum
    spec = specall[*,0]
    ; normalize it
    snorm = continuum_fit(x, spec)
    spec = temporary(spec)/temporary(snorm)

    ; Get vhelio, and wl coeffs for a file
    ; (vhelio is unnecessary, and wl coeffs are same for all)
    vhelio = fxpar(head, 'VHELIO', datatype=0.0d)
    wl0    = fxpar(head, 'CRVAL1', datatype=0.0d)
    wl1    = fxpar(head, 'CDELT1', datatype=0.0d)
    
    
    
    ; Make wl array
    ;logwl = wl0 + (x*wl1)
    ;log_dshift = alog10(1d0 + vhelio/bigc)

    ;logwl_shifted = logwl + log_dshift
    ;wl_grid = 10d0^(logwl)
    ;print, minmax(wl_grid)

    ; Save quantities
    vharr[i]  = vhelio
    wl0arr[i] = wl0    
    wl1arr[i] = wl1
    ;logwl_grid[*,i] = logwl
    ;log_shift[i] = log_dshift
    spectra[*,i] = spec
endfor

xx = rebin(x, nx, n)
; Make logarithmic wl scale and logarithmic wl shift
wl0 = rebin(reform(wl0arr, 1, n), nx, n)
wl1 = rebin(reform(wl1arr, 1, n), nx, n)
vhelio = rebin(reform(vharr, 1, n), nx, n)
logdop = alog10(1d0 + vhelio/bigc)
logdopneg = alog10(1d0 - vhelio/bigc)

logwl_grid = wl0 + (xx * wl1)
wl_grid = 10d0^logwl_grid
wl_shifted = 10d0^(logwl_grid + logdop)
wl_shiftedneg = 10d0^(logwl_grid + logdopneg)

minwl = max(wl_shifted[0,*]) > max(wl_shiftedneg[0,*])
maxwl = min(wl_shifted[-1,*]) < min(wl_shiftedneg[-1,*])

newgrid = minwl + ((maxwl-minwl)/nx)*x

; Spline all onto new grid

fir = 4990
las = 5400
!p.multi=0
window, 0, xsize = 1300, ysize = 1000
!p.multi=[0,1,3]
for i = 0, n-1 do begin
    
    ;spec0 = interpol(spectra[*,i], wl_grid[*,i], wl_shifted[*,i])
    ;spec  = interpol(spec0, wl_shifted[*,i], newgrid)
    spec  = interpol(spectra[*,i], wl_shifted[*,i], newgrid)
    
    if i eq 0 then begin
        plot, newgrid[fir:las], spec[fir:las], /xs, yr=[0.5,1.1], ps=3
        specsave = dblarr(n_elements(spec[fir:las]), n)
        
    endif else begin
        oplot, newgrid[fir:las], spec[fir:las], color = 200, ps=3
        
    endelse
        specsave[*,i] = spec[fir:las]
    
    wait, 0.15
    
endfor

for i = 0, n-1 do begin


    ;specneg0 = interpol(spectra[*,i], wl_grid[*,i], wl_shiftedneg[*,i])
    ;specneg = interpol(specneg0, wl_shiftedneg[*,i], newgrid)
    specneg = interpol(spectra[*,i], wl_shiftedneg[*,i], newgrid)
    if i eq 0 then begin
        plot, newgrid[fir:las], specneg[fir:las],/xs, yr=[0.5,1.1], ps=3
        specsaveneg = dblarr(n_elements(specneg[fir:las]), n)
    endif else begin
        oplot, newgrid[fir:las], specneg[fir:las], color = 200, ps=3
    endelse
    specsaveneg[*,i] = specneg[fir:las]
    wait, 0.15
endfor

for i = 0, n-1 do begin


    ;specneg0 = interpol(spectra[*,i], wl_grid[*,i], wl_shiftedneg[*,i])
    specns = interpol(spectra[*,i], wl_grid[*,i], newgrid)
    if i eq 0 then begin
        plot, newgrid[fir:las], specns[fir:las],/xs, yr=[0.5,1.1], ps=3
        specsavens = dblarr(n_elements(specns[fir:las]), n)
    endif else begin
        oplot, newgrid[fir:las], specns[fir:las], color = 200, ps=3
    endelse
    specsavens[*,i] = specns[fir:las]
    wait, 0.15
endfor



meanflux = mean(specsave, dimension=2)
meanfluxneg = mean(specsaveneg, dimension=2)
meanfluxns = mean(specsavens, dimension=2)
fluxerr = stddev(specsave, dimension=2)
fluxerrneg = stddev(specsaveneg, dimension=2)
fluxerrns = stddev(specsavens, dimension=2)

print, total(fluxerr)
print, total(fluxerrneg)
print, total(fluxerrns)

plot, newgrid[fir:las], meanflux, /xs, yr=[0.4, 1.2]
errplot, newgrid[fir:las], meanflux-fluxerr, meanflux+fluxerr

plot, newgrid[fir:las], meanfluxneg, /xs, yr=[0.4, 1.2]
errplot, newgrid[fir:las], meanfluxneg-fluxerrneg, meanfluxneg+fluxerrneg

plot, newgrid[fir:las], meanfluxns, /xs, yr=[0.4, 1.2]
errplot, newgrid[fir:las], meanfluxns-fluxerrns, meanfluxns+fluxerrns


!p.multi=0
stop

end
