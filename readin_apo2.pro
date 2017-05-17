;+
; NAME:
;
;
;
; PURPOSE:
;  To read APOGEE spectra into a structure for easy handling
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
;       Thu Mar 31 10:15:18 2016,
;       <stgilhool@iroquois.physics.upenn.edu>
; Changed continuum normalization to work on each chip independently      
;		
; Also made readin_apg, which should be used from now on, which saves
; the first 44 lines of the header as string arrays, and so doesn't
; need to use pointers.
; Actually, this one is better than readin_apg I think
;-

function readin_apo2, HDU=hdu, IRTF=irtf, NFILES=nfiles, FILES=files, FLUX=FLUX, DR12=DR12

if n_elements(nfiles) eq 0 and n_elements(files) eq 0 then begin
    nfiles=50
    filelist = lindgen(nfiles)
endif else if n_elements(nfiles) gt 0 and n_elements(files) gt 0 then begin
    print, "Ignoring NFILES, using FILES..."
    filelist = files
    nfiles = n_elements(filelist) 
endif else if n_elements(nfiles) gt 0 and n_elements(files) eq 0 then begin
    filelist = lindgen(nfiles)
endif else if n_elements(nfiles) eq 0 and n_elements(files) gt 0 then begin
    filelist = files
    nfiles = n_elements(filelist) 
endif


if n_elements(flux) eq 0 then flux = 0
if n_elements(irtf) eq 0 then irtf = 0
if n_elements(dr12) eq 0 then dr12 = 0

bigc = 299792.458 ;km/s
nx = 8575L

root_dir = '/home/stgilhool/APOGEE/'
data_path = root_dir + 'APOGEE_data/'

if irtf eq 1 then listfile = data_path + 'irtflist.txt' $
  else if dr12 then listfile = data_path + 'mdwarflist.txt' $
  else listfile = data_path + 'mdwarflist_dr13.txt'

; Readin file names
readcol, listfile, datafile, format="A", count=nfiles_tot

filenames = data_path + datafile

; Read in headers

header = ptrarr(nfiles, /allocate_heap)
output = ptrarr(nfiles, /allocate_heap)
spectra = dblarr(nx,nfiles)
errors = dblarr(nx,nfiles)

x = dindgen(nx)

;vharr = dblarr(nfiles)
;wl0arr = dblarr(nfiles)
;wl1arr = dblarr(nfiles)

;logwl_grid = dblarr(nx, nfiles)
;log_shift = dblarr(nfiles)

;for i = 0, nfiles-1 do begin
foreach fileidx, filelist, idx do begin
    ; Read in header for a file
    null = mrdfits(filenames[fileidx],0,head)
    ; save header
    *(header[idx]) = head
        
    ; Read in 1d spectrum array for a file
    specall = mrdfits(filenames[fileidx],1)
    errall = mrdfits(filenames[fileidx], 2)

    ; Save only the first (combined) spectrum
    spec = specall[*,0]
    err = errall[*,0]
    ; normalize it
    if flux eq 0 then begin
        ; Remove continuum
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

        ; Normalize and scale errors for each chip
        spec[bmin:bmax] = spec[bmin:bmax]/bnorm
        spec[gmin:gmax] = spec[gmin:gmax]/gnorm
        spec[rmin:rmax] = spec[rmin:rmax]/rnorm
        
        err[bmin:bmax] = err[bmin:bmax]/bnorm
        err[gmin:gmax] = err[gmin:gmax]/gnorm
        err[rmin:rmax] = err[rmin:rmax]/rnorm

        ; Zero the regions in between
        spec[0:bmin-1] = 0d0
        spec[bmax+1:gmin-1] = 0d0
        spec[gmax+1:rmin-1] = 0d0
        spec[rmax+1:-1] = 0d0

        err[0:bmin-1] = 0d0
        err[bmax+1:gmin-1] = 0d0
        err[gmax+1:rmin-1] = 0d0
        err[rmax+1:-1] = 0d0


    endif 

    ; Get vhelio, and wl coeffs for a file
    ; (vhelio is unnecessary, and wl coeffs are same for all)
    ;vhelio = fxpar(head, 'VHELIO', datatype=0.0d)
    wl0    = fxpar(head, 'CRVAL1', datatype=0.0d)
    wl1    = fxpar(head, 'CDELT1', datatype=0.0d)

    if n_elements(hdu) ne 0 then begin
        ; ADD check that hdu is valid
        *(output[idx]) = mrdfits(filenames[fileidx],hdu)
    endif
    
    ; Make wl array
    ;logwl = wl0 + (x*wl1)
    ;log_dshift = alog10(1d0 + vhelio/bigc)

    ;logwl_shifted = logwl + log_dshift
    ;wl_grid = 10d0^(logwl)
    ;print, minmax(wl_grid)

    ; Save quantities
    ;vharr[i]  = vhelio
    ;wl0arr[i] = wl0    
    ;wl1arr[i] = wl1
    ;logwl_grid[*,i] = logwl
    ;log_shift[i] = log_dshift
    spectra[*,idx] = spec
    errors[*,idx] = err
endforeach

; Make logarithmic wl scale and logarithmic wl shift
logwl_grid = wl0 + (x * wl1)
wl_grid = 10d0^logwl_grid

; Return requested info
outstr1 = {APG_IDX:0L, $
           HEADER:header[0], $
           SPECTRA:spectra[*,0], $
           ERRORS:errors[*,0], $
           LOGWL_GRID:logwl_grid, $
           WL_GRID:wl_grid, $
           OUTPUT:output[0] $
           }

outstr = replicate(outstr1, nfiles)

outstr.header = header
outstr.output = output
outstr.apg_idx = filelist

for s = 0, nfiles-1 do begin
    outstr[s].spectra = spectra[*,s]
    outstr[s].errors = errors[*,s]
    outstr[s].wl_grid = wl_grid
endfor
    
return, outstr

end
