function readin_spec
end

pro apo_regress, x_0, width


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


for i = 0, nfiles-1 do begin

    file_split = strsplit(files[i], '/', /extract)
    files_truncated[i] = file_split[-1]
endfor

; Readin irtf data
irtf = read_irtf_data()

n_irtf = n_elements(irtf)
irtf_files = strtrim(irtf.APG_file,2)
feh_afk = irtf.AFK
feh_afke= irtf.AFKE

    
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

badspecflag = intarr(n)

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

    print, "filenumber: ", i
    print, sdss_flagname('APOGEE_STARFLAG',fxpar(head,'ANDFLAG'))

    
    read, sto, prompt='0 to exclude, any key to store'
    if sto ne 0 then badspecflag[i] = 0 else badspecflag[i] = 1
    
    ;if i eq 0 then plot, spec, /xs else oplot, spec, color=200
    ;wait, 0.1

endfor


stop

end
