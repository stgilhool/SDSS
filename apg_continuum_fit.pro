; Do continuum fit on 3-chip APG spectra
; WARNING: other_vec is an input and output variable (it is altered
; and returned)
function apg_continuum_fit, spectrum, brange, grange, rrange, OTHER_VEC=other_vec

nx = n_elements(spectrum) 
if nx ne 8575L then message, "This function takes full spectrum (8575 pixels)"

; Define chip pixels
bmin = brange[0]
bmax = brange[1]
gmin = grange[0]
gmax = grange[1]
rmin = rrange[0]
rmax = rrange[1]

nblue  = bmax - bmin+1
ngreen = gmax - gmin+1
nred   = rmax - rmin+1

; Get continuum for each chip    
bnorm = continuum_fit(dindgen(nblue), spectrum[bmin:bmax])
gnorm = continuum_fit(dindgen(ngreen), spectrum[gmin:gmax])
rnorm = continuum_fit(dindgen(nred), spectrum[rmin:rmax])

; Normalize and scale errors for each chip
spec = spectrum
spec[bmin:bmax] = spec[bmin:bmax]/bnorm
spec[gmin:gmax] = spec[gmin:gmax]/gnorm
spec[rmin:rmax] = spec[rmin:rmax]/rnorm

; Zero the regions in between
spec[0:bmin-1] = 0d0
spec[bmax+1:gmin-1] = 0d0
spec[gmax+1:rmin-1] = 0d0
spec[rmax+1:-1] = 0d0
    
; Apply to other_vec if necessary
if n_elements(other_vec) eq nx then begin
    ; Normalize other_vec
    other_vec[bmin:bmax] = other_vec[bmin:bmax]/bnorm
    other_vec[gmin:gmax] = other_vec[gmin:gmax]/gnorm
    other_vec[rmin:rmax] = other_vec[rmin:rmax]/rnorm
    
    ; Zero the regions in between
    other_vec[0:bmin-1] = 0d0
    other_vec[bmax+1:gmin-1] = 0d0
    other_vec[gmax+1:rmin-1] = 0d0
    other_vec[rmax+1:-1] = 0d0
endif

return, spec

end
