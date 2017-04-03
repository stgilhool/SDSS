; Convert pixel position to wl (in Angstroms) 
; in APOGEE data
function apg_pix2wl, pix, log_pos=log_pos

; Log_pos, if set, will return the position of the input_wl
; on the logarithmic grid, which is technically correct.
; However, many of my functions and procedures linearly interpolate
; between the points of the wl grid.  It's a small error, so it's
; usually okay.  Leaving log_pos unset will return the correct
; pixel position for oversampled grids where I've interpolated
; linearly between points of the wl_grid
; In fact, probably due to floating point precision, there doesn't
; seem to be any difference in the output of this function
; (though the output of the apg_wl2pix function does depend of the
; setting of the log_pos keyword)
if n_elements(log_pos) eq 0 then log_pos = 0

; Constants
npix = 8575L

; Check input type
in_type = size(pix, /type)

case in_type of
    2: pixtype = 'int'
    3: pixtype = 'int'
    4: pixtype = 'float'
    5: pixtype = 'float'
    else: message, "Input must be integer or float"
endcase
        
; Check that pixel values are appropriate
if pix ge npix then message, "Pixel value too high"
if pix lt 0 then message, "Pixel value too low"

; Get wl info
pix_vec = dindgen(npix)

wl0 = 4.179d0
wl1 = 6d-6
log_wl = poly(pix_vec, [wl0,wl1])
wl_grid = 10d0^log_wl

; Find the wl at pixel position pix
pix = double(pix)
if log_pos then begin
    wlout_log = interpol(log_wl, pix_vec, pix)
    wlout = 10d0^wlout_log
endif else begin
    wlout = interpol(wl_grid, pix_vec, pix)
endelse

return, wlout

end
