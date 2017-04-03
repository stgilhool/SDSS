; Convert wl (in Angstroms) to pixel positions
; in APOGEE data
function apg_wl2pix, input_wl, round=round, log_pos=log_pos

if n_elements(round) eq 0 then round = 0 
; Log_pos, if set, will return the position of the input_wl
; on the logarithmic grid, which is technically correct.
; However, many of my functions and procedures linearly interpolate
; between the points of the wl grid.  It's a small error, so it's
; usually okay.  Leaving log_pos unset will return the correct
; pixel position for oversampled grids where I've interpolated
; linearly between points of the wl_grid
if n_elements(log_pos) eq 0 then log_pos = 0

; Get wl info
npix = 8575L
pix_vec = dindgen(npix)

wl0 = 4.179d0
wl1 = 6d-6
log_wl = poly(pix_vec, [wl0,wl1])
wl_grid = 10d0^log_wl

; Take log of input wl
log_input_wl = alog10(input_wl)


; Find where wl fits into wl_grid
; xpos is the pixel of the lower bound
xpos = value_locate(wl_grid, input_wl)
xpos_log = value_locate(log_wl, log_input_wl)

; Check that input is actually on the wl_grid
if xpos ge npix then message, "Input wl too high"
if xpos eq -1 then message, "Input wl too low"

; array of pixel values on either side of wl
pix_bound = [xpos, xpos+1]

; interpolate to get exact position
if log_pos then begin
    xpos_interp = interpol(pix_bound, log_wl[pix_bound], log_input_wl)
endif else begin
    xpos_interp = interpol(pix_bound, wl_grid[pix_bound], input_wl)
endelse

; if round is set, return nearest integer pixel value
xpos_interp_float = xpos_interp
if round then xpos_interp = round(xpos_interp_float)

return, xpos_interp

end
