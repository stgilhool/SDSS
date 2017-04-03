; readin apg headers as a structure
function get_apg_head

; file
hfile = '/home/stgilhool/APOGEE/apg_head_info.fits'

; read
out = mrdfits(hfile, 1)

return, out

end
