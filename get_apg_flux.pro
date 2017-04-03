; Readin the pre-made fits file with APG fluxes

function get_apg_flux

file = '/home/stgilhool/APOGEE/APOGEE_data/apg_flux.fits'

str = mrdfits(file, 1)

return, str

end
