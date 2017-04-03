; A tiny script to read in SDSS lsf coefficients from the file:
; '/home/stgilhool/APOGEE/APOGEE_data/lsf_coeffs.fits'
; It's a 28 coefficient x 1350 spectrum array, and the 28 lsf
; coefficients are an input to the function lsf_gh
; syntax - lsf = lsf_gh(xvector, xcentroid of lsf, lsf_coeffs)
function read_sdss_lsf_coeffs

lsf_file = '/home/stgilhool/APOGEE/APOGEE_data/lsf_coeffs.fits'

lsf_coeff_arr = mrdfits(lsf_file,0)

return, lsf_coeff_arr

end
