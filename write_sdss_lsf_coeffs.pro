; A one-time-use procedure to read in sdss lsf coefficients and then
; write them to a file, so that I don't have to read in all the files again
pro write_sdss_lsf_coeffs
napofiles = 1350
lsf_str = readin_apo(nfiles = napofiles, hdu = 8, /flux)

lsf_coeff_arr = dblarr(28, napofiles)

for i = 0, napofiles-1 do begin
    lsf_coeff_arr[0,i] = *(lsf_str[i].output)
endfor

outfile = '/home/stgilhool/APOGEE/APOGEE_data/lsf_coeffs.fits'

mwrfits, lsf_coeff_arr, outfile, /create

end
