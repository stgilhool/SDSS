; Quick routine to turn spec_results_dav_teff_vsini.txt into a fits
; file
pro vsini_lit_make_davison_fits

; Read in text file
davfil = '/home/stgilhool/APOGEE/vsini_literature/davison/spec_results_dav_teff_vsini.txt'
readcol, davfil, id1, id2, teff, vsini, e_vsini, format="A,A,L,F,F"

nlines = n_elements(vsini) 


outstr = []


for i = 0, nlines-1 do begin

    outstruct = {ID:strjoin([id1[i],id2[i]], ' '), $
                 TEFF:teff[i], $
                 VSINI:vsini[i], $
                 E_VSINI:e_vsini[i]}

    outstr = [outstr, outstruct]

endfor

mwrfits, outstr, '/home/stgilhool/APOGEE/vsini_literature/davison/davison_vsini.fits', /create

stop

end
    



