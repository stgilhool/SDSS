; Quick routine to turn vsini_newton_vsinionly.txt into a fits
; file
pro vsini_lit_make_newton_fits

; Read in text file
nwtfil = '/home/stgilhool/APOGEE/vsini_literature/vsini_newton_vsinionly.txt'
readcol, nwtfil, id, vsini, e_vsini, upper_limit, src, format="A,F,F,I,A"

nlines = n_elements(vsini) 


outstr = []


for i = 0, nlines-1 do begin

    

    outstruct = {ID:id[i], $
                 VSINI:vsini[i], $
                 E_VSINI:e_vsini[i], $
                 UPPER_LIMIT:upper_limit[i], $
                 SOURCE:src[i]}

    outstr = [outstr, outstruct]

endfor

mwrfits, outstr, '/home/stgilhool/APOGEE/vsini_literature/newton_vsini.fits', /create

stop

end
    



