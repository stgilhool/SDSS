; One-off routine to query SIMBAD to get RA and DEC for davison Vsini's

pro vsini_lit_add_radec_davison

; Read in Davison fits file
davfil = '/home/stgilhool/APOGEE/vsini_literature/davison/davison_vsini.fits'
outfile = '/home/stgilhool/APOGEE/vsini_literature/davison/davison_vsini2.fits'
dav = mrdfits(davfil,1)

nstars = n_elements(dav) 

; Initialize outputs
ra = []
dec = []
outstruct = []

; Loop through, get ra and dec, and add to structure
for i = 0, nstars-1 do begin

    ;Query Simbad to get ra and dec
    id = strtrim(dav[i].id,2)

    querysimbad, id, ra_i, dec_i, found=found, /cfa, /verbose

    if found eq 0 then message, "Star number "+strtrim(i,2)+" : "+id+" not found"

    ; If found, save stuff
    ra = [ra, ra_i]
    dec = [dec, dec_i]
    
    ; Update structure
    dstruct = dav[i]
    ostruct = create_struct(dstruct, 'RAJ2000', ra_i, 'DECJ2000', dec_i)

    outstruct = [outstruct, ostruct]

endfor

help, ra
help, dec

help, ostruct

stop
; Write the updated structure to disk
mwrfits, outstruct, outfile, /create

stop
end
