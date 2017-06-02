; Quick pro to readin davison table 11.1 (Projected Rotational
; Velocity Measurements from CSHELL)
pro readin_davison_tab11_1

; File path
outfile = '/home/stgilhool/APOGEE/vsini_literature/davison/tab11_1/dav_tab111.fits'
tabfile = '/home/stgilhool/APOGEE/vsini_literature/davison/tab11_1/table11_1.txt'
simfile = '/home/stgilhool/APOGEE/vsini_literature/davison/tab11_1/cshell_simbad.fits'

; Read the Davison File
readcol, tabfile, $
  ID_ALT1, $
  RA, $
  DEC, $
  VMAG, $
  M_V, $
  K_S, $
  V_K, $
  PARALLAX, $
  E_PARALLAX, $
  VSINI, $
  E_VSINI, $
  UPPER_LIMIT, $
  FORMAT = "A,A,A,F,F,F,F,F,F,F,F,B", $
  COMMENT = "#", $
  DELIMITER = "|", $
  COUNT=NSTARS

           
; Now read the file from SIMBAD output            
s = mrdfits(simfile,1)

; Loop through and make structure
outstruct = []

for i = 0, nstars-1 do begin

    ; Check that ID's match
    print, id_alt1[i], ' ', s[i].id_alt1
;     print, id_alt1[i], ' ', sim_id_alt1[i]
;     print, id_alt2[i], ' ', sim_id_alt2[i]
;    if (sim_id_2m[i] ne id_2m[i]) or $
;    (sim_id_alt1[i] ne id_alt1[i]) or $
;    (sim_id_alt2[i] ne id_alt2[i]) then $
;      message, "ID mismatch"
    dxstop
    ; Convert RA and DEC from sexigesimal to decimal
    ra_i = ten(ra[i])*360d0/24d0
    dec_i = ten(dec[i])

    sim_ra_i = s[i].ra
    sim_dec_i = s[i].dec

    ; Make structure
    struct_i = {ID_2M:s[i].id_2m, $
                ID_ALT1:ID_ALT1[i], $
                ID_ALT2:s[i].id_alt2, $
                RA:ra_i, $
                DEC:dec_i, $
                RA_SBD:sim_ra_i, $
                DEC_SBD:sim_dec_i, $
                UMAG_SBD:s[i].umag, $
                BMAG_SBD:s[i].BMAG, $
                VMAG_SBD:s[i].vmag, $
                RMAG_SBD:s[i].rmag, $
                IMAG_SBD:s[i].imag, $
                VMAG:VMAG[i], $
                M_V:M_V[i], $
                K_S:K_S[i], $
                V_K:V_K[i], $
                PARALLAX:PARALLAX[i], $
                E_PARALLAX:E_PARALLAX[i], $
                SP_TYPE:s[i].SP_TYPE, $
                VSINI:VSINI[i], $
                E_VSINI:E_VSINI[i], $
                REF_NUM_DAV:9, $
                REF_NAME:"Davison Thesis", $
                BIBCODE:"Davison Thesis", $
                UPPER_LIMIT:UPPER_LIMIT[i]}

    outstruct = [outstruct, struct_i]

endfor

mwrfits, outstruct, outfile, /create

stop

end
