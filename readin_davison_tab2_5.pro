; Quick pro to readin davison table 2.5 (Vsini from literature)
; WARNING: 2MASS ID's are missing final character in vdone_final.txt,
; but not in vdone_simbad_final.txt
pro readin_davison_tab2_5

; File path
outfile = '/home/stgilhool/APOGEE/vsini_literature/davison/tab2_5/dav_tab25.fits'
tabfile = '/home/stgilhool/APOGEE/vsini_literature/davison/tab2_5/vdone_final.txt'
simfile = '/home/stgilhool/APOGEE/vsini_literature/davison/tab2_5/vdone_simbad_final.txt'

; Read the Davison File
readcol, tabfile, $
  ID_2M, $
  ID_ALT1, $
  ID_ALT2, $
  RA, $
  DEC, $
  VMAG, $
  M_V, $
  K_S, $
  V_K, $
  PARALLAX, $
  E_PARALLAX, $
  VSINI, $
  REFERENCE, $
  UPPER_LIMIT, $
  FORMAT = "A,A,A,A,A,F,F,F,F,D,D,D,I,B", $
  SKIPLINE = 1, $
  DELIMITER = "|", $
  COUNT=NSTARS


; Here's a list of the references
bibcode =  ['1998A&A...331..581D', $
            '2015ApJ...801..106H', $
            '2009ApJ...704..975J', $
            '2012AJ....143...93R', $
            '2003ApJ...583..451M', $
            '2010ApJ...710..924R', $
            '2010AJ....139..504B', $
            '2014MNRAS.439.3094B', $
            'Davison Thesis', $
            '1986PASP...98.1233S', $
            '2015AJ....149..106D', $
            '2007A&A...467..259R', $
            '1998MNRAS.301.1031T', $
            '2013A&A...549A.109B', $
            '2007ApJ...656.1121R', $
            '2010A&A...514A..97L', $
            '2012ApJS..203...10T', $
            '2013AJ....146..154M']

refname = ['Delfosse et al. (1998)', $
           'Houdebine & Mullan (2015)', $
           'Jenkins et al. (2009)', $
           'Reiners et al. (2012)', $
           'Mohanty & Basri (2003)', $
           'Reiners & Basri (2010)', $
           'Browning et al. (2010)', $
           'Barnes et al. (2014)', $
           'Davison Thesis', $
           'Stauffer & Hartmann (1986)', $
           'Davison et al. (2015)', $
           'Reiners (2007)', $
           'Tinney & Reid (1998)', $
           'Bonfils et al. (2013a)', $
           'Reiners & Basri (2007)', $
           'L´opez-Santiago et al. (2010)', $
           'Tanner et al. (2012)', $
           'Mamajek et al. (2013)']
           
; Now read the file from SIMBAD output            
readcol, simfile, $
  SIM_IDX, $
  SIM_ID_2M, $
  SIM_ID_ALT1, $
  SIM_ID_ALT2, $
  SIM_RA, $
  SIM_DEC, $
  SIM_UMAG, $
  SIM_BMAG, $
  SIM_VMAG, $
  SIM_RMAG, $
  SIM_IMAG, $
  SP_TYPE, $
  TYPE, $
  NUM_BIB, $
  NOTE, $
  FORMAT = "I,A,A,A,A,A,F,F,F,F,F,A,A,I,I", $
  SKIPLINE = 1, $
  DELIMITER = "|", $
  /NAN, $
  COUNT=NSTARS_SIM




; Loop through and make structure
outstruct = []

for i = 0, nstars-1 do begin

    ; Check that ID's match
    ; print, id_2m[i], ' ', sim_id_2m[i]
;     print, id_alt1[i], ' ', sim_id_alt1[i]
;     print, id_alt2[i], ' ', sim_id_alt2[i]
;    if (sim_id_2m[i] ne id_2m[i]) or $
;    (sim_id_alt1[i] ne id_alt1[i]) or $
;    (sim_id_alt2[i] ne id_alt2[i]) then $
;      message, "ID mismatch"

    ; Convert RA and DEC from sexigesimal to decimal
    ra_i = ten(ra[i])*360d0/24d0
    dec_i = ten(dec[i])

    sim_ra_i = ten(sim_ra[i])*360d0/24d0
    sim_dec_i = ten(sim_dec[i])

    ; Look up what the reference is
    refnum = reference[i]-1

    struct_i = {ID_2M:SIM_ID_2M[i], $
                ID_ALT1:ID_ALT1[i], $
                ID_ALT2:ID_ALT2[i], $
                RA:ra_i, $
                DEC:dec_i, $
                RA_SBD:sim_ra_i, $
                DEC_SBD:sim_dec_i, $
                UMAG_SBD:SIM_UMAG[i], $
                BMAG_SBD:SIM_BMAG[i], $
                VMAG_SBD:SIM_VMAG[i], $
                RMAG_SBD:SIM_RMAG[i], $
                IMAG_SBD:SIM_IMAG[i], $
                VMAG:VMAG[i], $
                M_V:M_V[i], $
                K_S:K_S[i], $
                V_K:V_K[i], $
                PARALLAX:PARALLAX[i], $
                E_PARALLAX:E_PARALLAX[i], $
                SP_TYPE:SP_TYPE[i], $
                VSINI:VSINI[i], $
                E_VSINI:!values.f_nan, $
                REF_NUM_DAV:refnum, $
                REF_NAME:refname[refnum], $
                BIBCODE:bibcode[refnum], $
                UPPER_LIMIT:UPPER_LIMIT[i]}

    outstruct = [outstruct, struct_i]

endfor

mwrfits, outstruct, outfile, /create

stop

end
