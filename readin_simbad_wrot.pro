; Function to read in a parsed SIMBAD file after using
; simbad_parse_full.sed (which includes rot catalog)
function readin_simbad_wrot, filename

if n_params() eq 0 then message, "Must input filename"

if file_test(filename) eq 0 then message, "File does not exist"

; Now read the file from SIMBAD output            
readcol, filename, $
  IDX, $
  ID_ALT1, $
  ID_2M, $
  ID_ALT2, $
  TYPE, $
  RA, $
  DEC, $
  UMAG, $
  BMAG, $
  VMAG, $
  RMAG, $
  IMAG, $
  SP_TYPE, $
  NCITE, $
  NOTE, $
  NROT, $
  VSINI, $
  E_VSINI, $
  VSINI_MSG, $
  VSINI_Q, $
  VSINI_REF, $
  FORMAT = "I,A,A,A,A,A,A,F,F,F,F,F,A,I,I,I,F,F,A,A,A", $
  COMMENT = "#", $
  DELIMITER = "|", $
  /NAN, $
  COUNT=NSTARS


; Loop through, convert RA and DEC, and create structure
outstruct = []

for i = 0, nstars-1 do begin

    ; Convert RA and DEC from sexigesimal to decimal
    ra_i = ten(ra[i])*360d0/24d0
    dec_i = ten(dec[i])

    struct_i = {ID_2M:strtrim(ID_2M[i],2), $
                ID_ALT1:strtrim(ID_ALT1[i],2), $
                ID_ALT2:strtrim(ID_ALT2[i],2), $
                RA:ra_i, $
                DEC:dec_i, $
                UMAG:UMAG[i], $
                BMAG:BMAG[i], $
                VMAG:VMAG[i], $
                RMAG:RMAG[i], $
                IMAG:IMAG[i], $
                SP_TYPE:strtrim(SP_TYPE[i],2), $
                TYPE:TYPE[i], $
                NCITE:NCITE[i], $
                NOTE:NOTE[i], $
                NROT:NROT[i], $
                VSINI:VSINI[i], $
                E_VSINI:E_VSINI[i], $
                VSINI_MSG:VSINI_MSG[i], $
                VSINI_Q:VSINI_Q[i], $
                VSINI_REF:VSINI_REF[i]}

    outstruct = [outstruct, struct_i]

endfor

return, outstruct

end
