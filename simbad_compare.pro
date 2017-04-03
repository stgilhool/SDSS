; Compare the 2MASS names in a SIMBAD list to APG 2MASS names and find matches
pro simbad_compare, sim_match_idx, apg_match_idx, $
                    SIMBAD_FILE=simbad_file, APG_ID=apg_id, APG_STR=apg_str

;;; Default simbad file
if n_elements(simbad_file) eq 0 $
  then simbad_file = '/home/stgilhool/APOGEE/vsini_literature/davison/davison_vsini_done_simbad.txt'

;;; If no APG names are provided, read them in
if n_elements(apg_id) eq 0 then begin
    napgspec=1350
    ;;; Read in APOGEE headers
    apg_str = readin_apo(hdu = 0, nfiles = napgspec)
    apg_id = strarr(napgspec)
    ; Get objids for APG files
    for i = 0, napgspec-1 do begin
        head = *(apg_str[i].header)
        apg_id_temp = fxpar(head, 'OBJID')
        apg_id[i] = strtrim(apg_id_temp,2)
    endfor
endif

;;; Read results of SIMBAD conversion to 2MASS id's
readcol, simbad_file, fnum, orig_id, ids, delimiter = '|', format = 'L,A,A', count=nfil
sim_id = strarr(nfil)

sim_match_idx = []
apg_match_idx = []

; Loop through each plx file, manipulate id name, and
; look for a match in the APG files
for i = 0, nfil-1 do begin
    ; Change id from '2MASS Jxxxxxxxx+xxxxxxx' to '2Mxxxxxxxx+xxxxxxx'
    ; to match APG files
    sim_id_temp = mg_streplace(ids[i], '^2MASS J', '2M')
    ; Extract only the 2MASS id (get rid of alternate ids)
    sim_id_split = strsplit(sim_id_temp, /extract)
    sim_id_2m = sim_id_split[0]
    sim_id_2m = strtrim(sim_id_2m,2)
    sim_id[i] = sim_id_2m
    ; Compare with APG files
    apg_match_i = where(apg_id eq sim_id_2m, matchcnt)
    if matchcnt gt 1 then message, "More than 1 match!" $
      else if matchcnt eq 1 then begin
        sim_match_idx = [sim_match_idx, i]
        apg_match_idx = [apg_match_idx, apg_match_i]
    endif
        
endfor

; Check that the matches work
apg_matches = apg_id[apg_match_idx]
sim_matches = sim_id[sim_match_idx]

if array_equal(apg_matches, sim_matches) ne 1 then print, "WARNING: Arrays not equal"

;out_str = {ID:sim_id[sim_match_idx], $
;           SIM_MATCH_IDX:sim_match_idx, $
;           APG_MATCH_IDX:apg_match_idx
stop
return

end
