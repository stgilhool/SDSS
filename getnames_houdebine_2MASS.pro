pro getnames_houdebine_2MASS
; Match the 2MASS names from the Houdebine SIMBAD file
; to our APOGEE 2MASS names to determine if there are matches
houdpath = '/home/stgilhool/APOGEE/vsini_literature/houdebine/'
houdfile = houdpath + 'houdebine2MASS_readcol.txt'

apg_file = '/home/stgilhool/APOGEE/master_info5.fits'
apg_info = mrdfits(apg_file,1)
; Read in apg IDs
apg_2M = apg_info.objid

; Read in Houdebine IDs
readcol, houdfile, simnum, houd_2M, format = 'L,A'

nhoud = n_elements(houd_2M) 

simnum_vec = []
houd_2M_vec = []

; Loop through the Houdebine files and retrieve the APG info, if matched
for houdnum = 0, nhoud - 1 do begin
    houd_2M_i = strtrim(houd_2M[houdnum],2)
    houd_2M_fmt = mg_streplace(houd_2M_i, 'J', '2M')
    houd_idx = where(apg_2M eq houd_2M_fmt, match_chk)
    if match_chk eq 1 then begin
        ; Save match
        simnum_vec = [simnum_vec, simnum[houdnum]]
        houd_2M_vec = [houd_2M_vec, houd_2M_fmt]
        print, simnum[houdnum]
        print, houd_2M_fmt
        print, apg_info[houd_idx].vsini_lit
        print, apg_info[houd_idx].vsini_lit_src
        stop
    endif else if match_chk eq 0 then continue $
      else if match_chk gt 1 then message, "More than one match" 
endfor

stop
end
