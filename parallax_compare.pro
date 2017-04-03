; Match the APOGEE files with parallaxes from a big table in the
; literature
pro parallax_compare, plx_match_idx, apg_match_idx

; DESCRIPTION FROM DATA TABLE
; ================================================================================
; Byte-by-byte Description of file: datafile.txt
; --------------------------------------------------------------------------------
;    Bytes Format Units     Label          Explanations
; --------------------------------------------------------------------------------
;    1- 11 A11    ---       LSPM           LSPM North Catalog Designation
;   13- 20 A8     ---       GL             Gliese Catalog Number
;   22- 26 A5     ---       LHS            LHS Catalog Number
;   28- 32 I5     ---       NLTT           NLTT Catalog Number
;   34- 35 I2     ---       RAh            Right Ascension (hours)
;   37- 38 I2     ---       RAm            Right Ascension (minutes)
;   40- 44 F5.2   ---       RAs            Right Ascension (seconds)
;   46- 47 I2     deg       DEd            Declination (degrees)
;   49- 50 I2     ---       DEm            Declination (minutes)
;   52- 56 F5.2   ---       DEs            Declination (seconds)
;   58- 64 F7.4   ---       mu(RA)         Proper motion in RA * cos(Declination) in units of arcseconds / year
;   66- 72 F7.4   ---       mu(DEC)        Proper motion in declination in units of arcseconds / year
;   74- 82 F9.4   d         JD(start)      Julian epoch of the first data point in our fit
;   84- 92 F9.4   d         JD(end)        Julian epoch of the last data point in our fit
;   94- 99 F6.2   marcsec   plx(rel)       Relative parallax (milliarcseconds)
;  101-104 F4.2   marcsec   plx(abscorr)   Absolute parallax correction (milliarcseconds)
;  106-111 F6.2   marcsec   plx(abs)       Absolute parallax (milliarcseconds)
;  113-117 F5.2   marcsec e_plx(abs)       Absolute parallax error (milliarcseconds)
;  119-124 F6.3   mag       Jmag           2MASS J magnitude
;  126-131 F6.3   mag       Kmag           2MASS K magnitude
;  133-137 F5.3   solMass   Mass           Estimated mass (solar masses)
;  139-143 F5.3   solRad    Radius         Estimated radius (solar radii)
;  145-149 F5.3   solRad  e_Radius         Error in radius (solar radii)
;  151-155 I5     ---     o_plx(rel)       Number of data points in our fit
;  157-158 I2     ---       Num(refstars)  Number of reference stars used in our fit
;  160-164 A5     ---       tel            MEarth Telescope number the data was taken on


; Paths
par_path = '/home/stgilhool/APOGEE/parallax_literature/'
;par_filename = 'MEarth_Parallaxes.txt'
par_filename = 'parallax_data.fits'
par_file = par_path+par_filename
id2m_file = par_path + 'lspmto2mass.txt'

;;; Constants
napgspec = 1350

;;; Read in APOGEE headers
apg_id = strarr(napgspec)
hdrfile = '/home/stgilhool/APOGEE/tempheader.fits'
if file_test(hdrfile) then begin
    
    apg_info = mrdfits(hdrfile,1)
    apg_id_temp = apg_info.objid
    apg_id = strtrim(apg_id_temp,2)
    
endif else begin
    apg_str = readin_apo(hdu = 0, nfiles = napgspec)
    
    ; Get objids for APG files
    for i = 0, napgspec-1 do begin
        head = *(apg_str[i].header)
        apg_id_temp = fxpar(head, 'OBJID')
        apg_id[i] = strtrim(apg_id_temp,2)
    endfor
endelse


;;; Read in plx table
plx_str = mrdfits(par_file, 1)

; Make text file of names to input into SIMBAD
;lspm_name = plx_str.lspm
;lspm_name = 'LSPM'+lspm_name

;openw, lun, '/home/stgilhool/APOGEE/parallax_literature/lspm_namelist.txt', /get_lun
;foreach name, lspm_name do begin
;    printf, lun, name
;endforeach

;free_lun, lun
;stop

;;; Read results of SIMBAD conversion to 2MASS id's
readcol, id2m_file, fnum, lspm, ids, delimiter = '|', format = 'L,A,A', count=nplx
plx_id = strarr(nplx)

plx_match_idx = []
apg_match_idx = []

; Loop through each plx file, manipulate id name, and
; look for a match in the APG files
for i = 0, nplx-1 do begin
    ; Change id from '2MASS Jxxxxxxxx+xxxxxxx' to '2Mxxxxxxxx+xxxxxxx'
    ; to match APG files
    plx_id_temp = mg_streplace(ids[i], '^2MASS J', '2M')
    ; Extract only the 2MASS id (get rid of alternate ids)
    plx_id_split = strsplit(plx_id_temp, /extract)
    plx_id_2m = plx_id_split[0]
    plx_id_2m = strtrim(plx_id_2m,2)
    plx_id[i] = plx_id_2m
    ; Compare with APG files
    apg_match_i = where(apg_id eq plx_id_2m, matchcnt)
    if matchcnt gt 1 then message, "More than 1 match!" $
      else if matchcnt eq 1 then begin
        plx_match_idx = [plx_match_idx, i]
        apg_match_idx = [apg_match_idx, apg_match_i]
    endif
        
endfor

help, plx_match_idx
help, apg_match_idx
;stop
end


