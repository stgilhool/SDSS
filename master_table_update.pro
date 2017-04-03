; Take rfile_master.fits from vsini_fit, and add JHK & V mags,
; add proper motion, add avg RV, and parallax if available

pro master_table_update
;goto, startup
infile = '/home/stgilhool/APOGEE/vsini_results/initial_run/rfile_master.fits'

; Readin master file
instr = mrdfits(infile, 1)

napg = n_elements(instr) 

; Read in my FEH results
;feh_rfile = '/home/stgilhool/APOGEE/feh_results/feh_results_100.fits'
feh_rfile = '/home/stgilhool/APOGEE/feh_results/feh_results_redo.fits'
feh_steve = mrdfits(feh_rfile,0)


; Make new output structure
newstr = create_struct(['OBJID','J','H','K','V_SBD','V_LSPM','V_EST_LSPM','RA','DEC','rPM',$
                        'rPMra','rPMde','PM','PMra','PMde','RV_avg',$
                        'PLX_SBD','PLX_PTM','FEH_SG','FEH_IRTF'], $
                       '', $                       !
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan, $
                       !values.d_nan)

addstr = replicate(newstr, napg)

combstr = create_struct(instr[0], addstr[0])
outstr = replicate(combstr, napg)
; Readin APG headers
apgstr = readin_apo(hdu=0, nfiles = napg)

save, /variables, filename = 'master_tab_up_startup.sav'

startup:
restore, 'master_tab_up_startup.sav'

; Loop through and add to each record
for i = 0, napg-1 do begin
    
    
    ; Get APG header
    head = *(apgstr[i].header)
    
    ; Get JHK, RA and DEC
    objid = fxpar(head, 'OBJID')
    jmag = fxpar(head, 'J')
    hmag = fxpar(head, 'H')
    kmag = fxpar(head, 'K')
    ra = fxpar(head, 'RA')
    dec = fxpar(head, 'DEC')

    print, "Working on file: " + strtrim(i,2)
    print, "ID: " + objid

    ; Insert into add str
    addstr[i].OBJID = objid
    addstr[i].J = jmag
    addstr[i].H = hmag
    addstr[i].K = kmag
    addstr[i].RA = ra
    addstr[i].DEC = dec

    ; Add feh_sg to addstr
    addstr[i].feh_sg = feh_steve[i]

    ; Query SIMBAD for parallax and vmag, if available
    id_2mass = mg_streplace(objid, '2M', '2MASS J')
    querysimbad, id_2mass, ra_sbd, dec_sbd, id_sbd, found=found, /cfa, $
      parallax = plx, vmag = vmag
    if found then begin
        if n_elements(parallax) gt 0 then addstr[i].PLX_SBD = plx
        if n_elements(vmag) gt 0 then addstr[i].V_SBD = vmag
    endif

    ; Query VizieR for LSPM proper motions
    lspm = queryvizier('I/298/lspm_n', [ra,dec], 5, /allcolumns)
    ; Make sure we have the right match
    lspm_type = size(lspm, /type)
    if lspm_type ne 3 then begin
        ;if n_elements(lspm) gt 1 then begin
            ; Get just the right one
        id_2mass_lspm = mg_streplace(id_2mass, '2MASS J', '')
        lspm_match_idx = where(strtrim(lspm._2mass,2) eq $
                               strtrim(id_2mass_lspm,2), lspmcnt)
        if lspmcnt eq 1 then begin
            lspm = lspm[lspm_match_idx]

            addstr[i].rpm = lspm.rpm 
            addstr[i].rpmra = lspm.rpmra
            addstr[i].rpmde = lspm.rpmde
            addstr[i].pm = lspm.pm
            addstr[i].pmra = lspm.pmra
            addstr[i].pmde = lspm.pmde
            addstr[i].V_LSPM = lspm.vmag
            addstr[i].V_EST_LSPM = lspm.vemag
        endif else if lspmcnt gt 1 then $
          message, "More than 1 2MASS id match in LSPM search"
    endif

        


    ; Combine
    combstr = create_struct(instr[i], addstr[i])
    
    outstr[i] = combstr

endfor



mwrfits, outstr, '/home/stgilhool/APOGEE/master_table2.fits', /create

stop
end

    
