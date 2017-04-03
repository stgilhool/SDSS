; Redo the master table to include all the good stuff
pro nwt_vsini_update, outstr, nwt_idx

; Info for Newton's vsini's
nwt_id = ['J02170993+3526330', $
          'J02204625+0258375', $
          'J03205965+1854233', $
          'J03425325+2326495', $
          'J06000351+0242236', $
          'J07444018+0333089', $
          'J08294949+2646348', $
          'J08505062+5253462', $
          'J09301445+2630250', $
          'J09535523+2056460', $
          'J10163470+2751497', $
          'J10521423+0555098', $
          'J11474074+0015201', $
          'J12185939+1107338', $
          'J13003350+0541081', $
          'J13564148+4342587', $
          'J17195422+2630030', $
          'J18024624+3731048', $
          'J18452147+0711584', $
          'J19173151+2833147', $
          'J22245593+5200190', $
          'J03132299+0446293', $
          'J05011802+2237015', $
          'J06022918+4951561', $
          'J09002359+2150054', $
          'J11005043+1204108', $
          'J12265737+2700536', $
          'J16370146+3535456', $
          'J22081254+1036420', $
          'J23134727+2117294', $
          'J23354132+0611205' $
          ]
nwt_vsini = [ 28.2, $
              23.3, $
              8.0, $
              12.7, $
              5.8, $
              4.5, $
              8.1, $
              13.1, $
              6.7, $
              16.5, $
              4.0, $
              19.1, $
              5.6, $
              9.2, $
              16.8, $
              14.0, $
              1.9, $
              4.5, $
              16.1, $
              13.2, $
              4.5, $
              2.2, $
              8.8, $
              4.3, $
              20.0, $
              26.5, $
              13.5, $
              7.0, $
              18.6, $
              16.0, $
              9.8 $
              ]

nwt_upperlimit = [0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  1, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  1, $
                  0, $
                  0, $
                  1, $
                  1, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0, $
                  0 $
                  ]

nwt_vsini_err = [0.7, $
                 0.7, $
                 !values.f_nan, $
                 0.5, $
                 0.3, $
                 !values.f_nan, $
                 1.1, $
                 0.7, $
                 1.5, $
                 0.4, $
                 !values.f_nan, $
                 0.2, $
                 1.4, $
                 1.9, $
                 2.1, $
                 2.0, $
                 0.3, $
                 !values.f_nan, $
                 0.1, $
                 0.5, $
                 !values.f_nan, $
                 !values.f_nan, $
                 0.3, $
                 1.2, $
                 0.6, $
                 0.8, $
                 0.6, $
                 1.8, $
                 2.0, $
                 4.0, $
                 1.1 $
                 ]
                 
nnwt = n_elements(nwt_id) 

nwt_idx = []
for i = 0, nnwt-1 do begin
    ; Replace 'J' with '2M'
    nwt_2m_id = mg_streplace(nwt_id[i], 'J', '2M')
    nwt_2m_id = strtrim(nwt_2m_id,2)
    ; Match with APG id's
    nwt_match_idx = where(outstr.objid eq nwt_2m_id, nwtcnt)
    ; If match, store in outstr
    if nwtcnt eq 1 then begin
        nwt_idx=[nwt_idx,nwt_match_idx]
        outstr[nwt_match_idx].vsini_lit = nwt_vsini[i]
        outstr[nwt_match_idx].vsini_lit_src = 'NWT'
        ; Check if upper limit or not
        if nwt_upperlimit[i] eq 1 then begin
            outstr[nwt_match_idx].vsini_lit_hisig = 0d0
            outstr[nwt_match_idx].vsini_lit_losig = nwt_vsini[i]
        endif else begin
            outstr[nwt_match_idx].vsini_lit_hisig = nwt_vsini_err[i]
            outstr[nwt_match_idx].vsini_lit_losig = nwt_vsini_err[i]
        endelse
    endif else if nwtcnt gt 1 then message, "More than one NWT match"
endfor

return

end


pro rnr_vsini_update, outstr, rnr_idx

; Reiner's info
rnr_path = '/home/stgilhool/APOGEE/vsini_literature/reiners/'
rnr_file = rnr_path + 'vsini_reiners.fits'

rnr_str = mrdfits(rnr_file,1)
rnr_ra = rnr_str.raj2000
rnr_de = rnr_str.dej2000
rnr_ll = rnr_str.l_vsini
rnr_vsini = rnr_str.vsini

; APG info
apg_ra_arr = outstr.ra
apg_de_arr = outstr.dec

;See if there's a match in Reiner's data
srcor, rnr_ra, rnr_de, apg_ra_arr, apg_de_arr, 4000d0, out1, out2, /option, /spherical

; Get matches
r_vsini = rnr_str[out1].vsini

; Deal with errors and upper limits
r_hisig = replicate(3d0,n_elements(out1)) 
r_losig = replicate(3d0,n_elements(out1)) 
;upper limits
r_ll_idx = where(strtrim(rnr_str[out1].l_vsini,2) eq '<', rllcnt, complement=r_det_idx)
if rllcnt gt 0 then begin
    r_hisig[r_ll_idx] = 0.0
    r_losig[r_ll_idx] = r_vsini[r_ll_idx]
endif

; Store results
outstr[out2].vsini_lit = r_vsini
outstr[out2].vsini_lit_src = 'RNR'
outstr[out2].vsini_lit_hisig = r_hisig
outstr[out2].vsini_lit_losig = r_losig

rnr_idx = out2
return

end

pro dav_vsini_update, outstr, dav_idx

;;; These are the Davison matches
dav_id = ['LHS 1785', 'G108-021A', 'GJ 1105', 'GJ1182', 'LHS 3075', $
          'GJ 333.2B', 'GJ 363', 'GJ403']
dav_2M = ['2M05470907-0512106', $
          '2M06421118+0334527', $
          '2M07581269+4118134', $
          '2M14153253+0439312', $
          '2M15294392+4252498', $
          '2M09005033+0514293', $
          '2M09422327+5559015', $
          '2M10520440+1359509' $
          ]

dav_vsini = [4.5, $
             0.9, $
             2.0, $
             6.8, $
             2.5, $
             3.0, $
             3.0, $
             3.0 $
             ]

dav_vsini_hisig = [0.6, $
                   0.16, $
                   0.0, $
                   0.5, $
                   0.0, $
                   0.0, $
                   0.0, $
                   0.0 $
                  ]

dav_vsini_losig = [0.6, $
                   0.22, $
                   dav_vsini[2], $
                   0.5, $
                   dav_vsini[4], $
                   dav_vsini[5], $
                   dav_vsini[6], $
                   dav_vsini[7] $
                  ]

dav_idx=[]
; Loop through the Davison matches and retrieve the APG info
for i = 0, n_elements(dav_2M)-1 do begin
    dav_match_idx = where(outstr.objid eq dav_2M[i], davcnt)
    if davcnt ne 1 then message, "Perhaps there's a typo in your dav_2M"
    dav_idx = [dav_idx, dav_match_idx]
    ; Save result
    outstr[dav_match_idx].vsini_lit = dav_vsini[i]
    outstr[dav_match_idx].vsini_lit_src = 'DAV'
    outstr[dav_match_idx].vsini_lit_losig = dav_vsini_losig[i]
    outstr[dav_match_idx].vsini_lit_hisig = dav_vsini_hisig[i]
    
endfor

return

end

pro irtf_feh_update, outstr, irtf_idx



return
end

pro dtm_plx_update, outstr

;;; Add Dittman parallaxes (PLX_DTM)
par_path = '/home/stgilhool/APOGEE/parallax_literature/'
par_filename = 'parallax_data.fits'
par_file = par_path+par_filename
plx_str = mrdfits(par_file, 1)

; Get the matching indices between Dittman data set and APG
parallax_compare, plx_match_idx, apg_match_idx

; Store the correct entries
outstr[apg_match_idx].plx_dtm = plx_str.plx_abs[plx_match_idx]
outstr[apg_match_idx].eplx_dtm = plx_str.e_plx_abs[plx_match_idx]

return
end

pro master_table_final

; Read in master

masterfil = '/home/stgilhool/APOGEE/master_table2.fits'
masterstr = mrdfits(masterfil, 1)
napg = n_elements(masterstr) 

; Create Final structure
tags = ['OBJID', $
        'RA', $
        'DEC', $
        'RPM', $
        'RPMRA', $
        'RPMDE', $
        'PM', $
        'PMRA', $
        'PMDE', $
        'PLX_SBD', $
        'PLX_DTM', $
        'EPLX_DTM', $
        'JMAG', $
        'HMAG', $
        'KMAG', $
        'VMAG_SBD', $
        'VMAG_LSPM', $
        'VMAG_EST_LSPM', $
        'VMAG_IRTF', $
        'TEFF_APG', $
        'TEFF_VFIT', $
        'TEFF_IRTF', $
        'FEH_APG', $
        'FEH_IRTF', $
        'EFEH_IRTF', $
        'FEH_VFIT', $
        'FEH_RGRSS', $
        'VSINI_LIT', $
        'VSINI_LIT_LOSIG', $
        'VSINI_LIT_HISIG', $
        'VSINI_LIT_SRC', $
        'VSINI_VFIT', $
        'VSINI_VFIT_LOSIG', $
        'VSINI_VFIT_HISIG', $
        'RV_AVG', $
        'RV_AVG_IRTF', $
        'ROTATOR_IRTF']

; Initialize array of structures
outstr = []

; Re-make array of structures with desired tags
for i = 0, napg-1 do begin
    
    oldstr = masterstr[i]
    newstr = create_struct(tags, $
                           oldstr.objid, $
                           oldstr.ra, $
                           oldstr.dec, $
                           oldstr.rpm, $
                           oldstr.rpmra, $
                           oldstr.rpmde, $
                           oldstr.pm, $
                           oldstr.pmra, $
                           oldstr.pmde, $
                           oldstr.plx_sbd, $
                           !values.d_nan, $
                           !values.d_nan, $
                           oldstr.j, $
                           oldstr.h, $
                           oldstr.k, $
                           oldstr.v_sbd, $
                           oldstr.v_lspm, $
                           oldstr.v_est_lspm, $
                           !values.d_nan, $
                           oldstr.teff_apg, $
                           oldstr.teff_best, $
                           !values.d_nan, $
                           oldstr.feh_apg, $
                           oldstr.feh_irtf, $
                           !values.d_nan, $
                           oldstr.feh_best, $
                           oldstr.feh_sg, $
                           !values.d_nan, $
                           !values.d_nan, $
                           !values.d_nan, $
                           '', $
                           oldstr.vsini_best, $
                           !values.d_nan, $
                           !values.d_nan, $
                           !values.d_nan, $
                           !values.d_nan, $
                           !values.d_nan)
    ; Add to array of structures
    outstr = [outstr, newstr]

endfor
 
;;; Add Dittman PLX
dtm_plx_update, outstr
 
;;; Add literature vsini's
;NWT
nwt_vsini_update, outstr, nwt_idx
rnr_vsini_update, outstr, rnr_idx
dav_vsini_update, outstr, dav_idx

;;; Add irtf feh's, vmags, rv's, teff's and rotators
irtf_feh_update, outstr, irtf_idx
            
mwrfits, outstr, '/home/stgilhool/APOGEE/master_info.fits', /create
              
end
