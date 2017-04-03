pro houdebine_readin
; Match the 2MASS names from the Houdebine SIMBAD file
; to our APOGEE 2MASS names to determine if there are matches
houdpath = '/home/stgilhool/APOGEE/vsini_literature/houdebine/'
houdfile = houdpath + 'houdebine2MASS_pipedelimited.txt'

apg_file = '/home/stgilhool/APOGEE/master_info5.fits'
apg_info = mrdfits(apg_file,1)
; Read in apg IDs
apg_2M = apg_info.objid

; Read in Houdebine IDs
readcol, houdfile, simnum, houd_2M, houd_id, format = 'L,A,A', delimiter='|'

nhoud = n_elements(houd_2M) 

simnum_vec = []
houd_2M_vec = []
houd_id_vec = []
apg_idx_vec = []
; Loop through the Houdebine files and retrieve the APG info, if matched
for houdnum = 0, nhoud - 1 do begin
    houd_2M_i = strtrim(houd_2M[houdnum],2)
    houd_2M_fmt = mg_streplace(houd_2M_i, 'J', '2M')
    houd_id_i = strtrim(houd_id[houdnum],2)
    
    houd_2M_vec = [houd_2M_vec, houd_2M_fmt]
    houd_id_vec = [houd_id_vec, houd_id_i]

    houd_idx = where(apg_2M eq houd_2M_fmt, match_chk)
    if match_chk eq 1 then begin
        apg_idx_i = houd_idx[0]
        
    endif else apg_idx_i = !values.f_nan
    apg_idx_vec = [apg_idx_vec, apg_idx_i]

endfor


stop
;Read in table 5
houdtab5file = houdpath+'houd_tab5_readcol_pipe.txt'

readcol, houdtab5file, h5id, vmag5, R_I5, teff5, e_teff5, spt5, pi5, e_pi5, $
  Mv5, e_Mv5, rstar5, e_rstar5, vsini5, p_sini5, e_p5_lo, e_p5_hi, $
  feh5, feh_lit5, upper_limit5, delimiter = "|",$
  format = "A,F,F,L,L,A,F,F,F,F,F,F,F,F,F,F,F,F,L", /preserve_null
tab5_vec = replicate(5, n_elements(h5id)) 

;Read in table 4
houdtab4file = houdpath+'houd_tab4_readcol_pipe.txt'

readcol, houdtab4file, h4id, vmag4, R_I4, teff4, e_teff4, spt4, pi4, e_pi4, $
  Mv4, e_Mv4, rstar4, e_rstar4, vsini4, p_sini4, e_p4_lo, e_p4_hi, $
  feh4, feh_lit4, upper_limit4, delimiter = "|",$
  format = "A,F,F,L,L,A,F,F,F,F,F,F,F,F,F,F,F,F,L", /preserve_null
tab4_vec = replicate(4, n_elements(h4id)) 

;Read in table 3
houdtab3file = houdpath+'houd_tab3_readcol_pipe.txt'

readcol, houdtab3file, h3id, vmag3, R_I3, teff3, e_teff3, spt3, pi3, e_pi3, $
  Mv3, e_Mv3, rstar3, e_rstar3, vsini3, p_sini3, e_p3_lo, e_p3_hi, $
  feh3, upper_limit3, delimiter = "|",$
  format = "A,F,F,L,L,A,F,F,F,F,F,F,F,F,F,F,F,L", /preserve_null
tab3_vec = replicate(3, n_elements(h3id)) 

stop

; Output for each table
tab_vec = [tab3_vec,tab4_vec,tab5_vec]
id_vec = [h3id, h4id, h5id]
spt_vec = [spt3, spt4, spt5]
ul_vec = [upper_limit3,upper_limit4,upper_limit5]

;vmag
ntmp_idx = where(vmag3 eq -99999, nnan)
if nnan gt 0 then vmag3[ntmp_idx] = !values.f_nan
ntmp_idx = where(vmag4 eq -99999, nnan)
if nnan gt 0 then vmag4[ntmp_idx] = !values.f_nan
ntmp_idx = where(vmag5 eq -99999, nnan)
if nnan gt 0 then vmag5[ntmp_idx] = !values.f_nan
vmag_vec = [vmag3, vmag4, vmag5]

;teff
ntmp_idx = where(teff3 eq -99999, nnan)
if nnan gt 0 then teff3[ntmp_idx] = !values.f_nan
ntmp_idx = where(teff4 eq -99999, nnan)
if nnan gt 0 then teff4[ntmp_idx] = !values.f_nan
ntmp_idx = where(teff5 eq -99999, nnan)
if nnan gt 0 then teff5[ntmp_idx] = !values.f_nan
teff_vec = [teff3, teff4, teff5]

;vsini
ntmp_idx = where(vsini3 eq -99999, nnan)
if nnan gt 0 then vsini3[ntmp_idx] = !values.f_nan
ntmp_idx = where(vsini4 eq -99999, nnan)
if nnan gt 0 then vsini4[ntmp_idx] = !values.f_nan
ntmp_idx = where(vsini5 eq -99999, nnan)
if nnan gt 0 then vsini5[ntmp_idx] = !values.f_nan
vsini_vec = [vsini3, vsini4, vsini5]

outstr = {id:id_vec,spt:spt_vec,teff:teff_vec,vsini:vsini_vec,upper_limit:ul_vec,table:tab_vec}

help, outstr, /str
stop

mwrfits, outstr, 'houdebine_all.fits', /create

end
