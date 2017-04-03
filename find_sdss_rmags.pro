pro find_sdss_rmags


; Read in all the APG spectra
info_fil = '/home/stgilhool/APOGEE/master_info4.fits'

info = mrdfits(info_fil,1)

ra = info.ra
dec = info.dec

; Query vizier to get rmags
n_spec = 1350

out_struct = []
for spec = 0, n_spec-1 do begin
    
    sdss_arr = queryvizier('V/139/sdss9', [ra[spec],dec[spec]], 1d0/10d0, $
                            /allcolumns)
    
    ; check for no match
    null_result = size(sdss_arr, /type) eq 3
    
    if null_result then begin
        print, 'No matches for spec number: ' + strtrim(spec,2)
        rmag  = !values.d_nan
        ermag = !values.d_nan
        rrmag = !values.d_nan
        ;umag  = !values.d_nan
        ;eumag = !values.d_nan
    endif else begin
        ; make sure the sdss objs are stars
        star_idx = where((sdss_arr.cl eq 6) and $
                         (sdss_arr.gc eq 6), $
                         starcnt)
        if starcnt gt 0 then begin
            print, 'Found '+strtrim(starcnt,2)+' matches for spec number: '+strtrim(spec,2)
            sdss_arr = sdss_arr[star_idx]
            ; check v-k criterion
            kmag = info[spec].kmag
            vmag = sdss_arr.g_mag
            v_k = vmag - kmag
            v_k_criterion = (v_k ge 2.5) and (v_k le 10.)
            
            v_k_idx = where(v_k_criterion eq 1, vkcnt)
            if vkcnt gt 0 then begin
                sdss_arr = sdss_arr[v_k_idx]
                rad = sdss_arr._r

                ; check pm criterion
                apg_pmra  = info[spec].pmra*1d3
                apg_pmde  = info[spec].pmde*1d3
                sdss_pmra = sdss_arr.pmra
                sdss_pmde = sdss_arr.pmde
                
                if vkcnt eq 1 then begin
                    print, 'Final match found for spec number: '+strtrim(spec,2)
                    ;gmag  = double(sdss_arr.g_mag)
                    ;egmag = double(sdss_arr.e_g_mag)
                    ;umag  = double(sdss_arr.u_mag)
                    ;eumag = double(sdss_arr.e_u_mag)
                    rmag = double(sdss_arr.r_mag)
                    rrmag = double(sdss_arr.rmag)
                    ermag = double(sdss_arr.e_r_mag)
                endif else begin
                    ;print, 'Please choose which is the best match'
                    ;print, format='("APG_PMRA:  ", F6.1, " || ", T25, 100(F6.1, :, " | "))',apg_pmra, sdss_pmra
                    ;print, format='("APG_PMDE:  ", F6.1, " || ", T25, 100(F6.1, :, " | "))',apg_pmde, sdss_pmde
                    ;stop
                    print, 'Choosing the object with the least separation'
                    closest_idx = where(sdss_arr._r eq min(sdss_arr._r), closecnt)
                    if closecnt gt 1 then closest_idx = closest_idx[0]
                    sdss_arr = sdss_arr[closest_idx]
                    ;gmag  = double(sdss_arr.g_mag)
                    ;egmag = double(sdss_arr.e_g_mag)
                    ;umag  = double(sdss_arr.u_mag)
                    ;eumag = double(sdss_arr.e_u_mag)
                    rmag = double(sdss_arr.r_mag)
                    rrmag = double(sdss_arr.rmag)
                    ermag = double(sdss_arr.e_r_mag)
                endelse
            endif
        endif else begin
            print, 'No STAR matches for spec number: ' + strtrim(spec,2)
            rmag  = !values.d_nan
            ermag = !values.d_nan
            rrmag = !values.d_nan
            ;umag  = !values.d_nan
            ;eumag = !values.d_nan
        endelse
    endelse
    
;     if finite(gmag) then begin
;         print, strtrim(ra[spec],2) + ' | ' + strtrim(sdss_arr.raj2000,2)
;         print, strtrim(dec[spec],2) + ' | ' + strtrim(sdss_arr.dej2000,2)
;         print, "class (3-gal, 6-star) = " + strtrim(sdss_arr.cl,2)
;         help, sdss_arr.spcl
;         help, sdss_arr.subclass
;         help, sdss_arr.us 
;         help, sdss_arr.uc
;         help, sdss_arr.gs
;         help, sdss_arr.gc
;         stop
;     endif 
    print, format = '(T20, F, 3X, F, 3X, F)', rmag, rrmag, rmag-rrmag
    ; Add to master_info
    new_struct = create_struct(info[spec], 'R_MAG', rmag, 'E_R_MAG', ermag)
    out_struct = [out_struct, new_struct]
endfor
stop
mwrfits, out_struct, '/home/stgilhool/APOGEE/master_info5.fits', /create

end
