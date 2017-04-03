pro find_sdss_vmags


; Read in all the APG spectra
info_fil = '/home/stgilhool/APOGEE/master_info2.fits'

info = mrdfits(info_fil,1)

ra = info.ra
dec = info.dec

; Query vizier to get vmags
n_spec = 1350

out_struct = []
for spec = 0, n_spec-1 do begin
    
    sdss_arr = queryvizier('V/139/sdss9', [ra[spec],dec[spec]], 1d0/10d0, $
                            /allcolumns)
    
    if n_elements(sdss_arr) gt 1 then begin
        print, 'multiple matches for spec number: ' + strtrim(spec,2)
                
        ;print, "Which one is correct? Enter idx"
        ;stop
        star_idx = where(sdss_arr.cl eq 6, starcnt)
        if starcnt gt 0 then begin
            idx = star_idx[0]
            sdss_arr = sdss_arr[idx]
            gmag  = double(sdss_arr.g_mag)
            egmag = double(sdss_arr.e_g_mag)
            umag  = double(sdss_arr.u_mag)
            eumag = double(sdss_arr.e_u_mag)
        endif else begin
            print, 'no STAR matches for spec number: ' + strtrim(spec,2)
            gmag  = !values.d_nan
            egmag = !values.d_nan
            umag  = !values.d_nan
            eumag = !values.d_nan
        endelse
    endif else if n_elements(sdss_arr) eq 1 then begin
        if size(sdss_arr, /type) ne 3 then begin
            print, 'Exactly 1 match for spec number: ' + strtrim(spec,2)
            if sdss_arr.cl ne 6 then begin
                print, "But it's not a star"
                gmag  = !values.d_nan
                egmag = !values.d_nan
                umag  = !values.d_nan
                eumag = !values.d_nan
            endif else begin
                gmag  = double(sdss_arr.g_mag)
                egmag = double(sdss_arr.e_g_mag)
                umag  = double(sdss_arr.u_mag)
                eumag = double(sdss_arr.e_u_mag)
            endelse
        endif else begin
            print, 'no matches for spec number: ' + strtrim(spec,2)
            gmag  = !values.d_nan
            egmag = !values.d_nan
            umag  = !values.d_nan
            eumag = !values.d_nan
        endelse
    endif

    if finite(gmag) then begin
        print, strtrim(ra[spec],2) + ' | ' + strtrim(sdss_arr.raj2000,2)
        print, strtrim(dec[spec],2) + ' | ' + strtrim(sdss_arr.dej2000,2)
        print, "class (3-gal, 6-star) = " + strtrim(sdss_arr.cl,2)
        help, sdss_arr.spcl
        help, sdss_arr.subclass
        help, sdss_arr.us 
        help, sdss_arr.uc
        help, sdss_arr.gs
        help, sdss_arr.gc
    endif 

    ; Add to master_info
    new_struct = create_struct(info[spec], 'G_MAG', gmag, 'E_G_MAG', egmag, $
                               'U_MAG', umag, 'E_U_MAG', eumag)
    out_struct = [out_struct, new_struct]
endfor
stop
mwrfits, out_struct, '/home/stgilhool/APOGEE/master_info3.fits', /create

end
