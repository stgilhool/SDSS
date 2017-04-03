pro add_vmags_to_master_info


; Read in all the DTM parallax entries
info_fil = '/home/stgilhool/APOGEE/master_info.fits'

info = mrdfits(info_fil,1)

ra = info.ra
dec = info.dec

; Query vizier to get vmags
n_spec = 1350
vmag_apass = dblarr(n_spec)
evmag_apass = dblarr(n_spec)
out_struct = []
for spec = 0, n_spec-1 do begin
    
    apass_arr = queryvizier('II/336/apass9', [ra[spec],dec[spec]], 1d0/10d0, $
                            /allcolumns)
    if n_elements(apass_arr) gt 1 then begin
        print, 'multiple matches for spec number: ' + strtrim(spec,2)
        match_idx = where(apass_arr._r eq min(apass_arr._r))
    
        apass = apass_arr[match_idx]
        vmag_apass[spec] = apass.g_mag
        evmag_apass[spec] = apass.e_g_mag
    endif else if n_elements(apass_arr) eq 1 then begin
        if size(apass_arr, /type) ne 3 then begin
            print, 'Exactly 1 match for spec number: ' + strtrim(spec,2)
            vmag_apass[spec] = apass_arr.g_mag
            evmag_apass[spec] = apass_arr.e_g_mag
        endif else begin
            print, 'no matches for spec number: ' + strtrim(spec,2)
            vmag_apass[spec] = !values.d_nan
            evmag_apass[spec] = !values.d_nan
        endelse
    endif

    ; Add to master_info
    new_struct = create_struct(info[spec], 'VMAG_APASS', vmag_apass[spec], 'EVMAG_APASS', evmag_apass[spec])
    out_struct = [out_struct, new_struct]
endfor

mwrfits, out_struct, '/home/stgilhool/APOGEE/master_info2.fits', /create

end
