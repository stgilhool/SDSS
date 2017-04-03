pro find_tycho_vmags


; Read in all the APG spectra
;info_fil = '/home/stgilhool/APOGEE/master_info2.fits'
info_fil = '/home/stgilhool/APOGEE/master_info3.fits'

info = mrdfits(info_fil,1)

ra = info.ra
dec = info.dec

; Query vizier to get vmags
n_spec = 1350

out_struct = []
for spec = 0, n_spec-1 do begin
    
    tyc_arr = queryvizier('I/259/tyc2', [ra[spec],dec[spec]], 1d0/10d0, $
                            /allcolumns)
    
    ; check for no match
    null_result = size(tyc_arr, /type) eq 3
    
    if null_result then begin
        print, 'No matches for spec number: ' + strtrim(spec,2)
        vtmag  = !values.d_nan
        evtmag = !values.d_nan
        btmag  = !values.d_nan
        ebtmag = !values.d_nan
    endif else begin
        starcnt = n_elements(tyc_arr) 
        print, 'Found '+strtrim(starcnt,2)+' matches for spec number: '+strtrim(spec,2)
            
        ; check v-k criterion
        kmag = info[spec].kmag
        vmag = tyc_arr.vtmag
        v_k = vmag - kmag
        v_k_criterion = (v_k ge 2.5) and (v_k le 10.)
            
        v_k_idx = where(v_k_criterion eq 1, vkcnt)
        if vkcnt gt 0 then begin
            tyc_arr = tyc_arr[v_k_idx]
            rad = tyc_arr._r
            
            ; check pm criterion
            apg_pmra  = info[spec].pmra*1d3
            apg_pmde  = info[spec].pmde*1d3
            tyc_pmra = tyc_arr.pmra
            tyc_pmde = tyc_arr.pmde
            
            if vkcnt eq 1 then begin
                print, 'Final match found for spec number: '+strtrim(spec,2)
                vtmag  = double(tyc_arr.vtmag)
                evtmag = double(tyc_arr.e_vtmag)
                btmag  = double(tyc_arr.btmag)
                ebtmag = double(tyc_arr.e_btmag)
            endif else begin
                print, 'Please choose which is the best match'
                print, format='("APG_PMRA:  ", F6.1, " || ", T25, 100(F6.1, :, " | "))',apg_pmra, tyc_pmra
                print, format='("APG_PMDE:  ", F6.1, " || ", T25, 100(F6.1, :, " | "))',apg_pmde, tyc_pmde
                stop
                print, 'Choosing the object with the least separation'
                read, gdidx
                tyc_arr = tyc_arr[gdidx]
                ;closest_idx = where(tyc_arr._r eq min(tyc_arr._r), closecnt)
                ;if closecnt gt 1 then closest_idx = closest_idx[0]
                ;tyc_arr = tyc_arr[closest_idx]
                vtmag  = double(tyc_arr.vtmag)
                evtmag = double(tyc_arr.e_vtmag)
                btmag  = double(tyc_arr.btmag)
                ebtmag = double(tyc_arr.e_btmag)
            endelse
        endif else begin
            print, 'No matches within constraints for spec number: ' + strtrim(spec,2)
            vtmag  = !values.d_nan
            evtmag = !values.d_nan
            btmag  = !values.d_nan
            ebtmag = !values.d_nan
        endelse
    endelse

; Add to master_info
new_struct = create_struct(info[spec], 'VTMAG', vtmag, 'E_VTMAG', evtmag, $
                           'BTMAG', btmag, 'E_BTMAG', ebtmag)
out_struct = [out_struct, new_struct]
endfor
stop
mwrfits, out_struct, '/home/stgilhool/APOGEE/master_info4.fits', /create

end
