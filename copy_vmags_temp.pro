pro copy_vmags_temp


; Read in all the DTM parallax entries
info_fil = '/home/stgilhool/APOGEE/master_info.fits'
info_fil2 = '/home/stgilhool/APOGEE/master_info2.fits'

info = mrdfits(info_fil,1)
info2 = mrdfits(info_fil2,1)

ra = info.ra
dec = info.dec

n_spec = 1350
vmag_apass = info2.vmag_apass
evmag_apass = info2.evmag_apass
out_struct = []
for spec = 0, n_spec-1 do begin
    
    
    ; Add to master_info
    new_struct = create_struct(info[spec], 'VMAG_APASS', vmag_apass[spec], 'EVMAG_APASS', evmag_apass[spec])
    out_struct = [out_struct, new_struct]
endfor

mwrfits, out_struct, '/home/stgilhool/APOGEE/master_info2.fits', /create

end
