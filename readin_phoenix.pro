; Read in the PHOENIX files, interpolate, normalize
function readin_phoenix, teff, feh, logg, wl_grid, nonorm=nonorm

if n_elements(nonorm) eq 0 then nonorm=0
base_path = '/home/stgilhool/PHOENIX/'
df = -8d0 ;flux conversion factor for PHOENIX templates
outstruct = []
xvector = dindgen(n_elements(wl_grid))

foreach logg_i, logg, logg_idx do begin
    
    ;logg_path = base_path + 'logg' + strtrim(logg_i,2) + '/'
    logg_path = base_path + 'logg' + logg_i + '/'
    
    foreach teff_i, teff, teff_idx do begin
        
        teff_name = strtrim(teff_i,2)
        
        foreach feh_i, feh, feh_idx do begin
        
            file_name = teff_name + feh_i + '_result.fits'
            file = logg_path + file_name
            fcheck = file_test(file)
            if fcheck then begin
                phospec = mrdfits(file,1)
                wl = phospec.wave
                spec = phospec.flux
                realflux = 10d0^(spec+df)
                ; Interpolate onto wl_grid
                outspec = interpol(realflux, wl, wl_grid)
                ; Normalize
                if nonorm eq 0 then begin
                    norm = continuum_fit(xvector, outspec)
                    specnorm = outspec/norm
                endif else specnorm=outspec
                struct = {spec:specnorm, feh:feh_i, teff:teff_i, logg:logg_i}
                outstruct = [outstruct, struct]
            endif
        endforeach
    endforeach
endforeach
            
return, outstruct

end
