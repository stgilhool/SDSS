; Read in all PHOENIX templates, sort by wavelength, trim to APOGEE
; wavelengths, and write to fits

pro phoenix_reformat

;;; Read in PHOENIX
ppath = '/home/stgilhool/PHOENIX/'

pfeh = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
        '+0.3', '+0.5']

nteff = 20
teff_start = 20
pteff = lindgen(nteff)+teff_start
nfeh = n_elements(pfeh)
vsini_vec = dindgen(40)/2. + 8

cnt = 0
spectra = ptrarr(nteff, nfeh, /allocate_heap)
wl_grid = ptrarr(nteff, nfeh, /allocate_heap)

;;;Loop over teff
foreach teff, pteff, teff_idx do begin
        ;;;Loop over metallicity
    for feh = 0, nfeh-1 do begin

        pfile = ppath + 'logg5/lte0'+strtrim(teff,2)+'-5.0'+pfeh[feh]+'.BT-Settl.7'
            
        if file_test(pfile) eq 0 then begin
            *(spectra[teff_idx, feh]) = !values.d_nan
            *(wl_grid[teff_idx, feh]) = !values.d_nan
            print, "No spectrum number: "+strtrim(cnt,2)
            print, "with Teff = "+strtrim(teff*100L,2)+" K"
            print, "and FeH = "+pfeh[feh]
            print, ""
            cnt++
            continue
        endif
        df = -8d0
        ; Read in PHOENIX model
        readcol, pfile, pwl, c1, c2, c3, c4, format = "D,D,D,D,D"
        pws = sort(pwl)
        pwl = pwl[pws]
        
        c1 = c1[pws]
        c2 = c2[pws]
        c3 = c3[pws]
        c4 = c4[pws]
        
        xr_idx = where(pwl ge 1.4d4 and pwl le 1.7d4, nwl)
        print, "Number of APOGEE wl elements: "+strtrim(nwl,2)
        wl = pwl[xr_idx]
        f1 = c1[xr_idx]
        b1 = c2[xr_idx]
        c3 = c3[xr_idx]
        c4 = c4[xr_idx]
        
        *(spectra[teff_idx,feh]) = f1
        *(wl_grid[teff_idx,feh]) = wl
        
        print, "Saved spectrum number: "+strtrim(cnt,2)
        print, "with Teff = "+strtrim(teff*100L,2)+" K"
        print, "and FeH = "+pfeh[feh]
        print, ""
        cnt++
    endfor
endforeach

; Write a structure
outstr = {teff:pteff*100L, $
          feh:pfeh, $
          wl_grid:wl_grid, $
          spectra:spectra $
          }
outfile = '/home/stgilhool/PHOENIX/phoenix_cube.fits'
mwrfits, outstr, outfile, /create

stop
end
