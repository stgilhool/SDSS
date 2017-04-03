; Control file for vsini_grid
; Main program runs this to assign values to important variables
pro vsini_grid_info

common grid_info, bigc, oversamp, dlogwl, df, deltav, nplsf, npad, $
  apg_fnums, apg_type, apg_xr, $
  logg, feh_grid, teff_grid, vsini_grid

; Parameters that are unlikely to change
bigc = 2.99792458d5 ;km/s
oversamp = 7L
dlogwl= 6d-6
df= -8d0
deltav = ((10d0^(dlogwl/oversamp))-1d0)*bigc
nplsf = 141
npad = 5

; APOGEE info
apg_fnums = lindgen(1350) ;max 1350
apg_type = 'aspcapStar' ;aspcapStar or apStar
apg_xr = [2500,3200] ;must both be even or both be odd
apg_mask = 'none'

; Grid info
logg = 4.5
feh_grid = ['-4.0', '-3.5', '-3.0', '-2.5', '-2.0', '-1.5', '-1.0', '-0.5', '-0.0', $
            '+0.3', '+0.5']
teff_range = [2000,4200]
teff_step = 100
teff_grid = lindgen((teff_range[1]-teff_range[0])/teff_step + 1) * teff_step $
  + teff_range[0]

vsini_grid = [4d0,6d0,10d0,15d0,20d0,25d0,30d0,35d0,40d0,50d0,60d0,70d0,80d0,90d0,100d0]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;END DEFINITIONS.  BEGIN FILE MANAGEMENT;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

out_str = {bigc:bigc, $
          oversamp:oversamp, $
          dlogwl:dlogwl, $
          df:df, $
          deltav:deltav, $
          nplsf:nplsf, $
          npad:npad, $
          apg_fnums:apg_fnums, $
          apg_type:apg_type, $
          apg_mask:apg_mask, $
          apg_xr:apg_xr, $
          logg:logg, $
          feh_grid:feh_grid, $
          teff_grid:teff_grid, $
          vsini_grid:vsini_grid $
          }
          
outpath_base = '/home/stgilhool/APOGEE/vsini_results/'
runnum = 0 
terminate = 0

; Decide where to write structure of relevant variables, and write it
repeat begin
    
    outdir = 'run_'+strtrim(runnum,2)+'/'
    outpath = outpath_base + outdir
    outfile_name = outpath + 'run_'+strtrim(runnum,2)+'_info.fits'

    if file_test(outpath, /directory) then begin
        ; If the parent dir exists,
        ; Check if there's a parameter file in the directory
        if file_test(outfile_name) then begin
            ; If so, read it in and compare with current
            comp_str = mrdfits(outfile_name, 1)
            diff_list = compare_struct(comp_str, out_str)
            if n_elements(diff_list) eq 1 and diff_list[0].ndiff eq 0 then begin
                ; If it's the same, don't write to the directory, just stop
                print, "File already exists. Exiting..."
                terminate = 1
            endif else begin
                ; If there's a file and it's different, go on to next
                runnum++
            endelse
        endif else begin
            ; If there's no file, write the file and stop
            mwrfits, out_str, outfile_name, /create
            print, "Wrote file: "+outfile_name
            terminate = 1
        endelse
    endif else begin
        ; If there's no dir, make the directory and write the file
        file_mkdir, outpath
        print, "Created dir: " + outpath
        mwrfits, out_str, outfile_name, /create
        print, "Wrote file: "+outfile_name
        terminate = 1
    endelse

endrep until terminate eq 1

return

end
