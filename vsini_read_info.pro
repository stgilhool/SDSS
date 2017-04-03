; script to read in a particular info file, or the latest one

function vsini_read_info, RUN=run, RECENT=recent, SILENT=silent

path_base = '/home/stgilhool/APOGEE/vsini_results/'

if not keyword_set(silent) then silent = 0 else silent=1
; If recent isn't set
if not keyword_set(recent) then recent = 0
; But if run isn't set, automatically switch to recent
if n_elements(run) eq 0 then recent = 1 $
else begin
    ; If run is set, try to load the corresponding file
    recent = 0
    ;Check that run is appropriate type ADD
    ;Check that corresponding file exists
    input_dir = 'run_'+strtrim(run,2)+'/'
    input_path = path_base + input_dir
    infile = input_path + 'run_'+strtrim(run,2)+'_info.fits'
    
    if file_test(infile) then begin
        info_str = mrdfits(infile, 1)
    endif else begin
        ; If run file doesn't exist, give user the option to load most recent
        print, infile + " does not exist! Loading most recent file (y/n)?"
        repeat begin
            opt = get_kbrd()
            if opt eq 'n' or opt eq 'N' then begin
                print, "Exiting..."
                key_term = 1
                return, -1
            endif else if opt eq 'y' or opt eq 'Y' then begin
                print, "Loading most recent file"
                key_term = 1
                recent = 1
            endif else begin
                print, "Not a valid choice"
                key_term = 0
            endelse
        endrep until key_term eq 1
    endelse
endelse

if recent then begin
    ; Search for most recent info file and read it in
    ;;; CAUTION: this assumes there's at least one file
    run = 0
    terminate = 0
    repeat begin
        
        input_dir = 'run_'+strtrim(run,2)+'/'
        input_path = path_base + input_dir
        infile = input_path + 'run_'+strtrim(run,2)+'_info.fits'
                        
        if file_test(infile) then begin
            ; If the file exists, read it in and continue searching
            info_str = mrdfits(infile, 1)
            run++
        endif else begin
            ; If there's no file, stop iterating
            terminate = 1
        endelse
        
    endrep until terminate eq 1
endif    

if silent eq 0 then begin
    print, "====================================================================="
    print, "Read file " + infile
    print, "====================================================================="
    print, format='("Oversamp: ", T20, A)', strtrim(info_str.oversamp,2)
    print, format='("dlogwl: ", T20, A)', strtrim(info_str.dlogwl,2)
    print, format='("df: ", T20, A)', strtrim(info_str.df,2)
    print, format='("deltav: ", T20, A)', strtrim(info_str.deltav,2)
    print, format='("nplsf: ", T20, A)', strtrim(info_str.nplsf,2)
    print, format='("APG files: ", T20, A, " - ", A)', $
      strtrim(min(info_str.apg_fnums),2), $
      strtrim(max(info_str.apg_fnums),2)
    print, format='("APG type: ", T20, A)', strtrim(info_str.apg_type,2)
    print, format='("APG mask: ", T20, A)', strtrim(info_str.apg_mask,2)
    print, format='("APG pixels: ", T20, A, " - ", A)', $
      strtrim(info_str.apg_xr[0],2), $
      strtrim(info_str.apg_xr[1],2)
    print, format='("logg: ", T20, A)', strtrim(info_str.logg,2)
    print, format='("FeH grid pts: ", T20, 15(A, :, ", "))', $
      strtrim(info_str.feh_grid,2)
    print, format='("Teff grid range: ", T20, A, " - ", A)', $
      strtrim(min(info_str.teff_grid),2), $
      strtrim(max(info_str.teff_grid),2)
    ;print, format='("Teff grid pts: ", T20, 25(A, :, ", "))', $
      ;strtrim(info_str.teff_grid,2)
    print, format='("vsini grid pts: ", T20, 25(A5, :, ", "))', $
      strtrim(info_str.vsini_grid,2)
endif

return, info_str

end
