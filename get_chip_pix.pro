; Quick script to get the pixel ranges of the three APG chips

function get_chip_pix

outfile = '/home/stgilhool/APOGEE/APOGEE_data/chip_data.fits'

; If this has never been run before, run it
; Otherwise, just read the file which is the output of this function
if file_test(outfile) eq 0 then begin

    nspec=1350
    specinfo = readin_apo2(nfiles=nspec)
    
    ; Initialize
    bovermin_vec = lonarr(nspec)
    bovermax_vec = lonarr(nspec)
    govermin_vec = lonarr(nspec)
    govermax_vec = lonarr(nspec)
    rovermin_vec = lonarr(nspec)
    rovermax_vec = lonarr(nspec)
    
    for i = 0, nspec-1 do begin
        
        head = *specinfo[i].header
        
        bovermin_vec[i] = fxpar(head, 'BOVERMIN')
        bovermax_vec[i] = fxpar(head, 'BOVERMAX')
        govermin_vec[i] = fxpar(head, 'GOVERMIN')
        govermax_vec[i] = fxpar(head, 'GOVERMAX')
        rovermin_vec[i] = fxpar(head, 'ROVERMIN')
        rovermax_vec[i] = fxpar(head, 'ROVERMAX')
        
    endfor
    
    ; Global min max's
    bmin = max(bovermin_vec)
    bmax = min(bovermax_vec)
    
    gmin = max(govermin_vec)
    gmax = min(govermax_vec)
    
    rmin = max(rovermin_vec)
    rmax = min(rovermax_vec)
    
    ; Get wl's
    bmin_wl = apg_pix2wl(bmin, /log_pos)
    bmax_wl = apg_pix2wl(bmax, /log_pos)
    
    gmin_wl = apg_pix2wl(gmin, /log_pos)
    gmax_wl = apg_pix2wl(gmax, /log_pos)
    
    rmin_wl = apg_pix2wl(rmin, /log_pos)
    rmax_wl = apg_pix2wl(rmax, /log_pos)
    
    ; Output
    outstr = {brange_pix:[bmin,bmax], $
              grange_pix:[gmin,gmax], $
              rrange_pix:[rmin,rmax], $
              brange_wl:[bmin_wl,bmax_wl], $
              grange_wl:[gmin_wl,gmax_wl], $
              rrange_wl:[rmin_wl,rmax_wl]}

    ; Write the structure to disk
    mwrfits, outstr, outfile, /create

endif else begin

    outstr = mrdfits(outfile,1)

endelse

return, outstr

end
