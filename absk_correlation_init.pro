; Find features in APG spectra that correlate with parallax

pro absk_correlation_init

; Options
pixwt_or_glbwt = 'pixwt'
mask_flags_bit = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
apstar_or_aspcapstar = 'apstar'

mask_flags_dec = long(total(2L^mask_flags_bit))

; Readin PLX info
apg_info_file = '/home/stgilhool/APOGEE/master_info2.fits'
apg_info = mrdfits(apg_info_file, 1)

apg_info_id = apg_info.objid
apg_info_plx = apg_info.plx_dtm
apg_info_eplx = apg_info.eplx_dtm
apg_kmag = apg_info.kmag

;;; Readin spectra
apg_spec = get_apg_flux()

nspec = n_elements(apg_spec) 
npix = 8575L
pix_vec = dindgen(npix)
spec = dblarr(npix, nspec)
err  = dblarr(npix, nspec)
mask = dblarr(npix, nspec)
apg_spec_id = strarr(nspec)

outstruct = []

;;; Get spec, errs and masks, normalize, and apply masks
for specnum = 0, nspec-1 do begin
    
    ;Save objid
    apg_spec_id[specnum] = apg_spec[specnum].objid
    ; Check if it matches plx id
    if apg_spec_id[specnum] ne apg_info[specnum].objid $
      then message, "PLX structure doesn't match APG structure"
    
    ; Choose which type of weighting and get vectors
    case pixwt_or_glbwt of
        'pixwt': begin
            ; Get vectors for normalization
            current_flux = apg_spec[specnum].flux_pixwt
            current_ferr = apg_spec[specnum].e_flux_pixwt
        end
        'glbwt': begin
            ; Get vectors for normalization
            current_flux = apg_spec[specnum].flux_glbwt
            current_ferr = apg_spec[specnum].e_flux_glbwt
        end
    endcase

    ; Get other info for normalization
    brange = apg_spec[specnum].brange
    grange = apg_spec[specnum].grange
    rrange = apg_spec[specnum].rrange
    ; Normalize and Scale the errors as well
    current_spec = apg_continuum_fit(current_flux, brange, grange, rrange, $
                                     other_vec = current_ferr)
    
    ; Get the mask
    current_mask = apg_spec[specnum].mask_pixwt
    ; Determine which pixels are masked, given selected flags
    select_mask = current_mask and mask_flags_dec
    flg_idx = where(select_mask ne 0, nflg, complement=unflg_idx, $
                    ncomplement=nunflg)
    masked_mask = select_mask
    masked_mask[flg_idx] = 1
    ; Mask spec pixels
    ;masked_spec = interpol(current_spec[unflg_idx],
    ;unflg_idx, pix_vec)
    
    ; Adjust error pixels
    masked_err = current_ferr
    masked_err[flg_idx] = 1d8

    ; Save results
    ;spec[*,specnum] = current_spec
    ;err[*,specnum]  = masked_err
    ;mask[*,specnum] = masked_mask
    
    ; Add spectra, err, masks, etc to apg_info
    newstruct = create_struct(apg_info[specnum], 'SPEC',current_spec, $
                              'ERR',current_ferr, 'MASK',select_mask)
                              ;'FLG_IDX':flg_idx, 'NFLG':nflg, $
                              ;'GD_IDX':unflg_idx, 'NGD':nunflg)
    outstruct = [outstruct, newstruct]

endfor

; Write result
mwrfits, outstruct, '/home/stgilhool/APOGEE/absk_corr/absk_corr_init.fits', /create

end
