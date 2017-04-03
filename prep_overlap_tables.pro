pro prep_overlap_tables

;;; Paths

;aspcap for comparison
aspcap_file = '/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits'

; my stuff
sg_path = '/home/stgilhool/APOGEE/'
sg_file = sg_path + 'master_info6.fits'

; Readin my results
apg_sg_str = mrdfits(sg_file,1)

; Readin aspcap results
aspcap_str = mrdfits(aspcap_file,1)
;;;;;;;


;;;;;;;

; Get the indices for Vsini literature guys
lit_idx = where(finite(apg_sg_str.vsini_lit), nlit)
if nlit ne 0 then begin
    
    apg_sg = apg_sg_str[lit_idx]
    aspcap = aspcap_str[lit_idx]

endif else message, "Error in selecting literature overlaps"

;;; Readin apStar info
max_idx = max(lit_idx)

; Readin APG headers
apg_str = readin_apo(nfiles = max_idx+1, hdu=8)
; keep just the overlaps
spec_info = apg_str[lit_idx]
lsf_coeff = spec_info.output

final_struct = []
;;; Create each structure
for snum = 0, n_elements(lit_idx)-1 do begin

    ; get the structure from master_info
    minfo = apg_sg[snum]
    ; get the ASPCAP structure, and its teff and fe_h
    aspinfo = aspcap[snum]
    teff_asp = aspinfo.teff
    feh_asp = aspinfo.fe_h
    vsini_asp = aspinfo.vsini
    ; get the spectrum, etc.
    spec = spec_info[snum].spectra
    err = spec_info[snum].errors
    lsf = *lsf_coeff[snum]

    newstruct = create_struct(minfo, 'TEFF_ASPCAP', teff_asp, 'FEH_ASPCAP', feh_asp, 'VSINI_ASPCAP', vsini_asp, 'SPEC', spec, 'ESPEC', err, 'LSF_CO', lsf)

    final_struct = [final_struct, newstruct]

endfor

;stop
mwrfits, final_struct, '/home/stgilhool/APOGEE/overlap_table.fits', /create

stop
end
