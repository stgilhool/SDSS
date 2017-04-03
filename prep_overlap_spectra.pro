;+
; NAME: prep_overlap_spectra
;
;
;
; PURPOSE: Get needed info on Vsini literature overlaps and
;          prepare the APOGEE spectra for Q test by masking 
;          and smoothing them
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-



;;;; Util function to widen the pixel mask
function widen_mask, mask, width

wide_mask = smooth(mask, width, /edge_mirror)

mask_pix = where(wide_mask gt 0, nmask)

return_mask = wide_mask

if nmask gt 0 then return_mask[mask_pix] = 1
    
return, return_mask

end
;;;;

;;;;; Use the 2nd derivative as a cut for noisy pixels
function mask_by_deriv, spec
; First deriv
dspec = spec[1:*] - spec[0:-2]
; Second deriv
ddspec = dspec[1:*] - dspec[0:-2]

; Get stddev of 2nd deriv
dev = stddev(ddspec)

ddscat = abs(ddspec)

; Make mask

return, dev

end


;;;;; Use the 2nd derivative as a cut for noisy pixels
function mask_by_deriv2, spec
; First deriv
dspec = spec[1:*] - spec[0:-2]
; Second deriv
ddspec = dspec[1:*] - dspec[0:-2]

; Get stddev of 2nd deriv
dev = stddev(ddspec)

ddscat = abs(ddspec)

; Make mask
mask_pix = where(ddscat gt dev, nmask)

return, nmask

end


function mask_by_deriv3, spec
; First deriv
dspec = spec[1:*] - spec[0:-2]
; Second deriv
ddspec = dspec[1:*] - dspec[0:-2]

; Get stddev of 2nd deriv
dev = stddev(ddspec)
;dev = 0.1d0

ddscat = abs(ddspec)

; Make mask
mask_pix = where(ddscat gt 2*dev, nmask)

newspec = spec

; Mask the pixels
newspec[mask_pix+1] = !values.d_nan


return, newspec

end

function make_pixmask, aspmask

asp_flag_bits = [0,1,2,3,4,5,6,12,13,14]
asp_flag_dec = long(total(2L^asp_flag_bits))

mask_vec = asp_flag_dec and aspmask

return, mask_vec

end


function mask_by_pixmask, spec, aspmask

; Get mask vec
mask_vec = make_pixmask(aspmask)

; Get indices
mask_idx = where(mask_vec ne 0, nmsk)

; Mask
masked_spec = spec
if nmsk gt 0 then begin
    masked_spec[mask_idx] = !values.d_nan
endif

return, masked_spec

end

function mask_by_pixmask2, spec, aspmask, medspec

; Get mask vec
mask_vec = make_pixmask(aspmask)

; Get indices
mask_idx = where(mask_vec ne 0, nmsk)

; Mask
masked_spec = spec
if nmsk gt 0 then begin
    masked_spec[mask_idx] = medspec[mask_idx]
endif

return, masked_spec

end


pro prep_overlap_spectra, skip_load=skip_load

if n_elements(skip_load) eq 0 then skip_load = 1 
save_file_name = 'prep_overlap_spectra.sav'


; Gloabals
npix = 8575L
oversamp = 7
xvec = lindgen(npix)
xxvec = (dindgen(npix*oversamp)-(oversamp/2))/double(oversamp)
nxx = n_elements(xxvec)  

nplsf = 141
lsf_0 = xxvec[nxx/2]
lx = (dindgen(nplsf)-(nplsf/2))/oversamp + lsf_0


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

; Get rid of spectrum 22
lit0 = lit_idx[0:21]
lit1 = lit_idx[23:*]
lit_idx = [lit0,lit1]
nlit = n_elements(lit_idx) 

if nlit ne 0 then begin
    
    apg_sg = apg_sg_str[lit_idx]
    aspcap = aspcap_str[lit_idx]

endif else message, "Error in selecting literature overlaps"

if skip_load then goto, skip_load_mark

;;; Readin apStar info
max_idx = max(lit_idx)

; Readin APG headers
apg_str = readin_apo(nfiles = max_idx+1, hdu=8)
; keep just the overlaps
spec_info = apg_str[lit_idx]
lsf_coeff = spec_info.output
; delete
apg_str = 0

; Get mask (3)
mask_str = readin_apo(nfiles = max_idx+1, hdu=3)
temp_str = mask_str[lit_idx]

mask = temp_str.output

mask_arr = lonarr(npix,nlit)
for i = 0, nlit-1 do begin

    maski = *mask[i]
    maskii = maski[*,0]
    mask_arr[*,i] = maskii

endfor

;delete 
mask_str = 0
temp_str = 0

; Get sky (4)
sky_str = readin_apo(nfiles = max_idx+1, hdu=4)
temp_str = sky_str[lit_idx]

sky = temp_str.output

;delete 
sky_str = 0
temp_str = 0

; Get sky error (5)
esky_str = readin_apo(nfiles = max_idx+1, hdu=5)
temp_str = esky_str[lit_idx]

esky = temp_str.output
 
;delete
esky_str = 0
temp_str = 0

; Get tell (6)
; Get tell_err (7)


save, /variables, filename = save_file_name

skip_load_mark:
if skip_load then restore, save_file_name

aspcap = aspcap_str[lit_idx]
apg_sg = apg_sg_str[lit_idx]


wl = spec_info[0].wl_grid
ddwl = wl[1:-2]

wl_over = interpol(wl, xvec, xxvec)

specall = dindgen(npix, nlit)

medspec = median(spec_info.spectra, dimension=2)

maskall = lindgen(npix, nlit)
final_struct = []

for i = 0, n_elements(lit_idx)-1 do begin

    spec = spec_info[i].spectra
    pixmask = mask_arr[*,i]

    dev = mask_by_deriv(spec)
    nmask = mask_by_deriv2(spec)
    masked_spec = mask_by_deriv3(spec)

    pixmask_spec = mask_by_pixmask(spec,pixmask)
    ; make another vector replacing the flagged pixels with the median value
    pixmask2_spec = mask_by_pixmask2(spec,pixmask,medspec)

    spec_over = interpol(pixmask2_spec, wl, wl_over)

    
    print, "SNR: " + strtrim(aspcap[i].snr,2)
    print, "Dev: " + strtrim(dev,2)
    print, "Nmask: " + strtrim(nmask,2)

    ; Smooth spectrum with lsf
    ; make lsf
    lsf_coeffs = *lsf_coeff[i]

    lsf = lsf_gh(lx, lsf_0, lsf_coeffs)/oversamp

    convol_spec = convol(spec_over, lsf)
    
    final_spec = downsample_tophat(convol_spec, oversamp)

    specall[*,i] = final_spec

    ; Save mask
    indiv_mask = make_pixmask(pixmask)
    maskall[*,i] = indiv_mask




    ; Create each structure
    ; get the structure from master_info
    minfo = apg_sg[i]
    ; get the ASPCAP structure, and its teff and fe_h
    aspinfo = aspcap[i]
    teff_asp = aspinfo.teff
    feh_asp = aspinfo.fe_h
    vsini_asp = aspinfo.vsini
    ; get the spectrum, etc.
    
    err = spec_info[i].errors
    lsf = *lsf_coeff[i]

    newstruct = create_struct(minfo, 'TEFF_ASPCAP', teff_asp, 'FEH_ASPCAP', feh_asp, 'VSINI_ASPCAP', vsini_asp, 'SPEC', final_spec, 'ESPEC', err, 'MASK', indiv_mask, 'LSF_CO', lsf)

    final_struct = [final_struct, newstruct]

endfor

;stop
mwrfits, final_struct, '/home/stgilhool/APOGEE/overlap_table2.fits', /create

stop

end
