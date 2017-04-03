; Create a framework for choosing which spectra are bad by eye
pro apogee_vet

skip_load = 1
napo = 1350

; Master Info Table
mi_file = '/home/stgilhool/APOGEE/master_info6.fits'
mi_str = mrdfits(mi_file, 1)

; ASPCAP parameters
aspcap_str = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits',1)

; Readin APG spec
if skip_load then goto, skip_load_mark
; Read in the apg headers to get galactic latitude
apg_files = readin_apo(nfiles = napo, hdu=3)
apg_files_flux = readin_apo(nfiles = napo, /flux)

asp_files = readin_aspcap(nfiles = napo, hdu=4)

save, /variables, filename = 'apogee_data_all.sav'

skip_load_mark:
if skip_load then restore, 'apogee_data_all.sav'


; Loop through some spectra and look at relevant shit
window, 0, xs = 1600, ys=900
; for i = 0, 10 do begin

;     plot, asp_files[i].spec_model, xr=[400,2900], yr=[0,1.1], /ys, /xs
;     oplot, asp_files[i].spectra, color=200

;     dxstop

; endfor



;;; Loop through all spectra for vetting

; First, we'll look at the star-level bitmasks for bad stars
starflag = aspcap_str.starflag
andflag = aspcap_str.andflag
aspcapflag = aspcap_str.aspcapflag
aspcap_class = aspcap_str.aspcap_class
; Get the mask I care about for APG's STARFLAG
apg_starflag_bits = [0,3,4]
apg_starflag_vec = long(total(2L^apg_starflag_bits))
starflag_bad = starflag and apg_starflag_vec
; Get the mask I care about for ASP's ASPCAPFLAG
asp_starflag_bits = [23]
asp_starflag_vec = long(total(2L^asp_starflag_bits))
aspcapflag_bad = aspcapflag and asp_starflag_vec

; For the purposes of looking through subsets of the data, define
; which spectra to look at
prev_flag_vec = mrdfits('apogee_vet_flag.fits',1)
tot_flag_vec = prev_flag_vec + aspcapflag_bad
sg_only = where(tot_flag_vec eq 1)
asp_only = where(tot_flag_vec eq asp_starflag_vec, nasponly)
both = where(tot_flag_vec eq asp_starflag_vec+1)
; for now, look through the set we both flagged
orflag_vec = where(starflag ne 0, napgflag)
print, napgflag

;Read in the list of ASP BAD ONLY, which I re-flagged
aspokay_temp = mrdfits('apogee_vet_flag.fits',4)

asp_agree = where(aspokay_temp eq 1, naspagree)
asp_overrule = where(aspokay_temp eq 0, naspoverrule)

;speclist = sg_only
;speclist = asp_only
speclist = orflag_vec


; Now loop through
!p.multi = [0,1,3]
;eyeflag = lonarr(napo)
eyeflag = replicate(!values.f_nan, napo)

;i = 0
;while i lt napo do begin
foreach i, speclist, speclist_idx do begin

            
    ; Don't bother if it's a bad star
    ; if (starflag_bad[i] ne 0) or (aspcapflag_bad[i] ne 0) then begin
;         print, "Skipping spectrum number " + strtrim(i,2)
;         continue
;     endif
    
    ; See what is responsible for flagging as bad, if any
    in_sg_null = where(sg_only eq i, in_sg)
    in_both_null = where(both eq i, in_both)
    in_agree_null = where(asp_agree eq i, in_agree)
    in_overrule_null = where(asp_overrule eq i, in_overrule)

    flagtype_vec = [in_sg,in_both,in_agree,in_overrule]
    if total(flagtype_vec) gt 1 then message, "More than one flag type"
    
    if total(flagtype_vec) eq 0 then flagtype = 'None' else $
      if in_sg then flagtype = 'SG' else $
      if in_both then flagtype = 'BOTH' else $
      if in_agree then flagtype = 'AGREE' else $
      if in_overrule then flagtype = 'OVERRULE'

    ; Get the APG headers
    apghead_ptr = apg_files[i].header
    apghead = *(apghead_ptr)
    ; Get the pixel ranges for each chip
    bmin = fxpar(apghead, 'BMIN')
    bmax = fxpar(apghead, 'BMAX')
    gmin = fxpar(apghead, 'GMIN')
    gmax = fxpar(apghead, 'GMAX')
    rmin = fxpar(apghead, 'RMIN')
    rmax = fxpar(apghead, 'RMAX')
    
    ; Plot
    plot, asp_files[i].spec_model, xr=[bmin,bmax], yr=[0,1.1], /ys, /xs, $
      title = "File number: "+strtrim(i,2)+" | FlagType: "+flagtype, charsize = 1.5
    oplot, asp_files[i].spectra, color=200
    
    plot, asp_files[i].spec_model, xr=[gmin,gmax], yr=[0,1.1], /ys, /xs, $
      title = aspcap_str[i].starflags+' | '+aspcap_str[i].andflags, charsize = 1.5
    oplot, asp_files[i].spectra, color=200
    
    plot, asp_files[i].spec_model, xr=[rmin,rmax], yr=[0,1.1], /ys, /xs, $
      title = aspcap_str[i].aspcapflags, charsize = 1.5
    oplot, asp_files[i].spectra, color=200

    keep = ''
    read, keep, prompt="Keep? y/n/b: "
    if keep eq 'b' then begin
        ;i = i-1
        if speclist_idx ne 0 then speclist_idx = speclist_idx - 2
    endif else if keep eq 'y' then begin
        eyeflag[i] = 0 
        ;i = i+1
        
    endif else begin
        eyeflag[i] = 1
        ;i = i+1
    endelse
    
;endwhile
endforeach

mwrfits, eyeflag, 'apogee_vet_flag.fits'


stop
end
