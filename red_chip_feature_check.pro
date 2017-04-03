; Look at a particular section of the red chip to determine correlations
pro red_chip_feature_check

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

;;; Loop through all spectra for vetting
; Read-in gold sample flags
gold_sample = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/gold_sample.fits',0)
gs_idx = where(gold_sample eq 0, ngold, complement=bd_idx, ncomplement=nbad)
speclist = gs_idx

foreach i, speclist, speclist_idx do begin

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
    

    ; Smooth the region
    spec_chunk = (asp_files[i].spectra)[7340:7570]
    npix = n_elements(spec_chunk) 
    oversamp = 7
    x = dindgen(npix)
    xx = (dindgen(npix*oversamp) - (oversamp/2))/oversamp

    spec_interp = interpol(spec_chunk, x, xx)
    spec_smooth = gauss_smooth(spec_interp, 7)
    spec_smooth = gauss_smooth(spec_smooth, 7)

    ; Plot result of smoothing
    plot, xx, spec_interp
    oplot, xx, spec_smooth, color=200

    dxstop

;      Plot
;     plot, asp_files[i].spec_model, xr=[7340,7570], yr=[0,1.1], /ys, /xs, $
;       charsize = 1.5
;     oplot, asp_files[i].spectra, color=200

;     dxstop
    
;endwhile
endforeach


stop
end
