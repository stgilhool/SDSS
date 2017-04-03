; Basically a one-time use procedure to write out useful
; APG data and spectra to one fits file
pro write_apg_data

nx = 8575L

root_dir = '/home/stgilhool/APOGEE/'
data_path = root_dir + 'APOGEE_data/'

listfile = data_path + 'mdwarflist.txt'

; Readin file names
readcol, listfile, datafile, format="A", count=nfiles

files = data_path + datafile

; Initialize
x = dindgen(nx)
final_struct = []

for i = 0, nfiles-1 do begin
    ; Read in header for a file
    null = mrdfits(files[i],0,head)
    
    ; Get info from header
    objid  = fxpar(head, 'OBJID')
    ra     = fxpar(head, 'RA')
    dec    = fxpar(head, 'DEC')
    nvisits= fxpar(head, 'NVISITS')
    wl0    = fxpar(head, 'CRVAL1', datatype=0.0d)
    wl1    = fxpar(head, 'CDELT1', datatype=0.0d)
    snr    = fxpar(head, 'SNR')
    starflg= fxpar(head, 'STARFLAG')
    andflg = fxpar(head, 'ANDFLAG')

    bmin = fxpar(head, 'BOVERMIN')
    bmax = fxpar(head, 'BOVERMAX')
    gmin = fxpar(head, 'GOVERMIN')
    gmax = fxpar(head, 'GOVERMAX')
    rmin = fxpar(head, 'ROVERMIN')
    rmax = fxpar(head, 'ROVERMAX')

    ; Read in most data for a file
    fluxall = mrdfits(files[i],1)
    errall = mrdfits(files[i], 2)
    maskall = mrdfits(files[i], 3)
    skyall = mrdfits(files[i], 4)
    eskyall = mrdfits(files[i], 5)
    tellall = mrdfits(files[i], 6)
    etellall = mrdfits(files[i], 7)

    ; Save only the two combined spectra (or the only spec if 1 visit)
    if nvisits eq 1 then begin
        
        flux_pixwt = fluxall[*,0]
        flux_glbwt = fluxall[*,0]
        err_pixwt = errall[*,0]
        err_glbwt = errall[*,0]
        mask_pixwt = maskall[*,0]
        mask_glbwt = maskall[*,0]
        sky_pixwt = skyall[*,0]
        sky_glbwt = skyall[*,0]
        esky_pixwt = eskyall[*,0]
        esky_glbwt = eskyall[*,0]
        tell_pixwt = tellall[*,0]
        tell_glbwt = tellall[*,0]
        etell_pixwt = etellall[*,0]
        etell_glbwt = etellall[*,0]
    endif else if nvisits gt 1 then begin
        
        flux_pixwt = fluxall[*,0]
        flux_glbwt = fluxall[*,1]
        eflux_pixwt = errall[*,0]
        eflux_glbwt = errall[*,1]
        mask_pixwt = maskall[*,0]
        mask_glbwt = maskall[*,1]
        sky_pixwt = skyall[*,0]
        sky_glbwt = skyall[*,1]
        esky_pixwt = eskyall[*,0]
        esky_glbwt = eskyall[*,1]
        tell_pixwt = tellall[*,0]
        tell_glbwt = tellall[*,1]
        etell_pixwt = etellall[*,0]
        etell_glbwt = etellall[*,1]
    endif else message, "Nvisits problem da yo!"

    ; normalize fluxes
    ; Remove continuum
    nblue  = bmax - bmin+1
    ngreen = gmax - gmin+1
    nred   = rmax - rmin+1
    
    bnorm_pixwt = continuum_fit(dindgen(nblue), flux_pixwt[bmin:bmax])
    gnorm_pixwt = continuum_fit(dindgen(ngreen), flux_pixwt[gmin:gmax])
    rnorm_pixwt = continuum_fit(dindgen(nred), flux_pixwt[rmin:rmax])

    bnorm_glbwt = continuum_fit(dindgen(nblue), flux_glbwt[bmin:bmax])
    gnorm_glbwt = continuum_fit(dindgen(ngreen), flux_glbwt[gmin:gmax])
    rnorm_glbwt = continuum_fit(dindgen(nred), flux_glbwt[rmin:rmax])

    ; Normalize and scale errors for each chip
    spec_pixwt = flux_pixwt
    spec_pixwt[bmin:bmax] = flux_pixwt[bmin:bmax]/bnorm_pixwt
    spec_pixwt[gmin:gmax] = flux_pixwt[gmin:gmax]/gnorm_pixwt
    spec_pixwt[rmin:rmax] = flux_pixwt[rmin:rmax]/rnorm_pixwt
    
    spec_glbwt = flux_glbwt
    spec_glbwt[bmin:bmax] = flux_glbwt[bmin:bmax]/bnorm_glbwt
    spec_glbwt[gmin:gmax] = flux_glbwt[gmin:gmax]/gnorm_glbwt
    spec_glbwt[rmin:rmax] = flux_glbwt[rmin:rmax]/rnorm_glbwt

    espec_pixwt = eflux_pixwt
    espec_pixwt[bmin:bmax] = eflux_pixwt[bmin:bmax]/bnorm_pixwt
    espec_pixwt[gmin:gmax] = eflux_pixwt[gmin:gmax]/gnorm_pixwt
    espec_pixwt[rmin:rmax] = eflux_pixwt[rmin:rmax]/rnorm_pixwt

    espec_glbwt = eflux_glbwt
    espec_glbwt[bmin:bmax] = eflux_glbwt[bmin:bmax]/bnorm_glbwt	
    espec_glbwt[gmin:gmax] = eflux_glbwt[gmin:gmax]/gnorm_glbwt
    espec_glbwt[rmin:rmax] = eflux_glbwt[rmin:rmax]/rnorm_glbwt

    ; Zero the regions in between
    spec_pixwt[0:bmin-1] = 0d0
    spec_pixwt[bmax+1:gmin-1] = 0d0
    spec_pixwt[gmax+1:rmin-1] = 0d0
    spec_pixwt[rmax+1:-1] = 0d0
    
    espec_pixwt[0:bmin-1] = 0d0
    espec_pixwt[bmax+1:gmin-1] = 0d0
    espec_pixwt[gmax+1:rmin-1] = 0d0
    espec_pixwt[rmax+1:-1] = 0d0

    spec_glbwt[0:bmin-1] = 0d0	
    spec_glbwt[bmax+1:gmin-1] = 0d0
    spec_glbwt[gmax+1:rmin-1] = 0d0
    spec_glbwt[rmax+1:-1] = 0d0
    
    espec_glbwt[0:bmin-1] = 0d0
    espec_glbwt[bmax+1:gmin-1] = 0d0
    espec_glbwt[gmax+1:rmin-1] = 0d0
    espec_glbwt[rmax+1:-1] = 0d0

; Return requested info
    outstr = {OBJID:objid,$
              RA:ra,$
              DEC:dec,$
              NVISITS:nvisits,$
              SNR:snr,$
              STARFLAG:starflg,$
              ANDFLAG:andflg, $
              BRANGE:[bmin,bmax],$
              GRANGE:[gmin,gmax],$
              RRANGE:[rmin,rmax],$
              WL_COEFFS:[wl0,wl1],$
              FLUX_PIXWT:flux_pixwt,$
              FLUX_GLBWT:flux_glbwt,$
              E_FLUX_PIXWT:eflux_pixwt,$
              E_FLUX_GLBWT:eflux_glbwt,$
              ;SPEC_PIXWT:spec_pixwt,$
              ;SPEC_GLBWT:spec_glbwt,$
              ;E_SPEC_PIXWT:espec_pixwt,$
              ;E_SPEC_GLBWT:espec_glbwt,$
              MASK_PIXWT:mask_pixwt,$
              MASK_GLBWT:mask_glbwt};,$
              ;SKY_PIXWT:sky_pixwt,$
              ;SKY_GLBWT:sky_glbwt,$
              ;E_SKY_PIXWT:esky_pixwt,$
              ;E_SKY_GLBWT:esky_glbwt,$
              ;TELL_PIXWT:tell_pixwt,$
              ;TELL_GLBWT:tell_glbwt,$
              ;E_TELL_PIXWT:etell_pixwt,$
              ;E_TELL_GLBWT:etell_glbwt $
             ;}
    
    final_struct = [final_struct,outstr]
        
endfor
    
mwrfits, final_struct, '/home/stgilhool/APOGEE/APOGEE_data/apg_flux.fits', /create
;mwrfits, final_struct, '/home/stgilhool/APOGEE/APOGEE_data/apg_data.fits', /create

end
