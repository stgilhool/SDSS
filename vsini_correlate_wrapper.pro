pro vsini_correlate_wrapper, xmin, xmax, outstruct, files=files



display = 0
ps=0

;load rainbow color table with discrete color tags	      
loadct,39,/silent
setcolors,/system_variables,/silent,decomposed=0  ;requires setcolors.pro

;plotting parameters to set axes and fonts
!p.background=0    ;default background color
!p.color=255
!p.charsize=1.7		;text default size
!p.charthick=6		;text default thickness
!x.thick=5		;thicken x-axis 
!y.thick=5		;thicken y-axis
!p.font=1		;set default font
!p.thick=5		;set default plotting line thickness
!x.margin=[7,4]

;set psym=8 to circle by default
circind=findgen(40) * (!pi*2/39.)
defsysv,'!circ',transpose( [[cos(circind)],[sin(circind)]])  ;user symbol vertices
usersym,!circ(0,*),!circ(1,*),/fill  ;circle plot


bigc = 3d5
d_logwl = 6d-6
logwl0 = 4.179d0

if n_elements(files) eq 0 then begin
;files =
;[840,431,793,326,953,391,286,1076,1306,289,740,1142,60,43,622,479,523,867,520]
    files = [431, 793,326,397,289,1142,60,622]
;files = [289]
;teff_vec = [2600,2800,2800,2800,2900,3000,3000,3000,3100,3100,3100,3200,3200,3200,3300,3300,3400,3400,3400]
    ;vsini_known_vec =
    ;[19,14.5,56.5,32.25,10.75,11.25,9,12.75,15.25,25,12,26.25,7.25,16.5,4,4,4,4,4]
endif


mi = mrdfits('/home/stgilhool/APOGEE/master_info6.fits',1)
teff_vec = mi[files].teff_vfit
feh_vec = mi[files].feh_vfit
vsini_known_vec = mi[files].vsini_vfit
apg_info = readin_apo2(files=files, hdu=8)
;lsf_info = *apg_info[files].output
;dxstop, 'lsf_info'
nfiles = n_elements(files) 

; xmin=336L
; xmax=1500L
; xmin=1500L
; xmax=3100L

; Make wl grid
xpixels = lindgen(xmax-xmin+1)+xmin
logwl_grid = poly(xpixels,[logwl0,d_logwl])
wl_grid = 10d0^logwl_grid

;Oversample
oversamp = 11L
npix = n_elements(wl_grid) 
npix_over = npix*oversamp
x = dindgen(npix)
xx = (dindgen(npix_over)-(oversamp/2))/oversamp
wl_grid_over = interpol(wl_grid, x, xx)

xnorm = x-npix/2
xxnorm = xx-npix_over/(2*oversamp)

lag_over = lindgen(npix_over) - npix_over/2
vel_x_over = (10d0^(xxnorm*d_logwl)-1d0)*bigc
deltav = ((10d0^(d_logwl/oversamp))-1d0)*bigc

epsilon = 0.6d0
k1 = 0.60975d0 + 0.0639d0*epsilon + 0.0205*epsilon^2 + 0.021*epsilon^3

; for later use in making result go to zero quickly        
velx_zero_idx = where(abs(vel_x_over) ge 100,nfft,complement=mid_idx)
velx_decay_vec = abs(vel_x_over)-100
velx_decay_vec[mid_idx] = 0d0
velx_decay_vec = velx_decay_vec/50d0
decay_vec = exp(-1d0*velx_decay_vec)

; for later use.  the abscissa vector for the fft
sigvec_x = dindgen((npix_over - 1)/2) + 1
is_N_even = (npix_over MOD 2) EQ 0
if (is_N_even) then $
  sigvec = [0d0, sigvec_x, npix_over/2, -npix_over/2 + sigvec_x]/(npix_over*deltav) $
else $
  sigvec = [0d0, sigvec_x, -(npix_over/2 + 1) + sigvec_x]/(npix_over*deltav)

vsini_sigvec = k1/(sigvec)

vsini_recovered = []
outstruct = ptrarr(n_elements(files), /allocate_heap)

foreach apofile, files, fil_idx do begin

     
    ; Make APOGEE LSF
    ;lsf_sigma = 1.34d0
    npix_lsf = 35L * oversamp
    lsfx = (dindgen(npix_lsf)-(npix_lsf/2))/oversamp
    lsfco = *apg_info[fil_idx].output
    lsf0 = (xmax-xmin)/2 + xmin
    lsf = lsf_gh(lsfx+lsf0, lsf0, lsfco)
    lsf = lsf/total(lsf,/double)

    ; Data
    apo_spec = apg_info[fil_idx].spectra[xmin:xmax]
    apo_err = apg_info[fil_idx].errors[xmin:xmax]
    apo_spec_over_raw = interpol(apo_spec, x, xx)
    
    apo_mask_str = readin_apo2(files=[files[fil_idx]], hdu=3)
    apo_mask_arr = *apo_mask_str.output
    apo_mask = apo_mask_arr[xmin:xmax,0]
    ; Change mask to just the bits we want
    flag_bit = [0,1,2,3,4,5,6,12,13,14]
    flag_big = indgen(17)
    flag_decimal = long(total(2L^flag_bit))

    mask = flag_decimal and apo_mask
    mask_idx = where(mask gt 0, nmask)
    mask_cp = mask
    if nmask gt 0 then mask[mask_idx] = 1
    mask = byte(mask)

    ; Model parameters
    teff = teff_vec[fil_idx]
    vsini_guess = vsini_known_vec[fil_idx]
    ;feh = feh_vec[fil_idx]
    feh='-0.0'
    
    ; Load the PHOENIX template
    pho_spec_str = readin_phoenix([teff],[feh],['4.5'],wl_grid)
    pho_spec = pho_spec_str.spec
    pho_over_str = readin_phoenix([teff],[feh],['4.5'],wl_grid_over)
    pho_over = pho_over_str.spec

    vsini_vector = vsini_correlate(lsf=lsf, $
                                   display=display, $
                                   pho_spec=pho_spec, $
                                   pho_over=pho_over, $
                                   oversamp=oversamp, $
                                   apo_spec=apo_spec, $
                                   apo_err=apo_err, $
                                   mask=mask, $
                                   pix_x=x, $
                                   xx=xx, $
                                   vel_x_over=vel_x_over, $
                                   wl_grid=wl_grid, $
                                   lag_over=lag_over, $
                                   decay_vec=decay_vec, $
                                   sigvec=sigvec, $
                                   vsini_vfit=vsini_guess)

    sigpeak_str = vsini_vector

    vsini_vector = sigpeak_str.minima_vsini
    
    vs_tickvec = [60,30,25,20,15,10,8,5]	
    vs_tickval = (k1/vs_tickvec)*1e3

    plot, sigvec*1d3, sigvec*1d3, xr=[0,150], yr=[-15,-2], xs=8, /nodata
    oplot, sigpeak_str.sigvec*1d3, sigpeak_str.fftvec
    axis, xaxis=1, xticks=n_elements(vs_tickvec), $
      xtickn=string(vs_tickvec), xtickv=vs_tickval, charsize=1.5
    for i = 0, n_elements(sigpeak_str.minima_idx)-1 do begin
        oplot, replicate(sigpeak_str.minima_sigval[i]*1d3,2), [-100,100], linest=2, thick=1
    endfor
    
    ;read, mindx, prompt="Which one is the real minimum?"
    vsini_final = sigpeak_str.minima_vsini[0:2]

    ;vsini_final = sigpeak_str.minima_vsini[mindx]

    print, fil_idx
    print, "VSINI_VFIT = "+strtrim(vsini_guess,2)
    print, "VSINI_RECOVERED = "+strtrim(vsini_vector[0],2)
    print, vsini_vector

    ;vsini_recovered = [vsini_recovered, vsini_vector[0]]
    vsini_recovered = [[vsini_recovered], [vsini_final]]

    *outstruct[fil_idx] = sigpeak_str

endforeach

;plot, vsini_known_vec, vsini_recovered, ps=6
plot, [vsini_known_vec], [vsini_known_vec], /nodata
for i = 0, n_elements(vsini_known_vec) -1 do begin & $
    oplot, replicate(vsini_known_vec[i],3), vsini_recovered[*,i], ps=-8, thick=1 & $
endfor
oplot, lindgen(50), lindgen(50), linest=2

;stop

end
                                   
