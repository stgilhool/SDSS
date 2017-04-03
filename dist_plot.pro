; Make a plot of M_k vs V-K to compare with Cullen's

pro dist_plot

; Read in all the DTM parallax entries
info_fil = '/home/stgilhool/APOGEE/master_info.fits'

info_str = mrdfits(info_fil,1)

dtm_idx = where(finite(info_str.plx_dtm) eq 1, dtmcnt)

; Apparently, they all have Vmag's from simbad, somehow
info = info_str[dtm_idx]
ra = info.ra
dec = info.dec
vmag = info.vmag_sbd
vmag_est = info.vmag_est_lspm
plx_dtm = info.plx_dtm
eplx_dtm = info.eplx_dtm
kmag = info.kmag

; Calculate absolute magnitude
mag = kmag
plx = plx_dtm/1000d0 ;miliarcseconds -> arcseconds
eplx = eplx_dtm/1000d0

abs_mag = mag + 5d0 + 5d0*alog10(plx)
;eabs_mag = 5d0*alog10(eplx)
eabs_mag = 5d0*0.434*(eplx/plx)

; Query vizier to get vmags
n_spec = dtmcnt
vmag_apass = dblarr(n_spec)
evmag_apass = dblarr(n_spec)
for spec = 0, n_spec-1 do begin
    ktemp = kmag[spec]
    magcut = strtrim(2d0+ktemp,2)
    constraint = 'Vmag>='+magcut
    apass_arr = queryvizier('II/336/apass9', [ra[spec],dec[spec]], 1d0/10d0, $
                            /allcolumns)
    if n_elements(apass_arr) gt 1 then begin
        print, 'multiple matches for spec number: ' + strtrim(spec,2)
        match_idx = where(apass_arr._r eq min(apass_arr._r))
        ;match_idx = 0
        apass = apass_arr[match_idx]
        vmag_apass[spec] = apass.g_mag
        evmag_apass[spec] = apass.e_g_mag
    endif else if n_elements(apass_arr) eq 1 then begin
        if size(apass_arr, /type) ne 3 then begin
            print, 'Exactly 1 match for spec number: ' + strtrim(spec,2)
            vmag_apass[spec] = apass_arr.g_mag
            evmag_apass[spec] = apass_arr.e_g_mag
        endif else begin
            print, 'no matches for spec number: ' + strtrim(spec,2)
            vmag_apass[spec] = !values.d_nan
            evmag_apass[spec] = !values.d_nan
        endelse
    endif
endfor

; V-K
v_k = vmag_apass - kmag

; Get rid of nans
fin_idx = where(finite(vmag_apass) eq 1, fincnt)
ra = ra[fin_idx]
dec = dec[fin_idx]
abs_mag = abs_mag[fin_idx]
v_k = v_k[fin_idx]
eabs_mag = eabs_mag[fin_idx]
plx_out = plx[fin_idx]

; Plot
window, 0, xs=1500, ys=900
;plot, v_k, abs_mag, ps=6, yr = [max(abs_mag)+1, min(abs_mag)-1], xr =
;[3,9], $
plot, v_k, abs_mag, ps=6, yr = [12,3], xr = [3,9], $
  xtit = 'V-K', ytit = 'M_k', /xs, /ys
errplot, v_k, abs_mag+eabs_mag, abs_mag-eabs_mag

; Compare with Cullen
cinfo = mrdfits('/home/stgilhool/APOGEE/cullenparallax.fits',1)
cra = cinfo.ra[fin_idx]
cde = cinfo.dec[fin_idx]
cpa = cinfo.parallax[fin_idx]
cvmag = cinfo.vmag[fin_idx]
cmk = cinfo.mk[fin_idx]
svmag = vmag_apass[fin_idx]
for i = 0, fincnt - 1 do print, cra[i], ra[i], cde[i], dec[i], cvmag[i], svmag[i]

stop
openw, lun, 'vklist.txt', /get_lun
for i = 0, fincnt-1 do printf, lun, ra[i], dec[i], plx_out[i], v_k[i], abs_mag[i]

free_lun, lun
stop


end
