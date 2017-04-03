; plot M_K vs V-K for different vmags
pro show_vk_mk

; Readin
infile = '/home/stgilhool/APOGEE/master_info4.fits'
info = mrdfits(infile,1)

; Want to do M_K vs V-K for:
; VMAG_APASS
; VMAG_LSPM
; VMAG_EST_LSPM
; VMAG_SBD
; G_MAG (SDSS)
; VTMAG (TYCHO)

; Indices where we have both the V mag and the M_K
apass2_idx = where(finite(info.plx_dtm) and finite(info.vmag_apass), apass2cnt)
lspm2_idx = where(finite(info.plx_dtm) and finite(info.vmag_lspm), lspm2cnt)
lspm_est2_idx = where(finite(info.plx_dtm) and finite(info.vmag_est_lspm), lspmest2cnt)
sbd2_idx = where(finite(info.plx_dtm) and finite(info.vmag_sbd), sbd2cnt)
sdss2_idx = where(finite(info.plx_dtm) and finite(info.g_mag), sdss2cnt)
tycho2_idx = where(finite(info.plx_dtm) and finite(info.vtmag), tycho2cnt)
; Indices for just the V mags of each type
apass_idx = where(finite(info.vmag_apass), apasscnt)
lspm_idx = where(finite(info.vmag_lspm), lspmcnt)
lspm_est_idx = where(finite(info.vmag_est_lspm), lspmestcnt)
sbd_idx = where(finite(info.vmag_sbd), sbdcnt)
sdss_idx = where(finite(info.g_mag), sdsscnt)
tycho_idx = where(finite(info.vtmag), tychocnt)


; Get relevant vectors for all spectra
plx_vec = info.plx_dtm/1d3
e_plx_vec = info.eplx_dtm/1d3
k_vec = info.kmag
; Get their absolute k magnitudes
absk_vec = k_vec + 5d0 + 5d0*alog10(plx_vec)
eabsk_vec = 5d0*0.434*(e_plx_vec/plx_vec)

; Get V for each type
vmag_apass    = info[apass_idx].vmag_apass
vmag_lspm     = info[lspm_idx].vmag_lspm
vmag_est_lspm = info[lspm_est_idx].vmag_est_lspm
vmag_sbd      = info[sbd_idx].vmag_sbd
vmag_sdss     = info[sdss_idx].g_mag
vmag_tycho    = info[tycho_idx].vtmag

; Get V corresponding to PLX for each type
vmag2_apass    = info[apass2_idx].vmag_apass
vmag2_lspm     = info[lspm2_idx].vmag_lspm
vmag2_est_lspm = info[lspm_est2_idx].vmag_est_lspm
vmag2_sbd      = info[sbd2_idx].vmag_sbd
vmag2_sdss     = info[sdss2_idx].g_mag
vmag2_tycho    = info[tycho2_idx].vtmag

; Get V-K for each type
vk_apass    = vmag_apass - k_vec[apass_idx]
vk_lspm     = vmag_lspm - k_vec[lspm_idx]
vk_est_lspm = vmag_est_lspm - k_vec[lspm_est_idx]
vk_sbd      = vmag_sbd - k_vec[sbd_idx]
vk_sdss     = vmag_sdss - k_vec[sdss_idx]
vk_tycho    = vmag_tycho - k_vec[tycho_idx]

; Get V-K corresponding to PLX entries for each type
vk2_apass    = vmag2_apass - k_vec[apass2_idx]
vk2_lspm     = vmag2_lspm - k_vec[lspm2_idx]
vk2_est_lspm = vmag2_est_lspm - k_vec[lspm_est2_idx]
vk2_sbd      = vmag2_sbd - k_vec[sbd2_idx]
vk2_sdss     = vmag2_sdss - k_vec[sdss2_idx]
vk2_tycho    = vmag2_tycho - k_vec[tycho2_idx]

; Take just the absk's and eabsk's that we need
absk_apass    = absk_vec[apass2_idx]
absk_lspm     = absk_vec[lspm2_idx]
absk_est_lspm = absk_vec[lspm_est2_idx]
absk_sbd      = absk_vec[sbd2_idx]
absk_sdss     = absk_vec[sdss2_idx]
absk_tycho    = absk_vec[tycho2_idx]

; Make composite plot
window, 0, xs = 1700, ys = 900

plot, vk2_apass, absk_apass, title = "M_K vs V-K for all Vmag types", xtitle = "V-K", ytitle = "M_K", xr = [3,9], yr = [12,3],/xs,/ys, /nodata, charsize = 1.5
al_legend, ['APASS','LSPM','LSPM_EST','SBD','SDSS','TYCHO'], psym = [1,2,4,5,6,7], charsize = 1.5
oplot, vk2_apass, absk_apass, ps = 1
oplot, vk2_lspm, absk_lspm, ps = 2
oplot, vk2_est_lspm, absk_est_lspm, ps = 4
oplot, vk2_sbd, absk_sbd, ps = 5
oplot, vk2_sdss, absk_sdss, ps = 6  
oplot, vk2_tycho, absk_tycho, ps = 7


; Make multi plot
window, 1, xs = 1700, ys = 900
!p.multi = [0,3,2]

plot, vk2_apass, absk_apass, title = "M_K vs V-K for VMAG_APASS", xtitle = "V-K", ytitle = "M_K", xr = [3,9], yr = [12,3],ps=6,charsize=2,/xs,/ys
plot, vk2_lspm, absk_lspm, title = "M_K vs V-K for VMAG_LSPM", xtitle = "V-K", ytitle = "M_K", xr = [3,9], yr = [12,3],/xs,/ys, ps=6,charsize=2
plot, vk2_est_lspm, absk_est_lspm, title = "M_K vs V-K for VMAG_EST_LSPM", xtitle = "V-K", ytitle = "M_K", xr = [3,9], yr = [12,3],/xs,/ys, ps=6,charsize=2
plot, vk2_sbd, absk_sbd, title = "M_K vs V-K for VMAG_SBD", xtitle = "V-K", ytitle = "M_K", xr = [3,9], yr = [12,3],/xs,/ys, ps=6,charsize=2
plot, vk2_sdss, absk_sdss, title = "M_K vs V-K for SDSS U_MAG", xtitle = "V-K", ytitle = "M_K", xr = [3,9], yr = [12,3],/xs,/ys, ps=6,charsize=2
plot, vk2_tycho, absk_tycho, title = "M_K vs V-K for TYCHO VTMAG", xtitle = "V-K", ytitle = "M_K", xr = [3,9], yr = [12,3],/xs,/ys, ps=6,charsize=2


stop

end
