; Make black body plots for T=2500-4500
; 0.5-27 microns
pro bbody

max_wl = 2.7d5 ;27 microns in angstroms
min_wl = 5d3   ;0.5 microns in angstroms

wl_step_int = 1d0

n_steps_int = (max_wl-min_wl)/wl_step_int

wl_vec_int = dindgen(n_steps_int)*wl_step_int + min_wl

;window, 0, xs=1700, ys=1000
;!p.multi = [0,7,3]

teff = dindgen(21)*100 + 2500d0

int = []
flux = []

set_plot, 'ps'
device, filename = 'bblog.ps'
device, /color, bits=8
device, xs = 13, ys= 8, /inches
loadct, 13


for i = 0, 20 do begin
    
    int_i = planck(wl_vec_int, teff[i])
    
    int = [[int],[int_i]]

    ; Get the flux for each bin and sum up to 0.5um steps
    flux_i_hr = int_i * wl_step_int

    
    bin_size = 1d3 * wl_step_int

    n_bins = (max_wl-min_wl)/bin_size
    
    wl_vec = (dindgen(n_bins)*bin_size + (bin_size/2) +min_wl)/1d4 ; make it in microns
    
    
    flux_i_lr = dindgen(n_bins)
    for k = 0, n_bins-1 do begin
        
        flux_i_lr[k] = total(flux_i_hr[k*bin_size:(k+1)*bin_size - 1], /double)

    endfor
    
    flux = [[flux],[flux_i_lr]]

    if i eq 0 then plot, wl_vec, flux_i_lr, /ylog, xtit = "Wavelength(um)", $
      title = "Blackbody curves for Teff = 2500K-4500K in steps of 100K", $
      ytit = "Flux (erg/cm^2/s/um)", charsize=1.5, yr = [10d^4,10d^10] $
      else oplot, wl_vec, flux_i_lr, color = i*10

endfor



help, bb
help, n_bins
help, bin_size
help, flux_i_hr
help, flux_i_lr
help, flux
stop

end

    
