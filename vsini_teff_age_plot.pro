pro vsini_teff_age_plot

; Read in Vsini, Teff, and U,V,W (although we're just using W)
info_str = mrdfits('spacevelocity_wflags.fits',1)
uarr = info_str.u
uarr = -1d0*uarr ;Flip the sign due to silly convention
varr = info_str.v
warr = info_str.w
vsini = info_str.vsini
teff_rough = info_str.teff
p_membership = info_str.thickthin_prob

; Get a better teff
asp_str = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/allparams_dr13.fits',1)
teff_asp = asp_str.teff


; Read in gold-sample flags
flags = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/gold_sample.fits',0)
; Save flag indices
g_idx = where(flags eq 0, ngood, complement=b_idx, ncomplement=nbad)
napo = n_elements(flags) 
nstars = n_elements(g_idx) 

; Use the ASPCAP teff
teff_asp_gold = teff_asp[g_idx]
; Substitute the shitty teffs where asp in undefined
tnan_idx = where(teff_asp_gold le 0, ntnan)
teff = teff_asp_gold
teff[tnan_idx] = teff_rough[tnan_idx]


; Generate the age PDFS
stellar_age_estimator, age_str

age_vec = age_str.age
age_pdf_arr = age_str.prob

; Generate N samples for each star, based on their PDFs
nsamples = 1d5

age_realizations = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/vsini_teff_age_mcdata.fits',0)

; For each realization, bin by Age and Teff

tbinsize = 100
trangemin = 2900
trangemax = 4000

amin = 0
abinsize = 2

teff_hist = histogram(teff, binsize = tbinsize, min=trangemin, $
                      reverse_indices = ri_t)

; VALUE LOCATE Better way than histogram
tmin = min(teff)
tmax = max(teff)

tbinmin = floor(tmin/100.)*100.
tbinmax = ceil(tmax/100.)*100.

ntbins = (tbinmax-tbinmin)/100.

teff_bin_vec = lindgen(ntbins+1)*tbinsize + tbinmin ;[2900,4000]
teff_bin_ctr = lindgen(ntbins)*tbinsize + tbinmin + tbinsize/2

tbin_idx = value_locate(teff_bin_vec, teff)

nabins = 5
age_bin_vec = findgen(nabins+1)*abinsize + amin ;[0,10] 
age_bin_vec_vl = findgen(nabins)*abinsize + amin ;[0,10)
age_bin_ctr = lindgen(nabins)*abinsize + amin + abinsize/2

;tbin_idx_long = reform(rebin(tbin_idx,nstars,nsamples),nstars*nsamples,1)

abin_idx_long = []

;for trial = 0, nsamples-1 do begin
nsamples_test = 1000
for trial = 0, nsamples_test-1 do begin

    ; Get ages
    trial_age = reform(age_realizations[trial,*])
    
    ; Make histogram
    ; age_hist = histogram(trial_age, binsize = abinsize, min=amin, $
;                          reverse_indices = ri_a)

    ; VALUE LOCATE
    abin_idx = value_locate(age_bin_vec_vl, trial_age)
    
    abin_idx_long = [abin_idx_long, abin_idx]

    ;age_plot_vec = age_bin_ctr[abin_idx]
    ;teff_plot_vec = teff_bin_ctr[tbin_idx]

    ;test = plot3d(teff_plot_vec, age_plot_vec, vsini, 'o', /sym_filled, $
     ;             xr = [trangemin,trangemax])

    

endfor


window, 0, xs=1500, ys=900
set_plot, 'ps'
 device, filename = 'vsini_teff_age_plot.eps'
 device, /color, bits=8
 device, xs = 13, ys= 8, /inches
 loadct, 13



!p.multi = [0,1,5]

for age_bin = 0, nabins-1 do begin

    ; Find the indices in the whole trial where the age is in the
    ; age_bin bin
    abin_idx_idx = where(abin_idx_long eq age_bin, nage)
    ; Not checking that indices were found because they will always be
    ; there unless I change something

    ; Since the tbin_idx doesn't change by iteration, we can relate
    ; the position in the abin_idx_long to the other short vectors
    ; using the modulus (ie. if abin_idx_idx = nstars = the first index in
    ; the 2nd trial, longtoshort_idx = 0 = the first index in the short
    ; tbin_idx vector and the short vsini vector)
    longtoshort_idx = abin_idx_idx mod nstars
    tbin_idx_idx = tbin_idx[longtoshort_idx]
    ;tbin_idx_idx = tbin_idx_long[abin_idx_idx]
    teff_age_vec = teff_bin_ctr[tbin_idx_idx]
    vsini_age_vec = vsini[longtoshort_idx]

    rotfrac = dblarr(ntbins)
    rot_lims = dblarr(2,ntbins)

    for tstep = 0, ntbins-1 do begin
        
        tai_idx = where(tbin_idx_idx eq tstep, ntai)
        
        if ntai ne 0 then begin
            vsini_tai = vsini_age_vec[tai_idx]
            
            vrot_tai_idx = where(vsini_tai gt 5, nvrottai)

            rotfrac_tai = double(nvrottai)/double(ntai)

            ;rot_lims_tai = binomial_stat_errors(ntai, nvrottai, 95d0)
            
            
        endif else begin
            
            rotfrac_tai = 0
            ;rot_lims_tai = [0,1]
            
        endelse

        ; Save rotfrac and lims
        rotfrac[tstep] = rotfrac_tai
        ;rot_lims[*,tstep] = rot_lims_tai

    endfor

    help, rotfrac
    help, rotfrac_tai
    help, teff_bin_ctr
    dxstop
    

    plot, teff_bin_ctr, rotfrac, ps=6, symsize=0.5, xr = [trangemin,trangemax], $
      xtit = 'Teff', ytit = 'Vsini', tit = 'AGE : ' +strtrim(age_bin_vec[age_bin],2)+ $
      ' - ' +strtrim(age_bin_vec[age_bin+1],2)+ ' Gyrs', charsize = 1.5, charthick=2
    ;oploterror, teff_age_vec, rotfrac, rot_lims[1]-rotfrac, /hibar, ps=6
    ;oploterror, teff_age_vec, rotfrac, rotfrac-rot_lims[0], /lobar, ps=6
endfor

device, /close_file

stop
end
