; Procedure to estimate the age of stars based on Kumamoto, et al. 2017

pro stellar_age_estimator, output_structure

; Read in Vsini, Teff, and U,V,W (although we're just using W)
info_str = mrdfits('spacevelocity_wflags.fits',1)
uarr = info_str.u
uarr = -1d0*uarr ;Flip the sign due to silly convention
varr = info_str.v
warr = info_str.w
vsini = info_str.vsini
teff = info_str.teff
p_membership = info_str.thickthin_prob

; Read in gold-sample flags
flags = mrdfits('/home/stgilhool/APOGEE/APOGEE_data/gold_sample.fits',0)
; Save flag indices
;g_idx = where(flags eq 0, ngood, complement=b_idx, ncomplement=nbad)

napo = n_elements(flags) 
ngood = napo

; Make Age vector
max_age = 10 ;Gyr
min_age = 0.1 ;Gyr
d_age = 0.1
nsteps_age = round(((max_age-min_age)/d_age)) + 1

age = dindgen(nsteps_age) * d_age + min_age
; help, age
; print, minmax(age)
; dxstop

; Make W abscissa vector
max_w = floor(max(warr)+ 4*abs(max(warr)))
min_w = floor(min(warr)- 4*abs(min(warr)))
extremum = abs(max_w) > abs(min_w)
d_w = 1d0 ;km/s

nsteps_w = floor(2*extremum/d_w)
;help, nsteps_w
; Enforce oddness
nsteps_w = (nsteps_w mod 2 eq 0) ? nsteps_w+1 : nsteps_w
; help, nsteps_w
; dxstop

wvec = (dindgen(nsteps_w) - (long(nsteps_w)/2)) * d_w
; help, wvec
; print, minmax(wvec)
; dxstop

; Calculate z_dispersion at all times
z_disp = 7.36d0*age^0.59 ;km/s, with age in Gyr
u_disp = 19.4d0*age^0.37 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Next, we want to draw a bunch of samples from each age
; For a given star, we look at its W (vertical) velocity and ask at
; each stellar age, "What is the probability this star was drawn from
; this sample?"

; METHOD 1
; So the velocity dispersion ends up being the RMS of W, which is also
; the sigma for a Gaussian with that RMS.  I think this also means
; that the most likely age for a star will correspond to the age where
; the velocity dispersion matches the velocity.
; SO,
; 1) We generate a normalized gaussian at each age.  This is then the
; probability as a function of W velocity
; 2) We get that probability for each time step for each star
; 3) Look at that distribution and see what's there

; 1) We generate a normalized gaussian at each age.  This is then the
; probability as a function of W velocity

; MODIFIED PLAN
; 1) We skip generating the Gaussians, and just grab the probabilities
; for all stars at each time step by calculating the value of the
; Gaussians at the stars' velocities

prob_arr_z = []
prob_arr_r = []

foreach zdispersion, z_disp, disp_idx do begin

    
    rdispersion = u_disp[disp_idx]
    ; Make a gaussian
    ;pdf = gaussian_function(dispersion, width = nsteps_w, /normalize)
    
    ; Nevermind, let's directly calculate the probabilities by
    ; calculating the value of the Gaussian at the W velocity values
    ; for a normalized Gaussian with a given dispersion
    prob_vec_z = gaussian(warr,[1d0/(sqrt(2d0*!dpi)*zdispersion), 0.0, zdispersion])
    prob_vec_r = gaussian(uarr,[1d0/(sqrt(2d0*!dpi)*rdispersion), 0.0, rdispersion])

    ; Store the probabilties
    prob_arr_z = [[prob_arr_z],[prob_vec_z]]
    prob_arr_r = [[prob_arr_r],[prob_vec_r]]
    
    
endforeach

; 2) Each column should now represent the probabilities for a given
; star at each time step.  Let's look at them

!p.multi = [0,1,3]

prob_dist_arr = dblarr(nsteps_age, ngood)

foreach wvel, warr, w_idx do begin

    uvel = uarr[w_idx]

    ; Calculate the Z, R and combined probabilities
    prob_star_z = reform(prob_arr_z[w_idx,*])
    norm_factor_z = total(prob_star_z, /double)
    prob_star_norm_z = prob_star_z/norm_factor_z

    prob_star_r = reform(prob_arr_r[w_idx,*])
    norm_factor_r = total(prob_star_r, /double)
    prob_star_norm_r = prob_star_r/norm_factor_r

    combined_prob = prob_star_norm_z*prob_star_norm_r/total(prob_star_norm_z*prob_star_norm_r,/double)
    
    
    ; Save the star's pdf
    prob_dist_arr[*,w_idx] = combined_prob

    ; if w_idx ge 463 then begin
;         help, uvel
;         print, uvel
;         print, ''
     
    

;     ; Plot the distribution
;      plot, age, prob_star_norm_z, title = "W_LSR = "+strtrim(wvel,2)+ $
;        " | p(TD)/p(thin) = "+strtrim(p_membership[w_idx],2), charsize = 1.5, $
;        xtit = "Age (Gyr)", ytit = "Probability"

;      plot, age, prob_star_norm_r, title = "U_LSR = "+strtrim(uvel,2)+ $
;        " | p(TD)/p(thin) = "+strtrim(p_membership[w_idx],2), charsize = 1.5, $
;        xtit = "Age (Gyr)", ytit = "Probability"

;      plot, age, combined_prob, title = "Combined Prob", charsize = 1.5, $
;        xtit = "Age (Gyr)", ytit = "Probability"
    
;      stop
     
;      endif
        

endforeach

!p.multi=0
    
output_structure = {age:age, $
                    prob:prob_dist_arr}

end
