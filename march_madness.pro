; Randomly fill in the march madness bracket
; Assign an "upset percentage" to each game as a function of seeding
; difference, and then randomly draw to see if the underdog wins


function march_madness, seed_difference

ndiff = n_elements(seed_difference) 

; Assign a rough mapping between seed difference and upset percentage
difference_vec = lindgen(15)+1

upset_percentage_rough = 0.5d0 - (difference_vec-1d0)*0.4d0/max(difference_vec-1)

; For each game, perturb the upset percentage by a gaussian with 3%
; error, and then draw a winner

upset_percentage_input = interpol(upset_percentage_rough, difference_vec, seed_difference)

;upset_percentage_final = upset_percentage_input
upset_percentage_final = 0d0 > (upset_percentage_input * (1d0 + (randomn(seed,ndiff)*0.03d0))) < 1d0

;print, upset_percentage_final

; Draw random numbers, then check to see if they're greater than the
; upset percentage
random_draw = randomu(seed, ndiff)

; 1 means upset, 0 means no upset
upsets = random_draw le upset_percentage_final

null = where(upsets eq 1, nupsets, ncomplement=nwins)

;print, ''
;print, double(nupsets)/double(ndiff)

return, upsets

end
                    
