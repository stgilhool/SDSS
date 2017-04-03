pro transprob
; I could never work out how transit probability is equal to
; R_star/a.  It seemed to me that it should have some pi and stuff in
; there. I worked it out that it should be 2asin(R_star/a)/pi.  That
; matches my simulation perfectly, whereas R_star/a overestimates the
; transit probability a bit (and makes the biggest error at very small
; a (roughly 0.05au < a < 0.1 au).  So R_star/a is indeed an
; approximation, though I'm not sure exactly where it comes from.
; Maybe it's just a convenient relation that's reasonably close.

rstar = 7d5 ;km
au = 1.5d8 ;km

nsystems = 1d6

nsemi = 25
semimajor = (dindgen(nsemi)/100)*au+rstar 

transit_prob_sim = []

for i = 0, nsemi-1 do begin

    inclination = randomu(seed, nsystems)*!pi

    projected_height = semimajor[i] * cos(inclination)

    impact_parameter = abs(projected_height)/rstar

    transits = where(impact_parameter le 1, ntransits)

    transit_prob_i = ntransits/nsystems

    transit_prob_sim = [transit_prob_sim, transit_prob_i]

endfor

transit_prob_theory = rstar/semimajor

transit_prob_me = (2d0*asin(transit_prob_theory))/!pi

plot, semimajor, transit_prob_theory
oplot, semimajor, transit_prob_sim, color=200
dxstop
oplot, semimajor, transit_prob_me, color=99999


stop

end
