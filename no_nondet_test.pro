; Quick test to see what the probability of getting 0 non-detections
; is, as a function of N stars, if all stars are rotating quickly (50 km/s)

function sini_probability, sample
common pdf_info, age_pdf, age_vec
; Give the probability for each element in sample

sini_prob = sample/sqrt(1d0-sample^2)

return, sini_prob

end


pro no_nondet_test

