; Code to test if I can probe underlying v distribution
; given a distribution of i (cosi is uniform)
pro i_dist

; First, let's test that p(i) = sin(i) if cos(i)
; is uniform

; Make sin(i) distribution of i
i_vec = []

repeat begin

; create a bunch of random numbers between 0 and pi/2
num_gen = randomu(seed, 10e3) * (!pi/2.0)
; probability of each number is its sin
prob_num_gen = sin(num_gen)

; create more random numbers which will serve to check if
; we accept or reject the num_gen numbers
chk_gen = randomu(seed, 10e3)
; accept or reject
gd_num_gen = where(prob_num_gen ge chk_gen, gdcnt)

if gdcnt gt 0 then begin
    keepers = num_gen[gd_num_gen]
    ; add the keepers to the i_vec
    i_vec = [i_vec, keepers]
endif

endrep until n_elements(i_vec) ge 10e4

isrt = sort(i_vec)
i_vec_srt = i_vec[isrt]

;i_hist = histogram(i_vec_srt, binsize = 0.01, omin=om)
binsize = 0.01
i_hist = histogram(i_vec_srt, binsize = binsize, min=0)

i_hist = i_hist/total(i_hist*binsize)

i_vals = (dindgen(n_elements(i_hist))*binsize)

plot, i_vals, i_hist, ps=10;, yr=[0,0.015]
;oplot, i_vals, sin(i_vals)/total(sin(i_vals)*binsize), color=200
oplot, i_vals, sin(i_vals), color=200
stop
; Check that sini is distributed as tani

;;;;NOPE, sini is distributed as p(y) = y/sqrt(1-y^2)
; If y = sini, it's basically the same thing as tan, but
; the xscale changes such that we need the above form
sini_vec = sin(i_vec_srt)
;sini_hist = histogram(sini_vec, binsize = binsize, omin=om)
sini_hist = histogram(sini_vec, binsize = binsize, min=0)
sini_hist = sini_hist/total(sini_hist*binsize)
sini_vals = (dindgen(n_elements(sini_hist))*binsize)
plot, sini_vals, sini_hist, ps=10
;tani = tan(sini_vals)
tani = sini_vals/(sqrt(1d0-(sini_vals^2)))
oplot, sini_vals, tani, color=200
stop
; Another check
; cosi is uniformly distributed
cosi_vec = randomu(seed, 10e6)
i_acos_vec = acos(cosi_vec)
;i_acos_hist = histogram(i_acos_vec, binsize = binsize, omin=om)
i_acos_hist = histogram(i_acos_vec, binsize = binsize, min=0)

i_acos_hist = i_acos_hist/total(i_acos_hist*binsize)

i_acos_vals = (dindgen(n_elements(i_acos_hist))*binsize)

plot, i_vals, i_hist, ps = 10, xr = [0, max(i_vals) > max(i_acos_vals)]
oplot, i_acos_vals, i_acos_hist, color=200, ps = 10

stop

end
