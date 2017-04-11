function index_confidence, vector, confidence=confidence

if n_elements(confidence) eq 0 then confidence = 0.9

n_entries = n_elements(vector) 

n_inrange = ceil(n_entries * confidence)

sort_idx = sort(vector)

svec = vector(sort_idx)

; Run through all possible ranges and pick the smallest interval
n_possibilities = n_entries - n_inrange + 1

conf_interval = [0d0,0d0]


for trial = 0, n_possibilities-1 do begin

    range = svec[trial+n_inrange-1] - svec[trial]

    if trial eq 0 then begin 
        minrange = range
        conf_interval = [svec[trial], svec[trial+n_inrange-1]]
    endif
    
    if range lt minrange then begin
        minrange = range
        conf_interval = [svec[trial], svec[trial+n_inrange-1]]
    endif

    print, "Range is: "+strtrim(range,2)
    print, "Minrange is: "+strtrim(minrange,2)
    print, "Confidence bounds are: "+strtrim(conf_interval,2)
    print, ''

endfor

return, conf_interval

end
    
