function quick_total, current_cap
; Gets totals for the current state, plus if we dump out any one or
; two cups

c = current_cap
allowed_values = [c, total(c), total(c[0:1]), total(c[1:2]), (c[0]+c[2])]
gd_idx = where(allowed_values gt 0, ngd)
if ngd gt 0 then begin
    allowed_values = allowed_values[gd_idx]
endif else message, "i think there's a problem"

return, allowed_values
end


pro full_total, cap_hash
; Gets totals for current state and keeps track of the unique ones

curr = cap_hash['curr']
newtots = quick_total(curr)

oldhist = cap_hash['tot']
newhist = histogram(newtots, binsize = 1, omin=omi, omax=oma)

oldhist[omi:oma] = oldhist[omi:oma] + newhist
updated_hist = oldhist gt 0

; update hash
cap_hash['tot'] = updated_hist
cap_hash['ntot'] = total(updated_hist)

return
end


pro pour_cup, cap_hash, outcup, incup

c = cap_hash

; check if move is allowed (outcup has something, incup has avail)
outcurr = c['curr',outcup]
;outstat = outcurr gt 0

incurr = c['curr', incup]
inavail = c['cap', incup] - incurr
;instat = inavail gt 0

changeamt = outcurr < inavail

;if outstat and instat then begin
if changeamt gt 0 then begin
    
    outcurr = outcurr - changeamt
    incurr = incurr + changeamt

    if outcurr lt 0 or outcurr gt c['cap',outcup] then message, "problem with outcup"
    if incurr lt 0 or incurr gt c['cap', incup] then message, "problem with incup"

    ; update hash
    cap_hash['opstatus'] = 1
    cap_hash['curr',outcup] = outcurr
    cap_hash['curr',incup] = incurr

endif else begin
    cap_hash['opstatus'] = 0
endelse

return
end


pro empty_cup, cap_hash, outcup

cap = cap_hash['cap',outcup]
current = cap_hash['curr',outcup]

totals = cap_hash['tot']


if current ne 0 then begin
    current = 0
    ;update hash
    cap_hash['curr',outcup] = current
    
    cap_hash['opstatus'] = 1
endif else begin
    cap_hash['opstatus'] = 0
endelse

return
end


pro do_pour, cap_hash, opnum

case opnum of
    0: pour_cup, cap_hash, 0, 1
    1: pour_cup, cap_hash, 0, 2
    2: pour_cup, cap_hash, 1, 0
    3: pour_cup, cap_hash, 1, 2
    4: pour_cup, cap_hash, 2, 0
    5: pour_cup, cap_hash, 2, 1
    else: message, "That operation is not allowed"
endcase  

return
end


pro do_empty, cap_hash, opnum

case opnum of
    0: empty_cup, cap_hash, 0
    1: empty_cup, cap_hash, 1
    2: empty_cup, cap_hash, 2
    else: message, "That operation is not allowed"
endcase 

return
end



pro threeGlasses, total_cap

fullcheck = total_cap gt 0
nottoofullcheck = total_cap le 100

if total(fullcheck) ne 3 then message, "Each cup must hold at least 1 unit"
if total(nottoofullcheck) ne 3 then message, "Each cup must be less than 100 units"

; Get the easy totals
totals = quick_total(total_cap)
totals = totals[uniq(totals)]
tothist = histogram(totals, min=0, max=total(total_cap), binsize=1)
ntot = n_elements(totals) 

; Make hash to store data
cap_hash = hash('cap', total_cap, $
                'curr', total_cap, $
                'tot', tothist, $
                'ntot', ntot, $
                'opstatus', 1)

; Start doing operations
terminate = 0
trial_no_new_counter=0

while terminate eq 0 do begin
    
    ;refill
    cap_hash['curr'] = cap_hash['cap']
    refill = 0
    trial_total = cap_hash['ntot']

    while refill eq 0 do begin
        
        ;empty 1 cup to start
        emptystat = 0
        while emptystat eq 0 do begin
            
            ;do_empty, cap_hash, fix(randomu(seed,1)*3)
            curr_cap = cap_hash['curr']
            
            csort = sort(curr_cap)
            ccap_sort = curr_cap[csort]

            sort_pos_idx = where(ccap_sort gt 0, npos)
            ;if npos gt 0 then begin
                
            small_bin = sort(curr_cap[0])
            med_bin = sort(curr_cap[1])
            big_bin = sort(curr_cap[2])

            
            ;check status
            opstat = cap_hash['opstatus']
            if opstat then begin
                print, "Poured out a cup"
                ; if we've emptied the last cup, we have to start over
                if total(cap_hash['curr']) eq 0 then begin
                    refill = 1
                    print, 'Refilling...'
                    empty_another = 1
                    emptystat =1
                endif else begin
                    empty_another = 0
                    full_total, cap_hash
                    emptystat = 1
                endelse
            endif
        endwhile

        
        
        ; Now, pour until nothing happens
        no_new_counter = 0
        
        while empty_another eq 0 do begin
        
            ;remember how many totals we have now
            current_total = cap_hash['ntot']
            ; do a random pour operation
            do_pour, cap_hash, fix(randomu(seed,1)*6)
            ;check status
            opstat = cap_hash['opstatus']
            ; add new guys if legal op
            if opstat then begin
                full_total, cap_hash
        
                nnew = cap_hash['ntot']-current_total
                ; 10 consecutive no changes means we empty another cup
                if nnew eq 0 then no_new_counter++ else no_new_counter = 0
            endif
            
            print, "Total is: " + strtrim(cap_hash['ntot'],2)
            print, "No New Counter: " + strtrim(no_new_counter,2)
            print, ''

            if no_new_counter ge 25 then empty_another = 1
        endwhile

    endwhile

    new_trial_total = cap_hash['ntot']
    if new_trial_total eq trial_total then trial_no_new_counter++ $
      else trial_no_new_counter = trial_no_new_counter/2

    print, "Total is: " + strtrim(cap_hash['ntot'],2)
    print, "Trial No New Counter: " + strtrim(trial_no_new_counter,2)
    print, ''

    
    if trial_no_new_counter ge 1000 then terminate = 1
endwhile

final_total = cap_hash['ntot']
total_amts = where(cap_hash['tot'] gt 0, ntotcheck)

if final_total ne ntotcheck then print, "final thingies don't match"

print, final_total
print, total_amts

        

end
