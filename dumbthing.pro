pro dumbthing

n = lindgen(50)+1

p = (1/3.)*((2/3.)^n)

navg = total(n*p, /double)

cash_n = 1.5

print, navg
print, "you make " + strtrim(navg*1.5, 2) + " dollars on average"

window, 0
plot, n, p


stop

ntrials = 1d5

netcash = lonarr(ntrials)
rollsav = ptr_new(/allocate_heap)

for i = 0, ntrials-1 do begin
    terminate=0
    cash = 0
    while terminate ne 1 do begin
        roll = floor(randomu(seed, 1)*3.)
        *rollsav = [*rollsav, roll[0]]
        case roll[0] of
            0: terminate=1
            1: cash = cash + 1
            2: cash = cash + 2
            3: cash = cash + 2
        endcase
    
    endwhile
    
    netcash[i] = cash
endfor
        
h = histogram(netcash, reverse_indices=ri, omin=om, omax=omax)
window, 1
plot, h, ps = 10

hcheck = histogram(*rollsav)
window,2
plot, hcheck, ps=10

print, "Average from trials is: " + strtrim(mean(netcash),2)

stop

end
