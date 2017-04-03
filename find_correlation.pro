; Calculate the Pearson correlation coefficient for two vectors

function find_correlation, x, y

meanx = mean(x,/nan)
meany = mean(y,/nan)

diffx = x-meanx
diffy = y-meany

pcoeff = total(diffx*diffy, /double)/(sqrt(total(diffx^2, /double)) * sqrt(total(diffy^2, /double)))

return, pcoeff

end
