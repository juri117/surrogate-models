function res = rbf(x,knownX,lambdas,rbFunction)
  n = length(knownX);
  res = 0.;
  for i=1:n
    res += lambdas(i) * rbFunction(sqrt((x - knownX(i))**2));
  endfor
endfunction
