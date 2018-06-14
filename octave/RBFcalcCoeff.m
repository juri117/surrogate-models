function [b,coreFunction] = RBFcalcCoeff(rbfConstant,coords,values)
  [m n] = size(coords);

  coreFunction = @(x) exp(-(rbfConstant.*x).^2);   %%% define rbf core function

  A = zeros(m,m);     %%% init

  for i1=1:m          %%% loops to fill matrix A
    for i2=i1:m
      cummSum = 0;    %%% calc radius in n-dimensional space
      for k=1:n
        cummSum = cummSum + (coords(i1,k)-coords(i2,k)).^2;
      end
      radius = sqrt(cummSum);
      A(i1,i2) = coreFunction(radius);
      A(i2,i1) = A(i1,i2);  %%% save in matrix
    end
  end

  b=A\values;         %%% solve linear equation system