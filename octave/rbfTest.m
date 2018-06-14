rbfConstant = .01;
coords = [0., 2., 4., 6., 8., 10.]'
values = [0.95, 2.1432974268256815, 1.1371975046920717, 2.0745845018010733, 2.6433582466233814, -2.0940211108893703]'

[b,coreFunction] = RBFcalcCoeff(rbfConstant,coords,values)

n = 100;
res = zeros(n, 1);

xs = 0:0.01:10;

for i=1:length(xs)
  res(i) = rbf(xs(i), coords, b, coreFunction);
end

plot(xs, res);

