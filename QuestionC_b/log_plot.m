b = 200;
a = 10;
del = 10;
h = a:del:b;
nh = (b-a)/del + 1;
errs = zeros(1,nh);
Jac_diffs = zeros(1,nh);
iterations = zeros(1,nh);

for i=1:nh
    npts = del*i;
    npts
    [errs(1,i), Jac_diffs(1,i), iterations(1,i)] = solve_jacobi_iter(npts, npts);
end

delh = zeros(1,nh)

for i=1:nh
    delh(1,i) = 1/(h(i)+1);
end

l_delh = log(delh);
l_err = log(errs);
plot(l_delh, l_err)
xlabel('log h')
ylabel('log max-norm error')



    
