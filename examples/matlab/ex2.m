% CHECKS OF chiexp.m AND qfit.m 

addpath ../../lib/matlab

% we consider a two-point correlator of a free scalar field projected to
% zero momentum and we imagine to fit it over T time slices from zero
% Since we know the exact values of the parameters the chi^2 will be zero
% by definition; however we will see how different definitions of W lead to
% different expectation values of the chi^2
am=0.5;
T=16;
ncnfg=1000;

% analytical prediction with infinite Time extent
omega = acosh(0.5*am^2+1.);
Gt = @(x) exp(-abs(x)*omega)/(2*sinh(omega));

% covariance matrix
cov_theo = zeros(T,T);
for t1=0:T-1
    for t2=0:T-1
        cov_theo(t1+1,t2+1) = Gt(t1)*Gt(t2) + Gt(t1-t2)*Gt(0);
    end
end

fprintf('Covariance matrix eigenvalues\n');
disp(eig(cov_theo));

xx=[0:T-1];
yy=arrayfun(@(t) Gt(t),0:T-1);
dy=sqrt(arrayfun(@(i) cov_theo(i,i),1:T));

% fitted function and derivative
f = @(x,p) p(1)*exp(-p(2)*x);
df = @(x,p) [exp(-p(2)*x), -p(1)*x*exp(-p(2)*x)];

W=diag(1./dy.^2);

% definition of chi2
eps = @(p) arrayfun(@(i) f(xx(i),p) - yy(i),1:T);
chisquared = @(p) eps(p)*W*eps(p)';

% chi2 is exactly zero in this example
pars=[1/(2*sinh(omega)), omega];
chi2=chisquared(pars);
fprintf('amp  = %.4f \n',pars(1));
fprintf('mass = %.4f \n',pars(2));
fprintf('chi2 = %.4f / dof = %d \n',chi2,T-length(pars));

fprintf('\nW = diag of cov matrix \n');

% chiexp < dof due to W 
[ce,dce,nu,covest,ce2,dce2] = chiexp(xx,yy,W,pars,cov_theo,'df',df);
fprintf('chiexp = %.4f +- %.4f \n',ce,dce);
fprintf('chiexp (2nd method) = %.4f +- %.4f \n',ce2,dce2);
assert(ce < T-length(pars));

[Q,dQ,h] = qfit(chi2,nu,[],5000,'off');
fprintf('ChiExp from MC = %4.4f +- %.4f\n',mean(h),sqrt(var(h)/length(h)));

[U,S,V] = svd(cov_theo);
W = V*diag(1./diag(S))*U';

fprintf('\nW = inv of cov matrix \n');

% chiexp = dof thanks to W equal to cov
[ce,dce,nu,covest,ce2,dce2] = chiexp(xx,yy,W,pars,cov_theo,'df',df);
fprintf('chiexp = %.4f +- %.4f \n',ce,dce);
fprintf('chiexp (2nd method) = %.4f +- %.4f \n',ce2,dce2);
assert(abs(ce - T + length(pars)) < 1e-10);

[Q,dQ,h] = qfit(chi2,nu,[],5000,'off');
fprintf('ChiExp from MC = %4.4f +- %.4f\n',mean(h),sqrt(var(h)/length(h)));