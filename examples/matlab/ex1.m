% EXAMPLE OF derfit.m

addpath ../../lib/matlab

dat = importdata('../averages.dat')';

cc = dat(2,:);
xx = dat(1,:);
yy = dat(3,:);
dy = dat(4,:);

W=diag(1./dy.^2);
Ndat=length(xx);

func = @(x,p,c) (c==1)*p(1)*(1+p(4)*x) + (c==2)*p(2)*(1+p(4)*x) + (c==3)*p(3)*(1+p(4)*x);
dfunc = @(x,p,c) [(c==1)*(1+p(4)*x) ...
        (c==2)*(1+p(4)*x) ...
        (c==3)*(1+p(4)*x) ...
        (c==1)*p(1)*x + (c==2)*p(2)*x + (c==3)*p(3)*x];

% definition of chi2
eps = @(i,p) (func(xx(i),p,cc(i)) - yy(i))/dy(i);
chisquared = @(p) sum(arrayfun(@(i) eps(i,p)^2,1:Ndat));

% initial parameters
inip = [3 4 9 -1];

% minimization of chi2
[pars, chi2] = fminsearch(chisquared,inip)

% computation of derivatives
der = derfit(xx,yy,W,pars,'df',dfunc,'c',cc);

% computation of derivatives with test on gradient
der = derfit(xx,yy,W,pars,'df',dfunc,'f',func,'c',cc);

% error of parameters
dpars = (der.*repmat(dy',1,4)).^2;
dpars = sqrt(sum(dpars,1))

% plot
color={'r','b','k'};
xax=[0:1e-4:0.1]; Np=length(xax);
figure; hold on
for b=1:3
        errorbar((cc==b).*xx,(cc==b).*yy,(cc==b).*dy,['.' color{b}])
        errorbar(0,pars(b),dpars(b),'ok');
        plot(xax,arrayfun(@(i) func(xax(i),pars,b),1:Np),['--' color{b}])
end
hold off
title('Example of chiral extrapolation with combined fit')
