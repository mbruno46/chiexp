function [Q, dQ, h] = qfit(chi2,nu,ce,nmc,plotter)
%QFIT quality of fit
%   [Q,dQ] = QFIT(chi2,nu) computes the probability of the chi square
%   defined from the matrix nu, to be bigger than chi2, with a Monte Carlo
%   sampling process. 
%
%   [Q,dQ] = QFIT(chi2,nu,ce) computes the quality of fit using the formula
%   based on the eigenvalues of nu. ce must be a two-element array with the
%   expected chi square and its error, as obtained from chiexp. If empty
%   the Monte Carlo sample is used to evaluate Q and dQ.
%
%   [Q,dQ]= QFIT(chi2,nu,[],nmc,'on') computes the quality of fit using a
%   Monte Carlo sample of size nmc (default is 5000) and plots the
%   histogram of the probability distribution of the chi square defined
%   from nu, together with the value of chi2.
%
%   [Q,dQ,h] = QFIT(chi2,nu) returns the Monte Carlo history of the chi
%   square defined from nu in the output argument h.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2017-19 Mattia Bruno, Rainer Sommer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
    plotter='off';
    if nargin<4
        nmc=5000;
        if nargin<3
            ce=[];
            if nargin<2
                error('Not enough input arguments');
            end
        end
    end
end

% checks nu
[N,~] =size(nu);
if isreal(nu)
    if unique(abs(nu-nu')>1e-10)~=0
        error('nu is not symmetric')
    end
end

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if ~isOctave % rng not defined in Octave
    rng('shuffle','twister');  % change seeds
end
z=randn(N,nmc);

if ~isempty(ce)
    % computes eigenvalues and chops them
    eval=abs(eig(nu));
    eps=1e-14 * max(eval);
    eval=eval.*(eval>eps);

    chisq=eval'*(z.^2);
    Q=1.0-mean(chisq<chi2);

    x=chi2 - eval(2:end)'*(z(2:end,:).^2);
    x=x/eval(1);
    dQ = mean( (x>0).*exp(-x*0.5)*0.5./sqrt(abs(x)) );
    dQ = ce(2)/ce(1) * dQ;
else
    chisq=sum(z.*(nu*z) ,1);
    th=chisq<chi2; 
    Q=1.0-mean(th);
    dQ=sqrt(var(th)/nmc);
end

if nargout>2
    h=chisq;
end

if strcmp(plotter,'on')
    figure; hold on
    histogram(chisq,40,'Normalization','probability');
    yl=ylim;
    plot([chi2,chi2],[0,yl(2)])
    hold off
    title('Probability distribution')
end

