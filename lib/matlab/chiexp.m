function [ce,dce,nu,covest,ce2,dce2] = chiexp(X,Y,W,p,cov,varargin)
%CHIEXP expected chi square
%   [ce,dce] = CHIEXP(X,Y,W,p,cov,'field',value) computes the
%   expected chi square, following the derivation presented in <a
%   href="#">this Reference</a>. 
%   X and Y are arrays of length N; Y is assumed to contain the central 
%   values of the observable; W is the weight matrix, N-by-N, used in the 
%   fit; p is an array with the values of the parameters at the minimum.
%   if cov is a N-by-N matrix, the program takes it as the covariance 
%   matrix, assuming the user has previously computed it; if cov is a
%   M-by-N matrix, the program assumes that it contains the values of the
%   N observables Y on M different configurations.
%   The function returns the expectation value of the chi square and its
%   error (ce and dce).
%
%   [ce,dce] = CHIEXP(X,Y,W,p,cov,'df',df) computes the expectation value 
%   of the chi square using the analytical gradient of f provided by the 
%   user. df must be an array of function handles in the form df = @(x,p) 
%   [df/dp1 df/dp2 ...]
%
%   [ce,dce] = CHIEXP(X,Y,W,p,cov,'f',f) computes the expectation value of 
%   the chi square using the numerical gradient of f, which must be a 
%   function handle in the form f = @(x,p) f(x,p).
%
%   [ce,dce] = CHIEXP(X,Y,W,p,'df',df,'field',value) allows the 
%   user to provide additional arguments:
%     - 'Nrep': array with the number of configurations belonging to
%               each replica, e.g. Nrep=[N1, N2, N3 ...]
%     - 'Stau': value of the parameter Stau. Default is 1.5. If set to 0
%               absence of autocorrelations is assumed.
%     - 'plotter': flag to produce a plot of the expected chi square as a 
%               function of the window. Default is 'off'. Accepted values 
%               are 'on' and 'off'.
%     - 'dp': array with the errors of the parameters at the minimum. It is
%               used in the estimates of the numerical derivatives of f,
%               unless the gradient df is passed by the user. These errors
%               can be obtained from the function derfit.
%
%   Note that if cov is a N-by-N matrix 'Nrep', 'Stau' and 'plot' are ignored.
%
%   [ce,dce,nu,covest,ce2,dce2] = CHIEXP(X,Y,W,p,cov,'df',df) produces
%   additional output: the matrix nu=C^{1/2} W^{1/2} (1-P) W^{1/2} C^{1/2},
%   the covariance matrix covest defined from the same window determined 
%   for the expected chi square; if W is diagonal a second method to 
%   compute the expected chi square and its error is used and returned
%   as ce2 and dce2.
%
%   Examples:
%
%   See also <a href="https://mabruno.gitlab.io/chiexp">the project webpage</a>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2017-19 Mattia Bruno, Rainer Sommer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: - improve numerical derivative with more sophisticated function?
%       - if covest has neg eigenvalues return error, not warning
%       - in num. gradient check if dp*1e-4 is sensible

if nargin<7
    error('Not enough input arguments');
end

% checks X and Y
N=length(X);
if length(Y)~=N
    error(['Dimension mismatch: the number of elements in X and Y must be' ...
        ' the same']);
end

% checks W
s=size(W);
if s(1)==1
    if s(2)~=N
        error('Dimension mismatch between Y and W');
    else
        Wmat=diag(W);
    end
else
    if s(1)==s(2) && s(1)==N
        Wmat=W;
    else
        error('Dimension mismatch between Y and W');
    end
end

for k=1:2:length(varargin)
    if ~ischar(varargin{k})
        error('Additional input field passed with incorrect format');
    end
end

% checks p
Na=length(p);

% global varargin c and plot
c=[]; plotter='off'; dp=[];
for k=1:2:length(varargin)
    switch varargin{k}
        case 'c'
            c=varargin{k+1};
        case 'plotter'
            plotter=varargin{k+1};
        case 'dp'
            dp=varargin{k+1};
    end
end

% --------------------------------------------------------------

% varargin f vs df
df=[]; f=[]; 
for k=1:2:length(varargin)
    switch varargin{k}
        case 'f'
            f=varargin{k+1};
        case 'df'
            df=varargin{k+1};
    end
end

if isempty(df)
    if isempty(f)
        error('df (or f) is a required input');
    end
    if isempty(dp)
        EPS=1e-6;
        grad = numder(X,p,f,c,repmat(EPS,length(p)));
    else
        grad = numder(X,p,f,c,dp*1e-4);
    end 
else
    grad = zeros(N,Na);
    for i=1:N
        if isempty(c)
            grad(i,:) = df(X(:,i),p);
        else
            grad(i,:) = df(X(:,i),p,c(i));
        end            
    end
end

Wg = Wmat*grad;
H = grad'*Wmat*grad;
[U,S,V] = svd(H);
Hinv = V*inv(S)*U';

P=Wmat - Wg*Hinv*Wg';

% --------------------------------------------------------------

s=size(cov);
if s(1)==s(2) && s(1)==N
    ce=trace(P*cov);
    dce=0.0;
    if nargout>2
        covest=cov;
    end
    if nargout>4
        ce2=ce;
        dce2=dce;
    end
elseif s(2)==N
    ncnfg=s(1);
    yy = mean(cov,1);
    delta= cov - repmat(yy,ncnfg,1); % subtract mean, just in case

    Nrep=ncnfg; Stau=1.5;
    for k=1:2:length(varargin)
        switch varargin{k}
            case 'Nrep'
                Nrep=varargin{k+1};
                if sum(Nrep)~=ncnfg
                    error('Mismatch between Nrep and Y');
                end
            case 'Stau'
                Stau=varargin{k+1};
                if Stau<0
                    error('Unexpected value for Stau');
                end
        end
    end
    
    gam = gamma(delta,delta,Nrep);
    chiexp_t = compute_chiexp(gam,P,ncnfg);
    if Stau==0.0
        wopt=0;
    else
        wopt = find_Wopt(chiexp_t,Stau,ncnfg);
    end
    ce=chiexp_t(1) + 2*sum(chiexp_t(2:wopt+1));
    dce=2*sqrt(wopt/ncnfg)*ce;
        
    if strcmp(plotter,'on')
        figure; hold on
        plot([wopt,wopt],[0.0,0.8],'-r');
        plot([0:2*wopt-1],chiexp_t(1:2*wopt)/chiexp_t(1));
        hold off
        title(['Chiexp, Wopt = ', num2str(wopt)])
    end
    
    if nargout>2
        covest=gam(:,:,1) + 2*sum(gam(:,:,2:wopt+1),3);
        
        eval=eig(covest);
        if ~isempty(eval(eval<0))
            warning(['Estimated covariance matrix has negative eigenvalues'...
                ' with automatic window from chiexp']);
            decr=1;
            while ~isempty(eval(eval<0)) && (wopt-decr)>0
                covest=gam(:,:,1) + 2*sum(gam(:,:,2:wopt-decr+1),3);
                eval=eig(covest);
                decr=decr+1;
            end
            fprintf('Reducing summation window to %d, from %d \n', wopt-decr+1, wopt+1);
        end
        
        covest=covest/ncnfg;
    end
    
    if nargout>4
        if isdiag(Wmat)
            PP=Wg*Hinv*Wg';

            chiexp_t = compute_chiexp(gam,PP,ncnfg);
            wopt = find_Wopt(chiexp_t,Stau,ncnfg);
            ce2=chiexp_t(1) + 2*sum(chiexp_t(2:wopt+1));
            dce2=2*sqrt(wopt/ncnfg)*ce2;
            ce2=N-ce2;
        else
            ce2=ce;
            dce2=dce;
        end
    end
else
    error('Dimension mismatch between cov and Y');
end

% --------------------------------------------------------------

if nargout>2
    [U,S,V] = svd(covest);
    Csq = U*sqrt(S)*V';
    nu = Csq*P*Csq;
end

end


function chiexp_t = compute_chiexp(gam,mat,nc)
[~,~,N] = size(gam);
chiexp_t = zeros(1,N);
for t=1:N
    chiexp_t(t) = trace(mat*gam(:,:,t))/nc;
end
end


function gam = gamma(dat1,dat2,Nrep)
tmax = min(floor(min(Nrep)/2),1000);
[N, M] = size(dat1);
gam = zeros(M,M,tmax+1);
R = length(Nrep);

% sum over replica:
offset=1;
for r=1:R
    imax=offset-1+Nrep(r);
    aux1 = fft(dat1(offset:imax,:),[],1);
    aux2 = fft(dat2(offset:imax,:),[],1);
    for i=1:M
        for j=i:M
            tmp = ifft(aux1(:,i).*conj(aux2(:,j)));
            gam(i,j,:) = gam(i,j,:) + reshape(tmp(1:tmax+1),1,1,tmax+1);
        end
    end
    offset=offset+Nrep(r);
end

n = reshape(1./(N-R*(0:tmax)),1,1,tmax+1);
for i=1:M
    for j=i:M
        gam(i,j,:)=gam(i,j,:).*n;
        if i~=j
            gam(j,i,:)=gam(i,j,:);
        end
    end
end
end


function [Wopt, tmax] = find_Wopt(g,Stau,N)
Gint = 0;
flag=1;

for t=1:length(g)
    Gint=Gint+g(t+1)/g(1);
    if Gint <= 0, tauW=eps; else tauW = Stau/(log((Gint+1)/Gint)); end
    gW = exp(-t/tauW)-tauW/sqrt(t*N);
    if gW < 0                % this W is taken as optimal
      Wopt=t;
      tmax=min(length(g),2*t);
      flag=0;
      break;
    end
end

if flag
    warning(['Windowing condition failed up to W = ' num2str(length(g))]);
    Wopt=tmax;
end
end


function der = numder(xdata,pars,func,c,h)
[~, N]=size(xdata); Np=length(pars);
der = zeros(N,Np);
for j=1:Np
    pF=pars; pF(j)=pars(j)+h(j);
    pB=pars; pB(j)=pars(j)-h(j);
    if isempty(c)
        F = arrayfun(@(i) func(xdata(:,i),pF),1:N);
        B = arrayfun(@(i) func(xdata(:,i),pB),1:N);
    else
        F = arrayfun(@(i) func(xdata(:,i),pF,c(i)),1:N);
        B = arrayfun(@(i) func(xdata(:,i),pB,c(i)),1:N);
    end
    der(:,j) = [(F-B)./(2*h(j))]';
end
end

    