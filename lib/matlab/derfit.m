function D = derfit(X,Y,W,p,varargin)
%DERFIT derivatives of fit parameters
%   D = DERFIT(X,Y,W,p,'field',value) computes the derivatives of the 
%   parameters with respect to Y, at the minimum of the chi square. 
%   It can be used in the estimates of the error of the parameters using 
%   the chain rule. It returns a Na-by-N matrix, D, with the analytical 
%   derivatives of the parameters `der(i,a) = dp(a)/dy(i)` 
%   (Na is the number of parameters).
%   X and Y are arrays of length N; Y is assumed to contain the central 
%   values of the observable; W is the weight matrix, N-by-N, used in the 
%   fit; p is an array with the values of the parameters at the minimum.
%
%   D = DERFIT(X,Y,W,p,'df',df) computes the derivatives of the parameters
%   using the analytical gradient of f provided by the user. df must be an
%   array of function handles in the form `df = @(x,p) [df/dp1 df/dp2 ...]`
%
%   D = DERFIT(X,Y,W,p,'f',f) computes the derivatives of the parameters
%   using the numerical gradient of f. f must be a function handles in the
%   form `f = @(x,p) f(x,p)`.
%
%   D = DERFIT(X,Y,W,p,'df',df,'f',f) computes the derivatives of the 
%   parameters using the analytical gradient df. In this case f is used to
%   check (numerically) the gradient passed by the user.
%
%   D = DERFIT(X,Y,W,p,'df',df,'c',c) uses c, an array of length N which 
%   is passed to f and df, which in this case must obey the syntax 
%   `f = @(x,p,c) f(x,p,c)` (same for `df`). c is an additional parameter
%   useful in global fits.
%   
%   Examples:
%   >> X = [0,1,2,3...]
%   >> [Y,DY] = compute_central_values_and_errors()
%   >> W = diag(1./DY.^2)
%   >> f = @(x,p) p(1) + p(2)*x
%   >> df = @(x,p) [1, x]
%   >> D = DERFIT(X,Y,W,p,'df',df) % analytical
%   >> D = DERFIT(X,Y,W,p,'f',f) % numerical
%   >> D = DERFIT(X,Y,W,p,'df',df,'f',f) % numerical check
%
%   See also <a href="#">the project webpage</a>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2017-19 Mattia Bruno, Rainer Sommer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    error('Not enough input arguments');
end

N=length(X);
Na=length(p);

% checks W
s=size(W);
if s(1)~=s(2) && s(1)~=N
    error('Dimension mismatch between W and Y');
end
%if ~issymmetric(W)
%    error('W is not symmetric');
%end

% varargins
% DEFAULTS
df=[]; c=[]; f=[];

for k=1:2:length(varargin)
    if ~ischar(varargin{k})
        error('Additional input field passed with incorrect format');
    end
    switch varargin{k}
        case 'df'
            df=varargin{k+1};
        case 'f'
            f=varargin{k+1};
        case 'c'
            c=varargin{k+1};
            if length(c)~=N
                error(['Dimension mismatch, number of elements in Y and c' ...
                    'must be the same']);
            end
    end
end

    
if isempty(df)
    GRAD = numder(X,p,f,c,EPS);
else
    GRAD = zeros(N,Na);
    for i=1:N
        if isempty(c)
            GRAD(i,:) = df(X(:,i),p);        
        else
            GRAD(i,:) = df(X(:,i),p,c(i));
        end
    end
end

H = GRAD'*W*GRAD;
[U,S,V] = svd(H);
Hinv = V*inv(S)*U';

D=zeros(N,Na);

GG=W*GRAD;
for i=1:N
    D(i,:) = GG(i,:)*Hinv;
end

if ~isempty(df) && ~isempty(f)
    EPS=1e-6;
    
    % CHECK GRADIENT
    NUMDER = numder(X,p,f,c,EPS);
    compare = reshape(abs(GRAD - NUMDER),1,N*Na);
    disp(' ');
    disp(['Maximum difference between the provided gradient ' ...
        'and the numerical one']);
    [mm, idx] = max(compare);
    disp(['Step size = ' num2str(EPS) ' -> ' num2str(mm) ' at dfunc/dp(' ...
        num2str(round(idx/N)+1) ')']);
end
end


function der = numder(xdata,pars,func,c,h)
[~, N]=size(xdata); Np=length(pars);
der = zeros(N,Np);
for j=1:Np
    pF=pars; pF(j)=pars(j)+h;
    pB=pars; pB(j)=pars(j)-h;
    if isempty(c)
        F = arrayfun(@(i) func(xdata(:,i),pF),1:N);
        B = arrayfun(@(i) func(xdata(:,i),pB),1:N);
    else
        F = arrayfun(@(i) func(xdata(:,i),pF,c(i)),1:N);
        B = arrayfun(@(i) func(xdata(:,i),pB,c(i)),1:N);
    end
    der(:,j) = [(F-B)./(2*h)]';
end
end
