function [ MSWD,Wm,ErrWm ] = CorrWM( X,Y )
% Computes corrected weighted mean.

% Find input dimension
if nargin<2;
    [nl,nc]=size(X);
    if (nl==1 || nc==1);
        MSWD=NaN;
        Wm=X(1);
        ErrWm=X(2);
        fprintf('Warning, one sample input in CorrWm\n')
        return
    end
    
    if nl==2 && nc==2;
        fprintf('Warning CorrWm, 2 x 2 matrix, raw data is the default mode\n')
    end
    
    if nc==2;
        errX=X(:,2);
        X=X(:,1);
    elseif nl==2;
        errX=X(2,:);
        X=X(1,:);
    else
        MSWD=NaN;
        Wm=NaN;
        ErrWm=NaN;
        fprintf('No dimension of X equal to 2\n')
        return
    end
else
    X1=ones(1,length(X));
    X1(1:end)=X(1:end);
    X=X1;
    errX1=ones(1,length(Y));
    errX1(1:end)=Y(1:end);
    errX=errX1;
    if length(X)~=length(Y);
        MSWD=NaN;
        Wm=NaN;
        ErrWm=NaN;
        fprintf('length(X) different from length(Y)\n')
        return
    end
end

% Weighted mean
W=1./(errX.^2);
Wm=sum(W.*X)/sum(W);

% MSWD
MSWD=sum(((X-mean(X))./(errX)).^2)/(length(X)-1);

% Determination
Fact=1;
if MSWD>1;
    Fact=sqrt(MSWD);
end

% Error bar
ErrWm=Fact*sqrt(1/sum(W));

end

