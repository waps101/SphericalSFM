function [ R_est,t_est ] = SnP( x,X,method,niter,gamma )
%SNP Spherical-n-points using iteratively reweighted least squares
%
% Inputs:
% - X is 3 by n matrix of world points
% - x is 3 by n matrix of unit vector spherical image points
% - method is either 'soft', 'hard' or 'unconstrained' (default 'soft'). If
% method='soft' then the CVX toolbox is required.
% - if method is 'soft' or 'hard' niter is the number of iterations of
% iterative reweighting. If niter=1 no reweighting is applied. (default 1)
% - if method is 'soft' gamma is the weight for the dot product constraint
% (default 1)
%
% Outputs:
% R_est and t_est are the estimated rotation matrix and translation vector
%
% William Smith and Hao Guan 2017

% Set default parameters
if nargin<3
    method = 'hard';
end
if nargin<4
    niter = 1;
end
if nargin<5 && strcmp(method,'soft')
    gamma = 1;
end

npoints = size(x,2);

% Initialise weights to 1, i.e. first iteration is unweighted
weights = ones(npoints,1);

% Build the two matrices:
% 1. Evaluate negative dot product
C = -[x(1,:)'.*X(1,:)' x(1,:)'.*X(2,:)' x(1,:)'.*X(3,:)' ...
    x(2,:)'.*X(1,:)' x(2,:)'.*X(2,:)' x(2,:)'.*X(3,:)' ...
    x(3,:)'.*X(1,:)' x(3,:)'.*X(2,:)' x(3,:)'.*X(3,:)' ...
    x(1,:)' x(2,:)' x(3,:)'];
% 2. Evaluate cross product
A = [zeros(npoints,1)     zeros(npoints,1)      zeros(npoints,1)      ...
    -x(3,:)'.*X(1,:)'  -x(3,:)'.*X(2,:)'   -x(3,:)'.*X(3,:)'   ...
    x(2,:)'.*X(1,:)'   x(2,:)'.*X(2,:)'    x(2,:)'.*X(3,:)'    ...
    zeros(npoints,1)                  -x(3,:)'            x(2,:)';
    x(3,:)'.*X(1,:)'   x(3,:)'.*X(2,:)'    x(3,:)'.*X(3,:)'       ...
    zeros(npoints,1)     zeros(npoints,1)      zeros(npoints,1)    ...
    -x(1,:)'.*X(1,:)'  -x(1,:)'.*X(2,:)'   -x(1,:)'.*X(3,:)'    ...
    x(3,:)'            zeros(npoints,1)                  -x(1,:)';
    -x(2,:)'.*X(1,:)'  -x(2,:)'.*X(2,:)'   -x(2,:)'.*X(3,:)'      ...
    x(1,:)'.*X(1,:)'   x(1,:)'.*X(2,:)'    x(1,:)'.*X(3,:)'    ...
    zeros(npoints,1)     zeros(npoints,1)      zeros(npoints,1)     ...
    -x(2,:)'            x(1,:)'             zeros(npoints,1);];

% These are the linear equalities that are tried, corresponding to each
% element of row 1 of R being set to +/-1
ks = [1 1 2 2 3 3];
vals = [1 -1 1 -1 1 -1];

if strcmp(method,'hard')
    warning('off','optimlib:lsqlin:LinConstraints');
    options=optimset('Display','off');
    for iter = 1:niter
        A_scaled = A.*repmat([weights; weights; weights],[1 12]);
        maxdotp = -1e10;
        for i=1:6
            Aeq=zeros(1,12);
            Aeq(ks(i))=1;
            beq=vals(i);
            b = lsqlin(A_scaled,zeros(3*npoints,1),C,zeros(npoints,1),Aeq,beq,[],[],[],options);
            [ R_est,t_est ] = ClosestRot( b );
            dotp = ResidualError( x,X,R_est,t_est );
            if dotp>maxdotp
                maxdotp=dotp;
                R_best = R_est;
                t_best = t_est;
            end
        end
        R_est = R_best;
        t_est = t_best;
        x_pred = R_est*X+repmat(t_est,[1 npoints]);
        weights = 1./(sum(x_pred.^2,1))';
    end
    warning('on','optimlib:lsqlin:LinConstraints');
    
elseif strcmp(method,'soft')
    for iter = 1:niter
        A_scaled = A.*repmat([weights; weights; weights],[1 12]);
        
        maxdotp = -1e10;
        for i=1:6
            b = cvxsolve(A_scaled,C,ks(i),vals(i),gamma);
            [ R_est,t_est ] = ClosestRot( b );
            dotp = ResidualError( x,X,R_est,t_est );
            if dotp>maxdotp
                maxdotp=dotp;
                R_best = R_est;
                t_best = t_est;
            end
        end
        R_est = R_best;
        t_est = t_best;
        
        x_pred = R_est*X+repmat(t_est,[1 npoints]);
        weights = 1./(sum(x_pred.^2,1))';
        
    end
elseif strcmp(method,'unconstrained')
    [~,~,V] = svd(A);
    [ R_est,t_est ] = ClosestRot(V(:,12));
end

end

function b = cvxsolve(A,C,k,val,gamma)

cvx_begin quiet
variable b(12)
minimize( sum_square(A*b)  + gamma*sum_square_pos(C*b) )
subject to
b(k)==val
cvx_end

end

function dotp = ResidualError( x,X,R_est,t_est )

npoints = size(X,2);
x_pred = R_est*X+repmat(t_est,[1 npoints]);
norms = sqrt(sum(x_pred.^2,1));
x_pred(1,:)=x_pred(1,:)./norms;
x_pred(2,:)=x_pred(2,:)./norms;
x_pred(3,:)=x_pred(3,:)./norms;

dotp=sum(sum(x_pred.*x,1));

end

function [ R_est,t_est ] = ClosestRot( b )

R = reshape(b(1:9),3,3)';
T = b(10:12);

% Force R to be valid rotation matrix
[U,S,V] = svd(R);
R_est = U*V';

if det(R_est)<0
    R_est = -R_est;
    T = -T;
end

scales = R_est./R;
weights = abs(R_est);
weights = weights(~isinf(scales));
scales = scales(~isinf(scales));
weights = weights./sum(weights);
scale = sum(scales.*weights);

t_est = scale*T;

end