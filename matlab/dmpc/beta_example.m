%% example of Beta for one segment
% order of polynomial
n=2;
% dimension
dim = 3;

%delta is the matrix that transform the control points to polynomial
%coefficients
delta = zeros(n+1, n+1);
for k = 0:n
   for i = k:n
       delta(i+1, k+1) = (-1)^(i-k)*nchoosek(n, i)*nchoosek(i, k);
   end
end

% 2 forms of definitions for the transform
syms cx1 cx2 cx3 cy1 cy2 cy3 cz1 cz2 cz3 t

% second order Bernstein basis
basis = [1 t t^2];

%% Ctrl points vector base on sequence of X,Y,Z
% this form gives an easier to understand Beta matrix but control point
% sequence are abit weird
Beta1 = kron(eye(dim), delta);

% control points vector
c_vec1 = [cx1;cx2;cx3;cy1;cy2 ;cy3; cz1; cz2; cz3];

% Resultant polynomial expression B is a 3*1 vector with each entry is the
% polynomial expression for each axis
B1 = [basis, zeros(1,3), zeros(1,3);
        zeros(1,3), basis, zeros(1,3);
        zeros(1,3), zeros(1,3), basis] * Beta1 * c_vec1;


%% Ctrl points vector base on sequence of control points
% this form Beta matrix is mixed as well as the natural basis is mixed
% requires some matrix manupulation to get the Beta matrix
Beta2 = augment_array_ndim(delta, dim);

% control points vector
c_vec2 = [cx1; cy1;cz1; cx2; cy2; cz2;cx3; cy3;cz3];

% Resultant polynomial expression B is a 3*1 vector with each entry is the
% polynomial expression for each axis
B2 = [1 0 0 t 0 0 t^2 0 0;
         0 1 0 0 t 0 0 t^2 0;
         0 0 1 0 0 t 0 0 t^2]* Beta2*c_vec2;
     
 %% B1 is equal to B2 since only the control point vector different arrangement
