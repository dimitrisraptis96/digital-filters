% STEEPESTTDESCENT - Perform steepest descent to compute Wiener coefficients
%   
% SYNTAX
%
%   [W, Wt] = STEEPESTDESCENT( MU, N, R, P )
%
% INPUT
%
%   MU  Step value                              [scalar]
%   N   Iterations                              [scalar]
%   R   Auto-correlation matrix                 [m-by-m]
%   P   Cross-correlation vector                [m-vector]
%
% OUTPUT
%
%   W   Final Wiener coefficients               [m-vector]
%   WT  Evoluation/adaptation of Wiener coeff   [m-by-n]
%
% DESCRIPTION
%
%   STEEPESTDESCENT 
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      
%

function [ W, Wt ] = steepestDescent( mu, n, R, p )
    W = [0 ; 0];        % initial Wiener coefficients 
    Wt = zeros(2,n);
    Wt(:,1) = W;

    for k=2:n
         W = W + mu * (p - R*W);
         Wt(:,k) = W;
    end
end