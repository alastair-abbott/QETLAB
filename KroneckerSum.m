%%  KRONECKERSUM   Kronecker sum of two or more matrices
%   This function can be called in one of three ways:
%
%   KRN = KroneckerSum(A,B,C,...) returns the Kronecker sum A \oplus B
%   \oplus C \oplus ...
%
%   KRN = KroneckerSum(A,m) returns the Kronecker sum of A with itself m
%   times.
%
%   KRN = KroneckerSum(A), if A is a cell, returns the Kronecker sum of
%   all matrices contained within A.
%
%   URL: http://www.qetlab.com/KroneckerSum

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 8, 2026

function krn = KroneckerSum(A,varargin)

% If A is a cell, Kronecker sum together each element of the cell.
if(iscell(A))
    for j = 1:length(A)
        I{j} = speye(size(A{j},1));
    end

    krn = Tensor({A{1}, I{2:end}});
    for j = 2:length(A)
        krn = krn + Tensor({I{1:j-1}, A{j}, I{j+1:end}});
    end
    
% If two arguments were received and the second is a scalar, Kronecker sum A with itself that many times.
elseif(nargin == 2 && isscalar(varargin{1}))
    m = varargin{1};
    if(m == 1)
        krn = A;
    else
        I = speye(size(A,1));
        
        % Pre-compute these identity tensors to save time
        for j = 1:m-1
            IC{j} = Tensor(I,j);
        end
    
        krn = kron(A,IC{m-1});
        for j = 2:m-1
            krn = krn + Tensor(IC{j-1},A,IC{m-j});
        end
        krn = krn + kron(IC{m-1},A);
    end

% If two or more arguments were received and the second isn't a scalar, Kronecker sum together all arguments.
else
    krn = KroneckerSum({A,varargin{:}});
end