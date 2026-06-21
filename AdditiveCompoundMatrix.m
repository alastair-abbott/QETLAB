%%  ADDITIVECOMPOUNDMATRIX    Computes the r-th additive compound matrix of a given square matrix
%   This function has two required arguments:
%     A: an arbitrary matrix
%     r: a nonnegative integer
%
%   comp = AdditiveCompoundMatrix(A, r) returns the r-th additive compound matrix of the square matrix A
%
%   If r > n (the number of rows in A), then the r-th additive compound matrix of A is the 0x0 matrix.
%   Otherwise, the size of the result is (n choose r) x (n choose r).
%
%   URL: https://qetlab.com/AdditiveCompoundMatrix
%             
%   author: Nathaniel Johnston
%   package: QETLAB
%   last updated: January 8, 2026

function comp = AdditiveCompoundMatrix(A, r)
    m = size(A, 1);
    n = size(A, 2);
    if(m ~= n)
        error('AdditiveCompoundMatrix:InvalidDimensions', 'The matrix must be square.');
    end

    if r > n
        comp = [];
        return
    end

    P = AntisymmetricProjection(m,r,1);
    comp = P'*KroneckerSum(A,r)*P;
end
