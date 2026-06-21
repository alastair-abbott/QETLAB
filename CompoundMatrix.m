%%  COMPOUNDMATRIX    Computes the r-th (multiplicative) compound matrix of a given matrix
%   This function has two required arguments:
%     A: an arbitrary matrix
%     r: a nonnegative integer
%
%   comp = CompoundMatrix(A, r) returns the r-th compound matrix of A
%
%   If r > min{m, n}, then the r-th compound matrix of A is the 0x0 matrix.
%   Otherwise, the size of the result is (m choose r) x (n choose r).
%
%   URL: https://qetlab.com/CompoundMatrix
%             
%   author: Benjamin Talbot and Nathaniel Johnston
%   package: QETLAB
%   last updated: January 8, 2026

function comp = CompoundMatrix(A, r)
    m = size(A, 1);
    n = size(A, 2);
    if r > min(m, n)
        comp = [];
        return
    end

    % If A is symbolic, it is faster to compute its compound matrix via the
    % definition.
    if(isa(A,'sym')) % If A is symbolic, we have to make its compound matrix symbolic too.
        rows = nchoosek(1:m, r);
        cols = nchoosek(1:n, r);
        dims = [size(rows, 1), size(cols, 1)];
        comp = sym(zeros(dims));
        for i = 1:dims(1)
            for j = 1:dims(2)
                comp(i, j) = det(A(rows(i, :), cols(j, :)));
            end
        end

    % If A is not symbolic then it is quicker to just compute the compound
    % matrix via tensor products and projections.
    else
        P = AntisymmetricProjection(m,r,1);
        if(m == n)
            comp = P'*Tensor(A,r)*P;
        else
            Q = AntisymmetricProjection(n,r,1);
            comp = P'*Tensor(A,r)*Q;
        end
    end
end
