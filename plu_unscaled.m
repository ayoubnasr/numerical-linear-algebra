function [L, U, P] plu_unscaled(A, tol)

    n=length(A);
    P=eye(n);

    for k = 1 : n-1
        % Find the max pivot:
        [p, idx] = max(Abs(A(k:n, k)));

        % idx is counted starting from k, the start of the subvector k:n
        % So, we compute ik from idx
        ik = idx + k - 1

        if p > tol 
            % Note: If p < tol, then earlier elimination steps subtracted multiples of previous columns from this 
            % column until it vanished below the diagonal. This means our column is linearly dependent on previous columns.
            % This means the matrix is singular => Ax = b no longer has a unique solution.

            % Permute both A and P:
            A([k idx], :) = A([idx k], :);
            P([k idx], :) = A([idx k], :);

            % Perform Gauss Transform:
            % Get Multipliers:
            A(k+1:n, k) = A(k+1:n, k) / A(k,k);
            A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - A(k+1:n, k) * A(k, k+1:n);
    end

    U = triu(A);
    L = eye(n) + tril(A, -1)

