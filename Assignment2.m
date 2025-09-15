function x = naive_gauss(A)
    % Solve Ax = b using naive Gauss elimination without pivoting.
    %
    % Input:
    %   A - Augmented matrix of size n x (n+1), where the last column is b.
    %
    % Output:
    %   x - Solution vector of size n x 1.
    %
    % Throws:
    %   Error if a pivot is zero (indicating singular matrix).

    n = size(A, 1);
    
    % Forward elimination
    for k = 1:(n-1)
        if A(k,k) == 0
            error('Matrix is singular (zero pivot encountered at step %d)', k);
        end
        for i = (k+1):n
            factor = A(i,k) / A(k,k);
            A(i,k:(n+1)) = A(i,k:(n+1)) - factor * A(k,k:(n+1));
        end
    end
    
    % Back substitution
    x = zeros(n, 1);
    for i = n:-1:1
        if A(i,i) == 0
            error('Matrix is singular (zero pivot in back substitution at row %d)', i);
        end
        sum_ax = A(i,(i+1):n) * x((i+1):n);
        x(i) = (A(i,n+1) - sum_ax) / A(i,i);
    end
end

% Test script for the given system
A_test = [10, 2, -1, 27; -3, -5, 2, -61.5; 1, 1, 6, -21.5];
x = naive_gauss(A_test);
disp('Solution vector [x; y; z]:');
disp(x);