clc
clear

L = [  2            0            0            0
    -1          1.5            0            0
    0           -1          4/3            0
    0            0           -1         1.25];
U = [1         -0.5            0            0
    0            1         -2/3            0
    0            0            1        -0.75
    0            0            0            1];

A = L * U;

n = length(L);
x = zeros(n,1);
c = zeros(n,1);
d = zeros(n,1);

inverse  = zeros(n);


%% Get inv(L)
for k = 1 : n
    % It seems more natural to assign the rows of the identity matrix inside
    % the loop, since there is nothing special about c(1) 
    c(k) = 1;
    if (k > 1)
        c(k - 1) = 0;
    end
    
    % "i" needs to start from 1
    % Even if the inverse L is still a lower triagular matrix, it is more
    % straightforward and less error-prone to start from 1.
    for i = 1 : n
        sum = 0;
        
        for j = 1 : i - 1
            sum = sum + L(i, j) * d(j);
        end
        d(i) = (c(i) - sum) / L(i, i);        
    end
    
    L_inverse(:, k) = d;
    
    %% The code below has no change
    x(n) = d(n) / U(n, n);
    for i = n - 1 : -1 : 1
        sum = 0;
        for j = i + 1 : n
            sum = sum + U(i, j) * x(j);
        end
        x(i) = [d(i) - sum] / U(i, i);
    end
    inverse(:, k) = x;
end