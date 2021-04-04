function Ainv = inverse(A)
% Program to find Inverse of a Matrix A using  LU decomposition 
% Method for LU decomposition are Dolittle's and Crout's Algorithms
%==========================================================================
% INPUT:
%==========================================================================
% Coefficient Matrix, A (n-by-n)
%==========================================================================
% OUTPUT
%==========================================================================
% Inverse Matrix, Ainv
%==========================================================================
% Author : Arshad Afzal, India, Email: arshad.afzal@gmail.com 
%==========================================================================

fprintf('\n         ==================================================================================');
fprintf('\n                      Inverse of Matrix A using LU Decomposition');
fprintf('\n         ==================================================================================');
fprintf('\n         ----------------------------------------------------------------------------------');


%==========================================================================
[m,n] = size(A);
%==========================================================================
% Initialization
U = zeros(n,n); L = zeros(n,n);
B = zeros(n,n); Ainv = zeros(n,n);

choice = menu('ALGORITHM','DOLITTLE','CROUT');

if choice == 1
    

    for i = 1:n
        L(i,i) = 1;
    end

    %======================================================================
    % Dolittle's Algorithm
    %======================================================================

    for i = 1:n
        for j = i:n
            sum = 0;

            for k = 1:i-1
                sum = sum + L(i,k)*U(k,j);
            end

                U(i,j) = A(i,j)-sum;

             sum = 0;
             for k = 1:i-1
                sum = sum + L(j,k)*U(k,i);
             end
                L(j,i) = (A(j,i)-sum)/U(i,i);

         end  
    end

else
    
    
for i = 1:n
 U(i,i) = 1;
end

%==========================================================================
% Crout's Algorithm
%==========================================================================

for j = 1:n
    
    for i = j:n
        sum = 0;
        
        for k = 1:j-1
            sum = sum + L(i,k)*U(k,j);
        end
        
        L(i,j) = A(i,j)-sum;
            
        sum = 0;
        for k = 1:j-1
            sum = sum + L(j,k)*U(k,i);
        end
        
        for k = j+1:n
            U(j,k) = (A(j,k)-sum)/L(j,j);
        end
            
     end  
end 

end

%==========================================================================
% Main Program for calculating inverse of matrix
%==========================================================================
% Forward elimination, solve LB = I

b = eye(n);
for i = 1:m
    B(1,i) = b(1,i)/L(1,1);
    for k = 2:m
        sum = 0;
            for j = k-1:-1:1
              sum = sum + L(k,j)*B(j,i);
            end
        B(k,i) = (b(k,i)- sum)/L(k,k);
    end
end

% Backward substitution, solve U*Ainv = B
for i = 1:m
      Ainv(m,i) = B(m,i)/U(m,m);
      for k = m-1:-1:1
        sum = 0;
        for j = k+1:m
           sum = sum + U(k,j)*Ainv(j,i);
        end
        Ainv(k,i) = (B(k,i)- sum)/U(k,k);
      end
end

%==========================================================================
% Inverse of A
fprintf('\n  Inverse OF A :\n');
%==========================================================================
end
    