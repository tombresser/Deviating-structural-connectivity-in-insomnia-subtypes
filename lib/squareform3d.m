function M = squareform3d(Matrix)
% squareform function for matrices with 3rd dimension
% Behaves the same as squareform
% 3rd dimension is added as rows

    d = size(Matrix,3);
    
    if d==1
         z = size(Matrix,1);
         k = size(Matrix,2);
         n = (1+sqrt(1+8*k))/2;
         M = zeros(n,n,z);
         for i = 1:z
             M(:,:,i) = squareform(Matrix(i,:));
         end
    elseif d>1
        n = size(Matrix,1);
        k = n * (n - 1) / 2 ;
        M = zeros(d,k);

        for i = 1:d
            M(i,:) = squareform(Matrix(:,:,i));
        end
    end
end


%% testing
%A=ones(5,5);
%D=A-diag(diag(A));
%B=cat(3,A, A+1, A+2).*D;
%C=squareform3d(B)
%squareform3d(C)

