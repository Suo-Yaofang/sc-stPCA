function [weight_W,xi] = fpca_rebuild(H,n,L)

    [V, Lambda] = eig(H);

    eigenvalues = diag(Lambda);
    eigenvalues = real(eigenvalues);
    eigenvalues_prime = max(eigenvalues, 0);
 
    Lambda_prime = diag(eigenvalues_prime);
    
    H_prime = V * Lambda_prime * V';

    sqrt_eigenvalues = sqrt(eigenvalues_prime);

    sqrt_Lambda_prime = diag(sqrt_eigenvalues);
    
    M = V * sqrt_Lambda_prime;

    reconstruction_error = norm(H_prime - M * M', 'fro');
    fprintf('reconstruction_error (Frobenius): %e\n', reconstruction_error);

    [aa, eigvIdx]=sort(eigenvalues,'descend');
    V=V(:,eigvIdx);
   
    K=0;
    sum_xi=0;
    weight_W=zeros(n,L);
    KK=sum(aa > 0);

 if KK>0
    for ci=1:length(aa)
        if aa(ci)>0
            K=K+1;
            xi(K)=norm(V(:,ci)'*M);
            sum_xi=sum_xi+xi(K);
            WW(1:n,1:L,K)=reshape(V(:,ci),n,L);
            weight_W=WW(:,:,K)*xi(K)+weight_W;
        end
    end
 else 
     if aa(1)==0
        disp('There are no positive values, but there are zero values.')
        weight_W=reshape(V(:,length(aa)),n,L);
        xi=1;
        sum_xi=sum(xi);
      else disp('Only negative values.')
        weight_W=reshape(V(:,length(aa)),n,L);
        xi=1;
        sum_xi=sum(xi);
      end
 end
 
    weight_W=weight_W*(1/sum_xi);

end
 