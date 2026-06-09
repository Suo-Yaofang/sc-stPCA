function mu_min = compute_mu_min(X, L, use_iterative)

[n, m] = size(X);
if nargin < 2, L = 3; end          
if nargin < 3, use_iterative = []; end 

if m < 2
    error('Data matrix X must have at least 2 columns');
end
if L < 2
    error('Number of blocks L must be at least 2');
end

% Compute problem dimension
nL = n * L;
max_direct_size = 1e4; % Threshold for direct method (10,000)

% Check if X is full rank
if rank(X) < n
    warning('X is not full rank, results may be inaccurate');
end

P = X(:, 2:m);    % n x (m-1)
Q = X(:, 1:m-1);  % n x (m-1)

J_R = [speye(L-1); sparse(1, L-1)];  % L x (L-1)
J_S = [sparse(1, L-1); speye(L-1)];  % L x (L-1)

% Decide solution method
if isempty(use_iterative)
    use_iterative = (nL > max_direct_size);
end

if use_iterative
    fprintf('Using iterative method (nL = %d)...\n', nL);
    
    opts = struct();
    opts.tol = 1e-6;               
    opts.issym = true;             
    opts.isreal = true;           
    opts.maxit = 300;              
 
    try
        diag_blocks = cell(1, L);
        XXT = X * X';
        for i = 1:L
            diag_blocks{i} = XXT;
        end
        K_B_diag = blkdiag(diag_blocks{:});
        
        precond = ichol(K_B_diag, struct('type', 'ict', 'droptol', 1e-3));
        opts.M = @(x) precond \ (precond' \ x);
        fprintf('Using incomplete Cholesky preconditioner\n');
    catch ME
        fprintf('Preconditioner failed: %s\nUsing diagonal preconditioner\n', ME.message);
        diag_vals = repmat(sum(X.^2, 2), [L, 1]); 
        diag_vals = diag_vals + 1e-10 * max(diag_vals);
        opts.M = @(x) x ./ diag_vals;
    end
    
    [V, D] = eigs(@K_A_times_v, @K_B_times_v, nL, 1, 'smallestreal', opts);
    
else

    fprintf('Using direct method (nL = %d)...\n', nL);
   
    M = kron(J_R', P') - kron(J_S', Q');  % (m-1)(L-1) x nL

    K_A = M' * M;  % nL x nL
    
    %  K_B = I_L ⊗ (X X')
    XXT = X * X';
    if nL > 1000
        K_B = kron(speye(L), sparse(XXT));
    else
        K_B = kron(eye(L), XXT);
    end
    
    opts = struct();
    opts.tol = 1e-10;
    opts.issym = true;
    opts.isreal = true;
    [V, D] = eigs(K_A, K_B, 1, 'smallestreal');
end

function y = K_A_times_v(v)
        Z = reshape(v, [n, L]);
        
        R = Z * J_R; 
        S = Z * J_S; 
        
        diff = P' * R - Q' * S;
        
        temp1 = P * diff;
        temp2 = Q * diff;
        y1 = temp1 * J_R';
        y2 = temp2 * J_S';
        y = y1 - y2;
        y = y(:);  
end

 function y = K_B_times_v(v)
        V_mat = reshape(v, [n, L]);
        
        temp = X' * V_mat;  % m x L
        Y_mat = X * temp;   % n x L 
        
        y = Y_mat(:);
 end

mu_min = real(D(1, 1));
end

% rank(traindata)
% mu=compute_mu_min(traindata,10)
% lambda_min=1/(1+mu)
% a1=1-lambda_min