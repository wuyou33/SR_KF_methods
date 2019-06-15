% Square-root Kalman Filter
%      Type: Covariance filtering
%    Method: Cholesky-based implementation with upper triangular factors
%      From: One stage (condensed form), a priori form
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References: 
%   1. This is the upper triangular variant of Algorithm III.1 from
%      Park P., Kailath T. (1995) New square-root algorithms for Kalman 
%      filtering. IEEE Transactions on Automatic Control, 40(5), 895-899.
%      DOI: 10.1109/9.384225 
% ------------------------------------------------------------------- 
% Input:
%     matrices        - system matrices F,H,Q etc
%     initials_filter - initials x0,P0
%     measurements    - measurements (where y(t_k) is the k-th column)
% Output:
%     neg_LLF     - negative log LF
%     hatX        - filtered estimate (history) 
%     hatDP       - diag of the filtered error covariance (history)
% ------------------------------------------------------------------- 
function [neg_LLF,hatX,hatDP] = Riccati_KF_SR(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 

        [n,q]  = size(G);                % dimensions
        [m,n]  = size(H);                % dimensions
       N_total = size(measurements,2);   % number of measurements
          hatX = zeros(n,N_total+1);     % prelocate for efficiency
         hatDP = zeros(n,N_total+1);     % prelocate for efficiency

        if isdiag(Q), Q_sqrt = diag(sqrt(diag(Q))); else  Q_sqrt = chol(Q,'upper'); end; clear Q; % Cholesky factorization Q = A'*A
        if isdiag(R), R_sqrt = diag(sqrt(diag(R))); else  R_sqrt = chol(R,'upper'); end; clear R; 
        if isdiag(P), P_sqrt = diag(sqrt(diag(P))); else  P_sqrt = chol(P,'upper'); end;          % initials for the filter 
        
neg_LLF = 1/2*m*log(2*pi)*N_total;  % initial value for the neg Log LF
hatX(:,1) = X; hatDP(:,1) = diag(P); % save initials at the first entry
   for k = 1:N_total
             ek   = measurements(:,k) - H*X;   % residual
   % --- triangularize  ---          
        PreArray  = [R_sqrt,     zeros(m,n); 
                     P_sqrt*H',   P_sqrt*F';
                     zeros(q,m),  Q_sqrt*G';];
   [~,PostArray]  = qr(PreArray);

   % -- read-off the resulted quantities -------------------
        Rek_sqrt  = PostArray(1:m,1:m);
         bar_K_pk = PostArray(1:m,m+1:m+n)';      
           P_sqrt = PostArray(m+1:m+n,1+m:m+n);      
                X = F*X + bar_K_pk/Rek_sqrt'*ek;    
        
    neg_LLF = neg_LLF+log(det(Rek_sqrt))+1/2*ek'/Rek_sqrt/Rek_sqrt'*ek; 
    hatX(:,k+1) = X; hatDP(:,k+1) = diag(P_sqrt'*P_sqrt); 
 end;
end
