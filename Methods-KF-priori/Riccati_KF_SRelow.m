% ------------------------------------------------------------------- 
% extended Square-root Kalman Filter
%      Type: Covariance filtering
%    Method: Cholesky-based implementation with lower triangular factors
%      From: One stage (condensed form), a priori form, extended form
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References: 
%   1. This is the implementation of Algorithm III.2 from
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
function [neg_LLF,hatX,hatDP] = Riccati_KF_SRelow(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 

        [n,q]  = size(G);                % dimensions
        [m,n]  = size(H);                % dimensions
       N_total = size(measurements,2);   % number of measurements
          hatX = zeros(n,N_total+1);     % prelocate for efficiency
         hatDP = zeros(n,N_total+1);     % prelocate for efficiency

        if isdiag(Q), Q_sqrt = diag(sqrt(diag(Q))); else  Q_sqrt = chol(Q,'lower'); end; clear Q; % Cholesky factorization Q = A*A'
        if isdiag(R), R_sqrt = diag(sqrt(diag(R))); else  R_sqrt = chol(R,'lower'); end; clear R; 
        if isdiag(P), P_sqrt = diag(sqrt(diag(P))); else  P_sqrt = chol(P,'lower'); end;    
            PX = P_sqrt\X;               % initials for the filter 

neg_LLF = 1/2*m*log(2*pi)*N_total;  % initial value for the neg Log LF
hatX(:,1) = X; hatDP(:,1) = diag(P); % save initials at the first entry
   for k = 1:N_total
   % --- triangularize the first two (block) rows only ---          
        PreArray  = [R_sqrt,                        H*P_sqrt, zeros(m,q); 
                     zeros(n,m),                    F*P_sqrt, G*Q_sqrt;
                     -(R_sqrt\measurements(:,k))',  PX',      zeros(1,q);];
    [Oth,PostArray]  = qr(PreArray(1:end-1,:)');
         PostLastRow = Oth'*PreArray(end,:)'; 
         PostArray   = PostArray'; PostLastRow = PostLastRow';

   % -- read-off the resulted quantities -------------------
        Rek_sqrt  = PostArray(1:m,1:m);
         bar_K_pk = PostArray(m+1:m+n,1:m);      
           P_sqrt = PostArray(m+1:m+n,1+m:m+n);  
           bar_ek = -PostLastRow(1,1:m)';      
               PX = PostLastRow(1,m+1:m+n)';
                X = P_sqrt*PX;    
        
    neg_LLF = neg_LLF+log(det(Rek_sqrt))+1/2*bar_ek'*bar_ek; 
    hatX(:,k+1) = X; hatDP(:,k+1) = diag(P_sqrt*P_sqrt'); 
 end;
end
