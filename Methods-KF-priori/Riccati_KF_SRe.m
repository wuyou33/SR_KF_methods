% ------------------------------------------------------------------- 
% extended Square-root Kalman Filter
%      Type: Covariance filtering
%    Method: Cholesky-based implementation with upper triangular factors
%      From: One stage (condensed form), a priori form, extended form
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References: 
%   1. This is the upper variant of Algorithm III.2 from
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
function [neg_LLF,hatX,hatDP] = Riccati_KF_SRe(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 

        [n,q]  = size(G);                % dimensions
        [m,n]  = size(H);                % dimensions
       N_total = size(measurements,2);   % number of measurements
          hatX = zeros(n,N_total+1);     % prelocate for efficiency
         hatDP = zeros(n,N_total+1);     % prelocate for efficiency

        if isdiag(Q), Q_sqrt = diag(sqrt(diag(Q))); else  Q_sqrt = chol(Q,'upper'); end; clear Q; % Cholesky factorization Q = A*A'
        if isdiag(R), R_sqrt = diag(sqrt(diag(R))); else  R_sqrt = chol(R,'upper'); end; clear R; 
        if isdiag(P), P_sqrt = diag(sqrt(diag(P))); else  P_sqrt = chol(P,'upper'); end;    
            PX = P_sqrt'\X;              % initials for the filter 

neg_LLF = 1/2*m*log(2*pi)*N_total;   % initial value for the neg Log LF
hatX(:,1) = X; hatDP(:,1) = diag(P); % save initials at the first entry
   for k = 1:N_total
   % --- triangularize the first two (block) columns only ---          
    PreArray  = [R_sqrt,      zeros(m,n),   -(R_sqrt')\measurements(:,k); 
                 P_sqrt*H',   P_sqrt*F',     PX;
                 zeros(q,m),  Q_sqrt*G',     zeros(q,1)];

   [Orth,PostArray1] = qr(PreArray(:,1:m+n));  % Q is any orthogonal rotation that upper-triangularizes the first
                                               % two block columns of the pre-array
   PostArray2 = Orth'*PreArray(:,end);         % transform the lat row as well 
   PostArray  = [PostArray1, PostArray2];      % the full post-array

   % -- read-off the resulted quantities -------------------
        Rek_sqrt  = PostArray(1:m,1:m);
           P_sqrt = PostArray(m+1:m+n,m+1:m+n);  
           bar_ek = -PostArray(1:m,end);      
               PX = PostArray(m+1:m+n,end);
        
    neg_LLF = neg_LLF+log(det(Rek_sqrt))+1/2*bar_ek'*bar_ek; 
    hatX(:,k+1) = P_sqrt'*PX; hatDP(:,k+1) = diag(P_sqrt'*P_sqrt); 
 end;
end
