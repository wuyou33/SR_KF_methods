% ------------------------------------------------------------------- 
% Extended Square-root Kalman Filter 
%      Type: Covariance filtering
%    Method: Cholesky-based implementation with lower triangular factors
%      From: Two stages, a posteriori form
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References:
%   1. Park P., Kailath T. (1995) New square-root algorithms for Kalman 
%      filtering, IEEE Trans. Automat. Contr. 40(5):895-899.
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
function [neg_LLF,hatX,hatDP] = Riccati_KF_eSRlow(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 
          
        [m,n]  = size(H);                % dimensions
       N_total = size(measurements,2);   % number of measurements
          hatX = zeros(n,N_total+1);     % prelocate for efficiency
         hatDP = zeros(n,N_total+1);     % prelocate for efficiency

        if isdiag(Q), sqrtQ = diag(sqrt(diag(Q))); else  sqrtQ = chol(Q,'lower'); end; clear Q; % Cholesky factorization Q = A*A'
        if isdiag(R), sqrtR = diag(sqrt(diag(R))); else  sqrtR = chol(R,'lower'); end; clear R; 
        if isdiag(P), sqrtP = diag(sqrt(diag(P))); else  sqrtP = chol(P,'lower'); end;   
            PX = sqrtP\X;            % initials for the filter 
        
  neg_LLF = 1/2*m*log(2*pi)*N_total; % initial value for the neg Log LF
hatX(:,1) = X; hatDP(:,1) = diag(P); % save initials at the first entry
for k = 1:N_total                 
   [PX,sqrtP]                  = esrcf_predict(PX,sqrtP,F,G,sqrtQ); 
   [PX,sqrtP,norm_ek,sqrt_Rek] = esrcf_update(PX,sqrtP,measurements(:,k),H,sqrtR);
   
   neg_LLF = neg_LLF+log(det(sqrt_Rek))+1/2*(norm_ek')*norm_ek;
   hatX(:,k+1) = sqrtP*PX; hatDP(:,k+1) = diag(sqrtP*sqrtP'); % save estimates   
 end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Time update: a priori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PX,sqrtP] = esrcf_predict(PX,sqrtP,F,G,sqrtQ)
     [n,q]      = size(G);
     PreArray   = [F*sqrtP,  G*sqrtQ;
                   PX',       zeros(1,q);];
   % --- triangularize the first two (block) columns only ---          
   [Oth,PostArray] = qr(PreArray(1:n,:)');  
       PostLastRow = Oth'*PreArray(end,:)'; 
       PostArray   = PostArray'; PostLastRow = PostLastRow';

   % -- read-off the resulted quantities -------------------
    sqrtP        = PostArray(1:n,1:n);    % Predicted factor of P        
    PX           = PostLastRow(end,1:n)';  % Predicted PX element        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement update: a posteriori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PX,P_sqrt,bar_ek,Rek_sqrt] = esrcf_update(PX,P_sqrt,z,H,R_sqrt)
    [m,n]     = size(H);

   % --- triangularize the first two (block) rows only ---          
        PreArray  = [R_sqrt,                        H*P_sqrt; 
                     zeros(n,m),                    P_sqrt;
                     -(R_sqrt\z)',  PX';];
    [Oth,PostArray]  = qr(PreArray(1:end-1,:)');
         PostLastRow = Oth'*PreArray(end,:)'; 
         PostArray   = PostArray'; PostLastRow = PostLastRow';

   % -- read-off the resulted quantities -------------------
        Rek_sqrt  = PostArray(1:m,1:m);
         % bar_K_pk = PostArray(m+1:m+n,1:m);      
           P_sqrt = PostArray(m+1:m+n,1+m:m+n);  
           bar_ek = -PostLastRow(1,1:m)';      
               PX = PostLastRow(1,m+1:m+n)';
end
