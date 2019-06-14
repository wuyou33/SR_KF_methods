% ------------------------------------------------------------------- 
% SVD-based Maximum Correntropy Criterion Kalman Filter (MCC-KF)
%      Type: Covariance filtering
%    Method: SVD-based implementation
%      From: Two stages, a posteriori form, Cholesky decomposition for noise covariances
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References:
%   L. Wang, G. Libert, and P. Manneback (1992) "Kalman Filter Algorithm
%   based on Singular Value Decomposition", Proceedings of the 31st
%   Conference on Decision and Control. Tuczon, AZ, USA: IEEE, 1992,
%   pp. 1224–1229. See section "Algorithmic detailes". eq. (17), (22), (23) 
%   DOI: 10.1109/CDC.1992.371522 
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
function [neg_LLF,hatX,hatDP] = Riccati_KF_SVDSR(matrices,initials_filter,measurements)
       neg_LLF = NaN;                       % This method is not able to compute the log LF directly
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 
          
        [m,n]  = size(H);                % dimensions
       N_total = size(measurements,2);   % number of measurements
          hatX = zeros(n,N_total+1);     % prelocate for efficiency
         hatDP = zeros(n,N_total+1);     % prelocate for efficiency

        if isdiag(Q), Q_sqrt = diag(sqrt(diag(Q))); else  Q_sqrt = chol(Q,'upper'); end; clear Q; % Cholesky factorization Q = A'*A
        GQsqrt = G*Q_sqrt';                % compute once the new G*Q^{1/2}
        if isdiag(R), R_sqrt = diag(sqrt(diag(R))); else  R_sqrt = chol(R,'upper'); end; clear R; 
        [~,DP,QP] = svd(P); DPsqrt = DP.^(1/2); clear DP; 

hatX(:,1)  = X; hatDP(:,1) = diag(QP*DPsqrt^(2)*QP'); % save initials at the first entry
for k = 1:N_total      
   [X,DPsqrt,QP] = svd_sr_predict(X,DPsqrt,QP,F,GQsqrt); 
   [X,DPsqrt,QP] = svd_sr_update(X,DPsqrt,QP,measurements(:,k),H,R_sqrt);
 
   hatX(:,k+1) = X; hatDP(:,k+1) = diag(QP*DPsqrt^(2)*QP'); % save estimates  
 end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Time update: a priori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,DPsqrt,QP] = svd_sr_predict(X,DPsqrt,QP,F,G)
     [n,~]    = size(G);
     PreArray = [DPsqrt*QP'*F'; G';];
          
     [~,S,QP] = svd(PreArray); % Predicted SVD factors of P        
       DPsqrt = S(1:n,1:end);  % Predicted SVD factors of P     
            X = F*X;           % Predicted state estimate   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement update: a posteriori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,DPsqrt,QP] = svd_sr_update(X,DPsqrt,QP,z,H,R_sqrt)
      [~,n]     = size(H);
     PreArray   = [R_sqrt'\H*QP; inv(DPsqrt);];
     [~,Sm,Vm]  = svd(PreArray);            
         DPsqrt = inv(Sm(1:n,1:end));       
             QP = QP*Vm;                     % Filtered SVD factors of P                               
 
     Gain = (QP*DPsqrt*DPsqrt*QP')*H'/R_sqrt'/R_sqrt; % Gain
        X = X + Gain*(z-H*X);                           % Filtered state estimate        
end
     