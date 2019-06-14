% ------------------------------------------------------------------- 
% SVD-based Maximum Correntropy Criterion Kalman Filter (MCC-KF)
%      Type: Covariance filtering
%    Method: SVD-based implementation, "economy size"
%      From: Two stages, a posteriori form 
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References 
%   Kulikova M.V., Tsyganova J.V. (2017) "Improved discrete-time Kalman 
%   filtering within singular value decomposition". 
%   IET Control Theory & Applications, 11(15): 2412-2418
%   DOI:   10.1049/iet-cta.2016.1282 
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

function [neg_LLF,hatX,hatDP] = Riccati_KF_SVDe(matrices,initials_filter,measurements)
   [F,G,Qsys,H,R] = deal(matrices{:});       % get system matrices
          [X,P] = deal(initials_filter{:});  % initials for the filter 
   
        [m,n]  = size(H);                % dimensions
       N_total = size(measurements,2);   % number of measurements
          hatX = zeros(n,N_total+1);     % prelocate for efficiency
         hatDP = zeros(n,N_total+1);     % prelocate for efficiency

     [QQ,DQ] = svd(Qsys);      % SVD for process noise covariance
      Q_sqrt = QQ*DQ.^(1/2);   % compute once 
     [RQ,DR] = svd(R);         % SVD for measurement noise covariance
      R_sqrt = RQ*DR^(1/2);    % compute once
     [QP,DP] = svd(P); DPsqrt = DP.^(1/2); clear DP; 

neg_LLF = 1/2*m*log(2*pi)*N_total; % initial value for the neg Log LF
hatX(:,1)  = X; hatDP(:,1) = diag(QP*DPsqrt^(2)*QP'); % save initials at the first entry
for k = 1:N_total                       
   [X,DPsqrt,QP]              = svd_predict(X,DPsqrt,QP,F,G,Q_sqrt); 
   [X,DPsqrt,QP,ek,DRek,QRek] = svd_update(X,DPsqrt,QP,measurements(:,k),H,R_sqrt);
   
   neg_LLF = neg_LLF+1/2*log(det(DRek))+1/2*ek'*QRek/DRek*QRek'*ek; 
   hatX(:,k+1)= X; hatDP(:,k+1) = diag(QP*DPsqrt^(2)*QP'); % save estimates  
 end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Time update: a priori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,DPsqrt,QP] = svd_predict(X,DPsqrt,QP,F,G,Qsqrt)
      PreArray   = [F*QP*DPsqrt, G*Qsqrt];
     [QP,DPsqrt] = svd(PreArray,'econ'); 
               X = F*X;                  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement update: a posteriori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,DPsqrt,QP,residual,DRek,QRek] = svd_update(X,DPsqrt,QP,z,H,sqrtR)
     PreArray    = [H*QP*DPsqrt, sqrtR];
     [QRek,DRek_sqrt] = svd(PreArray,'econ');          % SVD factors of R_{e,k} 
           DRek = DRek_sqrt.^2;         
    
     residual  =  z-H*X;                               % residual
          Gain = QP*(DPsqrt^2)*QP'*H'*QRek/DRek*QRek'; % Kalman gain
             X = X + Gain*residual;                    % Filtered estimate 

       PreArray  = [(eye(size(QP,1)) - Gain*H)*QP*DPsqrt, Gain*sqrtR];  % SVD of P
     [QP,DPsqrt] = svd(PreArray,'econ'); 
end
