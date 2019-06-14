% ------------------------------------------------------------------- 
% classical Kalman Filter 
%      Type: Covariance filtering
%    Method: Conventional implementation
%      From: One stage, a priori form
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
%
% References:
%     D. Simon (2006) Optimal state estimation: Kalman, H-infinity, 
%     and nonlinear approaches, John Wiley & Sons.
% ------------------------------------------------------------------- 
% Input:
%     matrices        - system matrices F,H,Q etc
%     initials_filter - initials x0,P0
%     measurements    - measurements (where y(t_k) is the k-th column)
% Output:
%     neg_LLF     - negative log LF
%     predX       - a priori estimates (history) 
%     predDP      - diag of the predicted error covariance (history)
% ------------------------------------------------------------------- 
function [neg_LLF,predX,predDP] = Riccati_KF(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 
   
        [m,n]  = size(H);                 % dimensions
       N_total = size(measurements,2);    % number of measurements
         predX = zeros(n,N_total+1);      % prelocate for efficiency
        predDP = zeros(n,N_total+1);      % prelocate for efficiency
       neg_LLF = 1/2*m*log(2*pi)*N_total; % set initial value for the neg Log LF
 predX(:,1)  = X; predDP(:,1) = diag(P);  % save initials at the first entry
for k = 1:N_total               
    ek   = measurements(:,k) - H*X;   % residual
    Rek  = R + H*P*H';                % residual covariance matrix 
    K_pk = F*P*H'/Rek;                % Kalman Gain       
    
    X = F*X + K_pk*ek;                    % predicted state; formula (5.24) in the cited book; 
    P = F*P*F' + G*Q*G' - K_pk*Rek*K_pk'; % predicted error covariance; formula (5.26) in the cited book; 
        
    neg_LLF = neg_LLF+1/2*log(det(Rek))+1/2*ek'/Rek*ek; 
    predX(:,k+1) = X; predDP(:,k+1) = diag(P); 
 end;
end
