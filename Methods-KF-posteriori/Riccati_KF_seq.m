% ------------------------------------------------------------------- 
% Sequential Kalman Filter (component-wise measurement update)
%      Type: Covariance filtering
%    Method: Conventional implementation
%      From: Two stages, a posteriori form, sequential measurement update
% Recursion: Riccati-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References:
%   1. see Chapter 7 and Section 7.4.5 in 
%      Grewal, M.S., Andrews, A.P. (2015) Kalman filtering: theory and practice using MATLAB 
%      Prentice-Hall, New Jersey, 4th edn. 
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
function [neg_LLF,hatX,hatDP] = Riccati_KF_seq(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 
   
        [m,n]  = size(H);                 % dimensions
       N_total = size(measurements,2);    % number of measurements
         hatX  = zeros(n,N_total+1);      % prelocate for efficiency
         hatDP = zeros(n,N_total+1);      % prelocate for efficiency
       neg_LLF = 1/2*m*log(2*pi)*N_total; % set initial value for the neg Log LF
     hatX(:,1) = X; hatDP(:,1) = diag(P); % save initials at the first entry
fl=0; if isdiag(R)==0, [H,R,measurements,fl] = DecorrelateR(H,R,measurements); end; 

for k = 1:N_total              
              [X,P]  = kf_predict(X,P,F,G,Q); 
      for i=1:m
       [X,P,ei,alph] = kf_update(X,P,measurements(i,k),H(i,:),R(i,i)); 
       if fl==0, neg_LLF = neg_LLF+1/2*log(abs(alph))+1/2*(ei^2)/alph; % if measurement model is changed, then the LLF should be also
                 else neg_LLF = NaN; end;                              % transformed. In this case, we just set it to NaN, here.
      end;
      hatX(:,k+1)  = X; hatDP(:,k+1) = diag(P); 
 end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Time update: a priori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P] = kf_predict(X,P,F,G,Q)
    X = F*X;             % Predicted state estimate  
    P = F*P*F' + G*Q*G'; % Predicted error covariance 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement update: a posteriori estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,P,residual,cov_residual] = kf_update(X,P,z,H,R)
  residual     = z - H*X;            % residual
  cov_residual = R + H*P*H';         % residual covariance matrix 
  Kalman_gain  = P*H'/cov_residual;  % Gain matrix

  X = X + Kalman_gain*residual;      % Filtered state estimate
  P = P-Kalman_gain*H*P;             % Filtered error covariance
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Measurement de-correlation:  based on Cholesky d.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Hsys,Rsys,meas,flag_tr] = DecorrelateR(Hsys,Rsys,meas)
flag_tr = 0; 
if all(all(Rsys==diag(diag(Rsys))))==0 % if R is not diagonal 
   L = chol(Rsys,'lower');            % (L*L')
   Hsys = L\Hsys; 
   Rsys = eye(size(Rsys,2)); 
   meas = L\meas;
   flag_tr  = 1;
end;
end
