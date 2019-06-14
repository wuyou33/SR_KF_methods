% ------------------------------------------------------------------- 
% classical Kalman Filter 
%      Type: Covariance filtering
%    Method: Conventional implementation
%      From: One stage, a priori form
% Recursion: Chandrasekhar-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
%  1. The Algorithm is based on recursion (67) in the following paper (see eqs. (13'), (23), (24), (16)): 
%     Morf  M., Sidhu G. and  Kailath T. (1974) 
%     Some new algorithms for recursive estimation in constant, linear, discrete-time systems, 
%     IEEE Trans. Automat. Contr., vol. 19, no. 4, pp. 315-323, Aug. 1974.
%     DOI:   10.1109/TAC.1974.1100576 
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
function [neg_LLF,predX,predDP] = Chandrasekhar_KF3(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 
   
        [m,n]  = size(H);                 % dimensions
       N_total = size(measurements,2);    % number of measurements
         predX = zeros(n,N_total+1);      % prelocate for efficiency
        predDP = zeros(n,N_total+1);      % prelocate for efficiency
       neg_LLF = 1/2*m*log(2*pi)*N_total; % set initial value for the neg Log LF

predX(:,1) = X; predDP(:,1) = diag(P);    % save initials at the first entry
% --- Filter initials --------------------------------------------------------
invRek  = inv(R + H*P*H');                 % inverse of residual covariance 
    K_k = F*P*H';                          % initial gain
 Delta0 = F*P*F'+G*Q*G'-K_k*invRek*K_k'-P; % initial difference
 [L,M] = ldl(Delta0); alpha = rank(Delta0);       % initial decomposition deltaP = L*M*L' and displacemant rank alpha
  inds = find(diag(M)==0); M(inds,:) = []; M(:,inds) = []; L(:,inds)=[]; % cut the zero parts for low rank approximation
  if (alpha<size(M,1))&(isdiag(M)) % if alpha<size(M,1), then we should improve low-rank approximation, e.g. by deleting the same diag values 
     [~,Inds_no_duplicate] = intersect(diag(M),unique(diag(M)));  M = M(Inds_no_duplicate,Inds_no_duplicate); L = L(:,Inds_no_duplicate);
  end; 
% --- Filtering------ --------------------------------------------------------
invM = inv(M);  % this implementation updates the inverse of Rek and M matrices
for k = 1:N_total  
     ek = measurements(:,k) - H*X;  % residual
      X = F*X + K_k*invRek*ek;      % a priori estimate X
     neg_LLF = neg_LLF+1/2*log(1/det(invRek))+1/2*ek'*invRek*ek; 
     predX(:,k+1) = X; P = P + L*M*L';  predDP(:,k+1) = diag(P);            % save estimates

     %--- recursion for new L and invM ----------------
         L_next = (F-K_k*invRek*H)*L;
      invM_next = invM + L'*H'*invRek*H*L; M_next = inv(invM_next);         % inverse M of size alpha, once
         invRek = invRek - invRek*H*L*M_next*L'*H'*invRek; 
           K_k  = K_k + F*L*M*L'*H'; 
              L = L_next; invM = invM_next; M = M_next;                     % new iteration starts
end;
end
