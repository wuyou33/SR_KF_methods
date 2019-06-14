% ------------------------------------------------------------------- 
% classical Kalman Filter 
%      Type: Covariance filtering
%    Method: Conventional implementation
%      From: One stage, a priori form
% Recursion: Chandrasekhar-type underlying recursion
%   Authors: Maria Kulikova: maria dot kulikova at ist dot utl dot pt 
% ------------------------------------------------------------------- 
% References:
%  1. The Algorithm is based on recursion (67) in the following paper (see eqs. (13'), (26), (27), (16)): 
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
function [neg_LLF,predX,predDP] = Chandrasekhar_KF4(matrices,initials_filter,measurements)
   [F,G,Q,H,R] = deal(matrices{:});         % get system matrices
         [X,P] = deal(initials_filter{:});  % initials for the filter 
   
        [m,n]  = size(H);                 % dimensions
       N_total = size(measurements,2);    % number of measurements
         predX = zeros(n,N_total+1);      % prelocate for efficiency
        predDP = zeros(n,N_total+1);      % prelocate for efficiency
       neg_LLF = 1/2*m*log(2*pi)*N_total; % set initial value for the neg Log LF

predX(:,1) = X; predDP(:,1) = diag(P);    % save initials at the first entry
% --- Filter initials --------------------------------------------------------
Rek  = R + H*P*H';                       % initial residual covariance 
K_pk = F*P*H'/Rek;                       % initial gain
Delta0 = F*P*F'+G*Q*G'-K_pk*Rek*K_pk'-P; % initial difference
 [L,M] = ldl(Delta0); alpha = rank(Delta0); % initial decomposition deltaP = L*M*L' and displacemant rank alpha
  inds = find(diag(M)==0); M(inds,:) = []; M(:,inds) = []; L(:,inds)=[]; % cut the zero parts for low rank approximation
  if (alpha<size(M,1))&(isdiag(M)) % if alpha<size(M,1), then we should improve low-rank approximation, e.g. by deleting the same diag values 
     [~,Inds_no_duplicate] = intersect(diag(M),unique(diag(M)));  M = M(Inds_no_duplicate,Inds_no_duplicate); L = L(:,Inds_no_duplicate);
  end; 
% --- Filtering------ --------------------------------------------------------
invM = inv(M);   % this implementation updates the inverse of M 
for k = 1:N_total  
     ek = measurements(:,k) - H*X;     % residual
      X = F*X + K_pk*ek;               % a priori estimate X
     neg_LLF = neg_LLF+1/2*log(det(Rek))+1/2*ek'/Rek*ek; 
     predX(:,k+1) = X; P = P + L*M*L';  predDP(:,k+1) = diag(P);  

     %--- recursion for new L and invM ----------------
    Rek_next = Rek + H*L*M*L'*H'; 
   invM_next = invM + L'*H'/Rek*H*L;              % inverse old Rek of size m  
      L_next = (F - K_pk*H )*L;
       K_pk  = (K_pk*Rek + F*L*M*L'*H')/Rek_next; % inverse new Rek of size m
    % --- new cicle -------------- 
     L = L_next; Rek = Rek_next;  invM = invM_next;  M = inv(invM);  % inverse matrix of size alpha, once    
end;
end
