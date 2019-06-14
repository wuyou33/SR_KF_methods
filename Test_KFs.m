% ------------------------------------------------------------------- 
% Script for comparing various KF implementation methods.
% Authors: Maria Kulikova:     maria dot kulikova at ist dot utl dot pt 
% License: GNU General Public License Version 2 
% ------------------------------------------------------------------- 
clear all; close all; clc; warning off;

% ---- Parameters that you may change ------
rng('default')
MC_runs  = 10;                 % Number of Monte Carlo runs
N_total  = 100;                % Discrete-time instances
p = pwd; cd('Noises/'); noise_type = @noise_gauss; cd(p);  % Type of uncertainties 

% ----Load  Model to be examined ----
%p = pwd; cd('Models/'); [Fsys,Gsys,Qsys,Hsys,Rsys,P0,x0] = Model_satellite; cd(p); 
%p = pwd; cd('Models/'); [Fsys,Gsys,Qsys,Hsys,Rsys,P0,x0] = Model_resonator; cd(p); 
p = pwd; cd('Models/'); [Fsys,Gsys,Qsys,Hsys,Rsys,P0,x0] = Model_electrocardiogram; cd(p); 

%--- Filtering methods to be examined ----
p = pwd; cd('Methods-KF-posteriori/'); 
   handle_funs{1} = @Riccati_KF;         % Conventional
   handle_funs{2} = @Riccati_KF_seq;     % Sequential
   handle_funs{3} = @Riccati_KF_SR;      % Cholesky-based (upper triangular factors)
   handle_funs{4} = @Riccati_KF_SRe;     % extended Cholesky-based (upper triangular factors)
   handle_funs{5} = @Riccati_KF_SRlow;   % Cholesky-based (lower triangular factors)
   handle_funs{6} = @Riccati_KF_SRelow;  % extended Cholesky-based (lower triangular factors)
   handle_funs{7} = @Riccati_KF_SVD;     % SVD-based implementation (2017)
   handle_funs{8} = @Riccati_KF_SVDr;    % SVD-based method with improved robustness (2018)
   handle_funs{9} = @Riccati_KF_SVDe;    % "economy size" SVD-based method
   handle_funs{10} = @Riccati_KF_SVDSR;   % SVD-based with Cholesky decomposition for noise covariances (1992)
cd(p); 
start_alternative = size(handle_funs,2)+1;
p = pwd; cd('Methods-KF-priori/');      
   handle_funs{11} = @Riccati_KF;        % Conventional
   handle_funs{12} = @Riccati_KF_SR;     % Cholesky-based (upper triangular factors)
   handle_funs{13} = @Riccati_KF_SRe;    % extended Cholesky-based (upper triangular factors)
   handle_funs{14} = @Riccati_KF_SRlow;  % Cholesky-based (lower triangular factors)
   handle_funs{15} = @Riccati_KF_SRelow; % extended Cholesky-based (lower triangular factors)
   handle_funs{16} = @Chandrasekhar_KF1; % Conventional Chandrasekhar-based var 1
   handle_funs{17} = @Chandrasekhar_KF2; % Conventional Chandrasekhar-based var 2
   handle_funs{18} = @Chandrasekhar_KF3; % Conventional Chandrasekhar-based var 3
   handle_funs{19} = @Chandrasekhar_KF4; % Conventional Chandrasekhar-based var 4
cd(p); 

% --- Monte Carlo runs ----
 Number_Methods = size(handle_funs,2); % number of methods to be tested
 filters = cell(Number_Methods,1);     % prelocate for efficiency
 for exp_number = 1:MC_runs  
     fprintf(1,'Simulation #%d: \n',exp_number); 
     % ----Simulate system to get "true" state and measurements  ----
     [Times,Exact_StateVectorX,Measurements,~,~] = Simulate_Measurements(noise_type,{Fsys,Gsys,Qsys,Hsys,Rsys},{x0,P0},N_total);
     % ----Solve the inverse problem, i.e. perform filtering
     for i=1:Number_Methods;
       filters{i}.legend = func2str(handle_funs{i}); 
       tstart = tic;
       [LLF,hatX,hatDP] = feval(handle_funs{i},{Fsys,Gsys,Qsys,Hsys,Rsys},{x0,P0},Measurements);
       ElapsedTime = toc(tstart);
       filters{i}.hatX(:,:,exp_number)  = hatX;
       filters{i}.hatDP(:,:,exp_number) = hatDP;
       filters{i}.trueX(:,:,exp_number) = Exact_StateVectorX;
       filters{i}.yk(:,:,exp_number)    = Measurements;
       filters{i}.AE(:,:,exp_number)    = abs(hatX-Exact_StateVectorX); 
       filters{i}.Time(exp_number)      = ElapsedTime;
       filters{i}.neg_LLF(exp_number)   = LLF;
    end;
end;
% --- Estimation errors computation  ----
for i=1:Number_Methods;
   val1 = mean((filters{i}.AE).^2,3);             % mean by Monte Carlo runs 
   filters{i}.RMSE = sqrt(mean(val1,2));          % mean by time, i.e. RMSE in each component
end;

% --- Print the results  ----
fprintf(1,'\n \n ________ A POSTERIORI FORM (filtered X and P are computed with the related errors)_________ \n')
fprintf(1,'--------------------- \n'); fprintf(1,'  Filter Implementations:\t');       
for i=1:size(P0,1), fprintf(1,'RMSE_x(%d)\t',i); end; fprintf(1,'||RMSE||_2 \t av.CPU (s) \t neg LLF \n');
for i=1:start_alternative-1
 fprintf(1,'%d.%22s\t ',i,filters{i}.legend); 
 fprintf(1,'%8.4f\t',filters{i}.RMSE); fprintf(1,'%8.4f\t%8.4f\t%8.4f \n',norm(filters{i}.RMSE,2),mean(filters{i}.Time),filters{i}.neg_LLF(1));
end;
fprintf(1,'\n _________ A PRIORI FORM (predicted X and P are computed with the related errors)____________ \n')
fprintf(1,'--------------------- \n'); fprintf(1,'  Filter Implementations:\t');       
for i=1:size(P0,1), fprintf(1,'RMSE_x(%d)\t',i); end; fprintf(1,'||RMSE||_2 \t av.CPU (s) \t neg LLF \n');
for i=start_alternative:Number_Methods
 fprintf(1,'%d.%22s\t ',i,filters{i}.legend); 
 fprintf(1,'%8.4f\t',filters{i}.RMSE); fprintf(1,'%8.4f\t%8.4f\t%8.4f \n',norm(filters{i}.RMSE,2),mean(filters{i}.Time),filters{i}.neg_LLF(1));
end;

% ___ Illustrate, if you wish (ensure that all related entries coincide)
% Illustrate_XP(filters(1:start_alternative-1));               % Illustrate a posteriori values
% Illustrate_XP(filters(start_alternative:Number_Methods));    % Illustrate a priori values
