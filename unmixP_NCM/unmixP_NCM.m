function [P] = unmixP_NCM(X,E,cov,Parameters)
%% Input:
%        X - d-by-N HSI data, matrix
%        E - d-by-M endmember mean set, matrix
%        cov - d-by-d-by-M endmember covariance set, matrix
%  Output:
%        P - M-by-N proportion matrix
% Author: Alina Zare et al. Rewriten by Sheng Zou
% Department of Electrical and Computer Engineering, University of Florida
% 09/07/2017

%% Initialization of P
[d,N] = size(X);% d is the dimensionality and N is the number of pixels
M = size(E,2); % M is the nunmber of endmembers
P = 1/M.*ones(M,N); % intialize proportion values to be 1/M

%% Initialize all Likelihood and Prior Values
LogLikelihoodOld          = ComputeLogLikelihoodAll(X, E, P, cov, N); % compute the log likelihood

%%
for iteration = 2:Parameters.NumberIterations+1
    %% Metropolis-Hasting
    Y = randg(1, N, M) ; % generate some random values from gamma distribution
    v = sum(Y,2); %
    samples = (Y./v)'; % scale the random vector to sum-to-one, as the proposed proportion vector for each pixel
    LogLikelihoodNew      = ComputeLogLikelihoodAll(X, E, samples, cov, N); % compute the new log likelihood using the new proportion values
    Ratios                = exp(LogLikelihoodNew - LogLikelihoodOld); % always accept if the new log likelihood is larger than the old log likelihood, otherwise accept with a probability based on Ratios
    rands                 = rand(1,N);
    Vals                  = rands < Ratios;
    Vrep                  = repmat(Vals,M,1);
    P                     = samples.*Vrep  + P.*(1-Vrep); % Update the proportions
    
    LogLikelihoodOld = (1-Vals).*LogLikelihoodOld + (Vals).*LogLikelihoodNew; % Update the log likelihood values
    disp(strcat('iteration = ',num2str(iteration))) % display the index for current iteration
    disp(strcat('loglikelihood = ',num2str(sum(LogLikelihoodOld)))) % display the value of log likelihood
    
end


    function [LogLikelihoodAll] = ComputeLogLikelihoodAll(X, E, P, cov, N)
        %% compute log likelihood of all points
        
        term1 = zeros(1,N);
        term2 = zeros(1,N);
        
        for s = 1:size(E,2)
            statement(s)=isdiag(squeeze(cov(:,:,s))) && length(unique(diag(cov(:,:,s))))==1; % check if all the covariance matrices are diagonal and isotropic
        end
        
        if mean(statement) ==1 % all diagonal and isotropic
            for z = 1:size(P,1)
                a(z) = unique(diag(cov(:,:,z))); % the scalar on each covariance matrix
            end
            term1 = -.5*size(X,1)*log(a*(P.^2)); % compute the first term of the log likelihood equation (Ln of the determinant of the covariance)
            term2 = -.5*(sum((X-E*P).^2)./(a*(P.^2)));% compute the second term of the log likelihood equation (Ln of the exponential part)
        else % if the covariance matrices are full convariance, use a parallel computing method for faster computing 
            
            parfor t = 1:N
                P3D = zeros(1,1,size(E,2)); % Create a 1-by-1-by-M matrix
                P3D(:,:,1:size(E,2)) = P(:,t).^2; % Assign the proportion vector to P3D
                P3D_full = repmat(P3D,size(X,1),size(X,1)); % repmat the proportion to d-by-d-by-M
                term1(t) = -.5*log(det(sum(P3D_full.*cov,3)));% compute the first term of the log likelihood equation (Ln of the determinant of the covariance)
                term2(t) = -.5*(X(:,t)-E*P(:,t))'*(sum(P3D_full.*cov,3)\(X(:,t)-E*P(:,t)));% compute the second term of the log likelihood equation (Ln of the exponential part)
            end
            
        end
        
        LogLikelihoodAll = term1 + term2; % compute the whole term for log likelihood
    end



end