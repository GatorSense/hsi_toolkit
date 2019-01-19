function [E, IdxOfE, Xpca] = VCAbyDG(X,varargin)
%        [E, IdxOfE, Xpca] = VCAbyDG(X,'Endmembers',M,'SNR',r,'verbose',v)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A FUNCTION TO CALCULATE ENDMEMBERS USING Vertex Component Analysis (VCA)
% REFERENCE:
% José M. P. Nascimento and José M. B. Dias
% "Vertex Component Analysis: A Fast Algorithm to Unmix Hyperspectral Data"
% IEEE Trans. Geosci. Remote Sensing,
% April, 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% INPUTS:
%%%         X:          DATA MATRIX B WAVELENGTHS x N PIXELS
%%%         M:          NUMBER OF ENDMEMBERS TO ESTIMATE
%%% OUTPUTS:
%%%         E:          B x M MATRIX OF ESTIMATED ENDMEMBERS
%%%         IdxOfEinX:  INDICES OF ENDMEMBERS IN DATA X.
%%%                         (THE ENDMEMBERS ARE SELECTED FROM X)
%%%         XM:         X PROJECTED ONTO FIRST M COMPONENTS OF PCA
%
%%% OPTIONAL INPUTS:
%%% 'SNR'     r:   ESTIMATED SIGNAL TO NOISE RATIO IN dB
%%% 'verbose' v:   logical TOGGLE TO TURN DISPLAYS ON AND OFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Authors: José Nascimento (zen@isel.pt)
%%%          José Bioucas Dias (bioucas@lx.it.pt)
%%% Copyright (c)
%%% version: 2.1 (7-May-2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% MODIFED BY Darth Gader OCTOBER 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET DEFAULT AND PARSE INPUT PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose  = true; % default
InputSNR = 0;  % default this flag is zero,
% which means we estimate the SNR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking for input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim_in_par = length(varargin);
if (nargin - dim_in_par)~=1
    error('Wrong parameters');
elseif rem(dim_in_par,2) == 1
    error('Optional parameters should always go by pairs');
else
    for m = 1 : 2 : (dim_in_par-1)
        switch lower(varargin{m})
            
            case 'verbose'
                verbose   = varargin{m+1};
                
            case 'endmembers'
                M         = varargin{m+1};
                
            case 'snr'
                SNR      = varargin{m+1};
                InputSNR = 1;       % flag meaning that user gives SNR
            otherwise
                fprintf(1,'Unrecognized parameter:%s\n', varargin{m});
        end %switch
    end %for
end %if
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(X)
    error('there is no data');
else
    [B N]=size(X);
end

if (M<0 | M>B | rem(M,1)~=0),
    error('ENDMEMBER parameter must be integer between 1 and B');
end
%%
%%%%%%%%%%%%%%%%%%%%
%%% ESTIMATE SNR %%%
%%%%%%%%%%%%%%%%%%%%
if InputSNR==0,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% THE USER DID NOT INPUT AN SNR SO SNR IS CALCULATED HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALCULATE PCA AND PCT %%%
    MuX          = mean(X, 2);
    Xz           = X - repmat(MuX, [1, N]);
    SigmaX       = cov(X');
    [U,S,V]      = svds(SigmaX,M);
    Xpca         = U' * Xz;
    ProjComputed = true;
    
    %%% ESTIMATE SIGNAL TO NOISE %%%
    SNR = EstSNR(X,MuX,Xpca);
    
    %%% PRINT SNR %%%
    if verbose
        fprintf(1,'Estimated SNR  = %g[dB]\n',SNR);
    end
else
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% THE USER DID INPUT AN SNR SO NO SNR CALCUATION NEEDED
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if verbose
        fprintf(1,'Input SNR = %g[dB]\t',SNR); 
    end
end

%%% SET THRESHOLD TO DETERMINE IF NOISE LEVEL IS LOW OR HIGH %%%;
SNRThresh = 15 + 10*log10(M);

if verbose
    fprintf('SNRThresh= %f SNR= %f  Difference= %f\n', SNR, SNRThresh, SNRThresh-SNR)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF ESTIMATE SNR %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROJECTION.                                                %%%
%%% SELECT AND CALCULATE PROJECTION ONTO M DIMS                %%%
%%%                                                            %%%
%%% IF SNR IS LOW,IT IS ASSUMED THAT THERE IS STILL NOISE IN   %%%
%%% THE SIGNAL SO REDUCE DIM A BIT MORE.                       %%%
%%% ADD A CONSTANT VECTOR TO KEEP IN DIM M.                    %%%
%%%                                                            %%%
%%% IF SNR IS HIGH, PROJECT TO SIMPLEX IN DIM M-1 TO REDUCE    %%%
%%% EFFECTS OF VARIABLE ILLUMINATION                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SNR < SNRThresh
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN LOW SNR CASE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% PRINT MESSAGE %%%
    if verbose
        fprintf('Low SNR so Project onto Dimension M-1.\n');
    end
    
    %%% REDUCE SIZE OF PCT MATRIX TO M-1 %%%
    Dim    = M-1;
    MuX    = mean(X, 2);
    BigMuX = repmat(MuX, [1, N]);
    if ProjComputed
        U= U(:,1:Dim);
    else
        Xz      = X - BigMuX;
        SigmaX  = cov(X');
        [U,S,V] = svds(SigmaX,Dim);
        Xpca    =  U' * Xz;
    end
    
    %%% REDUCE DIMENSIONALITY IN PCA DOMAIN %%%
    XpcaReduced = Xpca(1:Dim,:);
    
    %%% RECONSTRUCT X "WITHOUT NOISE" BY PROJECTING BACK TO ORIGINAL SPACE %%%
    XNoiseFree =  U * XpcaReduced + BigMuX;
    
    %%% CONCATENATE  CONSTANT VEC = MAX NORM OF ALL DATA POINTS %%%
    BiggestNorm = sqrt(max(sum(XpcaReduced.^2,1)));
    YpcaReduced = [XpcaReduced ; BiggestNorm*ones(1,N)] ;
%%%%%%%%%%%%%%%%%%%%%%%%
%%% END LOW SNR CASE %%%
%%%%%%%%%%%%%%%%%%%%%%%%
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN HIGH SNR CASE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose
        fprintf(1,'High SNR so project onto dimension M\n',SNR); 
    end
    
    %%% CONTRUCT "PCA-LIKE" MATRIX BY DIAGONALIZING CORRELATION MATRIX     %%%
    %%% IF SUBTRACT THE MEAN, THEN MUPCA WILL BE 0, SO NO MEAN SUBTRACTION %%%
    [U,S,V] = svds(X*X'/N,M);
    
    %%% CALC PCT WITHOUT SUBTRACTING MEAN AND THEN RECONSTRUCT %%%
    XpcaReduced = U'*X;
    
    %%% RECONSTRUCT X VIA INVERSE XFORM ON REDUCED DIM PCA DATA       %%%
    %%% XXX PDG NOT SURE IF THIS IS CORRECT IN ORIGINAL CODE OR HERE. %%%
    %%% I THINK THE MEAN SHOULD BE ADDED LIKE IN LOW SNR CASE         %%%
    XNoiseFree =  U * XpcaReduced(1:M,:);      % again in dimension L (note that x_p has no null mean)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALCULATE NORMALIZED PROJECTION                           %%%
    %%% SEE LAST PARAGRAPH OF SECTION A AND FIGURE 4 IN VCA PAPER.%%%
    %%% MEAN OF THE PROJECT DATA IS USED AS u FROM THAT SECTION   %%%
    Mupca       = mean(XpcaReduced,2);
    Denom       = sum( XpcaReduced .* repmat(Mupca,[1 N]));
    YpcaReduced = XpcaReduced./ repmat( Denom ,[M 1]);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VCA ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%   INITIALIZE ARRAY OF INDICES OF ENDMEMBERS, IdxOfE
%%%   INITIALIZE MATRIX OF ENDMEMBERS IN PCA SPACE, Epca
%%%   DO M TIMES
%%%      PICK A RANDOM VECTOR w
%%%      USE w TO FIND f IN NULL SPACE OF Epca
%%%      PROJECT f ONTO ALL PCA DATA PTS
%%%      FIND INDEX OF PCA DATA PT WITH MAX PROJECTION
%%%      USE INDEX TO ADD PCA DATA PT TO Epca
%%%   END DO
%%%   
%%%   USE INDICES TO SELECT ENDMEMBERS IN SPECTRAL SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IdxOfE    = zeros(1,M);
Epca      = zeros(M,M);
Epca(M,1) = 1;
for m=1:M
    w                 = rand(M,1);
    f                 = w - Epca*pinv(Epca)*w;
    f                 = f / sqrt(sum(f.^2));
    v                 = f'*YpcaReduced;
    [v_max IdxOfE(m)] = max(abs(v));
    Epca(:,m)         = YpcaReduced(:,IdxOfE(m));
end

E = XNoiseFree(:,IdxOfE);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF VCA FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION TO ESTIMATE SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SNR = EstSNR(X,MuX,Xpca)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE ESTIMATED SNR IS EQUIVALENT TO THE RATIO OF 
%%% POWER OF SIGNAL TO POWER OF NOISE
%%% BUT IS COMPUTED AS
%%% POWER OF (SIGNAL + NOISE) TO POWER OF NOISE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[B N]=size(X); 
[M N]=size(Xpca);

%%% POWER OF SIGNAL + NOISE, SIGNAL, AND NOISE, RESP. %%%
Psn = sum(X(:).^2)/N;
Ps  = sum(Xpca(:).^2)/N + MuX'*MuX;
Pn  = Psn - Ps;
SNR = 10*log10(Ps/Pn);

%%% OLD SNR CODE AS IN VCA PAPER.  NOT MUCH DIFFERENT %%%
%%% BECAUSE (M/B) IS SMALL %%%
%%% SNR = 10*log10( (Ps - M/B*Psn)/(Psn - Ps) );
%%% fprintf('SNR= %8.4f   RatSNR= %8.4f\n', SNR, RatSNR);

return;
%%%%%%%%%%%%%%%
%%% THE END %%%
%%%%%%%%%%%%%%%