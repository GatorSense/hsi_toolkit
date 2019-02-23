function [SpecDist] = PlotSpectraDistribution(Spectra, WaveLengths,SampStuff, FigNum);
%function [SpecDist] = PlotSpectraDistribution(Spectra, WaveLengths,SampStuff, FigNum);
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A FUNCTION THAT CREATES & DISPLAYS SPECTRA AS A 2D HISTOGRAM   %%%
%%%    SPECTRA ARE ASSUMED REFLECTANCES OR EMISSIVITIES IN [0,1]   %%%
%%%    SPECTRA ARE MAPPED TO INTEGERS BETWEEN 0 AND 100 (OR < 100) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% INPUTS:
%%%   I1. Spectra IS A NUMBER SPECTRA x NUMBER BANDS...        %%%
%%%       ...ARRAY OF REFLECTANCES OR EMISSIVITIES             %%%
%%%   I2. WaveLengths IS A VECTOR OF THE SPECTRAL WAVELENGTHS  %%%
%%%   I3. SampStuff IS A VECTOR CONTAINING                     %%%
%%%          SampInt:       FRACTIONAL SIZE OF HISTOGRAM BINS  %%%
%%%          IntSampInt:    INT VERSION OF SampInt             %%%
%%%          IntTopReflect: INT VALUE OF MAX REF/EMIS BIN      %%%
%%%   I4. FigNum IS THE INDEX OF THE FIGURE TO USE FOR DISPLAY %%%
%%%          IF FigNum < 1, DO NOT DISPLAY ANYTHING            %%%
%%%
%%% OUTPUTS:
%%%   O1. SpecDist IS THE 2D HISTOGRAM                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                            %%%
%%% AUTHOR:      Darth Gader                                   %%%
%%% LAST UPDATE: 090218                                        %%%
%%%                                                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%% INITIALIZE PARAMETERS %%%
SampInt       = SampStuff(1);
IntSampInt    = SampStuff(2);
IntTopReflect = SampStuff(3);
SMOOTHSIZE    = [3,3];
NumSpec       = size(Spectra, 1);
NumWave       = size(Spectra, 2);
SpectDist     = zeros(IntTopReflect, NumWave);
assert(NumWave == length(WaveLengths), 'Wavelength sizes don''t match');

%%% MAP SPECTRA TO [0, 100] %%%
MappedSpectra = min(100, (Spectra*99)+1);
MappedSpectra = max(1, round(MappedSpectra/SampInt)*SampInt);

%%
%%% MAKE A HISTOGRAM FOR EACH WAVELENGTH %%%
for k = 1:NumWave
    SpecDist(:, k) = hist(MappedSpectra(:, k), 0:IntSampInt:IntTopReflect);
end

%%% SMOOTH BY TAKING A LOCAL MAX FOLLOWED BY A LOCAL AVERAGE %%%
SpecDist   = ordfilt2(SpecDist, 9, ones(3,3));
SpecDist   = conv2(SpecDist, (1/prod(SMOOTHSIZE))*ones(SMOOTHSIZE), 'same');

%%
%%% DISPLAY AS MESH %%%
if(FigNum > 0)
    XAxis      = WaveLengths;
    YAxis      = (0:IntSampInt:IntTopReflect)';
    figure(FigNum);
    mesh(XAxis, YAxis, SpecDist);
    title('Spectra Histogram')
    xlabel('Wavelength (nm)')
    ylabel('Reflectance')
end
%%% END OF FUNCTION %%%
%%%%%%%%%%%%%%%%%%%%%%%