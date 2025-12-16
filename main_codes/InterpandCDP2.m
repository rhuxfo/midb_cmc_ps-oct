function [CDP1, CDP2] = InterpandCDP2(bscans1,bscans2, Parameters)

MidLength = Parameters.alineLength/2;
PaddingLength1 = Parameters.InterpZeroPaddingFactor*Parameters.alineLength;
PaddedTransformedBscans1 = zeros(PaddingLength1-MidLength, size(bscans1,2));

%Channel 1
TransformedBscans1 = fft(bscans1);
TransformedBscans1 = TransformedBscans1(1:MidLength, :, :);
PDTB1 = vertcat(TransformedBscans1,PaddedTransformedBscans1);
ZeroPaddedBscans1 = real(ifft(PDTB1)) * Parameters.InterpZeroPaddingFactor;

%Channel 2
TransformedBscans2 = fft(bscans2);
TransformedBscans2 = TransformedBscans2(1:MidLength, :, :);
PDTB2 = vertcat(TransformedBscans2,PaddedTransformedBscans1);
ZeroPaddedBscans2 = real(ifft(PDTB2)) * Parameters.InterpZeroPaddingFactor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine bscan values for interpolated wavelengths
InterpolatedBscans1 = interp1(Parameters.Wavelengths_l, ZeroPaddedBscans1, Parameters.InterpolatedWavelengths,'linear','extrap'); % Interpolation
InterpolatedBscans2 = interp1(Parameters.Wavelengths_r, ZeroPaddedBscans2, Parameters.InterpolatedWavelengths,'linear','extrap'); % Interpolation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Dispersion compensation makes file sharper, add zeros to interpolated
    %buffer, and helps characterize coherence function
    if Parameters.dispersionComp == 1
        InterpolatedBscans1 = InterpolatedBscans1.*repmat(Parameters.PhaseCorrection1, [1, size(InterpolatedBscans1,2)]);
        InterpolatedBscans2 = InterpolatedBscans2.*repmat(Parameters.PhaseCorrection2, [1, size(InterpolatedBscans2,2)]);
    end
%Windowing
if Parameters.windowData == 1
    hammingWindow1 = repmat(hamming(Parameters.Start2-1), [1, 1]);
    hammingWindow2 = repmat(hamming(Parameters.Start2-1), [1, 1]);
    
    InterpolatedBscans1 = InterpolatedBscans1 .* hammingWindow1;
    InterpolatedBscans2 = InterpolatedBscans2 .* hammingWindow2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Calculate the complex depth profiles, adapted from Buffer2Jones.m
    AlineLength = (Parameters.OriginalLineLength*Parameters.CDPZeroPaddingFactor) / 2;
    
    CDP1 = fft(InterpolatedBscans1, Parameters.OriginalLineLength*Parameters.CDPZeroPaddingFactor);
    CDP2 = fft(InterpolatedBscans2, Parameters.OriginalLineLength*Parameters.CDPZeroPaddingFactor);
    
    CDP1(AlineLength+1:end, :, :) = []; %Cut fft in half due to symmetry of fft
    CDP1(1:(Parameters.AutoPeakCorrCut*Parameters.CDPZeroPaddingFactor), :, :) = []; %Remove DC offset
    
    CDP2(AlineLength+1:end, :, :) = []; %Cut fft in half due to symmetry of fft
    CDP2(1:(Parameters.AutoPeakCorrCut*Parameters.CDPZeroPaddingFactor), :, :) = []; %Remove DC offset
    
end %of function
