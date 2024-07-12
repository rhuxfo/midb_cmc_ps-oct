function [Wavelengths_l, Wavelengths_r,InterpolatedWavelengths, Ks]=InterpolateWavelengths3(Parameters)
% [Wavelengths_l, Wavelengths_r, InterpolatedWavelengths, Ks] = InterpolateWavelengths(Parameters)
% Parameters = [PaddingLength, OriginalLineLength1, Start1, OriginalLineLength2, Start2]
% Calculates 1024px-InterpolatedWavelengths vector
% Edited 07/2020 and renamed from interpolationwave.m

if length(Parameters) ~= 5
    clear Parameters
    Parameters = [4096 1024 1 1024 1025];
end

PaddingLength = Parameters(1);
OriginalLineLength1 = Parameters(2);
Start1 = Parameters(3);
OriginalLineLength2 = Parameters(4);
Start2 = Parameters(5);

% 060514: Pixel numbers and wavelength values to fit to a curve %Basler
% position_l= [60 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000]; %pixel numbers
% position_r= [1130 1168.5 1217 1265 1314 1362 1410 1459 1507 1556 1604 1653 1701.5 1750 1799 1847.5 1896 1946 1994 2043]; %pixel numbers
% wave_sample= [804.48 807.22 810.74 814.15 817.55 820.93 824.38 827.73 831.15 834.55 837.93 841.21 844.65 847.98 851.38...
%               854.68 857.96 861.31 864.64 867.95];

%05/25/2023 Octoplus
% position_l_white= [3 43 87 131 173.8 216.8 261.2 304.2 347.4 390.6 432.8 476.2 518.5 560.9 602.8 645.6 688.4 729.5 771.1 813.6 855.3 897.3 937.6]; %pixel numbers
% position_r_red= [1041 1082 1127 1172.2 1216.2 1260.3 1305.8 1350 1394.4 1439.1 1482.5 1527.4 1570.7 1614.5 1658.2 1702.1 1746.4 1789.1 1832.1 1876.1 1919.4 1962.8 2005]; %pixel numbers 
% wave_sample= [866.74 863.99 860.99 857.96 855.01 852.04 848.98 845.98 842.97 839.95 837.00 833.95 830.97 827.98 824.98 821.97 818.94 815.98 813.01 809.95 806.96 803.95 801.03];
% 
% position_l = position_l_white;
% position_r = position_r_red;

%11/13/2023 Octoplus
position_l_white= [1.2 58.2 130.5 202.3 275.1 348 420.7 494 567.1 641.4 715.2 788.6 863 938.5 1004]; %pixel numbers
position_r_red= [42 97 166.6 236 306.2 376.9 447 518 588.8 661.4 733.2 804.7 877.2 951 1014.8]; %pixel numbers 
wave_sample= [801.03 805.01 810.04 815.02 820.06 825.07 830.04 835.05 840.04 845.07 850.06 855.01 860.01 865.05 869.40];

position_l = position_l_white;
position_r = position_r_red+1024;

%Second order least-squares fitting to find the wavelengths of the whole spectrum
pcoeff1=polyfit(position_l,wave_sample,2); % 
x1=1:1024;
lamda_l=polyval(pcoeff1,x1);
pcoeff2=polyfit(position_r,wave_sample,2); % 
x2=1025:2048;
lamda_r=polyval(pcoeff2,x2);

W_l=lamda_l(Start1:Start1+OriginalLineLength1-1);
W_r=lamda_r(Start2-1024:Start2+OriginalLineLength2-1025);
xx1=linspace(Start1,Start1+OriginalLineLength1-1,PaddingLength);
Wavelengths_l=1e-9*interp1([Start1:Start1+OriginalLineLength1-1],W_l,xx1)';
xx2=linspace(Start2-1024,Start2+OriginalLineLength2-1025,PaddingLength);
Wavelengths_r=1e-9*interp1([Start2-1024:Start2+OriginalLineLength2-1025],W_r,xx2)';

clear x1 lamda_l x2 lamda_r W_l W_r xx1 xx2

minK = 2*pi / min([Wavelengths_l(end) Wavelengths_r(end)]);
maxK = 2*pi / max([Wavelengths_l(1) Wavelengths_r(1)]);

Ks=linspace(minK,maxK,OriginalLineLength1)';
InterpolatedWavelengths = (2*pi) ./ Ks; %Wavelengths we want; represents linear K-values
end
