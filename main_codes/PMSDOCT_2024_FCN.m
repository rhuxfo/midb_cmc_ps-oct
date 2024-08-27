function [status] = PMSDOCT_2024_FCN(P)
%% Parameter Initialization
%Setting up files
data_directory = P.dir; %data called from here
save_directory = P.Sdir; %Tile data saved here

data_filename = strcat(data_directory,P.baseN);
save_n = 'slice_';

%save folders
CompFolder = 'TComp'; % Full 3D data (Tiles)
EnFolder = 'Enface'; %Enface folder (Tiles)
StFolder = 'Stitched'; % Composite Slices folder (stitched tiles)

%contrast folders
c1 = 'CDP';
c2 = 'A1A2';
c3 = 'Orientation';
c4 = 'Retardance';
c5 = 'AbsoOri';
c6 = 'Cross';
c7 = 'Reflectivity';
folderNames = {c1,c2,c3,c4,c5,c6,c7};

%automatically create save folders if they do not exist
if P.autofolder ==1 
    folder3d = fullfile(P.Sdir,CompFolder);
    folderEnface = fullfile(P.Sdir,EnFolder);
    folderStitch = fullfile(P.Sdir,StFolder);
    
    if ~exist(folder3d,'dir')
        mkdir(folder3d);
        for i=1:length(folderNames)
            mkdir(fullfile(folder3d,folderNames{i}));
        end
    end
    if ~exist(folderEnface,'dir')
        mkdir(folderEnface);
        for i=1:length(folderNames)
            mkdir(fullfile(folderEnface,folderNames{i}));
        end
    end
    if ~exist(folderStitch,'dir')
        mkdir(folderStitch);
        for i=1:length(folderNames)
            mkdir(fullfile(folderStitch,folderNames{i}));
        end
    end
end

% Data size parameters
XTiles = P.XTiles;
YTiles = P.YTiles;
TileMtrx = zeros(XTiles,YTiles);
m=0;
for p=1:YTiles
    for q=1:XTiles
        TileMtrx(q,p)= 1+m;
        m = m+1;
    end
end

% Data size parameters
slice = P.Slices; %Number of slices / slice being analyzed
tilenum = P.tiles; % tile number in slice
scan = P.buffers; %Number scan in the slice

st = P.depthstart;
endc = P.depthcut;
ov = P.overlap;
Nthr = 6.5;
Fllip = P.Flip;
% Load variables
filePointer = fopen([data_filename, num2str(slice(1)),P.tileN,num2str(tilenum(1)),'_840_1.dat'], 'r', 'l');
headerStr = fgetl(filePointer); %Getting things from the labview
evalc(headerStr); %Evaluate header from labview to import variables
% Variables loaded: alineLength; alinePeriod; blineLength; buffersPerFile;
%                   clineLength; interAlineInterval; interscanInterval;
%                   stimDelay; stimDuration; stimEveryBScan; stimVoltage;
%                   XgalvoVoltageMax; XgalvoVoltageMin; YgalvoVoltageMax;
%                   YgalvoVoltageMin
fclose(filePointer);
Parameters.num_bscans = buffersPerFile;
Parameters.blineLength = blineLength;
Parameters.alineLength = alineLength;
Parameters.alines = (scan(end)*buffersPerFile);
linenum =1:blineLength;
Tile = zeros(endc,blineLength,Parameters.alines);

% Data processing actions
% 1 = do the action; 0 = skip the action
Parameters.dispersionComp = P.disper; %Dispersion compensation, optimizes shape of coherence peak
Parameters.windowData = P.wind;
Parameters.background = P.BGremoval;

% 1 = Calculate Contrast
calcCDP = 0;
calcCh1Ch2 = 1;
calcReflectivity = P.Flect;
calcRetardance = P.Retar;
calcCrossPolar = P.Cr;
calcOrientation = P.Orio;
calcAbsOrientation = P.AbOrio;
Enface = P.En;
dB = 1;

% 1 = Save
STileComp = P.TCsv;
SEnface = P.Ensv;
SStitch = P.Stsv;

%% Dispersion
if Parameters.dispersionComp==1 %need to make 1024
    %Software dispersion compensation creates phase correction vectors
    %using FT shifting properties
    dispcompfile1 = fullfile(P.DCf1);
    fid1 = fopen(dispcompfile1,'r');
    angledisp1 = fread(fid1,1024,'real*8');
    fclose(fid1);
    Parameters.PhaseCorrection1 = exp(-1i.*angledisp1);
    
    dispcompfile2 = fullfile(P.DCf2);
    fid2 = fopen(dispcompfile2,'r');
    angledisp2 = fread(fid2,1024,'real*8');
    fclose(fid2);
    Parameters.PhaseCorrection2 = exp(-1i.*angledisp2);
end
%%
% Calculate constant interpolation wavelengths in k-space
Parameters.InterpZeroPaddingFactor = 4; %Zero-padding factor for interpolation and for zero-padding the bscans
Parameters.CDPZeroPaddingFactor = 1; %Zero-padding factor in calculation of CDP to determine level of visualization
Parameters.OriginalLineLength = 1024;
InterpZeroPaddingLength = Parameters.InterpZeroPaddingFactor*Parameters.OriginalLineLength;
Start1 = 1;
Parameters.Start2 = Start1 + Parameters.OriginalLineLength; %1025

InterpolationParameters = [InterpZeroPaddingLength, Parameters.OriginalLineLength, Start1, Parameters.OriginalLineLength, Parameters.Start2];
[Parameters.Wavelengths_l, Parameters.Wavelengths_r, Parameters.InterpolatedWavelengths, ~] = InterpolateWavelengths3(InterpolationParameters); %left and right wavelengths refer to different polarization channels
Parameters.AutoPeakCorrCut = 10; %Cut low-frequency points to not see big dc offset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Data
for SliceInd=1:length(slice)
    if ~P.StitchOnly ==1       
        for TileInd = 1:length(tilenum)
            filename=strcat(data_filename,num2str(slice(SliceInd)),P.tileN,num2str(tilenum(TileInd)),'_840_');
            %Read in file, zero-pad, and interpolate
            fprintf('Processing tile %d  ...\n', tilenum(TileInd)); %Print out file being processed
            
            [Raw_1,~,Blines] = Read2024(filename,scan,Parameters);
            for Line = 1:length(linenum)
                               fprintf('Processing line %d ...\n',linenum(Line));
                b1 = Raw_1(:,:,linenum(Line));
                bscan1 = b1(1:Parameters.alineLength,:);
                bscan2 = b1(Parameters.alineLength+1:end,:);

                if Parameters.windowData == 1
                    hammingWindow1 = repmat(hamming(Parameters.Start2-1), [1, Parameters.blineLength]);
                    hammingWindow2 = repmat(hamming(Parameters.Start2-1), [1, Parameters.blineLength]);

                    bscan1 = bscan1 .* hammingWindow1;
                    bscan2 = bscan2 .* hammingWindow2;
                end

                [CDP2, CDP1] = InterpandCDP(bscan1, bscan2, Parameters);

                % if calcAbsOrientation == 1
                blin1 = Blines(1:Parameters.alineLength,linenum(Line));
                blin2 = Blines(Parameters.alineLength+1:end,linenum(Line));
                if Parameters.windowData == 1
                    blin1 = blin1 .* hammingWindow1;
                    blin2 = blin2 .* hammingWindow2;
                end
                [cdp2, cdp1] = InterpandCDP(blin1, blin2, Parameters);
                %end
                if Line ==1
                    ch1bg = abs(cdp1).^2;
                    ch1bg =  10*log10(ch1bg);
                    ch2bg = abs(cdp2).^2;
                    ch2bg =  10*log10(ch2bg);
                    CH1BG = mean(ch1bg,2);
                    CH2BG = mean(ch2bg,2);
                    lim1 = floor(CH1BG(150));
                    lim2 = floor(CH2BG(150));
                else
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                CDiv = CDP1./CDP2;
                MMM = CDP1.*conj(CDP2);

                CH1 = abs(CDP1).^2;
                CH2 = abs(CDP2).^2;
                reflectivity = CH1 + CH2;
                if dB == 1
                    CH1 = 10*log10(CH1);
                    CH2 = 10*log10(CH2);
                    reflectivity = 10*log10(reflectivity);
                end
                
                Retardance = (180/pi)*atan(abs(CDP2)./abs(CDP1));
                
                if calcOrientation == 1
                    ePi = exp(1i*pi);
                    Theta = (ePi./CDiv);
                    Theta2 = Theta./(abs(Theta));
                end
                if calcAbsOrientation == 1
                    a = 200;
                    tempAmp2 = abs(cdp2).^2;
                    [~,P2] = max(tempAmp2(a:a+100));
                    mmm = cdp1.*conj(cdp2);
                    px2 = P2+a-1;
                    RLO2 = mmm(px2);
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if calcCDP == 1
                    if Line ==1
                        Tile_CDP1 = Tile(st:endc,:,:);
                        Tile_CDP2 = Tile(st:endc,:,:);
                    end
                    Tile_CDP1(:,:,Line)= CDP1(st:endc,:);
                    Tile_CDP2(:,:,Line)= CDP2(st:endc,:);
                end
                if calcCh1Ch2 == 1
                    if Line ==1
                        Tile_ch1 = Tile(st:endc,:,:);
                        Tile_ch2 = Tile(st:endc,:,:);
                    end
                    Tile_ch1(:,:,Line)= CH1(st:endc,:);
                    Tile_ch2(:,:,Line)= CH2(st:endc,:);
                end
                if calcReflectivity ==1
                    if Line ==1
                        Tile_R = Tile(st:endc,:,:);
                    end
                    Tile_R(:,:,Line)= reflectivity(st:endc,:);
                end
                if calcRetardance == 1
                    if Line ==1
                        Tile_R2 = Tile(st:endc,:,:);
                    end
                    Tile_R2(:,:,Line)= Retardance(st:endc,:);
                end
                if calcCrossPolar == 1
                    if Line ==1
                        Tile_cross = Tile(st:endc,:,:);
                    end
                    Tile_cross(:,:,Line)= CH2(st:endc,:);
                end
                if calcOrientation == 1
                    if Line ==1
                        Tile_O = Tile(st:endc,:,:);
                    end
                    Tile_O(:,:,Line)= Theta2(st:endc,:);
                end
                if calcAbsOrientation == 1
                    if Line ==1
                        RLO2T = zeros(1,scan(end)*Parameters.num_bscans);
                        Tile_Om = Tile(st:endc,:,:);
                    end
                    RLO2T(1,Line) = RLO2;
                    Tile_Om(:,:,Line)= MMM(st:endc,:);
                end
            end
            clear Raw_1
            %% Calulate Contrast Enfaces
            if Enface ==1
                ch1Limit = lim1;
                ch2Limit = lim2;
                cut = 200;
                if calcReflectivity ==1
                    disp('Calculating Reflectivity Enface')
                    EnRef = CombomaskCross(Tile_R,ch1Limit+1,cut);
                end
                if calcCrossPolar == 1
                    disp('Calculating Cross Enface')
                    EnCr = CombomaskCross(Tile_cross,ch2Limit,cut);
                end
                if calcRetardance == 1
                    disp('Calculating Retardance Enface')
                    EnR = CombomaskR4(Tile_ch1,Tile_ch2,Tile_R2,ch1Limit,ch2Limit,cut);
                end
                if calcOrientation == 1
                    disp('Calculating Ori Enface')
                    EnO = Combomask4(Tile_ch1,Tile_ch2,Tile_O,ch1Limit,ch2Limit,cut);
                end
                if calcAbsOrientation == 1
                    disp('Calculating Abs Ori Enface')
                    EnO2 = Combomask4(Tile_ch1,Tile_ch2,Tile_Om,ch1Limit+Nthr,ch2Limit+Nthr,cut);
                    Off2 = RLO2T;
                    EnAO = EnO2./ Off2;
                end

            end
            %%
            if STileComp == 1
                disp('Saving 3D Tile')
                Save_base = strcat(save_directory,CompFolder);
                if calcCDP == 1
                    SaveF = fullfile(Save_base,c1);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_CDP1');
                    save(output_filename,'Tile_CDP1');
                    output_filename2 = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_CDP2');
                    save(output_filename2,'Tile_CDP2');
                end
                if calcCh1Ch2 == 1
                    SaveF = fullfile(Save_base,c2);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_CH1');
                    save(output_filename,'Tile_ch1');
                    
                end
                if calcReflectivity == 1
                    SaveF = fullfile(Save_base,c7);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_Reflect');
                    save(output_filename,'Tile_R');
                    
                end
                if calcCrossPolar == 1
                    SaveF = fullfile(Save_base,c6);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_Cross');
                    save(output_filename,'Tile_cross');
                    
                end
                if calcRetardance == 1
                    SaveF = fullfile(Save_base,c4);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_R');
                    save(output_filename,'Tile_R2');
                    
                end
                if calcOrientation == 1
                    SaveF = fullfile(Save_base,c3);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_Orien');
                    save(output_filename,'Tile_O');
                    
                end
                if calcAbsOrientation == 1
                    SaveF = fullfile(Save_base,c5);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_AbsOrien');
                    save(output_filename,'Tile_Om');
                    
                end
            end
            %%
            if SEnface ==1
                disp('Saving Enface')
                Save_base = strcat(save_directory,EnFolder);
                if calcReflectivity == 1
                     ON = strcat(save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnRef');
                    SaveFN = fullfile(Save_base,c7,ON);
                    save(SaveFN,'EnRef');
                    
                end
                if calcCrossPolar == 1
                     ON = strcat(save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnCr');
                    SaveFN = fullfile(Save_base,c6,ON);
                    save(SaveFN,'EnCr');
                end
                if calcRetardance == 1
                     ON = strcat(save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnR');
                    SaveFN = fullfile(Save_base,c4,ON);
                    save(SaveFN,'EnR');
                end
                if calcOrientation == 1
                     ON = strcat(save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnO');
                    SaveFN = fullfile(Save_base,c3,ON);
                    save(SaveFN,'EnO');
                end
                if calcAbsOrientation == 1
                    ON = strcat(save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnAO');
                    SaveFN = fullfile(Save_base,c5,ON);
                    save(SaveFN,'EnAO');
                    
                end
            end
            
        end %tile for loop
        status = 1;
    end
    if SStitch ==1
        disp('Stitching Tiles')
        Call_base = strcat(save_directory,EnFolder);
        Save_base = strcat(save_directory,StFolder);
        if calcCrossPolar ==1
            CallF = fullfile(Call_base,c6);
            SaveF = fullfile(Save_base,c6);
            [TEnCr]= MStitchFCN_mod(slice(SliceInd),6,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip);
        end
        if calcReflectivity == 1
            CallF = fullfile(Call_base,c7);
            SaveF = fullfile(Save_base,c7);
            [TEnRef]= MStitchFCN_mod(slice(SliceInd),7,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip);
        end
        if calcRetardance == 1
            CallF = fullfile(Call_base,c4);
            SaveF = fullfile(Save_base,c4);
            [TEnR]= MStitchFCN_mod(slice(SliceInd),4,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip);
        end
        if calcOrientation == 1
            CallF = fullfile(Call_base,c3);
            SaveF = fullfile(Save_base,c3);
            [TEnO]= MStitchFCN_mod(slice(SliceInd),3,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip);
        end
        if calcAbsOrientation == 1
            CallF = fullfile(Call_base,c5);
            SaveF = fullfile(Save_base,c5);
            [TEnAO]= MStitchFCN_mod(slice(SliceInd),5,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip);
        end
        status = 2;
    end
end %slice for loop
fprintf('Processing completed ...\n');
end
