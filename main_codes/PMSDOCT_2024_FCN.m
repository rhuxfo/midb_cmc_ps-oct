function [status] = PMSDOCT_2024_FCN(P)
%% Parameter Initialization
%Setting up files
data_directory = P.dir; %data called from here
save_directory = P.Sdir; %Tile data saved here

data_filename = strcat(data_directory,P.baseN);
save_n = 'slice_';

%save folders
CompFolder = 'TComp/'; % Full 3D data (Tiles)
EnFolder = 'Enface/'; %Enface folder (Tiles)
StFolder = 'Stitched/'; % Composite Slices folder (stitched tiles)

%contrast folders
c1 = 'CDP/';
c2 = 'A1A2/';
c3 = 'Orientation/';
c4 = 'Retardance/';
c5 = 'AbsoOri/';
c6 = 'Cross/';
c7 = 'Reflectivity/';
folderNames = {c1,c2,c3,c4,c5,c6,c7};

%automatically create save folders if they do not exist
if P.autofolder ==1 
    folder3d = strcat(P.Sdir,CompFolder);
    folderEnface = strcat(P.Sdir,EnFolder);
    folderStitch = strcat(P.Sdir,StFolder);
    
    if ~exist(folder3d,'dir')
        mkdir(folder3d);
        for i=1:length(folderNames)
            mkdir(strcat(folder3d,folderNames{i}));
        end
    end
    if ~exist(folderEnface,'dir')
        mkdir(folderEnface);
        for i=1:length(folderNames)
            mkdir(strcat(folderEnface,folderNames{i}));
        end
    end
    if ~exist(folderStitch,'dir')
        mkdir(folderStitch);
        for i=1:length(folderNames)
            mkdir(strcat(folderStitch,folderNames{i}));
        end
    end
end

status.Saved2 = 0;
status.Saved3 = 0;
status.Saved4 = 0;
status.Saved5 = 0;
status.Saved6 = 0;
status.Saved7 = 0;

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
lim1 = P.Ch_dB_limit;
lim2 = lim1;
px = P.Rline;
ov = P.overlap;

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
linenum =1:blineLength;
Tile = zeros(endc,blineLength,(scan(end)*buffersPerFile));

% Data processing actions
% 1 = do the action; 0 = skip the action
Parameters.dispersionComp = P.disper; %Dispersion compensation, optimizes shape of coherence peak
Parameters.windowData = P.wind;
Parameters.background = P.BGremoval;

% 1 = Calculate Contrast
calcCDP = P.CDP;
calcCh1Ch2 = P.CH12;
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
    dispcompfile1 = P.DCf1;
    fid1 = fopen(dispcompfile1,'r');
    angledisp1 = fread(fid1,1024,'real*8');
    fclose(fid1);
    Parameters.PhaseCorrection1 = exp(-1i.*angledisp1);
    
    dispcompfile2 = P.DCf2;
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
[Parameters.Wavelengths_l, Parameters.Wavelengths_r, Parameters.InterpolatedWavelengths, Ks] = InterpolateWavelengths3(InterpolationParameters); %left and right wavelengths refer to different polarization channels
Parameters.AutoPeakCorrCut = 10; %Cut low-frequency points to not see big dc offset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Data
for SliceInd=1:length(slice)
    if ~P.StitchOnly ==1
        
        for TileInd = 1:length(tilenum)
            filename=strcat(data_filename,num2str(slice(SliceInd)),P.tileN,num2str(tilenum(TileInd)),'_840_');
            %Read in file, zero-pad, and interpolate
            fprintf('Processing tile %d  ...\n', tilenum(TileInd)); %Print out file being processed
            
            [Raw_1,BG,Blines] = Read2024(filename,scan,Parameters);
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
                
                if calcAbsOrientation == 1
                    blin1 = Blines(1:Parameters.alineLength,linenum(Line));
                    blin2 = Blines(Parameters.alineLength+1:end,linenum(Line));
                    if Parameters.windowData == 1
                        blin1 = blin1 .* hammingWindow1;
                        blin2 = blin2 .* hammingWindow2;
                    end
                    [cdp2, cdp1] = InterpandCDP(blin1, blin2, Parameters);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                CDiv = CDP1./CDP2;
                
                CH1 = abs(CDP1).^2;
                CH2 = abs(CDP2).^2;
                reflectivity = CH1 + CH2;
                if dB == 1
                    CH1 = 10*log10(CH1);
                    CH2 = 10*log10(CH2);
                    reflectivity = 10*log10(reflectivity);
                end
                Retardance = (180/pi)*atan(abs(CDP2)./abs(CDP1));
                
                ePi = exp(1i*pi);
                Theta = (ePi./CDiv);
                Theta2 = Theta./(abs(Theta));
                indM1 = CH1>lim1;
                indM2 = CH2>lim2;
                indM3 = indM1+indM2;
                indM = indM3>1;
                
                
                if calcAbsOrientation == 1
                    cdiv = cdp1./cdp2;
                    tempTheta = (ePi./cdiv);
                    RLO = tempTheta(px)/(abs(tempTheta(px)));
                    %ThetaO = Theta2-RLO;
                end
                Theta3 = Theta2.*indM;
                
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
                        Tile_Ot = Tile(st:endc,:,:);
                    end
                    Tile_O(:,:,Line)= Theta2(st:endc,:);
                    Tile_Ot(:,:,Line)= Theta3(st:endc,:);
                end
                if calcAbsOrientation == 1
                    if Line ==1
                        %Tile_AO = Tile;
                        RLOT = zeros(1,scan(end)*Parameters.num_bscans);
                    end
                    %Tile_AO(:,:,Line)= ThetaO(1:endc,:);
                    RLOT(1,Line) = RLO;
                end
            end
            clear Raw_1
            %% Calulate Contrast Enfaces
            if Enface ==1
                ch1Limit = lim1;
                ch2Limit = lim2;
                cut = px-10;
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
                    EnR = CombomaskR4(Tile_ch1,Tile_cross,Tile_R2,ch1Limit,ch2Limit,cut);
                end
                if calcOrientation == 1
                    disp('Calculating Ori Enface')
                    %Tile_O2 = medfilt3(Tile_O);
                    EnO = Combomask4(Tile_ch1,Tile_ch2,Tile_O,ch1Limit,ch2Limit,cut);
                    EnO = 180/pi*(angle(EnO));
                    EnO = EnO/2;
                    if calcAbsOrientation == 1
                        disp('Calculating Abs Ori Enface')
                        Off = (180/pi)*(angle(RLOT)/2);
                        EnAO = (EnO - Off);
                        ind1 = EnAO > 90;
                        ind2 = EnAO < -90;
                        EnAO(ind1) = EnAO(ind1)-180;
                        EnAO(ind2) = EnAO(ind2)+180;
                    end
                end
                
            end
            %%
            if STileComp == 1
                disp('Saving 3D Tile')
                Save_base = strcat(save_directory,CompFolder);
                if calcCDP == 1
                    SaveF = strcat(Save_base,c1);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_CDP1');
                    save(output_filename,'Tile_CDP1');
                    output_filename2 = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_CDP2');
                    save(output_filename2,'Tile_CDP2');
                end
                if calcCh1Ch2 == 1
                    SaveF = strcat(Save_base,c2);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_CH1');
                    save(output_filename,'Tile_ch1');
                    
                end
                if calcReflectivity == 1
                    SaveF = strcat(Save_base,c7);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_Reflect');
                    save(output_filename,'Tile_R');
                    
                end
                if calcCrossPolar == 1
                    SaveF = strcat(Save_base,c6);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_Cross');
                    save(output_filename,'Tile_cross');
                    
                end
                if calcRetardance == 1
                    SaveF = strcat(Save_base,c4);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_R');
                    save(output_filename,'Tile_R2');
                    
                end
                if calcOrientation == 1
                    SaveF = strcat(Save_base,c3);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_Orien');
                    save(output_filename,'Tile_O');
                    
                end
                if calcAbsOrientation == 1
                    SaveF = strcat(Save_base,c5);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_AbsOrien');
                    save(output_filename,'Tile_AO');
                    
                end
            end
            %%
            if SEnface ==1
                disp('Saving Enface')
                Save_base = strcat(save_directory,EnFolder);
                if calcReflectivity == 1
                    SaveF = strcat(Save_base,c7);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnRef');
                    save(output_filename,'EnRef');
                    
                end
                if calcCrossPolar == 1
                    SaveF = strcat(Save_base,c6);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnCr');
                    save(output_filename,'EnCr');
                    
                end
                if calcRetardance == 1
                    SaveF = strcat(Save_base,c4);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnR');
                    save(output_filename,'EnR');
                    
                end
                if calcOrientation == 1
                    SaveF = strcat(Save_base,c3);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnO');
                    save(output_filename,'EnO');
                    
                end
                if calcAbsOrientation == 1
                    SaveF = strcat(Save_base,c5);
                    output_filename = strcat(SaveF,save_n,num2str(slice(SliceInd)),'_tile_',num2str(tilenum(TileInd)),'_EnAO');
                    save(output_filename,'EnAO');
                    
                end
            end
            
        end %tile for loop
    end
    if SStitch ==1
        disp('Stitching Tiles')
        Call_base = strcat(save_directory,EnFolder);
        Save_base = strcat(save_directory,StFolder);
        if calcCrossPolar ==1
            CallF = strcat(Call_base,c6);
            SaveF = strcat(Save_base,c6);
            [TEnCr]= MStitchFCN_mod(slice(SliceInd),6,SaveF,CallF,TileMtrx,alineLength,blineLength,P.overlap);
        end
        if calcReflectivity == 1
            CallF = strcat(Call_base,c7);
            SaveF = strcat(Save_base,c7);
            [TEnRef]= MStitchFCN_mod(slice(SliceInd),7,SaveF,CallF,TileMtrx,alineLength,blineLength,P.overlap);
        end
        if calcRetardance == 1
            CallF = strcat(Call_base,c4);
            SaveF = strcat(Save_base,c4);
            [TEnR]= MStitchFCN_mod(slice(SliceInd),4,SaveF,CallF,TileMtrx,alineLength,blineLength,P.overlap);
        end
        if calcOrientation == 1
            CallF = strcat(Call_base,c3);
            SaveF = strcat(Save_base,c3);
            [TEnO]= MStitchFCN_mod(slice(SliceInd),3,SaveF,CallF,TileMtrx,alineLength,blineLength,P.overlap);
        end
        if calcAbsOrientation == 1
            CallF = strcat(Call_base,c5);
            SaveF = strcat(Save_base,c5);
            [TEnAO]= MStitchFCN_mod(slice(SliceInd),5,SaveF,CallF,TileMtrx,alineLength,blineLength,P.overlap);
        end
    end
end %slice for loop
fprintf('Processing completed ...\n');
