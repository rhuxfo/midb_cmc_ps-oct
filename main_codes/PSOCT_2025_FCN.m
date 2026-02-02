function [status] = PSOCT_2025_FCN(P)
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
ImgFolder = 'jpegs'; % stitched image files (jpeg)

%contrast folders
c1 = 'CDP';
c4 = 'Retardance';
c5 = 'AbsoOri';
c6 = 'Cross';
c7 = 'Reflectivity';
folderNames = {c1,c4,c5,c6,c7};

%automatically create save folders if they do not exist
if P.autofolder ==1 
    folder3d = fullfile(P.Sdir,CompFolder);
    folderEnface = fullfile(P.Sdir,EnFolder);
    folderStitch = fullfile(P.Sdir,StFolder);
    folderImg = fullfile(P.Sdir,ImgFolder);
    
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
    if ~exist(folderImg,'dir')
        mkdir(folderImg);
        for i=1:length(folderNames)
            mkdir(fullfile(folderImg,folderNames{i}));
        end
    end
end

% Data size parameters
XTiles = P.XTiles;
YTiles = P.YTiles; 

TileMtrx= reshape(1:(XTiles*YTiles), [YTiles, XTiles])';

% Data size parameters
slice = P.Slices; %Number of slices / slice being analyzed
tilenum = P.tiles; % tile number in slice
scan = P.buffers; %Number scan in the slice

st = P.depthstart;
endc = P.depthcut;
ov = P.overlap;
Flip = P.Flip;
dBlimit = P.NoiseCut;
Ret_noise_level = P.NoiseCut;

% Load variables
Pname = strcat(data_filename, num2str(slice(1)),P.tileN,num2str(tilenum(1)),'_840_1.dat');
filePointer = fopen([data_filename, num2str(slice(1)),P.tileN,num2str(tilenum(1)),'_840_1.dat'], 'r', 'l');
fprintf('Calling');
disp(Pname)
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
%OT = zeros(255,104);
% Data processing actions
% 1 = do the action; 0 = skip the action
Parameters.dispersionComp = P.disper; %Dispersion compensation, optimizes shape of coherence peak
Parameters.windowData = P.wind;
Parameters.background = P.BGremoval;

% 1 = Calculate Contrast
calcCDP = 0;
calcReflectivity = P.Flect;
calcRetardance = P.Retar;
calcCrossPolar = P.Cr;
calcAbsOrientation = P.AbOrio;
Enface = P.En;
dB = 1;

% 1 = Save
STileComp = P.TCsv;
SEnface = P.Ensv;
SStitch = P.Stsv;
SImg = P.img;

if Parameters.background ==1
    filename = P.BG;
    load(filename);
    k = ones(20,20);
EnBG = BG2;
EnAO3 = convn(EnBG,k,'same')./convn(ones(size(EnBG)),k,'same');
end
%% Dispersion
if Parameters.dispersionComp==1 %need to make 1024
    %Software dispersion compensation creates phase correction vectors
    %using FT shifting properties
    dispcompfile1 = fullfile(P.DCf1);
    fid1 = fopen(dispcompfile1,'r');
    angledisp1 = fread(fid1,1024,'real*8');
    fclose(fid1);
    Parameters.PhaseCorrection2 = exp(-1i.*angledisp1);
    
    dispcompfile2 = fullfile(P.DCf2);
    fid2 = fopen(dispcompfile2,'r');
    angledisp2 = fread(fid2,1024,'real*8');
    fclose(fid2);
    Parameters.PhaseCorrection1 = exp(-1i.*angledisp2);
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
[Parameters.Wavelengths_l, Parameters.Wavelengths_r, Parameters.InterpolatedWavelengths, Ks] = InterpolateWavelengths4(InterpolationParameters); %left and right wavelengths refer to different polarization channels
Parameters.AutoPeakCorrCut = 10; %Cut low-frequency points to not see big dc offset
if Parameters.dispersionComp == 2
    sh = 40.5/(2*pi);
    zm = 2.35*10^-3;
    Arg = Ks*(sh)*(zm)/512;
    Parameters.PhaseCorrection1 = exp(-1i.*Arg);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Data
for SliceInd=1:length(slice)
    for TileInd = 1:length(tilenum)
        filename=strcat(data_directory,'Slice_',num2str(slice(SliceInd)),'_Tile_',num2str(tilenum(TileInd)),'_840_');
        fprintf('Processing tile %d  ...\n', tilenum(TileInd)); %Print out file being processed

        %Read in file, zero-pad, and interpolate
        [Raw,~,BG_lines] = Read2024(filename,scan,Parameters);
        
        for BLine = 1:length(linenum)
            fprintf('Processing line %d ...\n',linenum(BLine));
            b1 = Raw(:,:,linenum(BLine));
            Ch1 = b1(1:Parameters.alineLength,:);
            Ch2 = b1(Parameters.alineLength+1:end,:);
            
            %Interpolate and Compute complex depth profiles
            [CDP2, CDP1] = InterpandCDP2(Ch1, Ch2, Parameters);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Amp1 = abs(CDP1);
            Amp2 = abs(CDP2);
            
            %Reflectivity
            Reflectivity = (Amp1).^2 + (Amp2).^2;

            %if dB == 1
            Reflectivity = 10*log10(Reflectivity);
            Retardance = (Reflectivity-Ret_noise_level).*(Reflectivity>=Ret_noise_level).* exp(1i*atan(Amp2./Amp1));
            Amp1 = 10*log10(Amp1.^2);
            Amp2 = 10*log10(Amp2.^2); 
           % end

            %Retardance = Reflectivity.* exp(1i*atan(Amp2./Amp1));
            %Axis Orientation
            Weighted_DeltaPh = CDP2.*conj(CDP1);
            Weighted_DeltaPh_min = (Weighted_DeltaPh./(Amp1.*Amp2)).*(min(Amp1,Amp2).^2);%%%%

            if calcAbsOrientation == 1
             %Calibration line
                BG_ch1 = BG_lines(1:Parameters.alineLength,linenum(BLine));
                BG_ch2 = BG_lines(Parameters.alineLength+1:end,linenum(BLine));
                [BG_CDP2, BG_CDP1] = InterpandCDP2(BG_ch1, BG_ch2, Parameters);
                a = 180; %general location of peak
                Calib_Reflectivity = abs(BG_CDP1).^2 +abs(BG_CDP2).^2;
                [~,Calib_Loc] = max(Calib_Reflectivity(a:a+100));
                Calib_pixel = Calib_Loc+a-1;
                Calib_Weighted_DeltaPh = BG_CDP2.*conj(BG_CDP1);
                Calib_Orientation = Calib_Weighted_DeltaPh(Calib_pixel)./abs(Calib_Weighted_DeltaPh(Calib_pixel));
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if calcCDP == 1
                if BLine ==1
                    Tile_CDP1 = Tile(st:endc,:,:);
                    Tile_CDP2 = Tile(st:endc,:,:);
                end
                Tile_CDP1(:,:,BLine)= CDP1(st:endc,:);
                Tile_CDP2(:,:,BLine)= CDP2(st:endc,:);
            end

            if calcReflectivity ==1
                if BLine ==1
                    Tile_Ref = Tile(st:endc,:,:);
                end
                Tile_Ref(:,:,BLine)= Reflectivity(st:endc,:);
            end
            if calcRetardance == 1
                if BLine ==1
                    Tile_R = Tile(st:endc,:,:);
                end
                Tile_R(:,:,BLine)= Retardance(st:endc,:);
            end
            if calcCrossPolar == 1
                if BLine ==1
                    Tile_cross = Tile(st:endc,:,:);
                end
                Tile_cross(:,:,BLine)= Amp2(st:endc,:);
            end

            if calcAbsOrientation == 1
                if BLine ==1
                    Calib_Ori = zeros(1,scan(end)*Parameters.num_bscans);
                    %Tile_Ori = Tile(st:endc,:,:);
                    Tile_Ori_min = Tile(st:endc,:,:);
                end
                Calib_Ori(1,BLine) = Calib_Orientation;
                %Tile_Ori(:,:,BLine)= Weighted_DeltaPh(st:endc,:);
                Tile_Ori_min(:,:,BLine)= Weighted_DeltaPh_min(st:endc,:);
            end
        end
        clear Raw
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calulate Contrast Enfaces
        if Enface ==1
            cut = 185;
            if calcReflectivity ==1
                disp('Calculating Reflectivity Enface')
                EnRef = CombomaskCrossD(Tile_Ref,dBlimit+10,cut);
            end
            if calcCrossPolar == 1
                disp('Calculating Cross Enface')
                EnCr = CombomaskCrossD(Tile_cross,dBlimit,cut);
            end
            if calcRetardance == 1
                disp('Calculating Retardance Enface')
                %EnR = squeeze((180/pi)*angle(sum(Tile_R(1:cut,:,:))));
                EnR = squeeze(sum(Tile_R(1:cut,:,:)));
            end
            if calcAbsOrientation == 1
                disp('Calculating Abs Ori Enface')
                AO_DC_Offset=deg2rad(2*21); %Angle calculated based on the enface axis orientation (True-Read) 
                %Tile_Ori_Off = Tile_Ori.*conj(Calib_Ori).*exp(1i*AO_DC_Offset);
                Tile_Ori_Off_min = Tile_Ori_min.*conj(Calib_Ori).*exp(1i*AO_DC_Offset);
                %EnAO= squeeze((sum(Tile_Ori_Off(1:cut,:,:))));
                EnAO= squeeze((sum(Tile_Ori_Off_min(1:cut,:,:))));
                %Ori_test = squeeze(mean(Tile_Ori_min,2));
                %OT(:,tilenum(TileInd)) = squeeze(mean(Ori_test,2));
                % if tilenum(TileInd) == tilenum(end)
                %    TEnAOBG = MStitchFCN_mod_sub_out(EnAO,TileMtrx,Flip);
                % end
            end

        end
        %%
        if STileComp == 1
            disp('Saving 3D Tile')
            Save_base = strcat(save_directory,CompFolder);
            if calcCDP == 1
                ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_CDP1');
                SaveFN = fullfile(Save_base,c1,ON);
                save(SaveFN,'Tile_CDP1','-v7.3','-nocompression');
                ON2 = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_CDP2');
                SaveFN2 = fullfile(Save_base,c1,ON2);
                save(SaveFN2,'Tile_CDP2','-v7.3','-nocompression');
            end

            if calcReflectivity == 1
                ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_Reflect');
                SaveFN = fullfile(Save_base,c7,ON);
                save(SaveFN,'Tile_Ref','-v7.3','-nocompression');
            end

            if calcCrossPolar == 1
                ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_Cross');
                SaveFN = fullfile(Save_base,c6,ON);
                save(SaveFN,'Tile_cross','-v7.3','-nocompression');
            end

            if calcRetardance == 1
                ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_R');
                SaveFN = fullfile(Save_base,c4,ON);
                save(SaveFN,'Tile_R','-v7.3','-nocompression');
            end

            if calcAbsOrientation == 1
                ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_AbsOrien');
                SaveFN = fullfile(Save_base,c5,ON);
                save(SaveFN,'Tile_Ori','-v7.3','-nocompression');
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
        [TEnCr]= MStitchFCN_Vlad(slice(SliceInd),6,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,EnAO3);
    end
    if calcReflectivity == 1
        CallF = fullfile(Call_base,c7);
        SaveF = fullfile(Save_base,c7);
        [TEnRef]= MStitchFCN_Vlad(slice(SliceInd),7,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,EnAO3);
    end
    if calcRetardance == 1
        CallF = fullfile(Call_base,c4);
        SaveF = fullfile(Save_base,c4);
        [TEnR]= MStitchFCN_Vlad(slice(SliceInd),4,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,EnAO3);
    end
    if calcAbsOrientation == 1
        CallF = fullfile(Call_base,c5);
        SaveF = fullfile(Save_base,c5);
        [TEnAO]= MStitchFCN_Vlad(slice(SliceInd),5,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,EnAO3);
    end
    status = 2;
end
if SImg == 1
    Save_base = strcat(save_directory,ImgFolder);
    if calcReflectivity == 1
        ON = strcat(save_n,num2str(slice(SliceInd)),'_Ref.jpeg');
        SaveF = fullfile(Save_base,c7,ON);
        Temp = rescale(LimdB2D(75,57,TEnRef));
        Out = tiff23(Temp,SaveF,1,Flip);
    end
    if calcCrossPolar == 1
        ON = strcat(save_n,num2str(slice(SliceInd)),'_Cr.jpeg');
        SaveF = fullfile(Save_base,c6,ON);
        Temp = rescale(LimdB2D(60,5,TEnCr));
        Out = tiff23(Temp,SaveF,1,Flip);
    end
    if calcRetardance == 1
        ON = strcat(save_n,num2str(slice(SliceInd)),'_R.jpeg');
        SaveF = fullfile(Save_base,c4,ON);
        Temp = rescale(LimdB2D(60,5,TEnR));
        Out = tiff23(TEnR,SaveF,1,Flip);
    end
end %slice for loop
fprintf('Processing completed \n');
end




