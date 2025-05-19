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
ImgFolder = 'jpegs'; % stitched image files (jpeg)

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
Flip = P.Flip;
dBlimit = P.NoiseCut;

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
SImg = P.img;

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
                
                Tile_Reff= abs(CDP1)+ abs(CDP2);
                Retardance = Tile_Reff.* exp(1i*atan(abs(CDP2)./abs(CDP1)));
                %Retardance = (180/pi)*atan(abs(CDP2)./abs(CDP1));
                
                if calcOrientation == 1
                    ePi = exp(1i*pi);
                    Theta = (ePi./CDiv);
                    Theta2 = Theta./(abs(Theta));
                end
                if calcAbsOrientation == 1
                    a = 180;
                    tempAmp2 = abs(cdp2).^2;
                    [~,P2] = max(tempAmp2(a:a+100));
                    mmm = cdp1.*conj(cdp2);
                    px2 = P2+a-1;
                    RLO2 = mmm(px2)/abs(mmm(px2));
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
                ch1Limit = dBlimit;
                ch2Limit = dBlimit;
                cut = 160;
                if calcReflectivity ==1
                    disp('Calculating Reflectivity Enface')
                    EnRef = CombomaskCrossD(Tile_R,ch1Limit+1,cut);
                end
                if calcCrossPolar == 1
                    disp('Calculating Cross Enface')
                    EnCr = CombomaskCrossD(Tile_cross,ch2Limit,cut);
                end
                if calcRetardance == 1
                    disp('Calculating Retardance Enface')
                    EnR = squeeze((180/pi)*angle(sum(Tile_R2(1:cut,:,:))));
                end
                if calcOrientation == 1
                    disp('Calculating Ori Enface')
                    EnO = Combomask4(Tile_ch1,Tile_ch2,Tile_O,ch1Limit,ch2Limit,cut);
                end
                if calcAbsOrientation == 1
                    disp('Calculating Abs Ori Enface')
                    
                    Tile_Off = Tile_Om.*conj(RLO2T);
                    EnAO= squeeze((sum(Tile_Off(1:cut,:,:))));
                    if tilenum(TileInd) ==1
                        TEnAOBG = MStitchFCN_mod_sub_out(EnAO,TileMtrx,Flip);
                    end
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
                if calcCh1Ch2 == 1
                    ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_CH1');
                    SaveFN = fullfile(Save_base,c2,ON);
                    save(SaveFN,'Tile_ch1','-v7.3','-nocompression');
                end
                if calcReflectivity == 1
                    ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_Reflect');
                    SaveFN = fullfile(Save_base,c7,ON);
                    save(SaveFN,'Tile_R','-v7.3','-nocompression');
                    
                end
                if calcCrossPolar == 1
                    
                    ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_Cross');
                    SaveFN = fullfile(Save_base,c6,ON);
                    save(SaveFN,'Tile_cross','-v7.3','-nocompression');
                    
                end
                if calcRetardance == 1
                    ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_R');
                    SaveFN = fullfile(Save_base,c4,ON);
                    save(SaveFN,'Tile_R2','-v7.3','-nocompression');
                    
                end
                if calcOrientation == 1
                    ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_Orien');
                    SaveFN = fullfile(Save_base,c3,ON);
                    save(SaveFN,'Tile_O','-v7.3','-nocompression');
                    
                end
                if calcAbsOrientation == 1
                    ON = strcat(save_n,num2str(slice(SliceInd),'%03.f'),'_tile_',num2str(tilenum(TileInd),'%03.f'),'_AbsOrien');
                    SaveFN = fullfile(Save_base,c5,ON);
                    save(SaveFN,'Tile_Om','-v7.3','-nocompression');
                    
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
            [TEnCr]= MStitchFCN_mod2(slice(SliceInd),6,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,TEnAOBG);
        end
        if calcReflectivity == 1
            CallF = fullfile(Call_base,c7);
            SaveF = fullfile(Save_base,c7);
            [TEnRef]= MStitchFCN_mod2(slice(SliceInd),7,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,TEnAOBG);
        end
        if calcRetardance == 1
            CallF = fullfile(Call_base,c4);
            SaveF = fullfile(Save_base,c4);
            [TEnR]= MStitchFCN_mod2(slice(SliceInd),4,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,TEnAOBG);
        end
        if calcOrientation == 1
            CallF = fullfile(Call_base,c3);
            SaveF = fullfile(Save_base,c3);
            [TEnO]= MStitchFCN_mod2(slice(SliceInd),3,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,TEnAOBG);
        end
        if calcAbsOrientation == 1
            CallF = fullfile(Call_base,c5);
            SaveF = fullfile(Save_base,c5);
            [TEnAO]= MStitchFCN_mod2(slice(SliceInd),5,SaveF,CallF,TileMtrx,blineLength,Parameters.alines,ov,Flip,TEnAOBG);
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
    end
end %slice for loop
fprintf('Processing completed \n');
end
