%% Parameter Initialization
%Setting up files
data_directory ='C:\VisPSOCT\Recorded Data\Vis PS-OCT new\Vlad\'; %data called from here
save_directory = 'C:\VisPSOCT\Matlab Codes\Processed_Data\Vlad2\'; %data saved here
Parameter_directory = 'C:\VisPSOCT\Matlab Codes\Processed_Data\Parameters\';

%data_filename = strcat(data_directory,P.baseN);
save_n = 'Vlad_slice_';

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

% % Data size parameters
% XTiles = P.XTiles;
% YTiles = P.YTiles;
% TileMtrx = zeros(XTiles,YTiles);
% m=0;
% for p=1:YTiles
%     for q=1:XTiles
%         TileMtrx(q,p)= 1+m;
%         m = m+1;
%     end
% end

% % Stiching Parameters
% XTiles = 3;
% YTiles = 2;
% TileMtrx= reshape(1:(XTiles*YTiles), [YTiles, XTiles])';


% Data size parameters
slice = 147:147; %Number of slices / slice being analyzed
tilenum = 17:17; % tile number in slice
scan = 1:100; %Number scan in the slice

st = 20; %the pixel depth to start from in enface
endc = 650; %900 %the end pixel depth in enface %max:10
ov = 10;
Flip = 1;
dBlimit = 55; %noise floor

% Load variables
Parameters.num_bscans = 10; %Number of Buffers
Parameters.blineLength = 1000; %Number of A lines
Parameters.alineLength = 2048; %Number of pixels for one Channel
Parameters.alines = (scan(end)*Parameters.num_bscans); %Number of B lines
linenum =1:Parameters.blineLength;
Tile = zeros(endc,Parameters.blineLength,Parameters.alines);

% Data processing actions
% 1 = do the action; 0 = skip the action
Parameters.dispersionComp = 1; %Dispersion compensation, optimizes shape of coherence peak
Parameters.windowData = 1;
Parameters.background = 1;
Parameters.Shift=-0.6; %Shift of the the main peaks between channels
Parameters.CalibShift=-0.64; %Shift of the the Calib peak between channels

% 1 = Calculate Contrast
calcCDP = 1;
calcReflectivity = 1;
calcRetardance = 1;
calcCrossPolar = 1;
calcAbsOrientation = 1;
Enface = 1;
dB = 1;

% 1 = Save
STileComp = 0;
SEnface = 1;
SStitch = 0;
SImg = 0;

%% Dispersion
if Parameters.dispersionComp==1 %need to make 1024
    %Software dispersion compensation creates phase correction vectors using FT shifting properties
    %Main Dispersion
    dispcompfile1 = fullfile(Parameter_directory, 'DispComp_Oct25_Water1to4_NoShift_Windowed_Ch1.dat');
    fid1 = fopen(dispcompfile1,'r');
    angledisp1 = fread(fid1,2048,'real*8');
    fclose(fid1);
    Parameters.PhaseCorrection1 = exp(-1i.*angledisp1);

    dispcompfile2 = fullfile(Parameter_directory, 'DispComp_Oct25_Water1to4_NoShift_Windowed_Ch2.dat');
    fid2 = fopen(dispcompfile2,'r');
    angledisp2 = fread(fid2,2048,'real*8');
    fclose(fid2);
    Parameters.PhaseCorrection2 = exp(-1i.*angledisp2);

    %Calib Dispersion
    dispcompfile1 = fullfile(Parameter_directory, 'DispComp_Oct25_Calib_NoShift_Ch1.dat');
    fid1 = fopen(dispcompfile1,'r');
    angledisp1 = fread(fid1,2048,'real*8');
    fclose(fid1);
    Parameters.CalibCorrection1 = exp(-1i.*angledisp1);

    dispcompfile2 = fullfile(Parameter_directory, 'DispComp_Oct25_Calib_NoShift_Ch2.dat');
    fid2 = fopen(dispcompfile2,'r');
    angledisp2 = fread(fid2,2048,'real*8');
    fclose(fid2);
    Parameters.CalibCorrection2 = exp(-1i.*angledisp2);
end

%%
%Load Interpolation Parameters(constant interpolation wavelengths in k-space)
Parameters.InterpZeroPaddingFactor = 4; %Zero-padding factor for interpolation and for zero-padding the bscans
Parameters.OriginalLineLength = 2048;
Parameters.CDPZeroPaddingFactor = 1; %Zero-padding factor in calculation of CDP
Parameters.AutoPeakCorrCut = 0; %10 %Cut low-frequency points to not see big dc offset
Parameters.Start2=2048+1;

filePath = fullfile(Parameter_directory, 'Interp_parameters.mat');
load(filePath, 'Wavelengths_l','Wavelengths_r','InterpolatedWavelengths','Ks');
Parameters.Wavelengths_l=Wavelengths_l;
Parameters.Wavelengths_r= Wavelengths_r;
Parameters.InterpolatedWavelengths=InterpolatedWavelengths;
Parameters.Ks=Ks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Data
for SliceInd=1:length(slice)
    for TileInd = 1:length(tilenum)
        filename=strcat(data_directory,'Slice_',num2str(slice(SliceInd)),'_Tile_',num2str(tilenum(TileInd)),'_840_');
        fprintf('Processing tile %d  ...\n', tilenum(TileInd)); %Print out file being processed

        %Read in file, zero-pad, and interpolate
        [Raw,~,BG_lines] = VisPSOCT_Read2024(filename,scan,Parameters);
        
        for BLine = 1:length(linenum)
            fprintf('Processing line %d ...\n',linenum(BLine));
            b1 = Raw(:,:,linenum(BLine));
            Ch1 = b1(1:Parameters.alineLength,:);
            Ch2 = b1(Parameters.alineLength+1:end,:);
            
            %Interpolate and Compute complex depth profiles
            [CDP1, CDP2] = VisPSOCT_InterpandCDP(Ch1, Ch2, Parameters);

            %Calibration line
            BG_ch1 = BG_lines(1:Parameters.alineLength,linenum(BLine));
            BG_ch2 = BG_lines(Parameters.alineLength+1:end,linenum(BLine));
            [BG_CDP1, BG_CDP2] = VisPSOCT_InterpandCDP_Calib(BG_ch1, BG_ch2, Parameters);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Amp1 = abs(CDP1);
            Amp2 = abs(CDP2);
            %Reflectivity
            Reflectivity = (Amp1).^2 + (Amp2).^2;
            if dB == 1
                Reflectivity = 10*log10(Reflectivity);
            end

            %Retradance
            Retardance = Reflectivity.* exp(1i*atan(Amp2./Amp1));

            %Axis Orientation
            Weighted_DeltaPh = CDP2.*conj(CDP1);

            if calcAbsOrientation == 1
                a = 750; %general location of peak
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
                    Tile_Ori = Tile(st:endc,:,:);
                end
                Calib_Ori(1,BLine) = Calib_Orientation;
                Tile_Ori(:,:,BLine)= Weighted_DeltaPh(st:endc,:);
            end
        end
        clear Raw
        %% Calulate Contrast Enfaces
        if Enface ==1
            cut = endc-st;
            if calcReflectivity ==1
                disp('Calculating Reflectivity Enface')
                EnRef = CombomaskCrossD(Tile_Ref,dBlimit+1,cut);
            end
            if calcCrossPolar == 1
                disp('Calculating Cross Enface')
                EnCr = CombomaskCrossD(Tile_cross,dBlimit,cut);
            end
            if calcRetardance == 1
                disp('Calculating Retardance Enface')
                EnR = squeeze((180/pi)*angle(sum(Tile_R(1:cut,:,:))));
            end
            if calcAbsOrientation == 1
                disp('Calculating Abs Ori Enface')
                AO_DC_Offset=deg2rad(2*110); %Angle calculated based on the enface axis orientation (True-Read) 
                Tile_Ori_Off = Tile_Ori.*conj(Calib_Ori).*exp(1i*AO_DC_Offset);
                %Tile_Ori_Off = Tile_Ori .* reshape(conj(Calib_Ori), [1 1 1000]);
                EnAO= squeeze((sum(Tile_Ori_Off(1:cut,:,:))));
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
end %slice for loop
fprintf('Processing completed \n');

