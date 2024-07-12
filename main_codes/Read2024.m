function [Raw_Tile,BG,BGlines] = Read2024(filename_temp,scan,Parameters)

Raw_Tile = zeros(Parameters.alineLength*2, Parameters.blineLength, Parameters.num_bscans*100);

for ScanInd=1:length(scan)
    fileName=strcat(filename_temp,num2str(scan(ScanInd)),'.dat');
    %fprintf('Processing file %s ...\n', fileName); %Print out file being processed

    filePointer = fopen(fileName, 'r', 'b');
    headerStr = fgetl(filePointer);
    evalc(headerStr);
    Spectra_temp = fread(filePointer, 'int16');
    all_bscans = reshape(Spectra_temp,[Parameters.alineLength*2 Parameters.blineLength Parameters.num_bscans]); %3D array of 2D cross-sectional images, 2048
    clear Spectra_temp
    fclose(filePointer);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = scan(ScanInd)*10;
    b = a-9;
    Raw_Tile(:,:,b:a) = all_bscans;    
    BGlines(:,b:a) = mean(all_bscans,2);
end
BG = mean(Raw_Tile,[2 3]);
Raw_Tile =Raw_Tile-BG;
end