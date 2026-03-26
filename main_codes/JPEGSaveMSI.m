function Out1 = JPEGSaveMSI(Dir,Sdir,slices,Flip)

for i = 1:length(slices)
    f = slices(i);
    disp(f)
    fileNaRef = strcat('Slice_',num2str(f,'%03.f'),'_EnRef.mat');
    fileNaRef = fullfile(Dir,'Reflectivity',fileNaRef);
    fileNaR = strcat('Slice_',num2str(f,'%03.f'),'_EnR.mat');
    fileNaR = fullfile(Dir,'Retardance',fileNaR);
    fileNaCr = strcat('Slice_',num2str(f,'%03.f'),'_EnCr.mat');
    fileNaCr = fullfile(Dir,'Cross',fileNaCr);
    fileNaOri = strcat('Slice_',num2str(f,'%03.f'),'_EnAO.mat');
    fileNaOri = fullfile(Dir,'Orientation',fileNaOri);

    SaveNaRef = strcat('Slice_',num2str(f,'%03.f'),'_Ref.jpeg');
    SaveNaRef = fullfile(Sdir,'Reflectivity',SaveNaRef);
    SaveNaR = strcat('Slice_',num2str(f,'%03.f'),'_R.jpeg');
    SaveNaR = fullfile(Sdir,'Retardance',SaveNaR);
    SaveNaCr = strcat('Slice_',num2str(f,'%03.f'),'_Cr.jpeg');
    SaveNaCr = fullfile(Sdir,'Cross',SaveNaCr);
    SaveNaOri = strcat('Slice_',num2str(f,'%03.f'),'_Ori.jpeg');
    SaveNaOri = fullfile(Sdir,'Orientation',SaveNaOri);
    disp(SaveNaOri)

    load(fileNaRef)
    load(fileNaR)
    load(fileNaCr)
    load(fileNaOri)
    
    [N,E] = histcounts(TEnRef);
    [~,idx]= max(N(5:end));
    A = ceil(E(1,idx+5));
    
    clear N E

    if Flip == 1
        TEnRef = flip(TEnRef);
        TEnR = flip(TEnR);
        TEnCr = flip(TEnCr);
        TEnAO = flip(TEnAO);
    end
        TempRef = rescale(LimdB2D(A+12,A,TEnRef));
        TempCr = rescale(LimdB2D(A,55,TEnCr));
        TEnR = rad2deg(angle(TEnR));
        TEnR3 = LimdB2D(55,5,TEnR);
        TempR = rescale(TEnR3);
        
    % C = zeros(100);
    % for q = 1:100
    %     j = q/100;
    %     C(q) = sum(TEnR4 >j,'all');
    % end

    %C = rescale(C);
    %ind1 = C>.0005;
    %[~,Inx] = min(C(ind1));
    %ind2 = TEnR4>(Inx/100);
    %TEnR4(ind2) = 0;
    %Mt1 = TEnRef>A;
    TEnAO2 = (TEnAO.^-1)*(exp(1i*deg2rad(-270)));
    TempOri = Ori2RBG(TEnAO2,TempR);
    Out1 = tiff23(TempRef,SaveNaRef,1,Flip);
    Out2 = tiff23(TempCr,SaveNaCr,1,Flip);
    Out3 = tiff23(TempR,SaveNaR,1,Flip);
    imwrite(TempOri,SaveNaOri);
    
end
fprintf('Complete')
end















