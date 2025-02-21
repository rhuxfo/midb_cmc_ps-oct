function A = gifStack(CallDir,SaveDir,SliceNum)

Dir = CallDir;
Sdir = SaveDir;
slices = SliceNum;

for i = 1:length(slices)
    f = slices(i);
    fileNaRef = strcat(Dir,'Slice_',num2str(f),'_EnRef.mat');

    load(fileNaRef)

    Temp = rescale(LimdB2D(.55,.2,TEnRef));
    SaveNa = strcat(Sdir,'Slice_',num2str(f),'_Ref.jpeg');
    A = tiff23(Temp,SaveNa,1);

end

end
