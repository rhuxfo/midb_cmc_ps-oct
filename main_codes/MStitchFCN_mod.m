% contrast = 1,2,3,4,5,6, or 7
%1 : Chl
%2 : Ch2
%3 : Orientation
%4 : Retardance
%5 : Absolute Orientation
%6 : Cross Polarization
%7 : Reflectivity

function [EnStitch] = MStitchFCN_mod(slice,contrast,SaveFolder,Directory,TileMtrx,alines,blines,overlap)
Dir = Directory;
SaveDir = SaveFolder;

slicenum = slice;
n = contrast;

b1 = size(TileMtrx,1);
b2 = size(TileMtrx,2);
xp = 18;

X = blines;
Y = alines;
ov = round(((overlap/100)*blines));

A1 = (Y-ov);
A2 = (X-ov);
%c = A1/b1;
B1 = A1+ov;
B2 = A2+ov;
D1 = (X*b2);
D2 = (Y*b1);
d1 = (A2*b2)+ov;
d2 = (A1*b1)+ov;
Wms = BlendingMatrix(TileMtrx,ov,Y,X);

for s = 1:length(slicenum)
    if n ==1
        ch = 'CH1';
    elseif n ==2
        ch = 'CH2';
    elseif n ==3
        ch = 'EnO';
    elseif n ==4
        ch = 'EnR';
    elseif n == 5
        ch = 'EnAO';
    elseif n == 6
        ch = 'EnCr';
    elseif n == 7
        ch = 'EnRef';
    end

    %%
    for i = 1:b1
        for j = 1:b2
            filename = strcat(Dir,'slice_',num2str(slicenum(s)),'_Tile_',num2str(TileMtrx(i,j)),'_',ch);
            Sname = strcat(SaveDir,'Slice_',num2str(slicenum(s)),'_',ch);
            T = load(filename);

            if n == 3
                Tile = T.EnO;
            elseif n == 4
                Tile = T.EnR;
            elseif n == 5
                Tile = T.EnAO;
            elseif n == 6
                Tile = T.EnCr;
            elseif n == 7
                Tile = T.EnRef;
            end

            if i < 2
                ims{i,j} = Tile;
            else
                En = Tile;
                Epr = ims{i-1,j};
                En(1:xp+1,:) = Epr(end-ov:end-ov+xp,:);
                ims{i,j}= En;
            end
        end
    end

    % Merge tiles:
    ImR = zeros(D2,D1);
    for i1 = 1:b1
        for i2 = 1:b2
            ImR((A1*(i1-1)+(1:B1)),(A2*(i2-1)+(1:B2))) = ImR((A1*(i1-1)+(1:B1)),(A2*(i2-1)+(1:B2))) + (Wms{i1,i2}).*(ims{i1,i2});
        end
    end
    Stitched = ImR(1:d2,1:d1);
    if n == 3
        TEnO = rot90(Stitched);
        EnStitch = TEnO;
        save(Sname,"TEnO");
    elseif n == 4
        TEnR = rot90(Stitched);
        EnStitch = TEnR;
        save(Sname,"TEnR");
    elseif n == 5
        TEnAO = rot90(Stitched);
        EnStitch = TEnAO;
        save(Sname,"TEnAO");
    elseif n == 6
        TEnCr = rot90(Stitched);
        EnStitch = TEnCr;
        save(Sname,"TEnCr");
    elseif n == 7
        TEnRef = rot90(Stitched);
        EnStitch = TEnRef;
        save(Sname,"TEnRef");
    end
end
end