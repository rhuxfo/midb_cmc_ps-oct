%TEnAO3 = zeros(size(TEnAO,1),size(TEnAO,2),10);
%EnBG = sum(BGAvg,3)+sum(BGAvg,3)+sum(BGAvg,3)+sum(BGAvg,3);
function TEnAOBG = MStitchFCN_mod_sub_out(EnAO1, TMtrx,f)
k = ones(20,20);
blines = 1000;
overlap = 10;
alines = 1000;
%for w = 1:10
TileMtrx = TMtrx;
EnBG = EnAO1.*10;%sum(BGAvg,3).*(1+(w/4));
EnAO3 = convn(EnBG,k,'same')./convn(ones(size(EnBG)),k,'same');

%%
FLIP = f;

b1 = size(TileMtrx,1);
b2 = size(TileMtrx,2);

X = blines;
Y = alines;
ov = round(((overlap/100)*blines));

A1 = (Y-ov);
A2 = (X-ov);
D1 = (X+A2*b2);
D2 = (Y+A1*b1);
d1 = (A2*b2)+ov;
d2 = (A1*b1)+ov;
Wms = BlendingMatrix(TileMtrx,ov,Y,X);

    for i = 1:b1
        for j = 1:b2

            Tile = EnAO3;

            if FLIP == 1
                Tile = flip(Tile,2);
            else
            end
                ims{i,j} = Tile;
            
        end
    end
    % Merge tiles:
    ImR = zeros(D2,D1);
    for i1 = 1:b1
        for i2 = 1:b2
            ImR((A1*(i1-1)+(1:Y)),(A2*(i2-1)+(1:X))) = ImR((A1*(i1-1)+(1:Y)),(A2*(i2-1)+(1:X))) + (Wms{i1,i2}).*(ims{i1,i2});
        end
    end
    Stitched = ImR(1:d2,1:d1);
        TEnAOBG = rot90(Stitched);

% TEnAO2 = TEnAO - TEnAOBG;
% %TEnAO3(:,:,w) = TEnAO2;
% 
% %end
% 
% figure;
% %for v = 1:10
% imagesc(angle(TEnAO2)); colormap hsv; alpha(rescale(TEnR5)); set(gca,'color','black')
% %title(v)
% %pause(.1)
% %end

