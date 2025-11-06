function [Wms] = BlendingMatrix(TileMtrx,TileOverlap,THeight,TWidth)

b1 = size(TileMtrx,1);
b2 = size(TileMtrx,2);

X = TWidth;
Y = THeight;
ov = TileOverlap;

xcry = 1:ov;
ycry = 1/ov*xcry-(1/(2*ov));
yc = flip(ycry);
xc = transpose(yc);
xcf = flip(xc);

%Tile Type linear blending

% Surounded
W = ones(Y,X);
W(1:ov,:) = W(1:ov,:).*xcf;
W(:,1:ov) = W(:,1:ov).*ycry;
W(:,end-ov+1:end) = W(:,end-ov+1:end).*yc;
W(end-ov+1:end,:) = W(end-ov+1:end,:).*xc;
S = W;
% Upper Left Corner
W = ones(Y,X);
W(:,end-ov+1:end) = W(:,end-ov+1:end).*yc;
W(end-ov+1:end,:) = W(end-ov+1:end,:).*xc;
C1 = W;
% Upper edge
W = ones(Y,X);
W(:,1:ov) = W(:,1:ov).*ycry;
W(:,end-ov+1:end) = W(:,end-ov+1:end).*yc;
W(end-ov+1:end,:) = W(end-ov+1:end,:).*xc;
E1 = W;
%Upper Right Corner
W = ones(Y,X);
W(:,1:ov) = W(:,1:ov).*ycry;
W(end-ov+1:end,:) = W(end-ov+1:end,:).*xc;
C2 = W;
% Right Edge
W = ones(Y,X);
W(1:ov,:) = W(1:ov,:).*xcf;
W(:,1:ov) = W(:,1:ov).*ycry;
W(end-ov+1:end,:) = W(end-ov+1:end,:).*xc;
E2 = W;
% Lower Right Corner
W = ones(Y,X);
W(1:ov,:) = W(1:ov,:).*xcf;
W(:,1:ov) = W(:,1:ov).*ycry;
C3 = W;
% Lower edge
W = ones(Y,X);
W(1:ov,:) = W(1:ov,:).*xcf;
W(:,1:ov) = W(:,1:ov).*ycry;
W(:,end-ov+1:end) = W(:,end-ov+1:end).*yc;
E3 = W;
% Lower Left Corner
W = ones(Y,X);
W(1:ov,:) = W(1:ov,:).*xcf;
W(:,end-ov+1:end) = W(:,end-ov+1:end).*yc;
C4 = W;
% Left Edge
W = ones(Y,X);
W(1:ov,:) = W(1:ov,:).*xcf;
W(:,end-ov+1:end) = W(:,end-ov+1:end).*yc;
W(end-ov+1:end,:) = W(end-ov+1:end,:).*xc;
E4 = W;

for i = 1:b1
    for j= 1:b2
 if i == 1&&j == 1
        Wms{i,j} = C1;
 elseif i == 1&& j == b2
     Wms{i,j} = C2;
 elseif i == b1&&j == b2
     Wms{i,j} = C3;
 elseif i == b1&&j ==1
     Wms{i,j} = C4;

 elseif i == 1
     Wms{i,j} = E1;
 elseif j == 1
     Wms{i,j} = E4;
 elseif i == b1
     Wms{i,j} = E3;
 elseif j == b2
     Wms{i,j} = E2;
 else
     Wms{i,j} = S;
 end
    end
end

end


