function [EnO]= Combomask4(dB1,dB2,CD3,dB1limit,dB2limit,endp)

A = size(dB1,2);
B = size(dB1,3);
LdB1 = dB1(1:endp,:,:);
LdB2 = dB2(1:endp,:,:);

Tile_M1 = LdB1>dB1limit;
Tile_M2 = LdB2>dB2limit;
Tile_M3 = Tile_M1+Tile_M2;
Tile_indM = Tile_M3>1;

EnOReCal = zeros(A,B);

for x =1:A
    for y = 1:B
        GAngs = CD3(Tile_indM(:,x,y),x,y);
        TF = isempty(GAngs);
        if TF == 1
        else
            EnOReCal(x,y) = mean(GAngs);
        end
    end
end
EnO = EnOReCal;
end