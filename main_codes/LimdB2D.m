function [D]=LimdB2D(Upperlimit,Lowerlimit,dB)
UL = Upperlimit;
LL = Lowerlimit;
D = dB;

for j = 1:size(D,1)
    for i = 1:size(D,2)
        if D(j,i)>UL
            D(j,i)=UL;
        elseif D(j,i)<LL
            D(j,i) = LL;
        else
        end
    end
end
end
