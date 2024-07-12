function [En]= CombomaskCross(dB,dBlimit,endp)

A = size(dB,2);
B = size(dB,3);
LdB = dB(1:endp,:,:);
Tind = LdB>dBlimit;
MdB = LdB.*Tind;
EnOReCal = zeros(A,B);
for x =1:A
    for y = 1:B
        EnOReCal(x,y) = squeeze(mean(MdB(:,x,y)>0));
    end
end
En = EnOReCal;
end