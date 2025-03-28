function [En]= CombomaskCrossD(dB,dBlimit,endp)

A = size(dB,2);
B = size(dB,3);
LdB = dB(1:endp,:,:);
Tind = LdB>dBlimit;
MdB = LdB.*Tind;
EnOReCal = zeros(A,B)+dBlimit;
for x =1:A
    for y = 1:B
        el = nnz(MdB(:,x,y));
        EnOReCal(x,y) = squeeze(sum(MdB(:,x,y)))/el;
    end
end
En = EnOReCal;
end
