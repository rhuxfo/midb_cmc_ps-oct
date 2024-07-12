function [EnR]= CombomaskR4(dB1,dB2,Input,dB1limit,dB2limit,endp)

a = size(dB1,1);
b = size(dB1,2);
c = size(dB1,3);

L1 = dB1limit;
L2 = dB2limit;

MRe = zeros(a,b,c);
for i = 1:c
    for j = 1:b
        for k = 1:a
            if dB1(k,j,i) > L1||dB2(k,j,i) > L2
                    MRe(k,j,i) = Input(k,j,i);
            else

            end
        end
    end
end
%%

MRr = MRe(1:endp,:,:);
MRr2 = squeeze(mean(MRr,1));

EnR = medfilt2(MRr2); 
end