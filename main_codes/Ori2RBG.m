function TempOri = Ori2RBG(TEnAO,TEnR)

d1 = size(TEnR,1);
d2 = size(TEnR,2);
c = hsv(256);
disp(size(c,1))
Ang = angle(TEnAO)/2;
Ang2 = rescale(Ang,1,256);
Col = zeros(d1,d2,3);
disp(d1)
disp(d2)
for i = 1:d1
    for j = 1:d2 
    disp(i)
    disp(j)
        ind = round(Ang2(i,j)); 
        disp(ind)
        Col(i,j,:) = c(ind,:);
    end
end

TR0s = LimdB2D(.95,.05,rescale(TEnR));
TR0s = rescale(TR0s);

TempOri = Col.*TR0s;

end
