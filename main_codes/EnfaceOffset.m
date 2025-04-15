function Temp2 = EnfaceOffset(En,slice)
Temp2 = zeros(8500,11000);
slicenum = slice;

if slicenum < 20
    Offsetx = 2718;
    Offsety = 3600;
elseif slicenum > 19 && slicesnum < 33
    Offsetx = 1818;
    Offsety = 2700;
elseif slicesnum > 32 && slicesnum < 46
    Offsetx = 1818;
    Offsety = 1800;
elseif slicesnum > 45 && slicesnum < 70
    Offsetx = 918;
    Offsety = 1800;
elseif slicesnum > 69 && slicesnum < 90
    Offsetx = 918;
    Offsety = 900;
elseif slicesnum > 89 && slicesnum < 179
    Offsetx = 531;
    Offsety = 900;
elseif slicesnum > 178 && slicesnum < 202
    Offsetx = 531;
    Offsety = 0;
elseif slicesnum > 201 && slicesnum < 232
    Offsetx = 0;
    Offsety = 0;
elseif slicesnum > 231 && slicesnum < 270
    Offsetx = 0;
    Offsety = 900;
elseif slicesnum >269  && slicesnum < 285
    Offsetx = 531;
    Offsety = 900;
elseif slicesnum > 284 && slicesnum < 307
    Offsetx = 531;
    Offsety = 1800;
elseif slicesnum > 306 && slicesnum < 313
    Offsetx = 918;
    Offsety = 1800;
elseif slicesnum > 312 && slicesnum < 326
    Offsetx = 1818;
    Offsety = 2700;
elseif slicesnum > 325
    Offsetx = 2718;
    Offsety = 2700;
end

D1 = Offsetx+1;%size(Temp2,1) - (Offsetx + size(TEnR6,1));
D2 = Offsety+1;

% D1 = (size(Temp2,1)-size(TEnR6,1))/2;
% D2 = (size(Temp2,2)-size(TEnR6,2))/2;
Temp2(D1:(D1+size(En,1)-1),D2:(D2+size(En,2)-1)) = En;

end
