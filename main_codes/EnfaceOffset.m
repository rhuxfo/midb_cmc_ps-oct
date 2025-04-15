function Temp2 = EnfaceOffset(En,slicenum)
Temp2 = zeros(8500,11000);

if slicenum < 20
    Offsetx = 2718;
    Offsety = 3600;
elseif slicenum > 19 && slicenum < 33
    Offsetx = 1818;
    Offsety = 2700;
elseif slicenum > 32 && slicenum < 46
    Offsetx = 1818;
    Offsety = 1800;
elseif slicenum > 45 && slicenum < 70
    Offsetx = 918;
    Offsety = 1800;
elseif slicenum > 69 && slicenum < 90
    Offsetx = 918;
    Offsety = 900;
elseif slicenum > 89 && slicenum < 179
    Offsetx = 531;
    Offsety = 900;
elseif slicenum > 178 && slicenum < 202
    Offsetx = 531;
    Offsety = 0;
elseif slicenum > 201 && slicenum < 232
    Offsetx = 0;
    Offsety = 0;
elseif slicenum > 231 && slicenum < 270
    Offsetx = 0;
    Offsety = 900;
elseif slicenum >269  && slicenum < 285
    Offsetx = 531;
    Offsety = 900;
elseif slicenum > 284 && slicenum < 307
    Offsetx = 531;
    Offsety = 1800;
elseif slicenum > 306 && slicenum < 313
    Offsetx = 918;
    Offsety = 1800;
elseif slicenum > 312 && slicenum < 326
    Offsetx = 1818;
    Offsety = 2700;
elseif slicenum > 325
    Offsetx = 2718;
    Offsety = 2700;
end

D1 = Offsetx+1;%size(Temp2,1) - (Offsetx + size(En,1));
D2 = Offsety+1;

% D1 = (size(Temp2,1)-size(En,1))/2;
% D2 = (size(Temp2,2)-size(En,2))/2;
Temp2(D1:(D1+size(En,1)-1),D2:(D2+size(En,2)-1)) = En;

end
