
function A = tiff23(Image,Savename, newsize)
medflt=[2 2];
%%
Image1=rescale(Image);
img1=medfilt2(Image1(:,:),medflt);
   if newsize == 1
   else
img1 = imresize(img1,newsize,'nearest'); 
   end
imwrite(img1,Savename);

A = img1;

end