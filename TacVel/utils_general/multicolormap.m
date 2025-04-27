function finalfilename = multicolormap(drctry)
tiffz = getimages(drctry);

for i=1:length(tiffz)
[varnamez{i},imagedata{i}] = importfile(tiffz{i});
end
nimagez = length(imagedata);
% for ii=1:nimagez
%     figure(ii+1000), clf
%     image()

%%

%[542 655 768 881 994]
xx = 169:1743;
panelheight = 113;
y_offset = 542;
yy = zeros(nimagez,panelheight+1);
for ii=1:nimagez    
yy(ii,:) = y_offset+ ((ii-1)*panelheight:(ii)*panelheight);
end
yy = flipud(yy);

% forcetarget 1 : blue
% forcetarget 2 : red
% forcetarget 3 : green
% forcetarget 4 : magenta
%%
finalimage = imagedata{1,1}{1,1};

for j=2:nimagez
%     figure(254), clf
%     subplot(1,3,1)
%     imshow(finalimage)
%     drawnow
    altcolormapimage = imagedata{1,j}{1,1};
    finalimage(yy(j,:),xx,:)=altcolormapimage(yy(j,:),xx,:);
%     subplot(1,3,2)    
%     imshow(altcolormapimage)
%     subplot(1,3,3)    
%     imshow(difference_image)
%     drawnow
%     pause(1)
end
    
%
% figure(259), clf
% imagesc(finalimage)

finalfilename = tiffz{1};
finalfilename = [finalfilename(1:end-27),'FINAL.tiff'];
imwrite(finalimage,finalfilename)
winopen(finalfilename)
