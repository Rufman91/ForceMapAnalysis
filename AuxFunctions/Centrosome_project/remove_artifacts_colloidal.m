
%% BEGIN%%

NumberOfTheTipYouWantToUse = 2;
[Channel,Index] = E.CantileverTips{NumberOfTheTipYouWantToUse}.get_unprocessed_height_channel('Height (measured) (Trace)');
Based = AFMImage.subtract_line_fit_hist(Channel.Image,0.4);
Channel.Image = Based;
figure
imshow(Channel.Image,[])
DrawObject = drawfreehand;
% draw the mask in a way that EXCLUDES all the artifacts.
% So everything that should stay is covered.
% It's easier to draw this way. We invert the mask afterwards
Mask = ~DrawObject.createMask;
imshowpair(Channel.Image,Mask,'montage')
drawnow
pause(3)
CleanImage = Channel.Image;
CleanImage(Mask) = min(CleanImage,[],'all');
imshowpair(Channel.Image,CleanImage,'montage')


% Now we assign our image back to the place where the deconvolution
% algorithm will take it from. Be aware that this will overwrite a bit of original information.
% You can recover it as long as you have the 'Channel' variable in your workspace but otherwise its gone
%
% !!!CAUTION!!!
%
% Execute this step only if you are happy withcl the results from the clean up
E.CantileverTips{NumberOfTheTipYouWantToUse}.Channel(Index).Image = CleanImage;
% From here you can proceed the analysis as usual

%% END%%
