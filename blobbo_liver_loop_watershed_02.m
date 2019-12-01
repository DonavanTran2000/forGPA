%------------------------------------------------------------------------------------------------
% Demo to illustrate simple blob detection, measurement, and filtering.
% Requires the Image Processing Toolbox (IPT) because it demonstates some functions
% supplied by that toolbox, plus it uses the "coins" demo image supplied with that toolbox.
% If you have the IPT (you can check by typing ver on the command line), you should be able to
% run this demo code simply by copying and pasting this code into a new editor window,
% and then clicking the green "run" triangle on the toolbar.
% Running time = 7.5 seconds the first run and 2.5 seconds on subsequent runs.
% A similar Mathworks demo:
% http://www.mathworks.com/products/image/demos.html?file=/products/demos/shipping/images/ipexprops.html
% Code written and posted by ImageAnalyst, July 2009.
% Modified for MRI scans by C Rahul, February 2014.
%------------------------------------------------------------------------------------------------

% function BlobsDemo()
% echo on;
% Startup code.

tic; % Start timer.
clc; % Clear command window.
clear all; clc;% Get rid of variables from prior run of this m-file.
disp('Running BlobsDemo.m...'); % Message sent to command window.
workspace; % Show panel with all the variables.

% Change the current folder to the folder of this m-file.
% (The line of code below is from Brett Shoelson of The Mathworks.)
if(~isdeployed)
	cd(fileparts(which(mfilename)));
end



pcloud_out = [];

pcount = 0;

% %Load the 3D volume data
load('Liver01e_v3d_x1.mat');
[depth, width, height] = size(vox3d); %Extract the width, height, and depth


%Process the z-axis reconstruction, then combine the x- and y-'s for final
%results.

%Test the first 20 images of the z-axis, then we generalize it later.
for ii = 80:145
    
originalImage = vox3d(:, :, ii);

%Extract ROI (Region of Interest)
originalImage = originalImage(78:250, 172: 400);

% Image Smothing (Option 2)
%originalImage = adapthisteq(originalImage); %Use histogram equalization for contrast enhancement

% Image smoothing (Option 1) (acceptable)
% H = fspecial('average', 9);
% iprocessed = imfilter(originalImage, H, 'replicate');
% iprocessed = imcomplement(iprocessed);
% originalImage = iprocessed;




h = fspecial('sobel');
fd = double(originalImage); 
g = sqrt(imfilter(fd, h, 'replicate') .^2 + ...
         imfilter(fd, h, 'replicate') .^2);







%Plot all the images on the same frame
figure(1)
subplot(3, 3, 1);
imshow(originalImage); %Display the original image in gray scale (2-D), after filtered.

% Maximize the figure window.
%set(gcf, 'Position', get(0, 'ScreenSize'));
% Force it to display RIGHT NOW (otherwise it might not display until it's all done, unless you've stopped at a breakpoint.)
drawnow;
caption = sprintf('Original Ultrasound Doppler Image'); %For the first image
title(caption); 
axis square; % Make sure image is not artificially stretched because of screen's aspect ratio.


% %Extended minima transform of the image
im = imextendedmin(originalImage, 100);
fim = originalImage;
fim(im) = 175;

originalImage = fim;

% % Just for fun, let's get its histogram. 
% [pixelCount grayLevels] = imhist(originalImage);
% subplot(3, 3, 2);
% bar(pixelCount); title('Histogram of original image'); %Useful to know
% xlim([0 grayLevels(end)]); % Scale x axis manually.
subplot(3, 3, 2);
imshow(originalImage); %Display the original image in gray scale (2-D), after filtered.

% Maximize the figure window.
%set(gcf, 'Position', get(0, 'ScreenSize'));
% Force it to display RIGHT NOW (otherwise it might not display until it's all done, unless you've stopped at a breakpoint.)
drawnow;
caption = sprintf('Filtered Watershed Doppler Image'); %For the first image
title(caption); 
axis square; % Make sure image is not artificially stretched because of screen's aspect ratio.


% Filtering starts
sz = size(originalImage); %Extract the size (width and height) of the image
nr = sz(1); %width
nc = sz(2); %height

% Setting all pixels more than 180 and less than 60 to 0 to create a clear distinction
% between spine fluid and spine (3rd illustration image to put on the same
% frame). Maybe important~ 
for i=1:nr
    for j=1:nc
        %if originalImage(i,j) < 230 || originalImage(i,j) > 250
        if originalImage(i,j) < 100;
            originalImage(i,j) = 0;
        end
    end
end

%******************************IMPORTANT STEP******************************
% Threshold the image to get a binary image (only 0's and 1's) of class "logical."
% Option 1: use Ostuthresholding
% thresholdValue = graythresh(originalImage);
% binaryImage = imbinarize(originalImage, thresholdValue);

% Option 2: Set a defined threshold value
% thresholdValue = 110;
% %binaryImage = originalImage > thresholdValue; % Bright objects will be the chosen if you use >.
% binaryImage = originalImage < thresholdValue; % Dark objects will be the chosen if you use <.

% thresholdValue = graythresh(originalImage);
% binaryImage = imbinarize(originalImage, thresholdValue);

%thresholdValue = 140;
Lim = watershed(bwdist(im));
em = Lim == 0;

%The binary Image
binaryImage = imimposemin(g, im | em);

% Do a "hole fill" to get rid of any background pixels inside the blobs.
binaryImage = imbinarize(imfill(binaryImage, 'holes'));


%Skeletonize the image to convert the 2D pixel to only 1D pixel curve
%binaryImage = bwskel(binaryImage);


% Show the threshold as a vertical red bar on the histogram.
% hold on;
% maxYValue = ylim;
% hStemLines = stem(thresholdValue, maxYValue(2), 'r');
% %children = get(hStemLines, 'children');
% %set(children(2),'visible', 'off');
% 
% % Place a text label on the bar chart showing the threshold.
% annotationText = sprintf('Thresholded at %d gray levels', thresholdValue);
% 
% % For text(), the x and y need to be of the data class "double" so let's cast both to double.
% text(double(thresholdValue + 5), double(0.5 * maxYValue(2)), annotationText, 'FontSize', 10, 'Color', [0 .5 0]);
% text(double(thresholdValue - 70), double(0.94 * maxYValue(2)), 'Background', 'FontSize', 10, 'Color', [0 0 .5]);
% text(double(thresholdValue + 50), double(0.94 * maxYValue(2)), 'Foreground', 'FontSize', 10, 'Color', [0 0 .5]);

subplot(3, 3, 3); imagesc(originalImage); title('Binary Image, obtained after initial spine fluid segmentation'); axis square;
% Display the binary image.
subplot(3, 3, 4); imagesc(binaryImage); colormap(gray(256)); title('Binary Image, obtained by thresholding'); axis square;

labeledImage = bwlabel(binaryImage, 8);     % Label each blob so we can make measurements of it
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels

subplot(3, 3, 5);
imshow(labeledImage, []);
title('Labeled Image, from bwlabel()');
axis square;
subplot(3, 3, 6);
imagesc(coloredLabels);
axis square;
caption = sprintf('Pseudo colored labels, from label2rgb().\nBlobs are numbered from top to bottom, then from left to right.');
title(caption);

% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
blobMeasurements = regionprops(labeledImage, originalImage, 'all');
numberOfBlobs = size(blobMeasurements, 1);


%     disp("Watashi ga kitta")

    % bwboundaries() returns a cell array, where each cell contains the row/column coordinates for an object in the image.
    % Plot the borders of all the blobs on the original grayscale image using the coordinates returned by bwboundaries.
    subplot(3, 3, 7); 
    imagesc(originalImage);
    title('Outlines, from bwboundaries()'); axis square;
    hold on;
    [boundaries,L,N,A] = bwboundaries(binaryImage,8,'noholes');
    colors = ['b' 'g' 'r' 'c' 'm' 'y'];
    
    numberOfBoundaries = size(boundaries);

    poly_eq = zeros(numberOfBoundaries(1),3); % Setting the equations matrix

    for k = 1 : numberOfBoundaries

        boundary = boundaries{k};
        cidx = mod(k,length(colors))+1;
        plot(boundary(:,2), boundary(:,1),...
        colors(cidx),'LineWidth',2);
        %randomize text position for better visibility
        rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
        col = boundary(rndRow,2); row = boundary(rndRow,1);
        h = text(col+1, row-1, num2str(L(row,col)));
        set(h,'Color',colors(cidx),...
            'FontSize',6,'FontWeight','bold');
        % getting the polynomial equation for each blob
        num_of_points_in_blob = size(boundary);
        if(num_of_points_in_blob(1) >=4)
            eq = polyfit(boundary(:, 1),boundary(:, 2),2);
            poly_eq(k,:) = [eq(1,1) eq(1,2) eq(1,3)];
        end

    end

    %spy(A);
    hold on;


    fontSize = 6;	% Used to control size of "blob number" labels put atop the image.
    labelShiftX = -7;	% Used to align the labels in the centers of the coins.
    blobECD = zeros(1, numberOfBlobs);
    % Print header line in the command window.
    fprintf(1,'Blob #      Mean Intensity  Area   Perimeter    Centroid       Diameter\n');
    % Loop over all blobs printing their measurements to the command window.
    for k = 1 : numberOfBlobs           % Loop through all blobs.
        % Find the mean of each blob.  (R2008a has a better way where you can pass the original image
        % directly into regionprops.  The way below works for all versions including earlier versions.)
        thisBlobsPixels = blobMeasurements(k).PixelIdxList;  % Get list of pixels in current blob.
        meanGL = mean(originalImage(thisBlobsPixels)); % Find mean intensity (in original image!)
        meanGL2008a = blobMeasurements(k).MeanIntensity; % Mean again, but only for version >= R2008a
        blobArea = blobMeasurements(k).Area;		% Get area.
        blobPerimeter = blobMeasurements(k).Perimeter;		% Get perimeter.
        blobCentroid = blobMeasurements(k).Centroid;		% Get centroid.
        blobECD(k) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
        fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobECD(k));
        % Put the "blob number" labels on the "boundaries" grayscale image.
        text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', fontSize, 'FontWeight', 'Bold');
    end

    % Put the labels on the rgb labeled image also.
    subplot(3, 3, 5);
    for k = 1 : numberOfBlobs           % Loop through all blobs.
        blobCentroid = blobMeasurements(k).Centroid;		% Get centroid.
        text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', fontSize, 'FontWeight', 'Bold');
    end

    % FILTER
    % Find only those blobs
    % with a diameter between 150 and 220 and an area more than 110.
    allBlobDiameters = [blobMeasurements.EquivDiameter];
    allBlobAreas = [blobMeasurements.Area];
    % Get a list of the blobs that meet our criteria and we need to keep.
    allowableDiameterIndexes = (allBlobDiameters > 10) & (allBlobDiameters < 120);
    allowableAreaIndexes = (allBlobAreas > 200) & (allBlobAreas < 3000) ; % Take the small objects.
    keeperIndexes = find(allowableAreaIndexes & allowableDiameterIndexes);
    % Extract only those blobs that meet our criteria, and
    % eliminate those blobs that don't meet our criteria.
    % Note how we use ismember() to do this.
    keeperBlobsImage = ismember(labeledImage, keeperIndexes);
    % Re-label with only the keeper blobs kept.
    labeledDimeImage = bwlabel(keeperBlobsImage, 8);     % Label each blob so we can make measurements of it
    % Now we're done.  We have a labeled image of blobs that meet our specified criteria.
    subplot(3, 3, 8);
    imshow(labeledDimeImage, []);
    axis square;
    title('Blobs after area-filtering');



    kI = keeperIndexes(1)
    bb = boundaries{kI};
    nb = length(bb);
    z = ii.*(ones(nb,1));



    % Now use the keeper blob as a mask on the original image.
    % This will let us display the original image in the regions of the keeper blobs.
    maskedImageDime = originalImage; % Simply a copy at first.
    maskedImageDime(~keeperBlobsImage) = 0;  % Set all non-keeper pixels to zero.
    subplot(3, 3, 9);
    imshow(maskedImageDime);
    axis square;
    title('Only the 3 brightest dimes from the original image');

    % Now let's get the nickels (the larger coin type)
    keeperIndexes = find(allBlobAreas > 300);  % Take the larger objects.
    % Note how we use ismember to select the blobs that meet our criteria.
    nickelBinaryImage = ismember(labeledImage, keeperIndexes);
    maskedImageNickel = originalImage; % Simply a copy at first.
    maskedImageNickel(~nickelBinaryImage) = 0;  % Set all non-nickel pixels to zero.
    %subplot(3, 3, 9);
    %imshow(maskedImageNickel, []);
    %axis square;
    %title('Only the nickels from the original image');



    figure(2), hold on
    plot3(bb(:,1), bb(:,2), z, 'or')



    pcloud_temp = [bb z];

    pcloud_out = [pcloud_out; pcloud_temp];


    % Save each contour as a polygon variable for later analysis
    if nb > 0
        eval( ['poly' num2str(ii) ' = pcloud_temp;'] )
        pcount = pcount + 1;
        pnum(pcount) = ii;     % Keep track of planes with successful outlining
    end
    
    figure(2)
    axis equal
    grid on


    % Write Point Cloud
    bname = 'IVC'

    pfile = [bname '_pts.asc']

    fid0 = fopen(pfile, 'w');

    fprintf( fid0, '%3.6f %3.6f %3.6f\n', pcloud_out' );

    fclose(fid0);


    % Save Polygons
    save IVC_Poly_out pnum poly*


    elapsedTime = toc;
end




% Alert user that the demo is done and give them the option to save an image.
% message = sprintf('Finished running BlobsDemo.m.\n\nElapsed time = %.2f seconds.', elapsedTime);
% message = sprintf('%s\n\nCheck out the figure window for the images.\nCheck out the command window for the numerical results.', message);
% message = sprintf('%s\n\nDo you want to save the pseudo-colored image?', message);
% reply = questdlg(message, 'Save image?', 'Yes', 'No', 'No');
% % Note: reply will = '' for Upper right X, 'Yes' for Yes, and 'No' for No.
% if strcmpi(reply, 'Yes')
% 	% Ask user for a filename.
% 	FilterSpec = {'*.tif', 'TIFF images (*.tif)'; '*.*', 'All Files (*.*)'};
% 	DialogTitle = 'Save image file name';
% 	% Get the default filename.  Make sure it's in the folder where this m-file lives.
% 	% (If they run this file but the cd is another folder then pwd will show that folder, not this one.
% 	thisFile = mfilename('fullpath');
% 	[thisFolder, baseFileName, ext, version] = fileparts(thisFile);
% 	DefaultName = sprintf('%s/%s.tif', thisFolder, baseFileName);
% 	[fileName, specifiedFolder] = uiputfile(FilterSpec, DialogTitle, DefaultName);
% 	% Parse what they actually specified.
% 	[folder, baseFileName, ext, version] = fileparts(fileName);
% 	% Create the full filename, making sure it has a tif filename.
% 	fullImageFileName = fullfile(specifiedFolder, [baseFileName '.tif']);
% 	% Save the labeled image as a tif image.
% 	imwrite(uint8(coloredLabels), fullImageFileName);
% 	% Just for fun, read image back into the imtool utility to demonstrate that tool.
% 	tifimage = imread(fullImageFileName);
% 	imtool(tifimage, []);
% end
%
% message = sprintf('Would you like to crop out each coin to individual images?');
% reply = questdlg(message, 'Extract Individual Images?', 'Yes', 'No', 'Yes');
% % Note: reply will = '' for Upper right X, 'Yes' for Yes, and 'No' for No.
% if strcmpi(reply, 'Yes')
% 	figure;
% 	% Maximize the figure window.
% 	set(gcf, 'Position', get(0, 'ScreenSize'));
% 	for k = 1 : numberOfBlobs           % Loop through all blobs.
% 		% Find the bounding box of each blob.
% 		thisBlobsBoundingBox = blobMeasurements(k).BoundingBox;  % Get list of pixels in current blob.
% 		% Extract out this coin into it's own image.
% 		subImage = imcrop(originalImage, thisBlobsBoundingBox);
% 		% Display the image with informative caption.
% 		subplot(3, 4, k);
% 		imshow(subImage);
% 		caption = sprintf('Coin #%d\nDiameter = %.1f pixels\nArea = %d pixels', k, blobECD(k), blobMeasurements(k).Area);
% 		title(caption, 'FontSize', 14);
% 	end


%end