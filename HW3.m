
main();

function main()
% Clear command window.
clc;
% Delete all variables.
clear;
contrastAdjuster = 0.05;
threshold = 0.25;
kernelSize = 5;
neighbourhoodSize = 19;
countBirds('bird images/bird 1.jpg', contrastAdjuster, threshold, kernelSize, neighbourhoodSize);
countBirds('bird images/bird 2.jpg', contrastAdjuster, threshold, kernelSize, neighbourhoodSize);
countBirds('bird images/bird 3.bmp', contrastAdjuster, threshold, kernelSize, neighbourhoodSize);
return; % from main()
end

function countBirds(filePath, contrastAdjuster, threshold, kernelSize, neighbourhoodSize)
% Read the image on given filePath
img = readImages(filePath);
imwrite(img, strcat(filePath, '_1.jpg'));
% Convert image to grayscale
grayscale = convertToGrayscale(img);
figure, imshow(grayscale);
imwrite(grayscale, strcat(filePath, '_2.jpg'));
% Invert the pixel values of image
inv = invertImage(grayscale);
figure, imshow(inv);
imwrite(inv, strcat(filePath, '_3.jpg'));
% Adjust the contrast so more details are present
adjusted = adjustContrast(inv, contrastAdjuster);
figure, imshow(adjusted);
imwrite(adjusted, strcat(filePath, '_4.jpg'));
% Apply morphological operations erosion and dilation.
erosionImg = erosion(adjusted, neighbourhoodSize);
figure, imshow(erosionImg);
imwrite(erosionImg, strcat(filePath, '_5.jpg'));
background = dilation(erosionImg, neighbourhoodSize);
figure, imshow(background);
imwrite(background, strcat(filePath, '_6.jpg'));
% Get objects by removing the background
objects = adjusted - background;
figure, imshow(objects);
imwrite(objects, strcat(filePath, '_7.jpg'));
% Get binary image by a threshold
binaryImage = binarize(objects, threshold);
figure, imshow(binaryImage);
imwrite(binaryImage, strcat(filePath, '_8.jpg'));
% Create a new averaging kernel
kernel(1:kernelSize, 1:kernelSize) = 1/(kernelSize * kernelSize);
% Apply the averaging kernel
averagedBinaryImage = convolute(binaryImage, kernel);
% Apply Connected Component Analysis to find closed connected components
result = connectedComponentAnalysis(averagedBinaryImage);
% Get maximum values from result as object count
numberOfBirds = getMaxValue(result);
% Display the image
figure, imshow(averagedBinaryImage);
imwrite(averagedBinaryImage, strcat(filePath, '_9.jpg'));
% Display the values
disp(filePath);
disp(numberOfBirds);
return;
end

function img = readImages(filepath)
img = imread(filepath);
return;
end

function ret = convertToGrayscale(image)
channelsize = size(image,3);
imageSizeX = size(image,1);
imageSizeY = size(image,2);
ret = zeros(imageSizeX, imageSizeY);
for i = 1:imageSizeX
    for j = 1:imageSizeY
        pixel = [0, 0, 0];
        for k = 1:channelsize
            pixel(k) = image(i,j,k);
        end
        if (channelsize == 1)
            pixel(2) = pixel(1);
            pixel(3) = pixel(1);
        end
        avg = luminosityGrayscale(pixel(1), pixel(2), pixel(3));
        gs = double(avg) / double(255);
        ret(i,j) = gs;
    end
end
end

function val = luminosityGrayscale(r, g, b)
val = (0.21 * double(r)) + (0.72 * double(g)) + (0.07 * double(b));
end

% For given contrastAdjuster value, if pixel value is smaller set 0, if
% pixel value is greater than 1 - contrastAdjuster value set 1, else set
% normal pixel value.
function img = adjustContrast(image, contrastAdjuster)
imageSizeX = size(image,1);
imageSizeY = size(image,2);
img = zeros(imageSizeX, imageSizeY);
for i = 1:imageSizeX
    for j = 1:imageSizeY
        if (image(i,j) < contrastAdjuster)
            img(i,j) = 0;
        elseif (image(i,j) > (1 - contrastAdjuster))
            img(i,j) = 1;
        else
            img(i,j) = image(i,j);
        end
    end
end
return;
end

% Binarize the image for given threshold.
function img = binarize(image, threshold)
imageSizeX = size(image,1);
imageSizeY = size(image,2);
img = zeros(imageSizeX, imageSizeY);
for i = 1:imageSizeX
    for j = 1:imageSizeY
        if (image(i,j) < threshold)
            img(i,j) = 0;
        else
            img(i,j) = 1;
        end
    end
end
return;
end

% Invert the pixels
function img = invertImage(image)
imageSizeX = size(image,1);
imageSizeY = size(image,2);
img = zeros(imageSizeX, imageSizeY);
for i = 1:imageSizeX
    for j = 1:imageSizeY
        img(i,j) = 1 - image(i,j);
    end
end
return;
end

function ret = convolute(image, kernel)

kernelSizeX = size(kernel, 1);
kernelSizeY = size(kernel, 2);

imageSizeX = size(image,1);
imageSizeY = size(image,2);

flippedKernel = zeros(kernelSizeX, kernelSizeY);

for i = 1:kernelSizeX
    for j = 1:kernelSizeY
        flippedKernel(i, j) = kernel(kernelSizeX - i + 1, kernelSizeY - j + 1);
    end
end

imageX = 1;
imageY = 1;

ret = zeros(imageSizeX - (kernelSizeX - 1), imageSizeY - (kernelSizeX - 1));
for i = (1 + ((kernelSizeX - 1) / 2)):(imageSizeX - ((kernelSizeX - 1) / 2))
    for j = (1 + ((kernelSizeY - 1) / 2)):(imageSizeY - ((kernelSizeY - 1) / 2))
        for a = 1:kernelSizeX
            imageIndexerX = -((kernelSizeX - 1) / 2) + a - 1;
            for b = 1:kernelSizeY
                imageIndexerY = -((kernelSizeY - 1) / 2) + b - 1;
                ret(imageX,imageY) = ret(imageX,imageY) + image(i + imageIndexerX, j + imageIndexerY) * flippedKernel(a,b);
            end
        end
        imageY = imageY + 1;
    end
    imageX = imageX + 1;
    imageY = 1;
end
end

% Apply erosion for morphology, find the minimum for given neigbourhood
% size
function ret = erosion(image, neighbourhoodSize)
imageSizeX = size(image,1);
imageSizeY = size(image,2);
ret = zeros(imageSizeX,imageSizeY);

for i = 1:imageSizeX
    for j = 1:imageSizeY
        min = image(i,j);
        for x = -((neighbourhoodSize - 1) / 2) : ((neighbourhoodSize - 1) / 2)
            indexX = i + x;
            for y = -((neighbourhoodSize - 1) / 2) : ((neighbourhoodSize - 1) / 2)
                indexY = j + y;
                if ((indexX > 0) && (indexX < (imageSizeX + 1)) && (indexY > 0) && (indexY < (imageSizeY + 1)) && (image(indexX,indexY) < min))
                    min = image(indexX,indexY);
                end
            end
        end
        ret(i,j) = min;
    end
end
return;
end

% Apply dilation for morphology, find the maximum for given neigbourhood
% size
function ret = dilation(image, neighbourhoodSize)
imageSizeX = size(image,1);
imageSizeY = size(image,2);
ret = zeros(imageSizeX,imageSizeY);

for i = 1:imageSizeX
    for j = 1:imageSizeY
        max = image(i,j);
        for x = -((neighbourhoodSize - 1) / 2) : ((neighbourhoodSize - 1) / 2)
            indexX = i + x;
            for y = -((neighbourhoodSize - 1) / 2) : ((neighbourhoodSize - 1) / 2)
                indexY = j + y;
                if (indexX > 0 && indexX < imageSizeX + 1 && indexY > 0 && indexY < imageSizeY + 1 && image(indexX,indexY) > max)
                    max = image(indexX,indexY);
                end
            end
        end
        ret(i,j) = max;
    end
end
return;
end

% Find connected components and number them. Uses DFS
function returnImage = connectedComponentAnalysis(image)
visited = false(size(image));
imageSizeX = size(image,1);
imageSizeY = size(image,2);
returnImage = zeros(imageSizeX,imageSizeY);
counter = 1;

for row = 1 : imageSizeX
    for col = 1 : imageSizeY
        if image(row,col) == 0
            visited(row,col) = true;
        elseif visited(row,col)
            continue;
        else
            searchStack = [row col];
            while ~isempty(searchStack)
                loc = searchStack(1,:);
                searchStack(1,:) = [];
                if visited(loc(1),loc(2))
                    continue;
                end
                visited(loc(1),loc(2)) = true;
                returnImage(loc(1),loc(2)) = counter;
                [locs_y, locs_x] = meshgrid(loc(2)-1:loc(2)+1, loc(1)-1:loc(1)+1);
                locs_y = locs_y(:);
                locs_x = locs_x(:);
                out_of_bounds = locs_x < 1 | locs_x > imageSizeX | locs_y < 1 | locs_y > imageSizeY;
                
                locs_y(out_of_bounds) = [];
                locs_x(out_of_bounds) = [];
                is_visited = visited(sub2ind([imageSizeX imageSizeY], locs_x, locs_y));
                
                locs_y(is_visited) = [];
                locs_x(is_visited) = [];
                
                is_1 = image(sub2ind([imageSizeX imageSizeY], locs_x, locs_y));
                locs_y(~is_1) = [];
                locs_x(~is_1) = [];
                searchStack = [searchStack; [locs_x locs_y]];
            end
            counter = counter + 1;
        end
    end
end
return;
end

% Find max value in a matrix
function maxValue = getMaxValue(matrix)
sizeX = size(matrix,1);
sizeY = size(matrix,2);
maxValue = matrix(1,1);
for i=1:sizeX
    for j=1:sizeY
        if matrix(i,j) > maxValue
            maxValue = matrix(i,j);
        end
    end
end

return;
end