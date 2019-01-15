I = imread('images/Lenna.png');
channelSize = size(I,3);

hsize = 11;
sigma = 7;
gaussianKernel = fspecial('gaussian',hsize,sigma);

% Convert to greyscale
workingCopy = convertToGrayscale(I, channelSize);
% Noise reduction with gaussianfilter
workingCopy = convolute(workingCopy, gaussianKernel);
% Sobel Filters
Sx = [-1 0 1
      -2 0 2
      -1 0 1];
Sy = [1  2  1
      0  0  0
     -1 -2 -1];
 
Gx = convolute(workingCopy, (1/8) * Sx);
Gy = convolute(workingCopy, (1/8) * Sy);
slope = calculateSlope(Gx, Gy); %direction
magnitude = calculateMagnitude(Gx, Gy);%sum
supressedImage = nonmaximumSupression(slope, magnitude);
workingCopy = hysteresis(slope, magnitude, supressedImage);
imwrite(workingCopy,'test.jpg');
imshow(workingCopy);

function ret = hysteresis(slope, magnitude, workingCopy)
loThres = 100;
ret = zeros(size(slope,1), size(slope,2));
    isChanged = true;
    it = 0;
    while (isChanged)
        isChanged = false;
    end
    it = it + 1;
    for i = 1:size(slope, 1)
        for j = 1:size(slope, 2)
            if ((i < 3) ||i > size(slope, 1) - 2 || (j < 2) || (j > size(slope, 2)))
                currentDirection = magnitude(i, j);
                if (workingCopy(i,j) == 255)
                    workingCopy(i,j) = 64;
                    
                    if(currentDirection > 112.5 && currentDirection <= 157.5)
                        if(i > 1 && j > 1)
                            if ((loThres <= magnitude(i - 1, j - 1)) && (workingCopy(i - 1, j - 1) ~= 64) && (slope(i - 1, j - 1) > 112.5) && (slope(i - 1, j - 1) <= 157.5) && (magnitude(i-1,j-1) > magnitude(i, j-2)) && (magnitude(i-1,j-1) > magnitude(i-2, j)))
                                ret(i-1, j-1) = 255;
                                isChanged = true;
                            end
                        end
                        if ((i < size(slope, 1) - 1 ) && (j < size(slope, 2) - 1))
                            if ((loThres <= magnitude(i + 1, j + 1)) && (workingCopy(i + 1, j + 1) ~= 64) && (slope(i + 1, j + 1) > 112.5) && (slope(i + 1, j + 1) <= 157.5) && (magnitude(i-1,j-1) > magnitude(i, j+2)) && (magnitude(i-1,j-1) > magnitude(i+2, j)))
                                ret(i+1, j+1) = 255;
                                isChanged = true;
                            end
                        end
                    elseif(currentDirection > 67.5 && currentDirection <= 112.5)
                        if(i > 1)
                            if ((loThres <= magnitude(i - 1, j)) && (workingCopy(i - 1, j) ~= 64) && (slope(i - 1, j) > 67.5) && (slope(i - 1, j) <= 112.5) && (magnitude(i-1,j) > magnitude(i-1, j-1)) && (magnitude(i-1,j) > magnitude(i-1, j+1)))
                                ret(i-1, j-1) = 255;
                                isChanged = true;
                            end
                        end
                        if ((i < size(slope, 1) - 1 ))
                            if ((loThres <= magnitude(i + 1, j)) && (workingCopy(i + 1, j) ~= 64) && (slope(i + 1, j) > 67.5) && (slope(i + 1, j) <= 112.5) && (magnitude(i+1,j) > magnitude(i+1, j-1)) && (magnitude(i+1,j) > magnitude(i+1, j+1)))
                                ret(i+1, j+1) = 255;
                                isChanged = true;
                            end
                        end
                    elseif(currentDirection > 22.5 && currentDirection <= 67.5)
                        if((i < size(slope, 1) - 1) && j > 1)
                            if ((loThres <= magnitude(i + 1, j - 1)) && (workingCopy(i + 1, j - 1) ~= 64) && (slope(i + 1, j - 1) > 22.5) && (slope(i + 1, j - 1) <= 67.5) && (magnitude(i-1,j+1) > magnitude(i-2, j)) && (magnitude(i-1,j+1) > magnitude(i, j+2)))
                                ret(i-1, j-1) = 255;
                                isChanged = true;
                            end
                        end
                        if (i > 1 && (j < size(slope, 2) - 1))
                            if ((loThres <= magnitude(i- 1, j + 1)) && (workingCopy(i - 1, j + 1) ~= 64) && (slope(i - 1, j + 1) > 225.5) && (slope(i - 1, j + 1) <= 67.5) && (magnitude(i-1,j-1) > magnitude(i, j+2)) && (magnitude(i-1,j-1) > magnitude(i+2, j)))
                                ret(i+1, j+1) = 255;
                                isChanged = true;
                            end
                        end
                    end
                end
            end
        end
    end
end

function ret = nonmaximumSupression(slope, magnitude)
ret = zeros(size(slope,1), size(slope,2));
thresholdUp = 0.6;
  for i = 1:size(slope, 1)
        for j = 1:size(slope, 2)
            direction = ((atan(slope(i,j)) * 180) / pi);
            while (direction < 0)
                direction = direction + 180;
            end
            slope(i,j) = direction;
            if (magnitude(i,j) < thresholdUp)
                isEdge = true;
                if (direction > 112.5 && direction <= 157.5)
                    if  ((j > 1) && (i < size(slope, 1)) && magnitude(i,j) <= magnitude(i + 1, j - 1))
                        isEdge = false;
                    end
                    if  ((i > 1) && (j < size(slope, 2)) && magnitude(i,j) <= magnitude(i - 1, j + 1))
                        isEdge = false;
                    end
                elseif (direction > 67.5 && direction <= 112.5)
                    if  ((j > 1) && magnitude(i,j) < magnitude(i, j - 1))
                        isEdge = false;
                    end
                    if  ((j < size(slope, 2)) && magnitude(i,j) <= magnitude(i, j + 1))
                        isEdge = false;
                    end
                elseif  (direction > 22.5 && direction <= 67.5)
                    if  ((j > 1) && (i > 1) && magnitude(i,j) <= magnitude(i - 1, j - 1))
                        isEdge = false;
                    end
                    if  (( i< size(slope, 1)) && (j < size(slope, 2)) && magnitude(i,j) <= magnitude(i + 1, j + 1))
                        isEdge = false;
                    end
                else
                    if  ((i > 1) && magnitude(i,j) <= magnitude(i - 1, j))
                        isEdge = false;
                    end
                    if  ((i < size(slope, 2)) && magnitude(i,j) <= magnitude(i + 1, j))
                        isEdge = false;
                    end
                end

                if(isEdge)
                    ret(i, j) = 255;
                end
            end
        end
  end
end

function ret = calculateSlope(Gx, Gy)
ret = zeros(size(Gx,1), size(Gx,2));
    for i = 1:size(Gx, 1)
        for j = 1:size(Gx, 2)
        ret(i, j) = Gy(i, j) / Gx(i, j);
        end 
    end
end

function ret = calculateMagnitude(Gx, Gy)
ret = zeros(size(Gx,1), size(Gx,2));
    for i = 1:size(Gx, 1)
        for j = 1:size(Gx, 2)
        ret(i, j) = sqrt(Gx(i, j) * Gx(i, j)) + (Gy(i, j) * Gy(i, j));
        end 
    end
end

function ret = convolute(image, kernel)

kernelSizeX = size(kernel, 1);
kernelSizeY = size(kernel, 2);

imageSizeX = size(image,1);
imageSizeY = size(image,2);

flippedKernel = zeros(size(kernel, 1), size(kernel, 2));

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

function ret = convertToGrayscale(image, channelsize)
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

function val = averageGrayscale(r, g, b) 
    val = (r + g + b) / 3;
end 