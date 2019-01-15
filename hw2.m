main();

function main()
% Clear command window.
clc;
% Delete all variables.
clear;
% Define the color classes as an array.
classes = ["black" "blue" "green" "orange" "purple" "red" "white" "yellow"];
% Train on the training data.
classHistograms = training(classes);
% Test the algorithm on real data.
test(classes, classHistograms);
return; % from main()
end

function classHistograms = training(classes)
classHistograms = [];
% How many training files for each class.
classDataSize = 50;
% Choose a maximum value for histograms. Scale it from 0-255 to 0-127
histogramWidth = 128;
% Calculate the bin size for averaging.
spanSize = (256 / histogramWidth);
% For each color class.
for i = 1:length(classes)
    % Initialize the histograms array.
    classHistograms(:,:,i) = [zeros(histogramWidth,1) zeros(histogramWidth,1) zeros(histogramWidth,1)];
    % Check the directory for jpeg files.
    fileList = dir(strcat('customdataset/', classes(i),'/*.jpeg'));
    if (~isempty(fileList))
        % Initialize averaging histogram arrays.
        avgHistRed = zeros(histogramWidth, 1);
        avgHistGreen = zeros(histogramWidth, 1);
        avgHistBlue = zeros(histogramWidth, 1);
        corruptFiles = 0;
        % For each file
        for j = 1:classDataSize
            file = fileList(j);
            filePath = strcat(file.folder,'\', file.name);
            I = imread(filePath);
            % Scale each image to 350x350 px.
            I = imresize(I, [350 350]);
            channelSize = size(I,3);
            % Check if the image is RGB.
            if (channelSize == 3)
                chRed = I(:,:,1);
                chGreen = I(:,:,2);
                chBlue = I(:,:,3);
                % Calculate the histograms.
                [histRed, ~] = imhist(chRed);
                [histGreen, ~] = imhist(chGreen);
                [histBlue, ~] = imhist(chBlue);
                
                % Add histogram values.
                for k = 1:histogramWidth
                    for p = 1:spanSize
                        histogramIndex = ((k - 1) * spanSize) + p;
                        avgHistRed(k) = avgHistRed(k) + histRed(histogramIndex);
                        avgHistGreen(k) = avgHistGreen(k) + histGreen(histogramIndex);
                        avgHistBlue(k) = avgHistBlue(k) + histBlue(histogramIndex);
                    end
                    
                    %If bins are bigger than 1 average them.
                    if(spanSize > 1)
                        avgHistRed(k) = avgHistRed(k) / spanSize;
                        avgHistGreen(k) = avgHistGreen(k) / spanSize;
                        avgHistBlue(k) = avgHistBlue(k) / spanSize;
                    end
                end
                % #debug show class and file counter.
                % disp([classes(i), num2str(j),'.'])
            else
                % If image is now in RGB, increase corruptFiles counter.
                corruptFiles = corruptFiles + 1;
            end

        end 
        % Average the histograms.
        avgHistRed = avgHistRed / (classDataSize - corruptFiles);
        avgHistGreen = avgHistGreen / (classDataSize - corruptFiles);
        avgHistBlue = avgHistBlue / (classDataSize - corruptFiles);
        % Assign to the class histograms.
        classHistograms(:, 1, i) = avgHistRed;
        classHistograms(:, 2, i) = avgHistGreen;
        classHistograms(:, 3, i) = avgHistBlue;
    end
end

return;
end

function test(classes, classHistograms)

% Choose a maximum value for histograms. Scale it from 0-255 to 0-127
histogramWidth = 128;
% Calculate the bin size for averaging.
spanSize = (256 / histogramWidth);
% Check the directory for jpeg files.
fileList = dir(strcat('testdataset/*.jpeg'));
    if (~isempty(fileList))
        for fileindex = 1:length(fileList)
            file = fileList(fileindex);
            imageHistograms(:,:,1) = [zeros(histogramWidth,1) zeros(histogramWidth,1) zeros(histogramWidth,1)];
            I = imread(strcat(file.folder,'\', file.name));
            % Scale each image to 350x350 px.
            I = imresize(I, [350 350]);
            channelSize = size(I,3);
            % Check if the image is RGB.
            if (channelSize == 3)
                % Initialize averaging histogram arrays.
                avgHistRed = zeros(histogramWidth, 1);
                avgHistGreen = zeros(histogramWidth, 1);
                avgHistBlue = zeros(histogramWidth, 1);
                
                chRed = I(:,:,1);
                chGreen = I(:,:,2);
                chBlue = I(:,:,3);
                % Calculate the histograms.
                [histRed, ~] = imhist(chRed);
                [histGreen, ~] = imhist(chGreen);
                [histBlue, ~] = imhist(chBlue);
                % Add histogram values.
                for k = 1:histogramWidth
                    for p = 1:spanSize
                        histogramIndex = ((k - 1) * spanSize) + p;
                        avgHistRed(k) = avgHistRed(k) + histRed(histogramIndex);
                        avgHistGreen(k) = avgHistGreen(k) + histGreen(histogramIndex);
                        avgHistBlue(k) = avgHistBlue(k) + histBlue(histogramIndex);
                    end
                    %If bins are bigger than 1 average them.
                    if(spanSize > 1)
                        avgHistRed(k) = avgHistRed(k) / spanSize;
                        avgHistGreen(k) = avgHistGreen(k) / spanSize;
                        avgHistBlue(k) = avgHistBlue(k) / spanSize;
                    end
                end
                % Assign to the image histograms.
                imageHistograms(:, 1, 1) = avgHistRed;
                imageHistograms(:, 2, 1) = avgHistGreen;
                imageHistograms(:, 3, 1) = avgHistBlue;
            end
            % Initialize the distance array.
            distances = zeros(length(classes),1);
            for i = 1:length(classes)
                % Calculate the distance between classes and input image.
                distR = (imageHistograms(:, 1, 1) - classHistograms(:, 1, i)).^2;
                distG = (imageHistograms(:, 2, 1) - classHistograms(:, 2, i)).^2;
                distB = (imageHistograms(:, 3, 1) - classHistograms(:, 3, i)).^2;

                totalDistance = sqrt(sum(distR) + sum(distG) + sum(distB));
                distances(i,1) = totalDistance;
            end
            % Choose the minimum distance class.
            [~,index] = min(distances);
            disp([file.name, classes(index)]);
        end
    end
return;
end