% Start with a folder and get a list of all subfolders.
% Finds and prints names of all PNG, JPG, and TIF images in 
% that folder and all of its subfolders.
clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;

%% set up environment
home_dir = 'The\folder\directory\where\all\the\code\is\stored';
cd(home_dir);
addpath(genpath('scripts/')) % add subsidiary code to the home directory
series = 1;
timepoint = 1;
start_path = 'The\folder\directory\where\the\images\are\stored';


% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end
% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames)

% Process all image files in those folders.
for k = 1 : numberOfFolders
	% Get this folder and print it out.
	thisFolder = listOfFolderNames{k};
	fprintf('Processing folder %s\n', thisFolder);
    
	% Add on TIF files.
	filePattern = sprintf('%s/*.tif', thisFolder);
	baseFileNames = dir(filePattern);
	
    numberOfImageFiles = length(baseFileNames);
	% Now we have a list of all files in this folder.
	
	if numberOfImageFiles >= 1
		% Go through all those image files.
		for f = 1 : numberOfImageFiles
			fullFileName = fullfile(thisFolder, baseFileNames(f).name);
			fprintf('     Processing image file %s\n', fullFileName);
            hyb_reader = bfGetReader(fullFileName);
            num_hyb_channels = hyb_reader.getSizeC;

            hyb_fitc_channel = 1;
            hyb_cy3_channel = 2;
            hyb_TxRed_channel = 3;
            hyb_cy5_channel = 4;
            hyb_dapi_channel = 5;

            xlen(1) = hyb_reader.getSizeX; ylen(1) = hyb_reader.getSizeY; zlen(1) = hyb_reader.getSizeZ;
            fitc_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');
            cy3_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');
            TxRed_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');
            cy5_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');
            dapi_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');

            for z=1:zlen(1)
                fitc_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_fitc_channel,timepoint);
                cy3_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_cy3_channel,timepoint);
                TxRed_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_TxRed_channel,timepoint);
                cy5_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_cy5_channel,timepoint);
                dapi_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_dapi_channel,timepoint);
            end

            TxRed_subtract_stack{1}  = zeros(ylen(1),xlen(1),zlen(1),'uint16');
            
            % Subtract the bleedthrough from the Texas Red channel. Both
            % the slope and the intercept were obtained by fitting a line
            % through the data clound generated through plotting the pixel
            % intensities in the Texas Red and Cy3 channels.
            for z=1:zlen(1)
                 TxRed_subtract_stack{1}(:,:,z) = TxRed_stacks{1}(:,:,z)-56.6695-0.1804*cy3_stacks{1}(:,:,z);
            end

            TxRed_subtract_stack{1}(TxRed_subtract_stack{1}<0) = 0;


            comb_stack = cat(4, fitc_stacks{1}, cy3_stacks{1}, TxRed_subtract_stack{1}, cy5_stacks{1}, dapi_stacks{1});


            bfsave(comb_stack, fullfile('Directory\to\savei\the\Bleadthrough_subtracted\images\', baseFileNames(f).name), 'dimensionOrder', 'XYZCT');
		end
	else
		fprintf('     Folder %s has no image files in it.\n', thisFolder);
	end
end

