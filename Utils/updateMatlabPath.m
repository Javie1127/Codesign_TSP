function updateMatlabPath()

    currentFolder = pwd;
    addpath([currentFolder, '\Config']);
    addpath([currentFolder, '\Example_1']);
    addpath([currentFolder, '\Example_2']);
    addpath([currentFolder, '\Example_3']);
    addpath([currentFolder, '\Example_5']);
    addpath([currentFolder, '\Example_4']);
    addpath([currentFolder, '\Stats']);
    addpath([currentFolder, '\BCD_AP Algorithm']);
    addpath([currentFolder, '\Utils']);
end