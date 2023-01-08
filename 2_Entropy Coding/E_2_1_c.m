clc;clear;close all;

imageLena     = double(imread('lena.tif'));
imageSail     = double(imread('sail.tif'));
imageSmandril = double(imread('smandril.tif'));

pmfLena       = stats_marg(imageLena,0:255);
pmfSail       = stats_marg(imageSail,0:255);
pmfSmandril   = stats_marg(imageSmandril,0:255);

pixel_Lena=numel(imageLena);
pixel_Sail=numel(imageSail);
pixel_Smandril=numel(imageSmandril);

mergedPMF     = (pmfLena*pixel_Lena+pmfSail*pixel_Sail+pmfSmandril*pixel_Smandril)/(pixel_Lena+pixel_Sail+pixel_Smandril);


minCodeLengthLena     = min_code_length(mergedPMF,pmfLena);
minCodeLengthSail     = min_code_length(mergedPMF,pmfSail);
minCodeLengthSmandril = min_code_length(mergedPMF,pmfSmandril);

fprintf('--------------Using merged code table--------------\n');
fprintf('lena.tif      H = %.2f bit/pixel\n', minCodeLengthLena);
fprintf('sail.tif      H = %.2f bit/pixel\n', minCodeLengthSail);
fprintf('smandril.tif  H = %.2f bit/pixel\n', minCodeLengthSmandril);

%% Put all sub-functions which are called in your script here.
function pmf = stats_marg(image, range)
pmf=hist(image(:),range); % a vector
pmf=pmf/sum(pmf);
end

function H = min_code_length(pmf_table, pmf_image)
pmf_table_wo0=pmf_table(pmf_table>0);
H=-sum(pmf_image.*log2(pmf_table_wo0));

end