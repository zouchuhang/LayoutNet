
% Add local code directories to Matlab path
% addpaths;


imdir='../Images_resized/'; % directory with original images

workspcdir='../tempworkspace/'; % directory to save intermediate results
if ~exist(workspcdir,'dir')
    mkdir(workspcdir);
end

resdir='../Results/'; % This is where we will save final results using this demo script.
if ~exist(resdir,'dir')
    mkdir(resdir);
end



% You can run it on a single image as follows
imagename='indoor_0268.jpg';
[ boxlayout,surface_labels ] = getspatiallayout(imdir,imagename,workspcdir);


