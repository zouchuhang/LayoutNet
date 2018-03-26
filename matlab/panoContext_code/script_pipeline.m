add_path_init; % add folders
BEDROOMDATA; % define experiment parameters
config.bufname = './testfolder/'; % folder for saving intimediate result

aid = 22; % experiment data ID

compRoomHypot(aid); % get vanishing direction, normal map, and room layout
 
compObjectHypot(aid); % get object hypotheses

compDDSampling(aid); % data-driven sampling to get whole-room hypothese
 
compHolisticRanking( aid ); % holistic ranking with linear svm and visualize the top 1.
