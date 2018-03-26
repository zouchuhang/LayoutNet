function [depthorder reglabel regori rpo] = omap_depthorder(omap, vp, OMAP_FACTOR)
% function [reglabel regori rpo] = omap_depthorder(omap, vp, OMAP_FACTOR)
% % % % % % % depthorder{:,1} : smaller means closer to camera
% % % % % % % depthorder{:,2} : smaller means away from camera
% % % % % % % rpo(j,k) = 1 if region j is closer to camera than k (k is closer to vp than j)
% rpo{1}(j,k) = 1 if region j is above k
% rpo{2}(j,k) = 1 if region j is left of k
% rpo{3}(j,k) = 1 if region j is in front of k

[reglabel regori] = omap2region(omap);

omapsize = size(omap);
for scandir = 1:3
    scanlineendpts = vpscanlines(vp, scandir, omapsize, OMAP_FACTOR);
    rpo{scandir} = regionpartialorder(reglabel, regori, scanlineendpts);
    
%     figure;
%     for i = 1:size(scanlineendpts,1)
%         h = plot([scanlineendpts(i,1) scanlineendpts(i,3)],...
%             [scanlineendpts(i,2) scanlineendpts(i,4)],'x-');
%         axis([1 omapsize(2) 1 omapsize(1)]);
%         axis ij;
%         pause;
%         delete(h)
%     end
% 
%     depthorder{scandir,1} = getdepthorder(rpo{scandir});
% %     disp_omapdepth(reglabel, depthorder{scandir,1});
%     
%     depthorder{scandir,2} = getdepthorder(rpo{scandir}');
% %     disp_omapdepth(reglabel, depthorder{scandir,2});
depthorder = [];

end

%%
function depthorder = getdepthorder(rpo)
% smaller depthorder means closer to camera
% rpo = sparse(numreg, numreg)
% rpo(j,k) = 1 if region j is closer to camera than k (k is closer to vp than j)

% check if dag
if ~graphisdag(rpo)
    error('graph is not dag');
end

% figure; view(biograph(rpo));
% order = graphtopoorder(rpo);

% setup LP
numreg = size(rpo,1);
[row col] = find(rpo); % row(i) is closer to camera than col(i)
numrel = size(row,1);

% x(row(i)) + 1 < x(col(i))     % x is depth. smaller means closer to camera
% x(row(i)) - x(col(i)) < -1
% A(i,row(i)) = 1
% A(i,col(i)) = -1
% b(i) = -1
A = sparse(numrel,numreg);
b = -ones(numrel, 1);
for i = 1:numrel
    A(i,row(i)) = 1; 
    A(i,col(i)) = -1;
end

lb = ones(numreg, 1);
ub = lb + Inf;

Aeq = [];
beq = [];

f = ones(numreg, 1);

% run LP
depthorder = linprog(f,A,b,Aeq,beq,lb,ub);
depthorder = round(depthorder);



%%
function rpo = regionpartialorder(reglabel, regori, scanlineendpts)
% rpo = sparse(numreg, numreg)
% rpo(j,k) = 1 if region j is closer to camera than k (k is closer to vp than j)

numreg = length(regori);
rpo = sparse(numreg,numreg);

for i = 1:size(scanlineendpts,1)
    scanlinesamples = endpts2samples(scanlineendpts(i,:));
    
    regscan = reglabel(sub2ind(size(reglabel),scanlinesamples(:,2),scanlinesamples(:,1)));
    regscanuniq = setdiff(unique(regscan),0);
    if length(regscanuniq)>1
        for j = 1:length(regscanuniq)
            for k = j+1:length(regscanuniq)
                rj = regscanuniq(j);
                rk = regscanuniq(k);
                if regori(rj) ~= regori(rk) % if two regions have different orientation
                    minj = find(regscan==rj,1,'first');
                    maxj = find(regscan==rj,1,'last');
                    mink = find(regscan==rk,1,'first');
                    maxk = find(regscan==rk,1,'last');

                    if maxj<mink % j is closer to camera than k
                        rpo(rj,rk) = 1;
                    elseif maxk<minj % k is closer to camera than j
                        rpo(rk,rj) = 1;
                    end
                end
            end
        end
    end
end

%%
function samples = endpts2samples(endpts)
d = norm(endpts(1:2) - endpts(3:4));
samples = [linspace(endpts(1),endpts(3),d)' ...
    linspace(endpts(2),endpts(4),d)'];
samples = round(samples);

%%
function [reglabel regori] = omap2region(omap)
if length(size(omap)) ~= 3
    error('aoidsfoijfdsaoidsajfidsaj');
end

[reg1 num1] = bwlabel(omap(:,:,1));
[reg2 num2] = bwlabel(omap(:,:,2));
[reg3 num3] = bwlabel(omap(:,:,3));

% reglabel = reg1;
% reg2mask = find(reg2(:));
% reglabel(reg2mask) = reg2(reg2mask) + num1;
% reg3mask = find(reg3(:));
% reglabel(reg3mask) = reg3(reg3mask) + num1+num2;
regori = [repmat(1, num1, 1); ...
          repmat(2, num2, 1); ...
          repmat(3, num3, 1)];

regstat1 = regionprops(reg1, 'PixelIdxList', 'Area');
regstat2 = regionprops(reg2, 'PixelIdxList', 'Area');
regstat3 = regionprops(reg3, 'PixelIdxList', 'Area');
regstat = [regstat1; regstat2; regstat3];

MINAREA = 50;
largereg = [regstat.Area] > MINAREA;
regstat = regstat( largereg );

reglabel = zeros(size(omap,1), size(omap,2));
for i = 1:length(regstat)
    reglabel(regstat(i).PixelIdxList) = i;
end

regori = regori(largereg);

% 
% for i = 1:length(regstat)
%     if i<=length(regstat1)
%         regstat(i).dir = 1;
%     elseif i<=length(regstat1)+length(regstat2)
%         regstat(i).dir = 2;
%     else
%         regstat(i).dir = 3;
%     end
% end
% 
% figure; imshow(label2rgb(reg1))
% figure; imshow(label2rgb(reg2))
% figure; imshow(label2rgb(reg3))


%%
