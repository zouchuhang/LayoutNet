function ShowGTPolyg(img,gtPolyg,fignum)
%pfc={'r','g','b','k','w'};
 pfc={'r','r','r','r','r'};
 figure(fignum);
imshow(img,[]);hold on;



for f=1: numel(gtPolyg)
    if numel(gtPolyg{f})>0
      plot([gtPolyg{f}(:,1);gtPolyg{f}(1,1)],[gtPolyg{f}(:,2);gtPolyg{f}(1,2)],'LineWidth',4,...
            'Color',pfc{f});
    end
end

hold off;
 
   