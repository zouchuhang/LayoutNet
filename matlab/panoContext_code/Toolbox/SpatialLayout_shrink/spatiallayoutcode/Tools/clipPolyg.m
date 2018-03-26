function [Polyg]=clipPolyg(Polyg,h,w)

for f=1:numel(Polyg)
    
    if numel(Polyg{f})> 0
    
        
        Polyg{f}(:,1)=round(Polyg{f}(:,1));
        Polyg{f}(:,2)=round(Polyg{f}(:,2));

        
         tind=find(Polyg{f}(:,1) == w+1 | Polyg{f}(:,1)==w+2);
         Polyg{f}(tind,1)=w;
         
         tind=find(Polyg{f}(:,1) == 0);
         Polyg{f}(tind,1)=1;
         
         tind=find(Polyg{f}(:,2) == h+1 | Polyg{f}(:,2)==h+2);
         Polyg{f}(tind,2)=h;
         
         tind=find(Polyg{f}(:,2) == 0);
         Polyg{f}(tind,2)=1;
         
    end
end


return;