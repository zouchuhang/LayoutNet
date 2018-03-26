function [th1 th2]=gettheta_ind(xs,ys,vp,h,w,quantsiz)

if ~exist('quantsiz','var')
    quantsiz=500;
end

imcor=[0 h+1;0 0;w+1 0;w+1 h+1];
inim = vp(:,1)>=1 & vp(:,1)<=w & vp(:,2)>=1 & vp(:,2)<=h;
inimpts = xs>=1 & xs<=w & ys>=1 & ys<=h;

for vno=1:2
    if ~inim(vno)
        ll1=[vp(vno,1) imcor(1,1) vp(vno,2) imcor(1,2)];
        ll2=[vp(vno,1) imcor(2,1) vp(vno,2) imcor(2,2)];
        ll3=[vp(vno,1) imcor(3,1) vp(vno,2) imcor(3,2)];
        ll4=[vp(vno,1) imcor(4,1) vp(vno,2) imcor(4,2)];
        lla=[ll1;ll1;ll1;ll2;ll2;ll3];
        llb=[ll2;ll3;ll4;ll3;ll4;ll4];
        
        veca = lla(:,[2,4])-lla(:,[1,3]);
        vecb = llb(:,[2,4])-llb(:,[1,3]);
        norma=(sum(veca.*veca,2)).^.5;
        normb=(sum(vecb.*vecb,2)).^.5;
        %         theta = acosd(dot(veca,vecb,2)./norma./normb);
        theta = acosd(dot(veca./repmat(norma,[1 2]),vecb./repmat(normb,[1 2]),2));
        [vv ii]=max(theta);
        
        veca = veca(ii,:);
        vecb = vecb(ii,:);
        norma = norma(ii);
        normb = normb(ii);
        
        vecs = [xs(:)-vp(vno,1) ys(:)-vp(vno,2)];
        norms = (sum(vecs.*vecs,2)).^0.5;
        thetas{vno}=acosd((vecs./repmat(norms,[1 2])) * (veca'./norma));
        
        % if xs ys outside image
        vecan=[veca(2) -veca(1)];
        tt1=sum(vecs.*repmat(vecan,[length(xs) 1]),2);
        tt2=sum([w/2-vp(vno,1) h/2-vp(vno,2)].*vecan);
        tt3=1*(sign(tt1)==sign(tt2))-1*(sign(tt1)~=sign(tt2));
        
        thetas{vno}(~inimpts)=tt3(~inimpts).*thetas{vno}(~inimpts);
        theta_ind{vno} = thetas{vno}*quantsiz/theta(ii) + 1;
        
    else
        
        vecs = [xs(:)-vp(vno,1) ys(:)-vp(vno,2)];
        norms = (sum(vecs.*vecs,2)).^0.5;
        if vno==1
            veca=[vp(2,1)-vp(vno,1) vp(2,2)-vp(vno,2)];
        else
            veca=[vp(1,1)-vp(vno,1) vp(1,2)-vp(vno,2)];
        end
        
        th_ref = atan2(veca(2),veca(1));
        mat_ref = [cos(th_ref) sin(th_ref);-sin(th_ref) cos(th_ref)];
        vecs = mat_ref * vecs';
        thetas{vno} = atan2(vecs(2,:),vecs(1,:)) * 180/pi;
        tempinds = find(thetas{vno}<0);
        thetas{vno}(tempinds) = thetas{vno}(tempinds)+360;
        theta_ind{vno} = thetas{vno}*quantsiz/360 + 1;
        
    end
    theta_ind{vno} = round(theta_ind{vno});
end


th1=theta_ind{1};
th2=theta_ind{2};


return;