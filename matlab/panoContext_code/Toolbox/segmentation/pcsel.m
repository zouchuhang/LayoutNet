function [ colors, pc ] = pcsel( centers, ncolor )
%PCSEL: Patch Color Selector
%   [ colors, pc ] = pcsel( centers, ncolor ) returns patch colors and
%   patch color codes in colors and pc respectively according to input
%   arguments:
% centers: n x 2 array of n (x,y) coordinates of n patch centers.
% ncolor: number of colors to use
%
% See also k_means, plot
%
% By Yi Cao at Cranfield University on 27th October 2010
%
% Example:
%{ 
% Some random data
N=20000;
X = [randn(N,2)+ones(N,2); randn(N,2)-ones(N,2)];
% separated in 20 patches
[cidx,c] = k_means(X,20);
% use all 7 colors available in plot function
[colors, pc] = pcsel(c);
subplot(121)
for k=1:20,plot(X(cidx==k,1),X(cidx==k,2),pc{colors(k)}),hold on, end
hold off
axis('square')
% Demonstrate four color map theorem
[colors, pc] = pcsel(c,4);
subplot(122)
for k=1:20,plot(X(cidx==k,1),X(cidx==k,2),pc{colors(k)}),hold on, end
hold off
axis('square')
set(gcf,'Position',[100 100 560 280])
%}
if nargin < 2
    ncolor = 7;
end
ncolor = min(ncolor,7);
[~,ic] = sort(sum(abs(centers).^2,2));
n = size(centers,1);
colors(ic) = mod(0:n-1,ncolor)+1;
pc = {'b.','g.','r.','c.','m.','y.','k.'};
pc = pc(1:ncolor);
for k1=ncolor+1:n
    k0 = ic(1:k1-1);
    k = ic(k1);
    cd = sum(abs(centers(k0,:) - centers(k+zeros(k1-1,1),:)).^2,2);
    s = 0;
    cc = colors(k0);
    for c = 1:ncolor
        id = min(cd(cc==c));
        if id > s
            s = id;
            is = c;
        end
    end
    colors(k) = is;
    kc = find(cd'==s & cc==is, 1);
    kk = k0(kc);
    k0(kc) = k;
    k = kk;
    kd = sum(abs(centers(k0,:) - centers(k+zeros(numel(k0),1),:)).^2,2);
    ks = 0;
    cc = colors(k0);
    for c = 1:ncolor
        id = min(kd(cc==c));
        if id > ks
            ks = id;
            ik = c;
        end
    end
    if ks > s
        colors(k) = ik;
    end
end
end

