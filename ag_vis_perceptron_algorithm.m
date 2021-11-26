%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2021/11/26 Alexander German
% script visualizes the Perceptron Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mylim = [-1.5,1.5];
%sample randomly from the unit circle:
n=10;
r=sqrt(rand(n,1));
alpha=rand(n,1)*2*pi;
x = ones(n,3);
x(:,1) = r.*cos(alpha);
x(:,2) = r.*sin(alpha);
%random separator
d=2*rand-1;
beta=rand*2*pi;
w=[cos(beta),sin(beta),-d];
%labels
l = sign(w*x.');
colors = [[1,0,0];[0,0,0];[0,1,0]];%red (black) and green
%perceptron algorithm
wp=w*0;
wp_h=[];%history
ii_h=[];%history
m=1;
while m~=0
    m=0;
    for ii=1:1:n
        if wp*x(ii,:).'*l(ii)<=0
            wp=wp+x(ii,:)*l(ii);
            wp_h=[wp_h;wp];
            ii_h=[ii_h;ii];
            m=m+1;
        end
    end
end
%calculate separating lines
wp_h_d=-wp_h(:,3)./(wp_h(:,1).^2+wp_h(:,2).^2);%distance to origin
wp_h_y1 = wp_h_d.*wp_h(:,2)-(wp_h(:,1)./wp_h(:,2)).*(mylim(1)-wp_h_d.*wp_h(:,1));%left margin
wp_h_y2 = wp_h_d.*wp_h(:,2)-(wp_h(:,1)./wp_h(:,2)).*(mylim(2)-wp_h_d.*wp_h(:,1));%right margin
%% plot gif
filename = 'Perceptron_algorithm_animated.gif';
hold on
th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);
h = plot(xunit, yunit);
pbaspect([1 1 1])
xlim(mylim)
ylim(mylim)
scatter(x(:,1),x(:,2),[],colors(l+2,:),'filled');
for ii=1:1:numel(ii_h)
    pl=line([mylim(1) mylim(2)],[wp_h_y1(ii) wp_h_y2(ii)]);
    pl.Color = 'black';
    h1=quiver(wp_h_d(ii)*wp_h(ii,1),wp_h_d(ii)*wp_h(ii,2),wp_h(ii,1),wp_h(ii,2),0,'k');
    h2=scatter(x(ii_h(ii),1),x(ii_h(ii),2),[],'k','LineWidth',2);
    frame = getframe(gcf); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if ii == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
    pl.LineStyle = '--';
    delete(h1);
    delete(h2);
end
hold off
%plot 4x4
for kk=1:1:8
subplot(2,4,kk);
hold on
th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);
h = plot(xunit, yunit);
pbaspect([1 1 1])
xlim(mylim)
ylim(mylim)
scatter(x(:,1),x(:,2),[],colors(l+2,:),'filled');
for ii=1:1:kk-1
    pl=line([mylim(1) mylim(2)],[wp_h_y1(ii) wp_h_y2(ii)]);
    pl.Color = 'black';
    pl.LineStyle = '--';
end
pl=line([mylim(1) mylim(2)],[wp_h_y1(kk) wp_h_y2(kk)]);
pl.Color = 'black';
h1=quiver(wp_h_d(kk)*wp_h(kk,1),wp_h_d(kk)*wp_h(kk,2),wp_h(kk,1),wp_h(kk,2),0,'k');
h2=scatter(x(ii_h(kk),1),x(ii_h(kk),2),[],'k','LineWidth',2);
end
hold off
