q1 = imread('cube1.jpg');
q2 = imread('cube2.jpg');
run vl_setup.m

[f1, d1] = vl_sift( single(rgb2gray(q1)), 'PeakThresh', 1);
[f2, d2] = vl_sift( single(rgb2gray(q2)), 'PeakThresh', 1);

figure(14)
imagesc(q1);
hold on
vl_plotframe(f1);
hold off
axis equal

figure(15)
imagesc(q2);
vl_plotframe(f2);
hold off
axis equal


[matches ,scores] = vl_ubcmatch(d1,d2);

x1 = [f1(1,matches (1 ,:));f1(2,matches (1 ,:))];
x2 = [f2(1,matches (2 ,:));f2(2,matches (2 ,:))];

perm = randperm(size(matches ,2));
figure(16);
imagesc ([q1 q2]);
hold on;
plot([x1(1,perm (1:10)); x2(1,perm (1:10))+ size(q1 ,2)], ...
    [x1(2,perm (1:10));  x2(2,perm (1:10))] ,'-') ;

hold  off;

%% Ce 5 - Set up and solves DLT for trianhylation.

Xtriag = [];

for i = 1:length(x1)
    Mce5 = [cameras{1} -[x1(:,i);1]    zeros(3,1) ;
            cameras{2} zeros(3,1)       -[x2(:, i);1]];
    [U,S,V] = svd(Mce5);
    v = V(:,end);
    Xtriag = [Xtriag v(1:4)];
end

%% Project triangulated points onto images.
Xflat = pflat(Xtriag);
xproj1 = pflat(cameras{1}*Xflat);
xproj2 = pflat(cameras{2}*Xflat);



perm = randperm(size(matches ,2));
figure(16)
% imagesc ([q1 q2]);
hold on;
plot([xproj1(1,perm (1:10)); xproj2(1,perm (1:10))+ size(q1 ,2)], ...
    [xproj1(2,perm (1:10));  xproj2(2,perm (1:10))] ,'--') ;
hold  off;

good_points = (sqrt(sum((x1-xproj1 (1:2 ,:)).^2))  < 3 &   ...
    sqrt(sum((x2-xproj2 (1:2 ,:)).^2))  < 3);


Xgood = Xflat(:,good_points);


%%
figure(17)
plot3(Xgood(1,:), Xgood(2,:), Xgood(3,:), '.')
hold on
plot3([ Xmodel(1,startind );  Xmodel(1,endind )],...
    [Xmodel(2,startind );  Xmodel(2,endind )],...
    [Xmodel(3,startind );  Xmodel(3,endind)],'b-');
plotcams(cameras);
hold off
axis equal
%%

p1n = rq(P1)\P1;
p2n = rq(P2)\P2;
























