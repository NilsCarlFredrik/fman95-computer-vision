%% CE 1
load('compEx1data.mat')
%%
figure(1)
plot3(X(1,:), X(2,:), X(3,:), '.', 'Markersize', 1)
hold on
plotcams(P)
axis equal
hold off

%% compute and plot projection and original points from camera 1.
im = imread(imfiles{1});
visible = isfinite(x{1}(1,:));

x1 = pflat(P{1}*X);

figure(2)
imagesc(im)
hold on
plot(x{1}(1,visible), x{1}(2,visible), 'r*')
plot(x1(1,visible), x1(2,visible), 'b*', 'Markersize', 2)
axis equal
hold off

%% Compute and reconstructs using transformation T1
T1 = [1     0   0   0 ;
      0     4   0   0 ;
      0     0   1   0 ;
      1/10 1/10 0   1];
  
T1X = pflat(T1*X);

figure(3)
plot3(T1X(1,:), T1X(2,:), T1X(3,:), '.', 'Markersize', 1)
hold on
plotcams(P)
axis equal
hold off

%% Compute and reconstructs using transformation T2

T2 = [1     0   0   0 ;
      0     1   0   0 ;
      0     0   1   0 ;
      1/16 1/16 0   1];
  
T2X = pflat(T2*X);

figure(4)
plot3(T2X(1,:), T2X(2,:), T2X(3,:), '.', 'Markersize', 1)
hold on
plotcams(P)
axis equal
hold off

%% All reconstructions
figure(5)
plot3(X(1,:), X(2,:), X(3,:), '.', 'Markersize', 3)
hold on
plot3(T1X(1,:), T1X(2,:), T1X(3,:), '.', 'Markersize', 3)
plot3(T2X(1,:), T2X(2,:), T2X(3,:), '.', 'Markersize', 3)
plotcams(P)
axis equal
hold off

%% All projected

T1x1 = pflat(P{1}/T1*T1X);
T2x1 = pflat(P{1}/T2*T2X);


figure(6)
imagesc(im)
hold on
plot(x{1}(1,visible), x{1}(2,visible), 'b*')
plot(x1(1,visible), x1(2,visible), 'c*')
plot(T1x1(1,visible), T1x1(2,visible), 'r*')
plot(T2x1(1,visible), T2x1(2,visible), 'y.')


hold off

%% CE 2

K1 = rq(P{1}*T1)
K2 = rq(P{2}*T2)


%% CE 3 - Load and plot images
clear all
load('compEx3data.mat');
q1 = imread('cube1.JPG');
q2 = imread('cube2.JPG');

figure(7)
imagesc(q1)
hold on
plot(x{1}(1,:), x{1}(2,:), '*')
axis equal
hold off

figure(8)
imagesc(q2)
hold on
plot(x{2}(1,:), x{2}(2,:), '*')
axis equal
hold off

%% Compute xTilde and N for figure 1
mean1 = mean(x{1}(1:2,:),2);
std1 = std(x{1}(1:2,:),0,2);

N1 = [1/std1(1) 0           -mean1(1)/std1(1);
      0         1/std1(2)   -mean1(2)/std1(2);
      0         0           1               ];
 x1Tilde = N1*x{1};

%% Plot fig 1 osv
figure(9)
plot(x1Tilde(1,:), x1Tilde(2,:), '*')
hold on
plot([ x1Tilde(1,startind );  x1Tilde(1,endind )],...
    [x1Tilde(2,startind );  x1Tilde(2,endind )],...
    'b-');
hold off
axis equal

figure(10)
plot3(Xmodel(1,:), Xmodel(2,:), Xmodel(3,:), '*')
hold on
plot3([ Xmodel(1,startind );  Xmodel(1,endind )],...
    [Xmodel(2,startind );  Xmodel(2,endind )],...
    [Xmodel(3,startind );  Xmodel(3,endind)],'b-');
hold off
axis equal

%% Assemble M-matrix and compute v, norm(M*v) etc.

X=[Xmodel;ones(1, 37)];
M1 = [];

for i = 1:18
    M1(i*3 -2, 1:4) = X(1:4,i)';
    M1(i*3 -1, 5:8) = X(1:4,i)';
    M1(i*3 , 9:12) = X(1:4,i)';
    M1(i*3 -2:i*3, i+12) = -x1Tilde(1:3, i);
end

[U,S,V] = svd(M1);
v1 = V(:,end)
min(diag(S))
norm(M1*v1) %same regardless of which column?

P1Tilde = [v1(1:4)';v1(5:8)';v1(9:12)']
P1 = N1\P1Tilde

x1 = P1*X;
x1flat = pflat(x1);

%%
figure(7)
imagesc(q1)
hold on
plot(x{1}(1,:), x{1}(2,:), '*')
plot(x1flat(1,:),x1flat(2,:), 'y*')
axis equal
hold off

%% Compute xTilde and N for figure 
mean2 = mean(x{2}(1:2,:),2);
std2 = std(x{2}(1:2,:),0,2);

N2 = [1/std2(1) 0           -mean2(1)/std1(1);
      0         1/std2(2)   -mean2(2)/std1(2);
      0         0           1               ];
 
x2Tilde = N2*x{2};

%% Plot
figure(10)
plot(x2Tilde(1,:), x2Tilde(2,:), '*')
hold on
plot([ x2Tilde(1,startind );  x2Tilde(1,endind )],...
    [x2Tilde(2,startind );  x2Tilde(2,endind )],...
    'b-');
hold off
axis equal

figure(11)
plot3(Xmodel(1,:), Xmodel(2,:), Xmodel(3,:), '*')
hold on
plot3([ Xmodel(1,startind );  Xmodel(1,endind )],...
    [Xmodel(2,startind );  Xmodel(2,endind )],...
    [Xmodel(3,startind );  Xmodel(3,endind)],'b-');
hold off
axis equal

%% Assemble M-matrix and compute v, norm(M*v) etc.

X=[Xmodel;ones(1, 37)];
M2 = [];

for i = 1:18
    M2(i*3 -2, 1:4) = X(1:4,i)';
    M2(i*3 -1, 5:8) = X(1:4,i)';
    M2(i*3 , 9:12) = X(1:4,i)';
    M2(i*3 -2:i*3, i+12) = -x2Tilde(1:3, i);
end

[U,S,V] = svd(M2);
v2 = V(:,end);
min(diag(S))
norm(M2*v2) %same regardless of which column?

P2Tilde = [v2(1:4)';v2(5:8)';v2(9:12)']
P2 = N2\P2Tilde

x2 = P2*X;
x2flat = pflat(x2);

%%
figure(12)
imagesc(q1)
hold on
plot(x{1}(1,:), x{1}(2,:), '*')
plot(x1flat(1,:),x1flat(2,:), 'y*')
axis equal
hold off

%%
cameras = {P1,P2};

figure(13)
plot3(Xmodel(1,:), Xmodel(2,:), Xmodel(3,:), '*')
hold on
plot3([ Xmodel(1,startind );  Xmodel(1,endind )],...
    [Xmodel(2,startind );  Xmodel(2,endind )],...
    [Xmodel(3,startind );  Xmodel(3,endind)],'b-');
plotcams(cameras);
hold off
set(gca, 'Zdir', 'reverse')
axis equal

%%

[R1,Q1] = rq(P1);
[R2,Q2] = rq(P2);
K1ce3 = R1./R1(3,3)
K2ce3 = R2./R2(3,3)

%% CE 4

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
