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
v1 = V(:,end);
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
P1
[r,q] = rq(P1)

% Ambiguity??














