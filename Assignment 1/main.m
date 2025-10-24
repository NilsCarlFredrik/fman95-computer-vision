%
% READ ME: assignment1data need to be added to path to execute script
%
%% CE 1
load('compEx1.mat')

% Projects and plots points
flatmat2 = pflat(x2D);
figure(1)
plot(flatmat2(1,:), flatmat2(2,:),'x')
axis equal

flatmat3 = pflat(x3D);
figure(2)
plot3(flatmat3(1,:), flatmat3(2,:), flatmat3(3,:), 'x')
axis equal

%% CE 2

load('compEx2.mat')
im = imread('compEx2.JPG');

% Plot image, points and lines
figure(3)
imagesc(im)
colormap gray
hold on

plot(p1(1,:),p1(2,:),'*')
plot(p2(1,:),p2(2,:),'*')
plot(p3(1,:),p3(2,:),'*')

l1 = linsolve([p1' ; 0 0 1],[0;0;1]);
l2 = linsolve([p2' ; 0 0 1],[0;0;1]);
l3 = linsolve([p3' ; 0 0 1],[0;0;1]);

rital(l1)
rital(l2)
rital(l3)

% Not parallel in 3D since not in 2D

% Compute and plot intersection
intsec = linsolve([l2'; l3'; 0 0 1],[0;0;1]);
plot(intsec(1), intsec(2), '+', 'color', 'yellow')

% Distance
d = abs(l1(1)*intsec(1) + l1(2)*intsec(2) + l1(3))/sqrt(l1(1)^2+l1(2)^2)


hold off

%% CE 3
clear all
load('compEx3.mat');

H1 = [sqrt(3) -1 1 ; 1 sqrt(3) 1 ; 0 0 2];
H2 = [1 -2 1 ; 1 1 0 ; 0 0 1];
H3 = [1 1 0 ; 0 2 0 ; 0 0 1];
H4 = [sqrt(3) -1 1 ; 1 sqrt(3) 1 ; 1/4 1/2 2];

%Compute start and end points
start1 = pflat(H1*[startpoints ; ones(1,42)]);
end1 = pflat(H1*[endpoints ; ones(1,42)]);

start2 = pflat(H2*[startpoints ; ones(1,42)]);
end2 = pflat(H2*[endpoints ; ones(1,42)]);

start3 = pflat(H3*[startpoints ; ones(1,42)]);
end3 = pflat(H3*[endpoints ; ones(1,42)]);

start4 = pflat(H4*[startpoints ; ones(1,42)]);
end4 = pflat(H4*[endpoints ; ones(1,42)]);

% Plot each transformation together with original grid
figure(1)
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)], 'b-');
hold on
plot([start1(1,:); end1(1,:)], [start1(2,:); end1(2,:)], 'r-');
axis equal
hold off

figure(2)
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)], 'b-');
hold on
plot([start2(1,:); end2(1,:)], [start2(2,:); end2(2,:)], 'g-');
axis equal
hold off

figure(3)
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)], 'b-');
hold on
plot([start3(1,:); end3(1,:)], [start3(2,:); end3(2,:)], 'm-');
axis equal
hold off

figure(4)
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)], 'b-');
hold on
plot([start4(1,:); end4(1,:)], [start4(2,:); end4(2,:)], 'c-');
axis equal
hold off

%% CE 4

clear all
close all

load('compEx4.mat')

I1 = imread('compEx4im1.jpg');
I2 = imread('compEx4im2.jpg');

figure(1)
imagesc(I1)
colormap gray

figure(2)
imagesc(I2)
colormap gray

%%
% Computes cameracenters
camCent1 = pflat(null(P1))
camCent2 = pflat(null(P2))

% Compute principal axes
camDir1 = P1(3,1:3)
camDir2 = P2(3,1:3)

% Project data points
Uc = pflat(U);

% Plots 3D-model
figure(3)
plot3(Uc(1,:), Uc(2,:), Uc(3,:), '.', 'Markersize',2)

hold on
% Plot cameracenters and principial vectors
plot3(camCent1(1), camCent1(2), camCent1(3), '+', 'Markersize', 8)
quiver3(camCent1(1),camCent1(2),camCent1(3),camDir1(1),camDir1(2),camDir1(3),7)

plot3(camCent2(1), camCent2(2), camCent2(3), '+', 'Markersize', 8)
quiver3(camCent2(1),camCent2(2),camCent2(3),camDir2(1),camDir2(2),camDir2(3),7)

hold off

%%

% Plots original images with projection points from 3D-model.
proj1 = pflat(P1*U);
figure(4)
imagesc(I1)
colormap gray
hold on
plot(proj1(1,:), proj1(2,:), '.', 'Markersize', 2)


proj2 = pflat(P2*U);
figure(5)
imagesc(I2)
colormap gray
hold on
plot(proj2(1,:), proj2(2,:), '.', 'Markersize', 2)
























