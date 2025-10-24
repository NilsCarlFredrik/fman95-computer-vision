%% CE1 
% load images and setup vl
im_A = imread('a.jpg');
im_B = imread('b.jpg');
% im_A = imread('nils1.jpg');
% im_B = imread('nils2.jpg');
% im_A = imread('tavla1.jpg');
% im_B = imread('tavla2.jpg');

run vl_setup.m;

%% compute sift shit

[fA dA] = vl_sift( single(rgb2gray(im_A)) );
[fB dB] = vl_sift( single(rgb2gray(im_B)) );

%% plots 
figure(1)
imagesc(im_A);
hold on
vl_plotframe(fA);
hold off
axis equal

figure(2)
imagesc(im_B);
vl_plotframe(fB);
hold off
axis equal

%% compute matches
matches = vl_ubcmatch(dA,dB);

% xA = fA(1:2, matches(1,:));
% xB = fB(1:2, matches(1,:));

xA = [fA(1,matches (1 ,:)); fA(2,matches (1 ,:))];
xB = [fB(1,matches (2 ,:)); fB(2,matches (2 ,:))];

perm = randperm(size(matches ,2));

%% plot lines

figure(3);
imagesc ([im_A im_B]);
hold on;
plot([xA(1,perm (1:10)); xB(1,perm (1:10))+ size(im_A ,2)], ...
    [xA(2,perm (1:10));  xB(2,perm (1:10))] ,'-') ;
axis equal
hold  off;


%% compute H
cp = [];
for i=1:5000
    % Pick 4 random points
    rand = randperm(length(xA), 4);
    x = [xA(:,rand); ones(1,4)];
    y = [xB(:,rand); ones(1,4)];
    
    % Compute estimated H
    Me = zeros(12, 13);
    for j = 1:length(x)
        Me(3*j-2:3*j,:) = ...
         [x(:,j)' zeros(1,6) zeros(1,j-1) y(1,j) zeros(1,4-j);
         zeros(1,3) x(:,j)' zeros(1,3) zeros(1,j-1) y(2,j) zeros(1,4-j);
         zeros(1,6) x(:,j)' zeros(1,j-1) y(3,j) zeros(1,4-j)];
    end
    
    [U,S,V] = svd(Me);
    v = V(1:9, end);
    He = reshape(v, [3 3])';
    
    % Compute estimated points and find those with error < 5 pixels
    xBe = He*[xA ; ones(1, (length(xA)))];
    xBeFlat = pflat(xBe);
    err = sqrt((xB(1,:)-xBeFlat(1,:)).^2 + (xB(2,:)-xBeFlat(2,:)).^2);
    cps = find(err<5);
    
    % Save the best estimate of the loop
    if (length(cps) > length(cp))
        cp = cps;
        Hest = He;
    end
    
end

Hest = Hest./Hest(end, end);

length(cp)

%% 

% Transform image to comon coordinate system
Htform = projective2d(Hest');
Rout = imref2d ( size(im_A) ,[ -200 800] ,[ -400 600]);
% Rout = imref2d ( size(im_A) ,[ -1500 3800] ,[ -1400 4800]);
[Atransf] = imwarp(im_A,Htform ,'OutputView',Rout );
Idtform = projective2d(eye(3));
[Btransf] = imwarp(im_B, Idtform ,'OutputView',Rout );
AB = Btransf;
AB( Btransf < Atransf ) = Atransf( Btransf < Atransf );
figure(4)
imagesc( Rout . XWorldLimits , Rout . YWorldLimits ,AB );







%% Ce2

load("compEx2data.mat");
im1 =  imread("im1.jpg");
im2 =  imread("im2.jpg");
x1 = [x{1}; ones(1, length(x{1}))];
x2 = [x{2}; ones(1, length(x{1}))];

cp = [];
l1 = [];
l2 = [];
E = [];
F = [];

iter = 1000;
for i = 1:iter
    % get five random points
    randsel = randperm(length(x1), 5);
    x1r = x1(:, randsel);
    x2r = x2(:, randsel);
    
    % Calibrate with K
    x1rn = K\x1r;
    x2rn = K\x2r;
    
    % Calculate essential matrix estimation
    Ei = fivepoint_solver(x1rn, x2rn);
    
    for j = 1:length(Ei)
        
        %compute F
        Fe = K'\Ei{j}/K; 
        Fj = Fe./Fe(3,3);
        
        % compute lines
        l1j = Fj'*x2;
        l1j = l1j./sqrt(repmat(l1j(1,:).^2 + l1j(2,:).^2, [3 1])); %normalize
        
        l2j = Fj*x1;
        l2j = l2j./sqrt(repmat(l2j(1,:).^2 + l2j(2,:).^2, [3 1])); %normalize
        
        % compute distances to lines 
        dist1 = abs(sum(l1j.*x1));
        dist2 = abs(sum(l2j.*x2));
        
        % compute inliers
        cp1 = find(dist1 < 5);
        cp2 = find(dist2 < 5);
        
        % compute matching inliers for both images
        [cpj usch] = intersect(cp1, cp2);
        
        % save the best results
        if (length(cp) < length(cpj)) 
            cp = cpj;
            l1 = l1j;
            l2 = l2j;
            E = Ei{j};
            F = Fj;
        end
    end
end

disp('number of matching inliers:')
length(cp)

%% Plot 

figure(5)
imagesc(im1)
hold on
plot(x1(1,cp), x1(2,cp), 'y*', 'Markersize', 10) 
rital(l1(:,cp))
hold off
axis equal
            
figure(6)
imagesc(im2)
hold on
plot(x2(1,cp), x2(2,cp), 'y*', 'Markersize', 10) 
rital(l2(:,cp))
hold off
axis equal

%% compute cameras and triangulate

x1n = K\x1;
x2n = K\x2;

[U S V] = svd(E);

u3 = U(:,end);
W = [0 -1  0 ; 1 0 0 ; 0 0 1];

P1 = [eye(3) [0;0;0]];
% CAmera solutions
P2 = {[U*W*V' u3], [U*W*V' -u3], [U*W'*V' u3], [U*W'*V' -u3]};

in_front = 0;
bestP2 = 0;
bestX = 0;

for i=1:4
    
    % Triangulate
    X = [];
    for j = 1:length(x1n)
        Me = [P1, x1n(:,j) zeros(3,1) ; ...
               P2{i} zeros(3,1) x2n(:,j)];
        [Ut, St, Vt] = svd(Me);
        X(:,j) = pflat(Vt(1:4,end));
    end
    
    % Check number of points in front of camera
    temp = sum(P1(3,:)*X > 0 & P2{i}(3,:)*X > 0);
    if  temp >= in_front
        bestP2 = P2{i};
        bestX = X;
        in_front = temp;
        i
    end
end

P2f = K*bestP2;
X = bestX;
x1proj = pflat(K*P1*X);
x2proj = pflat(P2f*X);

%% plots projections of 3d and 3d model with cameras'

figure(7)
imagesc(im1)
hold on
plot(x1(1,cp), x1(2,cp), 'r*', 'Markersize', 5)
plot(x1proj(1,cp), x1proj(2,cp), 'b+', 'Markersize', 5)
hold off
axis equal

figure(8)
imagesc(im2)
hold on
plot(x2(1,:), x2(2,:), 'r*', 'Markersize', 5)
plot(x2proj(1,:), x2proj(2,:), 'b+', 'Markersize', 5)
hold off
axis equal
%%
figure(22)
plot3(X(1,cp), X(2,cp), X(3,cp), '.')
hold on
plotcams({P1, P2f})
hold off
axis equal


%% Plot histograms of errors

figure(9)
hist(sum(abs(x1(:,cp)-x1proj(:,cp))), 100);

figure(10)
hist(sum(abs(x2(:,cp)-x2proj(:,cp))), 100);

%% Compute rms errors
x1rms = mean(sqrt(sum((x1(:,cp)-x1proj(:,cp)).^2)))
x2rms = mean(sqrt(sum((x2(:,cp)-x2proj(:,cp)).^2)))



%% CE3

P = {K*P1 P2f};
U = X(:,cp);
u = {x1(:, cp) x2(:, cp)};

iter = 10;
gammak = 1;
rnorms = zeros(1, 10);
[r,J] = LinearizeReprojErr(P,U,u);
for i = 1:iter
    lastr = r;
    while 1
        deltav = -gammak*J'*r;
        [Pnew , Unew ] = update_solution(deltav ,P,U);
        [r,J] = LinearizeReprojErr(Pnew,Unew,u); 
        if norm(r) >= norm(lastr)
            gammak = .5* gammak;
        else
            P = Pnew;
            U = Unew;
            break;
        end
    end
    rnorms(i) = norm(r);
end
%%
figure(11)
plot(linspace(1, iter, iter), rnorms, '*-')

Pnew{2} = Pnew{2}./Pnew{2}(end, end);
x1e = pflat(Pnew{1}*Unew);
x2e = pflat(Pnew{2}*Unew);


figure(12)
imagesc(im1)
hold on
plot(x1(1,:), x1(2,:), 'b*', 'Markersize', 5) 
plot(x1e(1,:), x1e(2,:), 'r.', 'Markersize', 5)  
hold off
axis equal

figure(13)
imagesc(im2)
hold on
plot(x2(1,:), x2(2,:), 'b*', 'Markersize', 5) 
plot(x2e(1,:), x2e(2,:), 'r.', 'Markersize', 5) 
hold off
axis equal

% RMS error
disp('RMS x1')
mean(sqrt(sum((x1(:, cp)-x1e).^2)))
disp('RMS x2')
mean(sqrt(sum((x2(:, cp)-x2e).^2)))

%% CE4
P = {K*P1 P2f};
U = X(:,cp);
u = {x1(:, cp) x2(:, cp)};

iter = 10;
lambda = 1;
rnorms = zeros(1, 10);
[r,J] = LinearizeReprojErr(P,U,u);
C = J'*J+lambda*speye(size(J ,2));
c = J'*r;
deltav = -C\c;
[Pnew , Unew ] = update_solution(deltav ,P,U);
for i = 1:iter
    lastr = r;
    C = J'*J+lambda*speye(size(J ,2));
    c = J'*r;
    deltav = -C\c;
    [Pnew , Unew ] = update_solution(deltav ,Pnew,Unew);
    [r,J] = LinearizeReprojErr(Pnew,Unew,u);
    rnorms(i) = norm(r); 
end

%% plots

figure(14)
plot(linspace(1, iter, iter), rnorms,  '*-')
%%
x1e = pflat(Pnew{1}*Unew);
x2e = pflat(Pnew{2}*Unew);

figure(15)
imagesc(im1)
hold on
plot(x1(1,:), x1(2,:), 'b*', 'Markersize', 5) 
plot(x1e(1,:), x1e(2,:), 'r.', 'Markersize', 5) 
hold off
axis equal

figure(16)
imagesc(im2)
hold on
plot(x2(1,:), x2(2,:), 'b*', 'Markersize', 5) 
plot(x2e(1,:), x2e(2,:), 'r.', 'Markersize', 5) 
hold off
axis equal
%%
% RMS error
disp('RMS x1')
mean(sqrt(sum((x1(:, cp)-x1e).^2)))
disp('RMS x2')
mean(sqrt(sum((x2(:, cp)-x2e).^2)))
