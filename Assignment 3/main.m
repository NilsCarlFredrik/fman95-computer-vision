%% CE1

load('compEx1data.mat');

meanX = [mean(x{1}(1,:)), mean(x{2}(1,:))];
meanY = [mean(x{1}(2,:)), mean(x{2}(2,:))];

stdX = [std(x{1}(1,:)), std(x{2}(1,:))];
stdY = [std(x{1}(2,:)), std(x{2}(2,:))];

N1 = [1/stdX(1) 0           -meanX(1)/stdX(1) ; 
      0         1/stdY(1)   -meanY(1)/stdY(1) ; 
      0         0           1                 ];
  
N2 = [1/stdX(2) 0           -meanX(2)/stdX(2) ; 
      0         1/stdY(2)   -meanY(2)/stdY(2) ; 
      0         0           1                 ];


xn = {N1*x{1}; N2*x{2}};

% Create M
M=zeros(length(xn{1}),9);
for i = 1:length(xn{1})
    %xx = xn{2}(:,i)*xn{1}(:,i)';
    %M(i,:)=xx(:)';
    M(i,1:3) = xn{1}(1, i).*xn{2}(1:3, i);
    M(i,4:6) = xn{1}(2, i).*xn{2}(1:3, i);
    M(i,7:9) = xn{1}(3, i).*xn{2}(1:3, i);
end

% Svd 
[U,S,V] = svd(M);

% last column f_ij
v = V(:,end);

%assemble fundamental matrix
Fn = reshape(v, [3 3]);

det(Fn) %snall enough?
%figure(1)
%plot(diag(xn{2}'*Fn*xn{2}))

F = N2'*Fn*N1;
F = F./F(3,3)
l = F*x{1};
l = l./sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2, [3  1]));

%% plots

kronan1 = imread('kronan1.JPG');
kronan2 = imread('kronan2.JPG');

r = randi([1 length(x{1})], 1, 20);

figure(2)
imagesc(kronan2)
hold on
plot(x{2}(1,r), x{2}(2,r), 'y*', 'Markersize', 10) 
rital(l(:,r))
hold off
axis equal

figure(3)
hist(abs(sum(l.*x{2})), 100);
mean(abs(sum(l.*x{2})))

%% CE2

% Calculate P2 from F
[Ue2, Se2, Ve2] = svd(F');
e2 = Ve2(:, end);
% e2 = -e2;
% e2(3) = 0; %?????????
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
P1 = [eye(3, 3) zeros(3, 1)];
P2 = [e2x*F e2];

P1n = N1*P1;
P2n = N2*P2;

% Triangulation with DLT
Xa = ones(4, length(xn{1}));
for i = 1:length(xn{1})
    Me2 = [P1n(1:3,:) xn{1}(1:3,i) zeros(3,1);P2n(1:3,:) zeros(3, 1) xn{2}(1:3,i)];
    [Ue2, Se2, Ve2] = svd(Me2);
    ve2 = Ve2(:, end);
    Xa(1:4, i) = ve2(1:4);
end

% P1 = P1(:,[1 2 4 3]);
% P2 = P2(:,[1 2 4 3]);
Xa2 = pflat(Xa([1 2 4 3],:));
Xa = pflat(Xa);

x1p = pflat(P1*Xa);
x2p = pflat(P2*Xa);

%% Plots

figure(3)
imagesc(kronan2);
hold on;
axis equal;
randx = x{2}(:,r');
randl = l(:, r');
plot(x{2}(1,:), x{2}(2, :), 'r*');
plot(x2p(1, :), x2p(2, :), 'b+');


figure(4)
plot3(Xa2(1, :), Xa2(2, :), Xa2(3, :), '.');


%% CE3
load('compEx3data.mat')

% Normalize
xnorm = {K\x{1}, K\x{2}};

% Create M
Me3=zeros(length(xnorm{1}),9);
for i = 1:length(xnorm{1})
    %xx = xnorm{2}(:,i)*xnorm{1}(:,i)';
    %Me3(i,:)=xx(:)';
    Me3(i,1:3) = xnorm{1}(1, i).*xnorm{2}(1:3, i);
    Me3(i,4:6) = xnorm{1}(2, i).*xnorm{2}(1:3, i);
    Me3(i,7:9) = xnorm{1}(3, i).*xnorm{2}(1:3, i);
end

% solve with svd
[Ue3, Se3, Ve3] = svd(Me3);
ve3 = Ve3(:,end);

min(diag(Se3))
norm(Me3*ve3)


% create E
Eapprox = reshape(ve3, [3 3]);
[UE, SE, VE] = svd(Eapprox);
% 
% if det(UE*VE')>0 
%     E = UE*diag([1 1 0])*VE';
% else
%     VE = -VE;
%     E = UE*diag([1 1 0])*VE';
% end

if det(UE*VE') < 0
    VE = -VE;
end
E = UE*diag([1 1 0])*VE';

for i = 1:100
    xnorm{2}(:,i)'*E*xnorm{1}(:,i);
end

E=E./E(3,3);

Fe3 = K'\E/K; %?????+??????????????????????
Fe3 = Fe3/Fe3(3,3)
l2 = Fe3*x{1};
% l2 = l2./sqrt(repmat(l2(1 ,:).^2 + l2(2 ,:).^2 ,[3  1]));

%% Plot 

figure(5)
imagesc(kronan2)
hold on
plot(x{2}(1,r), x{2}(2,r), 'y*', 'Markersize', 10) 
rital(l2(:,r))
hold off
axis equal

figure(6)
hist(abs(sum(l2.*x{2})) ,100);


%% CE4

u3 = UE(:,end);
W = [0 -1  0 ; 1 0 0 ; 0 0 1];

% CAmera solutions
P2e4 = {[UE*W*VE' u3], [UE*W*VE' -u3], [UE*W'*VE' u3], [UE*W'*VE' -u3]};

in_front = 0;
bestP2 = 0;
bestX = 0;

for i=1:4
    
    % Triangulate
    X = [];
    for j = 1:length(xnorm{1})
        Me4 = [P1, xnorm{1}(:,j) zeros(3,1) ; ...
               P2e4{i} zeros(3,1) xnorm{2}(:,j)];
        [Ut, St, Vt] = svd(Me4);
        X(:,j) = pflat(Vt(1:4,end));
    end
    
    % Check number of points in front of camera
    temp = sum(P1(3,:)*X > 0 & P2e4{i}(3,:)*X > 0);
    if  temp >= in_front
        bestP2 = P2e4{i};
        bestX = X;
        in_front = temp;
    end
end


P2f = K*bestP2;
X = bestX;
Xflat = pflat(P2f*X);

%% plots
figure(7)
imagesc(kronan2)
hold on
plot(x{2}(1,:), x{2}(2,:), 'r*', 'Markersize', 5)
plot(Xflat(1,:), Xflat(2,:), 'b+', 'Markersize', 5)
hold off
axis equal

%%
figure(8)
plot3(X(1,:), X(2,:), X(3,:), '.')
hold on
plotcams({P1, P2f})
hold off
axis equal





















