U6 = 1/sqrt(2)*[1 -1 0 ; 1 1 0 ; 0 0 sqrt(2)];
u6 = [0 ; 0 ; 1];
V6 = [1 0 0 ; 0 0 -1 ; 0 1 0];
W6 = [0 -1  0 ; 1 0 0 ; 0 0 1];

P6 = {[U6*W6*V6' u6], [U6*W6*V6' -u6], [U6*W6'*V6' u6], [U6*W6'*V6' -u6]};

for i = 1:4
    figure(i)
    plot3([0],[0],[sqrt(2)*(-1)^(i+1)], '*')
    sqrt(2)*((-1)^i)
    hold on
    plotcams({P2e6{i} P2e6{5}})
    hold off
end