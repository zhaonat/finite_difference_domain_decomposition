%% create a triangle by specifying three points


point1 = [45,60];
point2 = [10,0]; 
point3 = [90,10];

TriangleCoords = [point1; point2; point3];
xmin = min(TriangleCoords(:,1));
xmax = max(TriangleCoords(:,1));
m1 = (point1(2)-point2(2))/(point1(1)-point2(1));
m2 = (point1(2)-point3(2))/(point1(1)-point3(1));
m3 = (point3(2)-point2(2))/(point3(1)-point2(1));
b1 = point1(2) -m1*point1(1);
b2 = point1(2) -m2*point1(1);
b3 = point2(2) -m3*point2(1);
line1 = [m1,b1];
line2 = [m2,b2];
line3 = [m3,b3];

x = linspace(1,100,100)
figure;
plot(m1*x+b1)
hold on;
plot(m2*x+b2)
plot(m3*x+b3)

Nx = 100; Ny=100;
N =[Nx,Ny];
M = ones(100,100);

[X,Y] = ind2sub(N,1:prod(N));





