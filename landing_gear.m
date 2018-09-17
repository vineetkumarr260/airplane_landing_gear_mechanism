clear all;
clc;

p1 = [6 0];
p2 = [-0.946699141100893 -4.247086478782448];
p3 = [2.428700058734363 -5.486475738094875];
p4 = [-5.055787223379352 -8.692224112360492];
p5 = [-4.242640687119286 -4.242640687119286];
p6 = [-10 -10];

M1 = [(p1(1)+p3(1))/2 (p1(2)+p3(2))/2];
M2 = [(p3(1)+p5(1))/2 (p3(2)+p5(2))/2];

S1 = -1*(p1(1)-p3(1))/(p1(2)-p3(2));
S2 = -1*(p3(1)-p5(1))/(p3(2)-p5(2));

C1 = M1(2)-S1*M1(1);
C2 = M2(2)-S2*M2(1);

A = [ S1 -1;
      S2 -1];
B = [ -C1;
      -C2];
X1 = (A\B)';

M3 = [(p2(1)+p4(1))/2 (p2(2)+p4(2))/2];
M4 = [(p4(1)+p6(1))/2 (p4(2)+p6(2))/2];

S3 = -1*(p2(1)-p4(1))/(p2(2)-p4(2));
S4 = -1*(p4(1)-p6(1))/(p4(2)-p6(2));

C3 = M3(2)-S3*M3(1);
C4 = M4(2)-S4*M4(1);

A = [ S3 -1;
      S4 -1];
B = [ -C3;
      -C4];
X2 = (A\B)';

r1 = pdist([X1(1),X1(2); X2(1),X2(2)],'euclidean');
r2 = pdist([X1(1) X1(2); p1(1) p1(2)],'euclidean');
r3 = pdist([p1(1),p1(2); p2(1),p2(2)],'euclidean');
r4 = pdist([X2(1) X2(2); p2(1) p2(2)],'euclidean');

theta_1 = 0;
theta_1fin = 0;
theta_2 = linspace(pi,45*(pi/180),200);
theta_3 = [];
theta_4 = [];
u_2 = [] ;

for i = 1:length(theta_2)
     c = (r3)^2-(r2)^2-(r4)^2-(r1)^2+2*r1*r2*cos(theta_2(i));
     a = 2*r4*((r1)-(r2*cos(theta_2(i))));
     b = -2*r2*r4*sin(theta_2(i));
     y = sqrt((a)^2+(b)^2-(c)^2);
     x = c+a ;
     u = (b+y)/(x);
     theta_4fin = 2*atan(u);
     cosang_3 = ((-r2*cos(theta_2(i)))+(r4*cos(theta_4fin))+r1)/(r3);
     sinang_3 = ((-r2*sin(theta_2(i)))+(r4*sin(theta_4fin)))/(r3) ;
     theta_3fin = atan2((sinang_3),(cosang_3)) ;
     theta_3 = [theta_3,theta_3fin];
     theta_4 = [theta_4,theta_4fin];
end

p0x = [X1(1) X2(1)];
p0y = [X1(2) X2(2)];
plot (p0x,p0y,'-oblack');
hold on;

p1x = [p1(1) p2(1)];
p1y = [p1(2) p2(2)];
plot(p1x,p1y,'r');
hold on;

p2x = [p3(1) p4(1)];
p2y = [p3(2) p4(2)];
plot(p2x,p2y,'g');
hold on;
   
p3x = [p5(1) p6(1)];
p3y = [p5(2) p6(2)];
plot(p3x, p3y,'m');
hold on;

ground=plot ([-15 10], [-18.3 -18.3],'-k');
hold on;
set(ground, 'LineWidth',4);

axis ([-15 10 -20 5]);
title('Airplane Landing Gear Mechanism');
xlabel('X [ft.]');
ylabel('Y [ft.]');
grid on;
hold on;
axis('equal');

for i = 1:length(theta_2)
    k2x = r2*cos(pi+theta_2(i))+X1(1);
    k2y = r2*sin(pi+theta_2(i))+X1(2);
    k4x = (r4)*cos(pi+theta_4(i))+X2(1);
    k4y = (r4)*sin(pi+theta_4(i))+X2(2);
    k5x = ((5+r4)*cos(pi+theta_4(i))+X2(1));
    k5y = ((5+r4)*sin(pi+theta_4(i))+X2(2));
    p2 = line([X1(1) k2x], [X1(2) k2y]);
    p3 = line([k2x k4x],[k2y k4y]);
    p4 = line([X2(1) k5x], [X2(2) k5y]);
    po2 = plot(k2x,k2y,'-oblue');
    po4 = plot(k4x,k4y,'-oblue');
    po5 = plot(k5x,k5y,'-ored');
    k5cen = [k5x k5y];
    circ1 = viscircles(k5cen,3);
    M(i) =getframe;
    delete(p2);
    delete(p3); 
    delete(p4); 
    delete(po2);
    delete(po4);
    delete(po5);
    delete(circ1);     
end

for i =length(theta_2):-1:1  
    k2x = r2*cos(pi+theta_2(i))+X1(1); 
    k2y = r2*sin(pi+theta_2(i))+X1(2); 
    k4x = r4*cos(pi+theta_4(i))+X2(1); 
    k4y = r4*sin(pi+theta_4(i))+X2(2); 
    k6x = (5+r4)*cos(pi+theta_4(i))+X2(1);
    k6y = (5+r4)*sin(pi+theta_4(i))+X2(2);
    p2 = line([X1(1) k2x], [X1(2) k2y]); 
    p3 = line([k2x k4x],[k2y k4y]); 
    p4 = line([X2(1) k6x], [X2(2) k6y]);
    po2 = plot(k2x,k2y,'-oblue');
    po4 = plot(k4x,k4y,'-oblue');
    po6 = plot(k6x,k6y,'-ored');
    k6cen = [k6x k6y]; 
    circ1=viscircles(k6cen,3);
    M(i) =getframe;
    if(i ~= 1)
        delete(p2); 
        delete(p3); 
        delete(p4); 
        delete(po2);
        delete(po4);
        delete(po6);
        delete(circ1);
    end
end
hold off;

figure 
axis equal
axis ([0.5 3.5 -1.5 2.5]) ;
xlabel('Theta 2[radians]'); 
ylabel('Theta 3[radians]'); 
title('Theta 3 vs Theta 2'); 
hold on;
grid on;
for i =1: length(theta_2)
   plot1 = plot(theta_2(i),theta_3(i),'-ro','MarkerSize',2);
   M(i) = getframe;
end 
hold off;

figure 
axis equal 
axis([0.5 3.5 -1.5 3.5]);
xlabel('Theta 2[radians]'); 
ylabel('Theta 4[radians]'); 
title('Theta 4 vs Theta 2');
hold on;
grid on;
for i = 1:length(theta_2)
    plot2 = plot(theta_2(i),theta_4(i),'-ro','MarkerSize',2);
    M(i) = getframe; 
end 
hold off;