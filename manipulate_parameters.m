%close all; clear all; clc
global theta0 r h kappa sigma l2;
global N; N=101;
global y_guess;

theta0=3*pi/4;
r=1.0;
h=10;
kappa = 1.0;
sigma=2.0;
l2=0.0;

L = 1;


solinit = bvpinit(linspace(0,1,101),@newguess,L);
sol = bvp4c(@mat4ode, @mat4bc, solinit);

xint = linspace(0,1);
Sxint = deval(sol,xint);

t = linspace(0,1,101);
y_guess = deval(sol,t);


%z = 0:0.1:3;
z = sol.y(5,:);
%r = sin(z);
r = sol.y(4,:);
t = linspace(0,2*pi,100);
[T,R] = meshgrid(t,r);
[~,Z] = meshgrid(t,z);
[X,Y] = pol2cart(T,R);
ccc=(sol.y(2,:) + sin(sol.y(1,:))/sol.y(4,:));
[~,C] = meshgrid(t,ccc);
%xx = [0.5*X(:); 0.75*X(:); X(:)];
%C = 1;
h = surf(X,Z,Y,C);
set(h, 'EdgeColor', 'none');
c = colorbar;
%c.Label.String = 'Mean Curvature';
%c.Label.FontSize=10;
c.FontSize = 20;
%h.FontSize = 20;
ax=gca;
%ax.XAxis.FontSize = 0.1;
ax.YAxis.FontSize = 20;
ax.ZAxis.FontSize = 20;
ax.XAxis.FontName = 'Times New Roman';
ax.YAxis.FontName = 'Times New Roman';
c.FontName = 'Times New Roman';
set(gca, 'XTick',[]);
%axis off;
box off;
grid off;

xlabel('')
ylabel('Z')
zlabel ('R')


function dydx = mat4ode(x,y,L) % equation being solved
global kappa sigma l2;
%sigma=1;
%l2=1;
dydx = [L*y(2)
        L*(    (-l2*cos(y(1)) + y(3)*sin(y(1)))/(kappa*y(4)) + (cos(y(1))/y(4))*(sin(y(1))/y(4) - y(2)))
        L*((kappa/2)*(y(2)*y(2) - (sin(y(1))*sin(y(1)))/(y(4)*y(4))) + sigma)
        L*cos(y(1))
        L*sin(y(1))];
end
%-------------------------------------------
function res = mat4bc(ya,yb,L) % boundary conditions
global theta0 r h;
%global r;
%r=1.0;
%h=2.0;
res = [ya(1) - theta0; yb(1) - pi + theta0; ya(4)-r; yb(4)-r; ya(5); yb(5)-h];
end
%-------------------------------------------

function v = newguess(q)
global N; 
global y_guess;
 
q = round(q*(N-1));
v = y_guess(:,q+1);
 
end
%...........................................
%..........................................
