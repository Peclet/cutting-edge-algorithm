function [f1,y_Ltan,f2,y_Rtan, L_slope, R_slope, xmin, y_xmin, y_xmax, x_ymax, ymax, y_mx, Pi, Qi, Si, Ti, xi, yi] = Geometry(X,Y,T1_X, T1_Y, T2_X, T2_Y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n=size(X,1);
xmin = min(X);
y_xmin = unique(Y(find(X==xmin)));
xmax=max(X);
y_xmax = unique(Y(find(X==xmax)));
y_mx = max(Y);
L1=zeros(n,2);
R1=zeros(n,2);
for i = 1:n
if X(i,1) <= T1_X
    L1(i,1) = X(i,1);
    L1(i,2) = unique(Y(find(X==X(i,1))));
elseif X(i,1) > T2_X
    R1(i,1) = X(i,1);
    R1(i,2) = unique(Y(find(X==X(i,1))));
end
end
Left_Data = L1(any(L1,2),:);
Right_Data = R1(any(R1,2),:);
h = mean(diff(Left_Data(:,1)));
dy = gradient(Left_Data(:,2), h);
L_slope = mean(dy);
h2 = mean(diff(Right_Data(:,1)));
dy2 = gradient(Right_Data(:,2), h);
R_slope = mean(dy2);
LConst = T1_Y-(L_slope*T1_X);
RConst = T2_Y - (R_slope*T2_X);
Bisect = mean(X);  %%%% Locates the X-axis Mid-Point
n = size(xmin:0.1:Bisect,2);
f1 = T1_X:0.1:(Bisect+(Bisect/1.5));
y_Ltan = (L_slope*f1')+LConst;
f2 = (Bisect-(Bisect/5)):0.005:T2_X;
y_Rtan = (R_slope(1,1)*f2') + RConst;  
%plot(f1, y_Ltan, 'b--');
%plot(f2, y_Rtan, 'b--');
[xi,yi] = polyxpoly(f1,y_Ltan,f2,y_Rtan);
[Pi,Qi] = polyxpoly(f2,y_Rtan,X,Y);
[Si,Ti] = polyxpoly(f1,y_Ltan,X,Y);
ymax = unique(max(Y));
x_ymax = unique(X(find(Y==ymax)));
mapshow(xi,yi,'DisplayType','point','Marker','o')
mapshow(Pi,Qi,'DisplayType','point','Marker','o')
mapshow(Si,Ti,'DisplayType','point','Marker','o')
plot([T2_X, xi],[T2_Y, yi],'b--');
plot([T1_X, xi],[T1_Y, yi],'b--');
plot(T2_X,T2_Y,'k*','LineWidth',4,'MarkerSize', 8)
plot(T1_X,T1_Y,'k*','LineWidth',4,'MarkerSize', 8)
plot(xi,yi,'r*','LineWidth',4,'MarkerSize', 8)
plot(x_ymax,ymax,'k*','LineWidth',4,'MarkerSize', 8)
end

