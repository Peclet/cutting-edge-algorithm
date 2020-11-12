clear all
clc
%%%% Reading the data generated from Zygo
DIR = 'C:\Users\doad224\Documents\Papers\Cutting Edge Paper\Code\K-Factor Test\K-Factor - Copy';
cd (DIR);
fid   = fopen('File.txt', 'r');
% fgetl(fid);
% fgetl(fid);
% fgetl(fid);
buffer =fread(fid,[1,Inf],'*char');
fclose(fid);
buffer = regexprep(buffer,'NO_DATA','Nan');
buffer = regexprep(buffer,'                     ','Nan ');
C = textscan(buffer, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','',...
    'MultipleDelimsAsOne', 1);
for i=1:10
    C{i} = fillmissing(C{i}, 'linear');
end
n=size(C{1},1);
Data = zeros(n,10);
for i = 1:10
    Data(:,i)=C{i};
end
% Data = textread('PVD1R.txt');
n = size(Data,1);
Output = zeros(6,15);





for z = 1:2:6


XY1 = Data(:,z:z+1);
XY1(:,1) = XY1(:,1).*1000;
X = XY1(:,1);

% FITT=fit( XY1(:,1), XY1(:,2)/1000, 'poly5', 'Normalize', 'on' );
% Y = FITT(X);
% plot(XY1(:,1),XY1(:,2)/1000)
% spl = spline(XY1(:,1),Y);
% fnplt(spl)
% hold on

Y = smooth(XY1(:,1),XY1(:,2),0.25,'rloess');
n=size(X,1);
% spl = spline(XY1(:,1),Y);
% tt = ppval(spl,X);

Fig = figure;
set(Fig,'defaultLegendAutoUpdate','off');
plot(XY1(:,1),(XY1(:,2)),'g-')
h1 = plot(X,Y,'r-', 'LineWidth',3, 'DisplayName', 'Smoothed Profile');
hold on
set(gca, 'FontName', 'Times New Roman')
xlabel('\it Distance (\mu m)', 'fontweight','bold','fontsize',12)
ylabel('\it Height (\mu m)','fontweight','bold','fontsize',12)

for j = 1
% Fig = figure;
% hold on
%  plot(XY1(:,1),XY1(:,2)/1000,'b.',X,Y,'r-')
% filename = strcat('Fig1_',num2str(j),'.jpg');
% saveas(Fig,filename);
% close(Fig);
% hold off

[TF,S1,S2] = ischange(Y,'linear','MaxNumChanges',3);
counter = 1;
for i = 1:((n-1))
    if S2(i+1,1) == S2(i,1)
    else
        idxr(counter,1)=i;
        counter = counter+1;
    end
end 

if size(idxr,1) == 3
TurnLoc = zeros(3,2);
TurnLoc(:,1) = XY1(idxr);
TurnLoc(:,2) = Y(idxr);
T1_X = TurnLoc(1,1);
T1_Y = TurnLoc(1,2);
T2_X = TurnLoc(3,1);
T2_Y = TurnLoc(3,2);
else
[TF,S1,S2] = ischange(Y,'linear','MaxNumChanges',4)
counter = 1;
for i = 1:((n-1))
    if S2(i+1,1) == S2(i,1)
    else
        idxr(counter,1)=i;
        counter = counter+1;
    end
end 
TurnLoc = zeros(4,2);
TurnLoc(:,1) = XY1(idxr);
TurnLoc(:,2) = Y(idxr);
T1_X = TurnLoc(2,1);
T1_Y = TurnLoc(2,2);
T2_X = TurnLoc(4,1);
T2_Y = TurnLoc(4,2);
end


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
n = size(xmin:0.01:Bisect,2);
f1 = T1_X:0.01:(Bisect+(Bisect/5));
y_Ltan = (L_slope*f1')+LConst;
f2 = (Bisect-(Bisect/5)):0.005:T2_X;
y_Rtan = (R_slope(1,1)*f2') + RConst;  
  plot(f1, y_Ltan, 'b--');
  plot(f2, y_Rtan, 'b--');

[xi,yi] = polyxpoly(f1,y_Ltan,f2,y_Rtan);
mapshow(xi,yi,'DisplayType','point','Marker','o')
[Pi,Qi] = polyxpoly(f2,y_Rtan,X,Y);
[Si,Ti] = polyxpoly(f1,y_Ltan,X,Y);
mapshow(Pi,Qi,'DisplayType','point','Marker','o')
mapshow(Si,Ti,'DisplayType','point','Marker','o')
plot([T2_X, xi],[T2_Y, yi],'b--');
plot([T1_X, xi],[T1_Y, yi],'b--');
legend(h1)


ymax = unique(max(Y));
x_ymax = unique(X(find(Y==ymax)));
plot(T2_X,T2_Y,'k*','LineWidth',4,'MarkerSize', 8)
plot(T1_X,T1_Y,'k*','LineWidth',4,'MarkerSize', 8)
plot(xi,yi,'r*','LineWidth',4,'MarkerSize', 8)
plot(x_ymax,ymax,'k*','LineWidth',4,'MarkerSize', 8)


Sisort = sort(Si);
Tisort = sort(Ti);
Pisort = sort(Pi);
Qisort = sort(Qi);
t1=isempty(Pi)
t2=isempty(Si)
if t1 == 1
Pisort = T2_X;
Qisort = T2_Y;
else
end

if t2 == 1
Sisort = T1_X;
Tisort = T1_Y;
else
end

% s_gamma = sqrt(((yi-Tisort)^2) + ((xi- Sisort)^2));
% s_alpha = sqrt(((yi-Qisort)^2) + ((xi- Pisort)^2));
s_alpha = sqrt(((yi-Tisort(end,1)).^2) + ((xi- Sisort(end,1)).^2));
s_gamma= sqrt(((yi-Qisort(end,1)).^2) + ((xi- Pisort(1,1)).^2));
K_Factor = s_gamma/s_alpha;

 p = polyval(f2, y_Rtan,1);
 G_L = atand(abs((L_slope-R_slope)/(1+(L_slope*R_slope))));
 w = (yi - y_xmin)/(xi-xmin);
% % % % 
 a=atand(p(1));


Dist_r = sqrt(((yi-ymax).^2)+((xi-x_ymax).^2)); %%% computes the shortest distance from the intersection of the lines representing the flank and face to the cutting edge profile
plot([xi, x_ymax],[yi, ymax],'r--'); %%% plot the shortest distance from the intersection of the lines representing the flank and face to the cutting edge profile

plot([xi, xi],[ymax, yi],'r--');
xlabel('Distance (um)')
ylabel('Height (um)')


slope2 = abs((yi-ymax)/(xi-x_ymax));
f = 90-atand(slope2);
Ang_bet = (atan((T1_Y-yi)/(T1_X-xi))-atan((T2_Y-yi)/(T2_X-xi)))*180/pi;
Cir_Rad = Dist_r/(sqrt(1+(1/(tand(Ang_bet/2))))-1);
Output(z,1) = K_Factor;
Output(z,2) = s_alpha;
Output(z,3) = s_gamma;
Output(z,4) = f;
Output(z,5) = Dist_r;
Output(z,6) = Ang_bet;
Output(z,7) = Cir_Rad;


annotation('textbox', [0.62, 0.74, 0.1, 0.1], 'String', "Edge Radius = " + Cir_Rad,'String', "Edge Radius = " + Cir_Rad)
filename = strcat('Fig_',num2str(z),'.jpg');
saveas(Fig,filename);
% close(Fig);
hold off

end
end
Output = Output(any(Output,2),:);
Output(z+1,:) = mean(Output);


