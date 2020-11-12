clear all
clc

%%% This Block imports the files and arranges the data in a format the algorithm can work with
File = readtable('C:\Users\doad224\Documents\Papers\Cutting Edge Paper\Code\K-Factor Test\K-Factor - Copy\Tab40D.csv');
File = table2array(File);
N_rows =  size(File,1);
N_Slices = 40;
Split = N_rows/N_Slices;
Z = zeros(Split,N_Slices*2);
X = zeros(Split,N_Slices);
Y = zeros(Split,N_Slices);
for i=0:(N_Slices-1)
    X(:,(i+1)) = File(((Split*i+1):(Split*(i+1))), 1);  %%Extracts the X-Values accross the Initial Dataset, using the specified Slice
    Y(:,(i+1)) = File(((Split*i+1):(Split*(i+1))), 2);  %%Extracts the X-Values accross the Initial Dataset, using the specified Slice 
end
plot(X(:,1),Y(:,1), 'k')

Y = fillmissing(Y,'Linear',1,'EndValues','none');
for i = 1:(N_Slices*2)
    if mod(i,2)==1
        Z(:,i) = X(:,(i+1)/2);
    else
      Z(:,i) = Y(:,i/2);  
    end
end
Data=Z;
%%%  End of File Import Block 


for z = 1:2:(N_Slices*2)
XY1 = Data(:,z:z+1);
XY1(:,1) = XY1(:,1);
X = XY1(:,1);
Y = smooth(XY1(:,1),XY1(:,2),0.35,'rloess');
n=size(X,1);

Fig = figure;
set(Fig,'defaultLegendAutoUpdate','off');
plot(XY1(:,1),(XY1(:,2)),'g-')
h1 = plot(X,Y,'r-', 'LineWidth',3, 'DisplayName', 'Smoothed Profile');
hold on
set(gca, 'FontName', 'Times New Roman')
xlabel('\it Distance (\mu m)', 'fontweight','bold','fontsize',12)
ylabel('\it Height (\mu m)','fontweight','bold','fontsize',12)
[TF,S1,S2] = ischange(Y,'linear','MaxNumChanges',3);
[T1_X,T1_Y,T2_X, T2_Y]= TurningPoints(XY1,Y,S2,n);
Left_height = T1_Y;
Right_height = T2_Y;
Max_height = max(Y);
Min_height = min(Y);
Proceed = PROCEED(Max_height, Min_height,Left_height,Right_height);

if Proceed == 1
[f1,y_Ltan,f2,y_Rtan, L_slope, R_slope, xmin, y_xmin, y_xmax, x_ymax, ymax, y_mx, Pi, Qi, Si, Ti, xi, yi] = Geometry(X,Y,T1_X, T1_Y, T2_X, T2_Y);  
legend(h1)
Sisort = sort(Si);
Tisort = sort(Ti);
Pisort = sort(Pi);
Qisort = sort(Qi);
t1=isempty(Pi);
t2=isempty(Si);
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
if Proceed == 1
filename = strcat('Fig_',num2str(z),'.jpg');
saveas(Fig,filename);
else
end
% close(Fig);
hold off
else
end
end
savefig(Tayo)
%%Outputs
Output = Output(any(Output,2),:);
Output(end+1,:)= mean(Output(:,:));














Deviation = std(Output(:,7)); %Computes the deviation on the Edgeradius
Mean = mean(Output(:,7));
SixDev = (max(Output(:,7)-Mean) /3); %Computes the Six deviation on the Edgeradius
UppBound = Mean + (3*SixDev);
LowBound = Mean - (3*SixDev);

for i = 12:-1:1
    if (Output(i,7)>= LowBound) && (Output(i,7)<= UppBound)
    
    else
        Output(i,:)=[];
    end
end

        


