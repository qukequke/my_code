%----------------------------------------------------------------
% This code can be used to generate Figure 6-21 thru23
%----------------------------------------------------------------
% This file requires the following files to be present in the same
% directory:
%
% CoutUssFletcher.mat 

clear all
close all
clc
 
%---Radar parameters------------------------------------------------
pulses = 128;             % # no of pulses          
burst = 128;              % # no of bursts 
c = 3.0e8;                % speed of EM wave [m/s]
f0 = 9e9;                 % Starting frequency of SFR radar system [Hz]
bw = 125e6;               % Frequency bandwidth [Hz]
% bw = 250e6; 
T1 = (pulses-1)/bw;       % Pulse duration [s]
PRF = 35e3;               % Pulse repetition frequency [Hz]
T2 = 1/PRF;               % Pulse repetition interval [s]
 
%---target parameters------------------------------------------------
theta0 = 0;               % Initial angle of target's wrt target [degree]
% w = 1.2;                  % Angular velocity [degree/s] 
Vr = 30;                 % radial velocity of EM wave [m/s]
ar = 0.16;                % radial accelation of EM wave [m/s^2]
R0 = 16e3;                 % target's initial distance from radar [m]
dr = c/(2*bw);            % range resolution [m]  
% W = w*pi/180;             % Angular velocity [rad/s] 
W = 0.02;

%% 自己加 调整参数
% scatter_flag = false; %是否加散射点强度
scatter_flag = false;
gather_data = false;
m = 2000;  %收集数据个数
filename = '+_50_100_V10';
%% 
data = zeros(140, 140, 2, m);
label = zeros(140, 140, m);
%---load the coordinates of the scattering centers on the fighter------
load duo6
if gather_data == false
    m = 1
end
for iq = 1:m;
% Xc = round(randn(1,60)*15);
% Yc = round(randn(1,60)*15);

% Xc = round(unifrnd(-50, 50, 1, 100))
% Yc = round(unifrnd(-50, 50, 1, 100))

    new_mat = zeros(140, 220);
    scatter_power = abs(randn(1,length(Xc)));
    for iii = 1:1:length(Xc);
%         new_mat(Yc(iii)+70, Xc(iii)+110) = scatter_power(iii)
        if scatter_flag == true
            new_mat(Yc(iii)+70, Xc(iii)+110) = scatter_power(iii)
        else
            new_mat(Yc(iii)+70, Xc(iii)+110) = 1
        end
    end
%     figure;
%     imagesc(new_mat)
    new_mat = flipud(new_mat);
%      new_mat = fliplr(flipud(new_mat))
    new_mat = new_mat(:, 40:179);
%     imagesc(fliplr(flipud(new_mat)))
%     colormap gray;
%     figure();
%     plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0])
%     set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Normal');
%     axis([-60 60 -60 60])
%     xlabel('X [m]'); ylabel('Y [m]');
    
%     figure();
%     plot(Xc+90,Yc+70,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0])
%     set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Normal');
%     axis([0 160 0 140])
%     xlabel('X [m]'); ylabel('Y [m]');
    
    %---Scattering centers in cylindirical coordinates------------------------
    [theta,r] = cart2pol(Xc,Yc);
    % [theta,r] = cart2pol(x,y);

    Theta = theta+theta0*0.017455329; %add initial angle

    i = 1:pulses*burst;
    T = T1/2+2*R0/c+(i-1)*T2;%calculate time vector
    Rvr = Vr*T+(0.5*ar)*(T.^2);%Range Displacement due to radial vel. & acc. 
    Tetw = W*T;% Rotational Displacement due to angular vel. 

    i = 1:pulses;
    df = (i-1)*1/T1; % Frequency incrementation between pulses
    k = (4*pi*(f0+df))/c; 
    k_fac = ones(burst,1)*k; 

    %------Calculate backscattered E-field-------------------------------------
%             Es(burst,pulses)=0.0;
            Es = zeros(burst, pulses);
%             for scat=1:1:1;    
            for scat=1:1:length(Xc);  
                arg = (Tetw - theta(scat));
                rngterm = R0 + Rvr - r(scat)*sin(arg);
                range = reshape(rngterm,pulses,burst);
                range = range.';
                phase = k_fac.* range;
                if scatter_flag == true
                    Ess = scatter_power(scat)*exp(-j*phase);
                else
                    Ess = exp(-j*phase);
                end
%                 Ess = scatter_power(scat)*exp(j*phase);
%                 Ess = exp(j*phase);
                Es = Es+Ess;
            end
            Es = Es.';

            
%% ---这里是最小熵算法-------------------------------------
% JTF Representation of range cell  
EsMp = reshape(Es,1,pulses*burst);
S = spectrogram(EsMp,128,64,120);
[a,b] = size(S);

% h = figure;
% matplot2((1:a),(1:b),abs(S),50);
% set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
% colormap(1-gray);
% xlabel('time pulses');
% ylabel('frequency index');
% title('Spectrogram');

%Prepare time and frequency vectors
f = (f0+df);% frequency vector
T = reshape(T,pulses,burst); %prepare time matrix
F = f.'*ones(1,burst); %prepare frequency matrix
 
% Searching the motion parameters via min. entropy method  
syc=1;
V = 0:0.2:30;
A = 0:0.02:2;
m = 0; 
for Vest = V;
    m = m+1;
    n = 0;
    for iv = A;
        n = n+1;
        VI(syc,1:2) = [Vest,iv];
        S = exp((j*4*pi*F/c).*(Vest*T+(0.5*iv)*(T.^2)));
        Scheck = Es.*S;
        ISAR = abs(fftshift(fft2((Scheck))));
        SumU = sum(sum(ISAR));
        I = (ISAR/SumU);
        Emat = I.*log10(I);
        EI(m,n) = -(sum(sum(Emat)));
        syc = syc+1;
    end    
end
 
[dummy,mm] = min(min(EI.')); %Find index for estimated velocity
[dummy,nn] = min(min(EI));   %Find index for estimated acceleration
%---Figure 8_10 ---------------------------------------------
h =surfc(A,V,EI);
colormap(gray)
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight', 'Bold'); 
ylabel('Translational velocity [m/s]'); 
xlabel('Translational acceleration [m/s^2]');
zlabel ('Entropy value')
 
% Form the mathing phase for compensation
Sconj = exp((j*4*pi*F/c).*(V(mm)*T+(0.5*A(nn)*(T.^2))));
% Compansate
S_Duz = Es.*Sconj;

ISAR3 = fftshift(fft2((S_Duz), 140, 220));
    if Vr == 10
        ISAR3=circshift(ISAR3,[0 -50]); %v = 10的时候 用这个平移一下
    end;
    ISAR3 = ISAR3/255;
    ISAR3 = ISAR3(:, 40:179);
    ISAR3 = flipud(ISAR3)
ISAR3 = circshift(ISAR3, [0, -28]);
figure;
colormap gray;
imagesc(abs(ISAR3))
title('最小熵结果');
%% ---------------OSS --------------------------------------------
d=rand(128,128);
figure
colormap(gray(256))
% imagesc(255-d);
xlabel('距离维单元数'); 
ylabel('方位维单元数');
% f=load('Es_v5_a0.06.mat');%zhendong1
f=Es;

an=angle(f);
f=abs(f);
g=20*log10(f);
figure(100)
imagesc(f);
xlabel('距离维单元数'); 
ylabel('方位维单元数');
f=circshift(f,[3 19]);%有噪声时选取[67 42]
g=20*log10(f);
% figure(99)
% imagesc(f);
% xlabel('距离维单元数'); 
% ylabel('方位维单元数');
% angle=angle.angle;
% angle=fftshift(angle);
% angle=circshift(angle,[0 0]);
% d=ones(564,850);
% an=load('00.mat');
% an=an.an;
%  an=circshift(an,[-13 -43]);
supp(1)=60;
supp(2)=60;
for iter=180:10:180
    for beta=0.7:0.01:0.7
      [mask RfacR RfacF  RFD]= OSS_samplecode3(f,supp,iter, beta,0,d,0,an); 
% [mask RfacR RfacF  RFD]= OSS_samplecode3(f,supp,iter, beta,1,d,0,an); 
% [S,dbs]5 = USample2D(RFD,100,16);                      %升采样并显示二维切面
% save('OSS.mat','RFD')
 figure;
%  G=20*log10(abs(RFD.')+1e-6);
%  gm=max(max(G));
%  gn=gm-14;%显示动态范围40dB
%  G=255/(gm-gn)*(G-gn).*(G>gn);
%  imagesc(G)
% colormap(gray)
% imagesc(255-RFD.')
 G=20*log10(abs(RFD)+1e-6);
 gm=max(max(G));
 gn=gm-14;%显示动态范围40dB
 G=255/(gm-gn)*(G-gn).*(G>gn);
 colormap gray;
 imagesc(G)
 title('OSS');
%  title([iter, beta,supp(1),supp(2)]);
%   imagesc(255-G
% colormap(gray)
% figure
%         mesh(double(RFD));
%         view(45,50);
%         axis tight
 toc;
% figure
% imagesc(abs(fftshift(fft2(G))));
  
    end
end





%% 

    %---Figure 6.23---------------------------------------------------
    ISAR = fftshift(fft2((Es), 140, 220));
    if Vr == 10
        ISAR=circshift(ISAR,[0 -50]); %v = 10的时候 用这个平移一下
    end;
    ISAR = ISAR/255;
    ISAR = ISAR(:, 40:179);
    ISAR = flipud(ISAR)
data(:, :, 1, iq) = real(ISAR);
data(:, :, 2, iq) = imag(ISAR);

label(:, :, iq) = new_mat;
end

figure;
%     imagesc(new_mat)
% new_mat = new_mat(:, 40:179);
% new_mat = fliplr(flipud(new_mat))
imagesc(new_mat)
colormap gray;
figure();
plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor',[1 0 0])
set(gca,'FontName', 'Arial', 'FontSize',14,'FontWeight', 'Normal');
axis([-60 60 -60 60])
xlabel('X [m]'); ylabel('Y [m]');
figure;
imagesc(abs(ISAR));
colormap gray;

if gather_data == true
    save(filename, 'data', 'label');
end

%     figure(jj+2);
%     
%     colormap(gray)
%     imagesc(1-ISAR);
%     axis xy;
%     data(:, :,jj)=1-(abs(fftshift(fft2((Es), 120, 120))))/255;
%     label(:, :,jj)=1-dddd
% % end
% save('mat_save/data_train', 'data');
% save('mat_save/label_train', 'label');
