% Five_level_diff_split
%Calculate Five level different splitting, ACZeeman Energy with different detuning

clear;
close all;

freqrange = [18:.001:22]*1e6;
B=20/.6998;
% freqrange = [1:.001:3]*1e6;
% B=2/.6998;


%SI Unit
% h_bar=2*pi;
h= 6.62606896*10^-34;            % joules? 
h_bar = h/(2*pi);            % joules? 
m_87 = 1.44316060*10^-25; % kg
g = 9.8 ; %m/s/s



[E_J, Eham_Hz, CGg]  = eight_E_cg_mFham_func( B );
% Eham_Hz = Eham_Hz;
% F=2
w0 = Eham_Hz(5)-Eham_Hz(4);             %E(2,-1)-E(2,-2)
delta0 = Eham_Hz(6)-Eham_Hz(5)-w0;       %E(2,0)-E(2,-1)-basicsplitting
delta1 = Eham_Hz(7)-Eham_Hz(5)-2*w0;       %E(2,1)-E(2,-1)-basicsplitting
delta2 = Eham_Hz(8)-Eham_Hz(5)-3*w0;       %E(2,2)-E(2,-1)-basicsplitting
% [basicsplitting, delta3,delta4,delta5] = Breit_Rabi_Five_Level( B );
% F=1
w1 = Eham_Hz(2)-Eham_Hz(1);             %E(1,-1)-E(1,0)
delta11 = Eham_Hz(3)-Eham_Hz(2)-w1;       %E(1,1)-E(1,0)-basicsplitting

%% old calculation: brait-Rabi
% B = B/1e4;                      %T
% [Ebr, Ebr_Hz] = Breit_Rabi_Eight_Level_array( B );
% Ebr_Hz = squeeze(Ebr_Hz); % we calculate in frequencies later
% % F=2
% w0 = Ebr_Hz(2,2)-Ebr_Hz(2,1);             %E(2,-1)-E(2,-2)
% delta0 = Ebr_Hz(2,3)-Ebr_Hz(2,2)-w0;       %E(2,0)-E(2,-1)-basicsplitting
% delta1 = Ebr_Hz(2,4)-Ebr_Hz(2,2)-2*w0;       %E(2,1)-E(2,-1)-basicsplitting
% delta2 = Ebr_Hz(2,5)-Ebr_Hz(2,2)-3*w0;       %E(2,2)-E(2,-1)-basicsplitting
% % [basicsplitting, delta3,delta4,delta5] = Breit_Rabi_Five_Level( B );
% % F=1
% w1 = Ebr_Hz(1,3)-Ebr_Hz(1,4);             %E(1,-1)-E(1,0)
% delta11 = Ebr_Hz(1,2)-Ebr_Hz(1,3)-w1;       %E(1,1)-E(1,0)-basicsplitting

% resonances = [w0+delta2-delta1, w0+delta1-delta0, w0+delta0, w0]
             % |2,+2> - |2,+1>                  |2,+1> - |2,0>             |2,0> - |2,-1>        |2,-1> - |2,-2> 
             
%% calculations %%
% dRabidR = -Rabi0/ R;%  dRabidR = -dRabi / dR;
% dRabiFactor = .12; % like charlie
% dRabiFactor*Rabi0;
% newRabi = Rabi0-dRabi;
Rabi0 = 5e5;        % power -> 1/2 MHz rabi strength
Rabitiny = 1e-12;
dRabi = 5e3 ;   % 1 khz / unit
dR = .000001;         % um / unit
newRabi = Rabi0-dRabi; 

% F=2
Rabiab = Rabi0*CGg(5,4);          %Hz
Rabibc = Rabi0*CGg(3,4);          %Hz
Rabicd = Rabi0*CGg(3,2);          %Hz
Rabide = Rabi0*CGg(1,2);                    %Hz

tRabiab = Rabitiny*Rabiab;          %Hz
tRabibc = Rabitiny*Rabibc;          %Hz
tRabicd = Rabitiny*Rabicd;          %Hz
tRabide = Rabitiny*Rabide;                    %Hz

dRabiab = newRabi*CGg(5,4);          %Hz
dRabibc = newRabi*CGg(3,4);          %Hz
dRabicd = newRabi*CGg(3,2);          %Hz
dRabide = newRabi*CGg(1,2);                    %Hz

% F=1
Rabifg = Rabi0*CGg(6,7);          %Hz
Rabigh = Rabi0*CGg(7,8);          %Hz

tRabifg = Rabitiny*Rabifg ;          %Hz
tRabigh = Rabitiny*Rabigh ;          %Hz

dRabifg = newRabi*CGg(6,7);          %Hz
dRabigh = newRabi*CGg(7,8);          %Hz

%% preallocate
numd = length(freqrange);
E2 = zeros(5,numd); E2low = zeros(5,numd); F2 = zeros(5,numd);x2 = zeros(1,numd); 
E1 = zeros(3,numd); E1low = zeros(3,numd); F1 = zeros(3,numd);x1 = zeros(1,numd);
V2 = zeros(5,5,numd); V1= zeros(3,3,numd); 

%% asd 

for countf=1:length(freqrange)
    wrf = freqrange(countf);
    
H2tot = [ 2*(wrf-w0)    Rabiab/2        0           0                       0;...
      Rabiab/2          (wrf-w0)        Rabibc/2    0                           0;...
      0                 Rabibc/2        delta0      Rabicd/2            0 ;...
      0                 0               Rabicd/2    -(wrf-w0)+delta1    Rabide/2;...
      0                 0               0           Rabide/2            -2*(wrf-w0)+delta2
]*h;

H2low = [ 2*(wrf-w0)   tRabiab/2           0           0                   0;...
    tRabiab/2          (wrf-w0)        tRabibc/2   0                   0;...
    0                  tRabibc/2      delta0      tRabicd/2           0 ;...
    0                  0              tRabicd/2   -(wrf-w0)+delta1    tRabide/2;...
    0                  0              0           tRabide/2           -2*(wrf-w0)+delta2
]*h;

H2d = [ 2*(wrf-w0)   dRabiab/2                0           0                           0;...
      dRabiab/2             (wrf-w0)        dRabibc/2    0                           0;...
      0                     dRabibc/2               delta0      dRabicd/2                    0;...
      0                     0                       dRabicd/2   -(wrf-w0)+delta1     dRabide/2;...
      0                     0                       0           dRabide/2                    -2*(wrf-w0)+delta2
]*h;

H1tot = [   (wrf-w1)        Rabifg/2            0   ;...
            Rabifg/2                0                   Rabigh/2;...
            0                       Rabigh/2            delta11-(wrf-w1)  ;...
]*h;

H1low = [   (wrf-w1)         tRabifg/2           0   ;...
            tRabifg/2                0                   tRabigh/2;...
            0                       tRabigh/2            delta11-(wrf-w1)  ;...
]*h;

H1d = [   (wrf-w1)          dRabifg/2            0   ;...
            dRabifg/2                0                   dRabigh/2;...
            0                       dRabigh/2            delta11-(wrf-w1)  ;...
]*h;


[V1tot,D1tot] = eig(H1tot);
[~,D1low] = eig(H1low);
[~,D1d] = eig(H1d);

[V2tot,D2tot] = eig(H2tot);
[~,D2low] = eig(H2low);
[~,D2d] = eig(H2d);


E2(:,countf) = diag(D2tot)/h/1e6;
E2low(:,countf) = diag(D2low)/h/1e6;
F2(:,countf) = -diag(D2tot-D2d)/dR/m_87/g;
V2(:,:,countf) = V2tot; 

E1(:,countf) = diag(D1tot)/h/1e6;
E1low(:,countf) = diag(D1low)/h/1e6;
F1(:,countf) = -diag(D1tot-D1d)/dR/m_87/g;
V1(:,:,countf) = V1tot; 

% x2(countf)=((wrf-w1)+w0)*10^-6;
% x1(countf)=((wrf-w1)+w1)*10^-6;
% x1(countf)= wrf;
% x2(countf)= wrf;
end
x1= freqrange/1e6;
x2= freqrange/1e6;


%% plotting %%
% color2 = ['rmgcb'];
% % color1 = ['mgc'];
% colors = [[0,.75,.75];...%c
%           [0,1,0];...%g
%           [1,0,1];...%m
%           [1,0,0];...%r
%           [1,0,1];...%m
%           [0,1,0];...%g
%           [0,.75,.75];...%c
%           [0,0,1]];    %b

% F=2
% figure('position', [560   260   829   688]);  % energy large rabi
% eight_mFHamBRplotter

%% constants
h    = 6.62606896e-34;  %Js
hbar = h/2/pi;          %Js
Ahfs = h*3.417341305452145e9;% 5s1/2
% muB  = 9.27400915e-24      %J/T
muB  = h*1.399624604e6;      %J / G
gs =  2.0023193043622;
gi = -0.0009951414;
% gl =  0.99999369;
% gj =  2.00233113;

txsz = 12;

%% input parameters
B = [0:.01:41];   % Teslas
% B = Bingauss;   % this calculation is in gauss. 

%% states
magI = 3/2; vecI = magI:-1:-magI; densq = (2*magI+1);
magS = 1/2; vecS = magS:-1:-magS;
magF(2) = magI+magS;
magF(1) = magI-magS;

%% Hyperfine 
IdS(2) = .5*(magF(2)*(magF(2)+1) -magI*(magI+1) -magS*(magS+1));
IdS(1) = .5*(magF(1)*(magF(1)+1) -magI*(magI+1) -magS*(magS+1));
Hhfs(2) = Ahfs * IdS(2);
Hhfs(1) = Ahfs * IdS(1);
% (Hhfs2-Hhfs1)/h;      % 6.834...GHz
Hhfs = diag([Hhfs(2),Hhfs(2),Hhfs(2),Hhfs(2),Hhfs(2), Hhfs(1),Hhfs(1),Hhfs(1)]);
% Hhfs = [Hhfs2,Hhfs2,Hhfs2,Hhfs2,Hhfs2, Hhfs1,Hhfs1,Hhfs1]';


%% magnetic
% construct B in the S, I basis
% col 1 = i, col 2 = s 
SIb = [[vecI,vecI];[magS,magS,magS,magS],-[magS,magS,magS,magS]]';
Hisd = zeros(size(SIb,1),1);
for s=1:1:size(SIb,1)
    Hisd(s) = gi*SIb(s,1) + gs*SIb(s,2);
end

%% Spin matrices
%just organizing thoughts here, geting locations:
% bra state (trans ) ket state : +=- sig + / pi / sig -
            % bra state /
%           %   +2      +1        0      -1      -2 |    +1      0     -1   % /ket state
% CGcoef=[[      0,  21+22,       0,      0,      0,  11+22,      0,      0];...%+2
%         [  22+21,      0,   20+21,      0,      0,  11=21,  10+21,      0];...%+1
%         [      0,  21+20,       0, 2m1+20,      0,  11-20,  10=20, 1m1+20];...%0
%         [      0,      0,  20+2m1,      0,2m2+2m1,      0, 10-2m1,1m1=2m1];...%-1
%         [      0,      0,       0,2m1+2m2,      0,      0,      0,1m1-2m2];...%-2------
%         [  22+11,  21=11,   20-11,      0,      0,      0,  10-11,      0];...%+1
%         [      0,  21+10,   20=10, 2m1-10,      0,  11-10,      0, 1m1-10];...%0  
%         [      0,      0,  20+1m1,2m1=1m1,2m2-1m1,      0, 10-1m1,      0]   ];%-1

% sigma plus : flip S spin only down ket to up bra
Sp = [  sqrt((2+ 2)*(2- 1));...         % 2 2<->2 1 RF
        sqrt((2+ 1)*(2- 0));...         % 2 1<->2 0 RF
        sqrt((2+ 0)*(2--1));...         % 2 0<->2-1 RF
        sqrt((2+-1)*(2--2));...         % 2-1<->2-2 RF
        sqrt((2+ 2)*(2+ 1));...         % 2 2<->1 1 MW
        sqrt((2+ 1)*(2+ 0));...         % 2 1<->2 0 MW
        sqrt((2+ 0)*(2+-1))  ]./densq./2;  % 2 0<->2-1 MW
                                    % half from +- def
          % +2  +1      0    -1    -2   +1     0    -1
SpM = [ [    0,Sp(1),    0,    0,    0,Sp(5),    0,    0];...%+2
        [Sp(1),    0,Sp(2),    0,    0,    0,Sp(6),    0];...%+1
        [    0,Sp(2),    0,Sp(3),    0,    0,    0,Sp(7)];...%0
        [    0,    0,Sp(3),    0,Sp(4),    0,    0,    0];...%-1
        [    0,    0,    0,Sp(4),    0,    0,    0,    0];...%-2
        [Sp(5),    0,    0,    0,    0,    0,    0,    0];...%+1
        [    0,Sp(6),    0,    0,    0,    0,    0,    0];...%0
        [    0,    0,Sp(7),    0,    0,    0,    0,    0]   ];%-1
% sigma minus : flip S spin only up ket to down bra
Sm = [  -sqrt((2-+1)*(2+ 0));...         % 1 0<->1 1 RF
        -sqrt((2-+0)*(2-1));...         % 1-1<->1 0 RF
        -sqrt((2--2)*(2--1));...         % 2-2<->1-1 MW
        -sqrt((2--1)*(2- 0));...         % 2-1<->1 0 MW
        -sqrt((2- 0)*(2- 1))  ]./densq./2;  % 2 0<->1 1 MW
                                    % half from +- def
          % +2    +1     0    -1    -2 |  +1     0    -1
SmM = [ [    0,    0,    0,    0,    0,    0,    0,    0];...%+2
        [    0,    0,    0,    0,    0,    0,    0,    0];...%+1
        [    0,    0,    0,    0,    0,Sm(5),    0,    0];...%0
        [    0,    0,    0,    0,    0,    0,Sm(4),    0];...%-1
        [    0,    0,    0,    0,    0,    0,    0,Sm(3)];...%-2------
        [    0,    0,Sm(5),    0,    0,    0,Sm(1),    0];...%+1
        [    0,    0,    0,Sm(4),    0,Sm(1),    0,Sm(2)];...%0
        [    0,    0,    0,    0,Sm(3),    0,Sm(2),    0]   ];%-1
% pi transition: no flip S, but get m_s factor from op. on kets
% basically two terms per term
       % - spin up half      	  spin down half       
Sz = [  -sqrt((2 +1)*(2- 1)) + -sqrt((2-+1)*(2+ 1)) ;...            % 2 1<->1 1 
        -sqrt((2 +0)*(2- 0)) + -sqrt((2- 0)*(2+ 0)) ;...            % 2 0<->1 0 
        -sqrt((2+-1)*(2--1)) + -sqrt((2--1)*(2+-1))   ]./2./densq;  % 2-1<->1-1 
                            %  - from spin down       1/2 from m_s
          % +2    +1     0    -1    -2 |  +1     0    -1
SzM = [ [    0,    0,    0,    0,    0,    0,    0,    0];...%+2
        [    0,    0,    0,    0,    0,Sz(1),    0,    0];...%+1
        [    0,    0,    0,    0,    0,    0,Sz(2),    0];...%0
        [    0,    0,    0,    0,    0,    0,    0,Sz(3)];...%-1
        [    0,    0,    0,    0,    0,    0,    0,    0];...%-2------
        [    0,Sz(1),    0,    0,    0,    0,    0,    0];...%+1
        [    0,    0,Sz(2),    0,    0,    0,    0,    0];...%0
        [    0,    0,    0,Sz(3),    0,    0,    0,    0]   ];%-1
    
% sigma plus : flip I spin only -1/2 ket to +1/2 bra
            % spin ups          spin downs
Ip = [  sqrt((2+ 2)*(2+ 1))+ sqrt((2-+2)*(2- 1));...            % 2 2<->2 1 RF
        sqrt((2+ 1)*(2+ 0))+ sqrt((2- 1)*(2- 0));...            % 2 1<->2 0 RF
        sqrt((2+ 0)*(2+-1))+ sqrt((2- 0)*(2--1));...            % 2 0<->2-1 RF
        sqrt((2+-1)*(2+-2))+ sqrt((2--1)*(2--2));...            % 2-1<->2-2 RF
       -sqrt((2+ 2)*(2- 1))+ sqrt((2- 2)*(2+ 1));...            % 2 2<->1 1 MW
       -sqrt((2+ 1)*(2- 0))+ sqrt((2- 1)*(2+ 0));...            % 2 1<->1 0 MW
       -sqrt((2+ 0)*(2--1))+ sqrt((2- 0)*(2+-1))  ]./densq./2;  % 2 0<->1-1 MW
                                                    % half from +- def
          % +2  +1      0    -1    -2   +1     0    -1
IpM = [ [    0,Ip(1),    0,    0,    0,Ip(5),    0,    0];...%+2
        [Ip(1),    0,Ip(2),    0,    0,    0,Ip(6),    0];...%+1
        [    0,Ip(2),    0,Ip(3),    0,    0,    0,Ip(7)];...%0
        [    0,    0,Ip(3),    0,Ip(4),    0,    0,    0];...%-1
        [    0,    0,    0,Ip(4),    0,    0,    0,    0];...%-2
        [Ip(5),    0,    0,    0,    0,    0,    0,    0];...%+1
        [    0,Ip(6),    0,    0,    0,    0,    0,    0];...%0
        [    0,    0,Ip(7),    0,    0,    0,    0,    0]   ];%-1
% sigma plus : flip I spin only -1/2 ket to +1/2 bra
            % spin ups          spin downs
Im = [ -sqrt((2+ 0)*(2- 1))+ sqrt((2- 0)*(2+ 0));...            % 1 0<->1 1 RF
       -sqrt((2+-1)*(2- 0))+ sqrt((2--1)*(2+ 0));...            % 1-1<->1 0 RF
       -sqrt((2+-2)*(2--1))+ sqrt((2--2)*(2+-1));...            % 2-2<->1-1 MW
       -sqrt((2+-1)*(2- 0))+ sqrt((2--1)*(2+ 0));...            % 2-1<->1 0 MW
       -sqrt((2+ 0)*(2-+1))+ sqrt((2- 0)*(2+ 1))  ]./densq./2;  % 2 0<->1 1 MW
                                                    % half from +- def
          % +2    +1     0    -1    -2 |  +1     0    -1
ImM = [ [    0,    0,    0,    0,    0,    0,    0,    0];...%+2
        [    0,    0,    0,    0,    0,    0,    0,    0];...%+1
        [    0,    0,    0,    0,    0,Im(5),    0,    0];...%0
        [    0,    0,    0,    0,    0,    0,Im(4),    0];...%-1
        [    0,    0,    0,    0,    0,    0,    0,Im(3)];...%-2------
        [    0,    0,Im(5),    0,    0,    0,Im(1),    0];...%+1
        [    0,    0,    0,Im(4),    0,Im(1),    0,Im(2)];...%0
        [    0,    0,    0,    0,Im(3),    0,Im(2),    0]   ];%-1
% pi transition: no flip S, but get m_s factor from op. on kets
% basically two terms per term
            % spin up half    m_I    spin down half       +-1/2 from Sz op
Iz = [  -sqrt((2+ 1)*(2- 1))*( 1/2) + sqrt((2-+1)*(2+ 1))*( 3/2) ;...       % 2 1<->1 1 
        -sqrt((2+ 0)*(2- 0))*(-1/2) + sqrt((2- 0)*(2+ 0))*( 1/2) ;...       % 2 0<->1 0 
        -sqrt((2+-1)*(2--1))*(-3/2) + sqrt((2--1)*(2+-1))*(-1/2) ]./densq;  % 2-1<->1-1 
          % +2    +1     0    -1    -2 |  +1     0    -1
IzM = [ [    0,    0,    0,    0,    0,    0,    0,    0];...%+2
        [    0,    0,    0,    0,    0,Iz(1),    0,    0];...%+1
        [    0,    0,    0,    0,    0,    0,Iz(2),    0];...%0
        [    0,    0,    0,    0,    0,    0,    0,Iz(3)];...%-1
        [    0,    0,    0,    0,    0,    0,    0,    0];...%-2------
        [    0,Iz(1),    0,    0,    0,    0,    0,    0];...%+1
        [    0,    0,Iz(2),    0,    0,    0,    0,    0];...%0
        [    0,    0,    0,Iz(3),    0,    0,    0,    0]   ];%-1

%% H setup
% HB = (muB/hbar)*(gs*ms + gi*mi)*B
HisB = zeros(size(SIb,1),size(B,2));
Hfdiag = zeros(8,8,size(B,2));
Hoffup = zeros(8,8,size(B,2));
Hoffdn = zeros(8,8,size(B,2));
Hftot = zeros(8,8,size(B,2)); %#ok<PREALL>
H_MHz = zeros(8,size(B,2));
E_Hz = zeros(8,size(B,2));
E_Joules = zeros(8,size(B,2));
V= zeros(8,8,size(B,2));
Vr= zeros(8,8,size(B,2));
D= zeros(8,8,size(B,2));
W= zeros(8,8,size(B,2));
CGs = zeros(8,8,size(B,2));
CGi = zeros(8,8,size(B,2));
HisB(:,:) = muB*repmat(Hisd,1,size(B,2)).*repmat(B,8,1);

% self-B terms from I-S, in "F" basis
Hfdiag(1,1,:) = (2+2)/densq*HisB(1,:) + 0;
Hfdiag(2,2,:) = (2+1)/densq*HisB(2,:) + (2-1)/densq*HisB(5,:);
Hfdiag(3,3,:) = (2+0)/densq*HisB(3,:) + (2-0)/densq*HisB(6,:);
Hfdiag(4,4,:) = (2-1)/densq*HisB(4,:) + (2+1)/densq*HisB(7,:);
Hfdiag(5,5,:) = 0                     + (2+2)/densq*HisB(8,:);
Hfdiag(6,6,:) = (2-1)/densq*HisB(2,:) + (2+1)/densq*HisB(5,:);
Hfdiag(7,7,:) = (2+0)/densq*HisB(3,:) + (2+0)/densq*HisB(6,:);
Hfdiag(8,8,:) = (2+1)/densq*HisB(4,:) + (2-1)/densq*HisB(7,:);

%off- diagonal terms mixing F=1/F=2 from spin up : only +1/0/-1
Hoffup(2,6,:)=-sqrt((2+1)*(2-1))/densq*HisB(2,:);
Hoffup(3,7,:)=-sqrt((2+0)*(2-0))/densq*HisB(3,:);
Hoffup(4,8,:)=-sqrt((2-1)*(2+1))/densq*HisB(4,:);
Hoffup(6,2,:)=-sqrt((2+1)*(2-1))/densq*HisB(2,:);
Hoffup(7,3,:)=-sqrt((2+0)*(2-0))/densq*HisB(3,:);
Hoffup(8,4,:)=-sqrt((2-1)*(2+1))/densq*HisB(4,:);

%off- diagonal terms mixing F=1/F=2 from spin down : only +1/0/-1
Hoffdn(2,6,:)=sqrt((2-1)*(2+1))/densq*HisB(5,:);
Hoffdn(3,7,:)=sqrt((2+0)*(2-0))/densq*HisB(6,:);
Hoffdn(4,8,:)=sqrt((2+1)*(2-1))/densq*HisB(7,:);
Hoffdn(6,2,:)=sqrt((2-1)*(2+1))/densq*HisB(5,:);
Hoffdn(7,3,:)=sqrt((2+0)*(2-0))/densq*HisB(6,:);
Hoffdn(8,4,:)=sqrt((2+1)*(2-1))/densq*HisB(7,:);

Hftot = Hfdiag + Hoffup + Hoffdn + repmat(Hhfs,1,1,size(B,2));

% B- eigenvalue loop
for b = 1:1:size(B,2)
[~,D(:,:,b)] = eig(Hftot(:,:,b));        % eig can only work on 2-D matrices
H_MHz(:,b) = diag(D(:,:,b))/h/1e6;
E_Hz(:,b) = diag(D(:,:,b))/h;
E_Joules(:,b) = diag(D(:,:,b));
end

%% plotting
figure('position', [200 100 1300 600], 'renderer', 'painter');
tsz = 16;
topup = .02;
eachx = .175;
eachy = .39;
buffx = (1-4*eachx)/5;
buffy = (1-2*eachy)/3;
state2 = {'|+ +\rangle';'| + \rangle';'| 0 \rangle';'| - \rangle';'|- -\rangle';};
state1 = {'| +'' \rangle';'| 0'' \rangle';'| -'' \rangle';};
bare2 = {'|-2,+2\omega\rangle';'|-1,+\omega\rangle';'|0,+0\omega\rangle';'|+1,-\omega\rangle';'|+2,-2\omega\rangle';};
bare1 = {'|+1,+\omega\rangle';'|0,+0\omega\rangle';'|-1,-\omega\rangle'};

xoff = .01;


states = {'|1,+1\rangle';'|1, 0\rangle';'|1,-1\rangle';'|2,-2\rangle';'|2,-1\rangle';'|2, 0\rangle';'|2,+1\rangle';'|2,+2\rangle';};
% colors = 'cgmrmgcb'; 
colors = [[0,.75,1];...%c
          [0,.875,0];...%g
          [1,0,1];...%m
          [1,0,0];...%r
          [1,0,1];...%m
          [0,.875,0];...%g
          [0,.75,1];...%c
          [0,0,1]];    %b
% g(1)=subplot(2,4,1);hold on;
g(1)=axes('position', [1*buffx+0*eachx,2*buffy+1*eachy+topup ,eachx,eachy] );hold on;
for s = 4:1:6
plot(B,H_MHz(s,:)-H_MHz(6,1),  'color', colors(s,:), 'linewidth', 1); 
text(B(end),(H_MHz(s,end)-H_MHz(6,1))/1.2, states{s}, 'color', colors(s,:),'verticalalignment', 'bottom','horizontalalignment', 'right', 'fontsize', txsz );
end
for s = 7:1:8
plot(B,H_MHz(s,:)-H_MHz(6,1),  'color', colors(s,:), 'linewidth', 1); 
text(B(end),(H_MHz(s,end)-H_MHz(6,1))/1.2, states{s}, 'color', colors(s,:),'verticalalignment', 'top','horizontalalignment', 'right', 'fontsize', txsz );
end
plot(20/(1.399624/2)*[1,1],[min(min(H_MHz(4,:)-H_MHz(6,1))),max(max(H_MHz(8,:)-H_MHz(6,1)))], 'k--');
% plot(8.5/(1.399624/2)*[1,1],[min(min(H_MHz(4,:)-H_MHz(6,1))),max(max(H_MHz(8,:)-H_MHz(6,1)))], 'k--');
title('F=2 dc Zeeman', 'fontsize', tsz );ylabel('\Delta{}E_{DCZ} /h, MHz');
grid on; axis tight; box on;%xlabel('Gauss');
xticks([0:10:50]);yticks([-40:20:40]);
text(1,-47, '(a2)', 'color', 'k', 'fontsize', 20 );

% g(5)=subplot(2,4,5);hold on;
g(5)=axes('position', [1*buffx+0*eachx,1*buffy+0*eachy,eachx,eachy] );hold on;
for s = 1
    plot(B,H_MHz(s,:)-H_MHz(2,1),  'color', colors(s,:), 'linewidth', 1); 
    text(B(end),H_MHz(s,end)-H_MHz(2,1), states{s}, 'color', colors(s,:),'verticalalignment', 'top','horizontalalignment', 'right', 'fontsize', txsz );
end
for s = 2:1:3
    plot(B,H_MHz(s,:)-H_MHz(2,1),  'color', colors(s,:), 'linewidth', 1); 
    text(B(end),H_MHz(s,end)-H_MHz(2,1), states{s}, 'color', colors(s,:),'verticalalignment', 'bottom','horizontalalignment', 'right', 'fontsize', txsz );
end
plot(20/(1.399624/2)*[1,1],[min(min(H_MHz(4,:)-H_MHz(6,1))),max(max(H_MHz(8,:)-H_MHz(6,1)))], 'k--');
% plot(8.5/(1.399624/2)*[1,1],[min(min(H_MHz(4,:)-H_MHz(6,1))),max(max(H_MHz(8,:)-H_MHz(6,1)))], 'k--');
title('F=1 dc Zeeman', 'fontsize', tsz );ylabel('\Delta{}E_{DCZ} /h, MHz'); xlabel('Gauss');
grid on;axis tight;box on;
xticks([0:10:50]);yticks([-40:20:40]);
text(1,49, '(a1)', 'color', 'k', 'fontsize', 20 );

% linkaxes([g], 'y');

smallxticks = [19.5:.25:20.5];
bigxticks = [18:.5:22];
smallyticks = [-1.5:.5:1];
diff21 = diff(Eham_Hz); off21 = (diff21(6)-diff21(7))/1e6;

right = find(x2==max(smallxticks));

% g(2)=subplot(2,4,2);hold on; 
g(2)=axes('position', [2*buffx+1*eachx,2*buffy+1*eachy+topup ,eachx,eachy] );hold on;
for s=1:1:3
    plot(x2,E2low(s,:)+off21, '-', 'color', .6*ones(1,3)); 
    text(xoff+x2(right),E2low(s,right)+off21, bare2{6-s}, 'color', .3*ones(1,3),'verticalalignment', 'middle','horizontalalignment', 'left', 'fontsize', txsz );
end
for s=1:1:3
    plot(x2,E2(s,:)+off21,'color', colors(3+s,:), 'linewidth', 1); 
end

for s=4:1:5
    plot(x2,E2low(s,:)+off21, 'color', .6*ones(1,3)); 
    text(xoff+x2(right),E2low(s,right)+off21, bare2{6-s}, 'color', .3*ones(1,3),'verticalalignment', 'middle','horizontalalignment', 'left','fontsize', txsz );
end
for s=4:1:5
    plot(x2,E2(s,:)+off21, 'color',colors(3+s,:), 'linewidth', 1); 
end
label2pos = [2010,2008,2253,2210,2010];
text(x2(label2pos(1)),(E2(5,label2pos(1))+off21)+.05, state2{1}, 'color', colors(8,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );
text(x2(label2pos(2)),(E2(4,label2pos(2))+off21)+.02, state2{2}, 'color', colors(7,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );
text(x2(label2pos(3)),(E2(3,label2pos(3))+off21), state2{3}, 'color', colors(6,:),'verticalalignment', 'top','horizontalalignment', 'left', 'fontsize', txsz );
text(x2(label2pos(4)),(E2(2,label2pos(4))+off21)+.01, state2{4}, 'color', colors(5,:),'verticalalignment', 'top','horizontalalignment', 'center', 'fontsize', txsz );
text(x2(label2pos(5)),(E2(1,label2pos(5))+off21)-.02, state2{5}, 'color', colors(4,:),'verticalalignment', 'top','horizontalalignment', 'center', 'fontsize', txsz );


hold off; grid on;box on;
% xlabel('Applied RF (MHz)');
ylabel('E_{ACZ} /h, MHz');
title('Eigenstate Energies F=2', 'fontsize', tsz );
xticks(smallxticks );yticks(smallyticks );
% ylim([-1.5,1]); 
xlim([19.5,20.5]);
text(20.03,-1.25, '(b2)', 'color', 'k', 'fontsize', 20 );


% g(3) = subplot(2,4,3);hold on ;
g(3)=axes('position', [3*buffx+2*eachx,2*buffy+1*eachy+topup ,eachx,eachy] );hold on;
for s=1:1:5
    plot(x2,(E2(s,:)-E2low(s,:))*1e3, 'color', colors(3+s,:), 'linewidth', 1); 
end
label2pos = [2000,2080,2053,2115,2170];
text(x2(label2pos(1)),(E2(1,label2pos(1))-E2low(1,label2pos(1)))*1e3-20, state2{5}, 'color', colors(4,:),'verticalalignment', 'top','horizontalalignment', 'center', 'fontsize', txsz );
text(x2(label2pos(2)),(E2(2,label2pos(2))-E2low(2,label2pos(2)))*1e3, state2{4}, 'color', colors(5,:),'verticalalignment', 'top','horizontalalignment', 'left', 'fontsize', txsz );
text(x2(label2pos(3)),(E2(3,label2pos(3))-E2low(3,label2pos(3)))*1e3, state2{3}, 'color', colors(6,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );
text(x2(label2pos(4)),(E2(4,label2pos(4))-E2low(4,label2pos(4)))*1e3, state2{2}, 'color', colors(7,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );
text(x2(label2pos(5)),(E2(5,label2pos(5))-E2low(5,label2pos(5)))*1e3, state2{1}, 'color', colors(8,:),'verticalalignment', 'bottom','horizontalalignment', 'left', 'fontsize', txsz );

hold off;grid on;box on;
xticks(bigxticks );
% xlabel('Applied RF (MHz)');
ylabel('\Delta{}E_{ACZ} (kHz)');
title('AC Zeeman Energy F=2', 'fontsize', tsz );
xticks(smallxticks );%yticks(smallyticks );
ylim([-150,250]);
xlim([19.5,20.5]);
text(19.54,-100, '(c2)', 'color', 'k', 'fontsize', 20);

% g(4) = subplot(2,4,4); hold on;
g(4)=axes('position', [4*buffx+3*eachx,2*buffy+1*eachy+topup ,eachx,eachy] );hold on;
for s=1:1:5 %        eigen vector from init index 1 (mf+2) into other states
%     plot(x2, squeeze(V2(s,5,:).^2), 'color', colors(3+s,:), 'linewidth', 1);
    plot(x2, F2(s,:), 'color', colors(3+s,:), 'linewidth', 1);
end
% label2pos = [2540,2600,2010,1400,1420];
% text(x2(label2pos(1))+.1,squeeze(V2(1,5,label2pos(1)).^2), states{4}, 'color', colors(4,:),'verticalalignment', 'middle','horizontalalignment', 'left', 'fontsize', txsz );
% text(x2(label2pos(2)),squeeze(V2(2,5,label2pos(2)).^2), states{5}, 'color', colors(5,:),'verticalalignment', 'bottom','horizontalalignment', 'left', 'fontsize', txsz );
% text(x2(label2pos(3)),squeeze(V2(3,5,label2pos(3)).^2)+.1,    states{6}, 'color', colors(6,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );
% text(x2(label2pos(4)),squeeze(V2(4,5,label2pos(4)).^2), states{7}, 'color', colors(7,:),'verticalalignment', 'bottom','horizontalalignment', 'right', 'fontsize', txsz );
% text(x2(label2pos(5)),squeeze(V2(5,5,label2pos(5)).^2)-.025, states{8}, 'color', colors(8,:),'verticalalignment', 'middle','horizontalalignment', 'right', 'fontsize', txsz );


hold off;grid on;box on;
xticks(bigxticks );
% xlabel('Applied RF (MHz)');
ylabel('Force / mg');
title('ACZ Force F=2', 'fontsize', tsz );
% set(get(gca, 'title'), 'position', [20,1,0]);
% xticks(smallxticks );%yticks(smallyticks );
ylim([-1.5,1]);
xlim([19,21]);
text(19.05,-1.2, '(d2)', 'color', 'k', 'fontsize', 20);

% F=1
% g(6)=subplot(2,4,6); hold on;
g(6)=axes('position', [2*buffx+1*eachx,1*buffy+0*eachy,eachx,eachy] );hold on;
hold on; 
for s=1
    plot(x1,E1low(s,:), 'color', .6*ones(1,3)); 
    text(xoff+x1(right),E1low(s,right), bare1{4-s}, 'color', .3*ones(1,3),'verticalalignment', 'middle','horizontalalignment', 'left','fontsize', txsz );
end
for s=1
    plot(x1,E1(s,:),  'color', colors(4-s,:), 'linewidth', 1); 
end
for s=2:1:3
    plot(x1,E1low(s,:), 'color', .6*ones(1,3)); 
    text(xoff+x1(right),E1low(s,right), bare1{4-s}, 'color', .3*ones(1,3),'verticalalignment', 'middle','horizontalalignment', 'left', 'fontsize', txsz );
end
for s=2:1:3
    plot(x1,E1(s,:), 'color',  colors(4-s,:), 'linewidth', 1); 
end
label1pos = [2090,2350,2090];
text(x1(label1pos(1)),E1(1,label1pos(1)), state1{3}, 'color', colors(3,:),'verticalalignment', 'top','horizontalalignment', 'center', 'fontsize', txsz );
text(x1(label1pos(2)),E1(2,label1pos(2))+.02, state1{2}, 'color', colors(2,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );
text(x1(label1pos(3)),E1(3,label1pos(3))+.05, state1{1}, 'color', colors(1,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );

hold off;  grid on;box on;
xticks(smallxticks );yticks(smallyticks );
xlabel('Applied RF (MHz)');
ylabel('E_{ACZ} /h, MHz');
title('Eigenstate Energies F=1', 'fontsize', tsz );
% linkaxes([g(1) g(3)], 'xy');
ylim([-1.5,1]); xlim([19.5,20.5]);
text(20.03,.75, '(b1)', 'color', 'k', 'fontsize', 20 );

% g(7) = subplot(2,4,7);hold on ;
g(7)=axes('position', [3*buffx+2*eachx,1*buffy+0*eachy,eachx,eachy] );hold on;
for s=1:1:3
    plot(x1,(E1(s,:)-E1low(s,:))*1e3, 'color', colors(4-s,:), 'linewidth', 1); 
end
label1pos = [2080,2080,2090];
text(x1(label1pos(1)),(E1(1,label1pos(1))-E1low(1,label1pos(1)))*1e3-15, state1{3}, 'color', colors(3,:),'verticalalignment', 'top','horizontalalignment', 'center', 'fontsize', txsz );
text(x1(label1pos(2)),(E1(2,label1pos(2))-E1low(2,label1pos(2)))*1e3-10, state1{2}, 'color', colors(2,:),'verticalalignment', 'top','horizontalalignment', 'center', 'fontsize', txsz );
text(x1(label1pos(3)),(E1(3,label1pos(3))-E1low(3,label1pos(3)))*1e3+20, state1{1}, 'color', colors(1,:),'verticalalignment', 'bottom','horizontalalignment', 'center', 'fontsize', txsz );

hold off;grid on;box on;
xticks(bigxticks );
xlabel('Applied RF (MHz)');
ylabel('\Delta{}E_{ACZ} (kHz)');
title('AC Zeeman Energy F=1', 'fontsize', tsz );
xticks(smallxticks );%yticks(smallyticks );
ylim([-150,250]);
xlim([19.5,20.5]);
text(19.54,210, '(c1)', 'color', 'k', 'fontsize', 20 );



% g(8) = subplot(2,4,8); hold on;
g(8)=axes('position', [4*buffx+3*eachx,1*buffy+0*eachy,eachx,eachy] );hold on;
for s=1:1:3 %        eigen vector from init index 1 (mf+2) into other states
%     plot(x1, squeeze(V1(s,3,:).^2), 'color', colors(s,:), 'linewidth', 1);
    plot(x1, F1(s,:), 'color', colors(4-s,:), 'linewidth', 1);
end
% 
% label1pos = [1969,2150,2179];
% text(x1(label1pos(1))-.05,squeeze(V1(3,3,label1pos(1)).^2), states{3}, 'color', colors(3,:),'verticalalignment', 'middle','horizontalalignment', 'right', 'fontsize', txsz );
% text(x1(label1pos(2))+.1,squeeze(V1(2,3,label1pos(2)).^2), states{2}, 'color', colors(2,:),'verticalalignment', 'middle','horizontalalignment', 'left', 'fontsize', txsz );
% text(x1(label1pos(3))+.05,squeeze(V1(1,3,label1pos(3)).^2), states{1}, 'color', colors(1,:),'verticalalignment', 'middle','horizontalalignment', 'left', 'fontsize', txsz );


hold off;grid on;box on;
xticks(bigxticks );
xlabel('Applied RF (MHz)');
ylabel('Force / mg');
title('ACZ Force F=1', 'fontsize', tsz );
% set(get(gca, 'title'), 'position', [20,1,0]);
% xticks(smallxticks );%yticks(smallyticks );
ylim([-1.5,1]);
xlim([19,21]);
text(19.05,.75, '(d1)', 'color', 'k', 'fontsize', 20 );


set(g([3,4,7,8]), 'YAxisLocation','right');
% set(g,'FontSize',12)
% linkaxes([g], 'x');
% xlim([19.1,20.9]);


% linkaxes([g(2) g(4)], 'xy');

% ylim([-2,2]);
% linkaxes(h, 'xy');
% 
% figure('position', [367         410        1340         484] ,'renderer', 'painter');
% 
% 
% for p = 1:1:5
%     subplot(2,5,6-p); hold on;
%     for s=1:1:5 %        eigen vector from init index 1 (mf+2) into other states
%         plot(x2, squeeze(V2(s,p,:).^2), 'color', colors(3+s,:), 'linewidth', 1);
%     end
%     if p==5
%         xlabel('Applied RF (MHz)');
%         ylabel('Population Ratio');
%     end
%     grid on; box on; 
%     xlim([19,21]);    ylim([0,1]);
%     xticks([19:.5:21]);
%     title(state2{6-p}, 'color', colors(3+p,:));
% end
% 
% for p = 1:1:3
%     subplot(2,5,6+p); hold on;
%     for s=1:1:3 %        eigen vector from init index 1 (mf+2) into other states
%         plot(x1, squeeze(V1(s,4-p,:).^2), 'color', colors(s,:), 'linewidth', 1);
%     end
%     grid on; box on; 
%     xlim([19,21]);    ylim([0,1]);
%     xticks([19:.5:21]);
%     title(state1{p}, 'color', colors(p,:));
% end
% 
% 
% 



