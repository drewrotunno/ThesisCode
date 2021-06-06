clear; close all;
% hamiltonian in the I/S basis to diagonalize better than the F basis. 
% There's so many signs in here, there may be an error, but I did my best. 

%% constants
h    = 6.62606896e-34;  %Js
hbar = h/2/pi;          %Js
Ahfs = h*3.417341305452145e9;% 5s1/2
% muB  = 9.27400915e-24      %J/T
muB  = h*1.399624604e6;      %J / G
gs =  2.0023193043622;
% gl =  0.99999369;
% gj =  2.00233113;
gi = -0.0009951414;

%% parameters
% B = [0:1:10000]*1e-4;   % Teslas
B = [.001,12,.1:.1:.9,1:1:10000];   % Gausses
diffBval = 52;         %check diff at XX Gauss
Bloc = find(B==diffBval);

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
H_GHz = zeros(8,size(B,2));
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
[V(:,:,b),D(:,:,b)] = eig(Hftot(:,:,b));        % eig can only work on 2-D matrices
% re-order into +2+10-1-2/+10-1 order (from low-to-high Energies) 
% Vr(:,:,b) = [V(:,8,b),V(:,7,b),V(:,6,b),V(:,5,b),V(:,4,b),V(:,1,b),V(:,2,b),V(:,3,b)];
%render convenient signs
Vr(:,:,b)=[ V(:,8,b)* sign(V(1,8,b)),...
            V(:,7,b)* sign(V(2,7,b)),...
            V(:,6,b)* sign(V(3,6,b)),...
            V(:,5,b)* sign(V(4,5,b)),...
            V(:,4,b)* sign(V(5,4,b)),...
            V(:,1,b)*-sign(V(6,1,b)),...
            V(:,2,b)*-sign(V(7,2,b)),...
            V(:,3,b)*-sign(V(8,3,b))];
% Vr2(:,:,b) = Vr(:,:,b).^2;
% Vr(Vr<1e-6) = 0
H_GHz(:,b) = diag(D(:,:,b))/h/1e9;
% I think for these vector-calculation, should use the 2-d mats in the loop 
% Not obvious how to do this operation across and indep. of B dimension
% sigma plus RF
CGs(1,2,b) = (Vr(:,1,b)'*SpM*Vr(:,2,b)); 
CGs(2,3,b) = (Vr(:,2,b)'*SpM*Vr(:,3,b)); 
CGs(3,4,b) = (Vr(:,3,b)'*SpM*Vr(:,4,b)); 
CGs(4,5,b) = (Vr(:,4,b)'*SpM*Vr(:,5,b)); 
% sigma plus microwave
CGs(1,6,b) = (Vr(:,1,b)'*SpM*Vr(:,6,b));
CGs(2,7,b) = (Vr(:,2,b)'*SpM*Vr(:,7,b));
CGs(3,8,b) = (Vr(:,3,b)'*SpM*Vr(:,8,b));
% sigma minus RF
CGs(7,8,b) = (Vr(:,7,b)'*SmM*Vr(:,8,b));
CGs(6,7,b) = (Vr(:,6,b)'*SmM*Vr(:,7,b));
% sigma minus Microwave
CGs(3,6,b) = (Vr(:,3,b)'*SmM*Vr(:,6,b));
CGs(4,7,b) = (Vr(:,4,b)'*SmM*Vr(:,7,b));
CGs(5,8,b) = (Vr(:,5,b)'*SmM*Vr(:,8,b));
% pi trans Microwave
CGs(2,6,b) = (Vr(:,2,b)'*SzM*Vr(:,6,b));
CGs(3,7,b) = (Vr(:,3,b)'*SzM*Vr(:,7,b));
CGs(4,8,b) = (Vr(:,4,b)'*SzM*Vr(:,8,b));
% copy to the transpose
CGs(:,:,b) = CGs(:,:,b) +CGs(:,:,b)';

% sigma plus RF
CGi(1,2,b) = (Vr(:,1,b)'*IpM*Vr(:,2,b));
CGi(2,3,b) = (Vr(:,2,b)'*IpM*Vr(:,3,b));
CGi(3,4,b) = (Vr(:,3,b)'*IpM*Vr(:,4,b));
CGi(4,5,b) = (Vr(:,4,b)'*IpM*Vr(:,5,b));
% sigma plus microwave
CGi(1,6,b) = (Vr(:,1,b)'*IpM*Vr(:,6,b));
CGi(2,7,b) = (Vr(:,2,b)'*IpM*Vr(:,7,b));
CGi(3,8,b) = (Vr(:,3,b)'*IpM*Vr(:,8,b));
% sigma minus RF
CGi(7,8,b) = (Vr(:,7,b)'*ImM*Vr(:,8,b));
CGi(6,7,b) = (Vr(:,6,b)'*ImM*Vr(:,7,b));
% sigma minus Microwave
CGi(3,6,b) = (Vr(:,3,b)'*ImM*Vr(:,6,b));
CGi(4,7,b) = (Vr(:,4,b)'*ImM*Vr(:,7,b));
CGi(5,8,b) = (Vr(:,5,b)'*ImM*Vr(:,8,b));
% pi trans Microwave
CGi(2,6,b) = (Vr(:,2,b)'*IzM*Vr(:,6,b));
CGi(3,7,b) = (Vr(:,3,b)'*IzM*Vr(:,7,b));
CGi(4,8,b) = (Vr(:,4,b)'*IzM*Vr(:,8,b));
% copy to the transpose
CGi(:,:,b) = CGi(:,:,b) +CGi(:,:,b)';
end

CGsg = CGs.*gs;
CGig = CGi.*gi;
CGg = CGig+CGsg; 


%% plotting
figure; colors = 'mkgrgkmb'; hold on;
for s = 1:1:8
%     plot(B,H_GHz(s,:), [colors(s),'x']);
    plot(B,H_GHz(s,:), colors(s));
end
ylabel('Energy, GHz'); xlabel('Gauss');

figure('position', [90         178        1603         800], 'renderer', 'painter');

% RF with g's, combined
subplot(3,5,1); hold on; title('RF F=2 g_iI_+ + g_sS_+'); grid on;
plot(B, squeeze(CGg(1,2,:)), 'b');plot(B, squeeze(CGg(2,1,:)), '--m');
plot(B, squeeze(CGg(2,3,:)), 'm');plot(B, squeeze(CGg(3,2,:)), '--k');
plot(B, squeeze(CGg(3,4,:)), 'k');plot(B, squeeze(CGg(4,3,:)), '--g');
plot(B, squeeze(CGg(4,5,:)), 'g');plot(B, squeeze(CGg(5,4,:)), '--r');
% RF S+ only
subplot(3,5,6); hold on; title('RF F=2 S_+'); grid on;
plot(B, squeeze(CGs(1,2,:)), 'b');plot(B, squeeze(CGs(2,1,:)), '--m');
plot(B, squeeze(CGs(2,3,:)), 'm');plot(B, squeeze(CGs(3,2,:)), '--k');
plot(B, squeeze(CGs(3,4,:)), 'k');plot(B, squeeze(CGs(4,3,:)), '--g');
plot(B, squeeze(CGs(4,5,:)), 'g');plot(B, squeeze(CGs(5,4,:)), '--r');
% RF I+ only
subplot(3,5,11); hold on; title('RF F=2 I_+'); grid on;
plot(B, squeeze(CGi(1,2,:)), 'b');plot(B, squeeze(CGi(2,1,:)), '--m');
plot(B, squeeze(CGi(2,3,:)), 'm');plot(B, squeeze(CGi(3,2,:)), '--k');
plot(B, squeeze(CGi(3,4,:)), 'k');plot(B, squeeze(CGi(4,3,:)), '--g');
plot(B, squeeze(CGi(4,5,:)), 'g');plot(B, squeeze(CGi(5,4,:)), '--r');

% RF with g's, sigma -
subplot(3,5,2); hold on; title('RF F=1 g_iI_- + g_sS_-'); grid on;
plot(B, squeeze(CGg(6,7,:)), 'k');plot(B, squeeze(CGg(7,6,:)), '--g');
plot(B, squeeze(CGg(7,8,:)), 'm');plot(B, squeeze(CGg(8,7,:)), '--k');
% rf S- only
subplot(3,5,7); hold on; title('RF F=1 S_-'); grid on;
plot(B, squeeze(CGs(6,7,:)), 'k');plot(B, squeeze(CGs(7,6,:)), '--g');
plot(B, squeeze(CGs(7,8,:)), 'm');plot(B, squeeze(CGs(8,7,:)), '--k');
% RF I- only
subplot(3,5,12); hold on; title('RF F=1 I_-'); grid on;
plot(B, squeeze(CGi(6,7,:)), 'k');plot(B, squeeze(CGi(7,6,:)), '--g');
plot(B, squeeze(CGi(7,8,:)), 'm');plot(B, squeeze(CGi(8,7,:)), '--k');

%Microwave with g's, combined sigma + transitions
subplot(3,5,3); hold on; title('\mu{}W g_iI_+ + g_sS_+'); grid on;
plot(B, squeeze(CGg(1,6,:)), 'b');plot(B, squeeze(CGg(6,1,:)), '--m');
plot(B, squeeze(CGg(2,7,:)), 'm');plot(B, squeeze(CGg(7,2,:)), '--k');
plot(B, squeeze(CGg(3,8,:)), 'k');plot(B, squeeze(CGg(8,3,:)), '--g');
% microwave S+ only
subplot(3,5,8); hold on; title('\mu{}W S_+'); grid on;
plot(B, squeeze(CGs(1,6,:)), 'b');plot(B, squeeze(CGs(6,1,:)), '--m');
plot(B, squeeze(CGs(2,7,:)), 'm');plot(B, squeeze(CGs(7,2,:)), '--k');
plot(B, squeeze(CGs(3,8,:)), 'k');plot(B, squeeze(CGs(8,3,:)), '--g');
% microwave I+ only
subplot(3,5,13);hold on; title('\mu{}W I_+'); grid on;
plot(B, squeeze(CGi(1,6,:)), 'b');plot(B, squeeze(CGi(6,1,:)), '--m');
plot(B, squeeze(CGi(2,7,:)), 'm');plot(B, squeeze(CGi(7,2,:)), '--k');
plot(B, squeeze(CGi(3,8,:)), 'k');plot(B, squeeze(CGi(8,3,:)), '--g');

%Microwave with g's, combined PI
subplot(3,5,4); hold on; title('\mu{}W g_iI_z + g_sS_z'); grid on;
plot(B, squeeze(CGg(2,6,:)), 'm');plot(B, squeeze(CGg(6,2,:)), ':w');
plot(B, squeeze(CGg(3,7,:)), 'k');plot(B, squeeze(CGg(7,3,:)), ':w');
plot(B, squeeze(CGg(4,8,:)), 'g');plot(B, squeeze(CGg(8,4,:)), ':w');
% microwave Sz only
subplot(3,5,9); hold on; title('\mu{}W S_z'); grid on;
plot(B, squeeze(CGs(2,6,:)), 'm');plot(B, squeeze(CGs(6,2,:)), ':w');
plot(B, squeeze(CGs(3,7,:)), 'k');plot(B, squeeze(CGs(7,3,:)), ':w');
plot(B, squeeze(CGs(4,8,:)), 'g');plot(B, squeeze(CGs(8,4,:)), ':w');
% microwave Iz only
subplot(3,5,14);hold on; title('\mu{}W I_z'); grid on;
plot(B, squeeze(CGi(2,6,:)), 'm');plot(B, squeeze(CGi(6,2,:)), ':w');
plot(B, squeeze(CGi(3,7,:)), 'k');plot(B, squeeze(CGi(7,3,:)), ':w');
plot(B, squeeze(CGi(4,8,:)), 'g');plot(B, squeeze(CGi(8,4,:)), ':w');

%Microwave with g's, combined sigma - transitions
subplot(3,5,5); hold on; title('\mu{}W g_iI_- + g_sS_-'); grid on;
plot(B, squeeze(CGg(3,6,:)), 'k');plot(B, squeeze(CGg(6,3,:)), '--m');
plot(B, squeeze(CGg(4,7,:)), 'g');plot(B, squeeze(CGg(7,4,:)), '--k');
plot(B, squeeze(CGg(5,8,:)), 'r');plot(B, squeeze(CGg(8,5,:)), '--g');
% microwave S+ only
subplot(3,5,10); hold on; title('\mu{}W S_-'); grid on;
plot(B, squeeze(CGs(3,6,:)), 'k');plot(B, squeeze(CGs(6,3,:)), '--m');
plot(B, squeeze(CGs(4,7,:)), 'g');plot(B, squeeze(CGs(7,4,:)), '--k');
plot(B, squeeze(CGs(5,8,:)), 'r');plot(B, squeeze(CGs(8,5,:)), '--g');
% microwave I+ only
subplot(3,5,15);hold on; title('\mu{}W I_-'); grid on;
plot(B, squeeze(CGi(3,6,:)), 'k');plot(B, squeeze(CGi(6,3,:)), '--m');
plot(B, squeeze(CGi(4,7,:)), 'g');plot(B, squeeze(CGi(7,4,:)), '--k');
plot(B, squeeze(CGi(5,8,:)), 'r');plot(B, squeeze(CGi(8,5,:)), '--g');


%% find diff from square roots at 12 Gauss
% RF outer
difBg(1)= (CGg(2,1,Bloc) - sqrt(2/8))./sqrt(2/8);
difBg(2)= (CGg(3,2,Bloc) - sqrt(3/8))./sqrt(3/8);
difBg(3)= (CGg(4,3,Bloc) - sqrt(3/8))./sqrt(3/8);
difBg(4)= (CGg(5,4,Bloc) - sqrt(2/8))./sqrt(2/8);
difBg.*100

difBs(1)= (2*CGs(2,1,Bloc) - sqrt(2/8))./sqrt(2/8);
difBs(2)= (2*CGs(3,2,Bloc) - sqrt(3/8))./sqrt(3/8);
difBs(3)= (2*CGs(4,3,Bloc) - sqrt(3/8))./sqrt(3/8);
difBs(4)= (2*CGs(5,4,Bloc) - sqrt(2/8))./sqrt(2/8);
difBs.*100

difBd(1)= (2*CGs(2,1,Bloc) - CGg(2,1,Bloc))./CGg(2,1,Bloc);
difBd(2)= (2*CGs(3,2,Bloc) - CGg(3,2,Bloc))./CGg(3,2,Bloc);
difBd(3)= (2*CGs(4,3,Bloc) - CGg(4,3,Bloc))./CGg(4,3,Bloc);
difBd(4)= (2*CGs(5,4,Bloc) - CGg(5,4,Bloc))./CGg(5,4,Bloc);
difBd.*100

