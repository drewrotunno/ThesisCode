% Based on Silvester's book example code. 
% Jm = -i * w * mu * sigma * height^2 / 2pi * sum(Jn * log Dmn) + sigma*E
% solve for  (U+i*Z) J = Volts by inverting inv(U+i*Z) * Volts = J

% much wider than tall [===============]
% center (x) symmetry +   [Ncol=====1|1=====Ncol] ...
clear all; close all;
tic

%% settings
freq = [10, 50, 100, 250, 500, 1000,2500,5000 ];     % few kHz 
% freq = [2.4:.1:2.8]*1e3;
width = 80e-3;              %aluminum sheet
height = .65e-3;             %
Nrow = 2;                   %<-----set for precision
% detheight = .65e-3;
detheight = .0025;      % 1/r fit
% detheight = .00025;      % fudge

%% animation
doanimation = 1;
% doanimation = 0;
phasestep = 5; 
repeats = 3;    
delay = .1;

%% constants 
Dmm = .4470491559036622;        %geometric mean of a square from itself
g = 2.5e7;                %  conductivity
% g = 1.6e7;        %fudge
E0 = 1;             % voltage across, Volts. 
mu0 = 4*pi*1e-7;

%% creation
dy = height/Nrow;
dx = dy;                    %must be squares
Ncol = floor(width/2/dy);     %if too high, cut width or Nrow
N = Nrow*Ncol
% if(N>2200);disp('no thank you');pause;end %takes way too long. 
%       cancel if casual, comment if hardcore

% xplot = linspace(-width/2,width/2,Ncol);        %NOT the same as square center positions. 
% yplot = linspace(-height/2,height/2,Nrow);      %same. Use xn, yn values for that
xplot = [dx/2:dx:(Ncol-1/2)*dx];
yplot = [dy/2:dy:(Nrow-1/2)*dy];
fullxplot = [-fliplr(xplot),xplot];

Bsamplex = round(-60e-3:.5e-3:60e-3,4);
dety = detheight+height;
Bx = cell(numel(freq),1);

halfwidth = width/2;
halfheight = height/2;

Z = zeros(N);
X = zeros(N);
D = zeros(N);
% ynmat = zeros(numel(Bsamplex),N);
% ymmat = zeros(numel(Bsamplex),N);
Dother = zeros(N);
ydist = zeros(Nrow,numel(freq));
xdist = zeros(numel(freq),Ncol*2);
curnt = cell(numel(freq),1);
currentforB= cell(numel(freq),1);
normcurrent= cell(numel(freq),1);
secondnorm = zeros(numel(freq));
totalcurrent = zeros(numel(freq),1);
ymax = 0;


%% matrix population D,X

for f=1:1:numel(freq)
    % angl = zeros(10);
    const = freq(f)*mu0*g*dy*dx;

    for n=1:1:N             % n is main, m is others
            yn = (Nrow - 1/2 - floor((n-1)/Ncol) )*dy;
            xn = -(1/2 + mod(n-1,Ncol))*dx;
        for m=1:1:N
              ym = (Nrow - 1/2 - floor((m-1)/Ncol) )*dy;
              xm = -(1/2 + mod(m-1,Ncol))*dx;
            if(n==m)
                X(n,n) = const*log(Dmm*dx* 2*xn );
            else
                X(m,n) = const*log(hypot(xn-xm,yn-ym)*hypot(xn+xm,yn-ym));
            end
        end
    end
    
%%  Calculation

    Z = eye(N)+1i.*X;
        if(abs(imag(Z(1,1)))<1e-6)     %convergence check...takes forever! 
             disp(Z(1:3,1:3));pause;
        end
    volts = g.*ones(N,1);
    thiscurrent = Z\volts;

    curnt{f} = [fliplr(reshape(thiscurrent, Ncol,Nrow)'),reshape(thiscurrent, Ncol,Nrow)'];
    totalcurrent(f) = sum(sum(abs(curnt{f})));
    normcurrent{f} = curnt{f};%./totalcurrent(f);
    currentforB{f} = thiscurrent;%./sum(abs(thiscurrent));
%     currentforB{f} = currentforB{f}./mean(mean(normcurrent{1}));
%     secondnorm = mean(mean(normcurrent{1}));
%     normcurrent{f} = normcurrent{f}./secondnorm;
    
    ydist(:,f) = mean(abs(normcurrent{f}),2);
    xdist(f,:) = mean((normcurrent{f}),1);
%     ymean(f) = mean(ydist(:,f));
%     xmean(f) = mean(xdist(f,:));
    yabs(:,f) = mean(abs(normcurrent{f}),2);
    xabs(f,:) = mean(abs(normcurrent{f}),1);

%     zmax = max(zmax, max(max(abs(currentforB{f}))));
end
toc


% Bfield calculation
for f=1:1:numel(freq)
    Bx{f} = zeros(1,numel(Bsamplex));
    for d=1:1:numel(Bsamplex)
        for n = 1:1:N
            yn = (Nrow - 1/2 - floor((n-1)/Ncol) )*dy;
            xn = -(1/2 + mod(n-1,Ncol))*dx;
%                 ynmat(d,n) = yn;
%                 xnmat(d,n) = xn;
            Bx{f}(d) =  Bx{f}(d)+ ...
                (mu0*currentforB{f}(n)/4/pi)*...
                ((dety-yn)/hypot(dety-yn, Bsamplex(d)-xn)^2 ...
                +(dety-yn)/hypot(dety-yn, Bsamplex(d)+xn)^2);
%                 ((dety-yn)/hypot(dety-yn, Bsamplex(d)-xn) ...
%                 +(dety-yn)/hypot(dety-yn, Bsamplex(d)+xn));
        end
    end
%     Bx{f} = Bx{f}./(mean(abs(Bx{f}(find(Bsamplex==-.04):find(Bsamplex==.04)))));
    ymax = max(ymax, max(max(abs(Bx{f}))));
end



% 


% 
figure( 'position', [     9    49   800   900]) ;
% color = ['brmgcky']; colormap('jet');
   subplot(2,2,1); hold on; title('amplitude, current');grid on;
   subplot(2,2,3); hold on;title('phase, current');grid on;
   subplot(2,2,2); hold on; title('amplitude, magnetic');grid on;
   subplot(2,2,4); hold on;title('phase, magnetic');grid on;
   
for f=1:1:numel(freq)
   subplot(2,2,1);
   plot(fullxplot, xabs(f,:));
      
   subplot(2,2,3);
   plot(fullxplot, rad2deg(unwrap(angle(xdist(f,:)))));
   
   subplot(2,2,2);
   plot(Bsamplex, abs(Bx{f}));
      
   subplot(2,2,4);
   plot(Bsamplex, -rad2deg(unwrap(angle(Bx{f}))));
     
end


if(doanimation)

%% animation
figure('position', [860    70   600   900]);
yJmax = max(max(abs(xdist)));
   subplot(2,1,1); hold on; title('Magnetic field ');xlabel('transverse position (m)');ylabel('');grid on;
   subplot(2,1,2); hold on;title('Current density');xlabel('transverse position (m)');ylabel('');grid on;
for phase = deg2rad(0:phasestep:360*repeats)
        subplot(2,1,1);  cla;
        subplot(2,1,2);  cla;
    for f=1:1:numel(freq)
        subplot(2,1,1);  hold on;
        plot(Bsamplex, real(Bx{f}*exp(-1i*phase)));
        xlim([min(Bsamplex) max(Bsamplex)]);
        ylim([-ymax ymax]);
        
        subplot(2,1,2);  hold on;
        plot(fullxplot, real(xdist(f,:)*exp(1i*phase)));
        xlim([min(fullxplot) max(fullxplot)]);
        ylim([-yJmax yJmax]);
    end
    pause(delay)
end
end
