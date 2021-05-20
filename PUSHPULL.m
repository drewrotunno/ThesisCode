% function plotwitherrorbars_play()
%% compare the theory and experiment for push pull RF ACZ data

close all;
clear all;
 
%% options
whichforce = 'y';
% whichforce = 'x';
% whichforce = 'xy';

% load('newGradRabi.mat');
% load('newGradRabibigERROR.mat');        %fudged in a 3 micron systematic error
% load('gradrabi20dec18');
% load('gradrabi07dec18');
% load('gradrabi07dec18v2');

powerfactor = sqrt(10^.134);                    % between 5vva (rabi meas.) and 10vva (push/pull)

% % from fits:
Rabi = 3.707623580923737e+05;                       % day of the push/pulls, re-analyzed 20Mar2019
errorRabi = 9.862788845343422e+02;                  %  weighted AVG of the two 'high' sets
% Rabi = 3.707603626818930e+05;                       % day of the push/pulls, re-averaged 26Mar2019
% errorRabi = 9.857921359053023e+02;                  %  weighted AVG of the two 'high' sets

Rabi = powerfactor*Rabi;
errorRabi = powerfactor*errorRabi;

gradlowfreq = 3.525896671828356e+05;                % Feb 21 data, re-analyzed on 20Mar2019
gradlowfreqerr = 1.178903498091165e+03;             %  "
gradhighfreq = 3.744615410762737e+05;               %  " 
gradhighfreqerr = 9.438510242136450e+02;            %  " 
[gradmidfreq,gradmidfreqerr] = weightedAVG([gradlowfreq,gradhighfreq],[gradlowfreqerr,gradhighfreqerr]);



% gradfact = powerfactor;
% gradfact = Rabi / gradlowfreq;        %better?
% gradfacterr = gradfact*hypot(errorRabi/Rabi,gradlowfreqerr/gradlowfreq);
% gradfact = Rabi / gradhighfreq;         %proper?
% gradfacterr = hypot(errorRabi/Rabi,gradlowfreqerr/gradhighfreq);
gradfact = Rabi / gradmidfreq;         %in the middle?
gradfacterr = hypot(errorRabi/Rabi,gradmidfreqerr/gradmidfreq);

% avgGradRabi = [3296347218.24157];                   %  weighted AVG of 4 -> power scaled
% avgGradRabiError = [146215965.287278];              %  and 2 at high power. all combos of high/low locations
avgGradRabi = [3.327787908770734e9]*gradfact;                    %  line fit of 6 pts
avgGradRabiError = [.55e9]*gradfact;                 %  

avgGradRabiError = hypot(avgGradRabiError,gradfacterr);

avgdR = 1e-6;           % any small number will work


%% take in data, sort
load('DECpushpulldata');   

% column 1 is frequency
% column 2 is x pixel
% column 3 is x pixel error
% column 4 is y pixel
% column 5 is y pixel error
% cells: { -2, -1, 0, 1, 2}
adata{5} = sortrows(p2adata, 1);
adata{4}  = sortrows(p1adata, 1);
adata{3}  = sortrows(z0adata, 1);
adata{2}  = sortrows(m1adata, 1);
adata{1}  = sortrows(m2adata, 1);

bdata{5} = sortrows(p2bdata, 1);
bdata{4}  = sortrows(p1bdata, 1);
bdata{3}  = sortrows(z0bdata, 1);
bdata{2}  = sortrows(m1bdata, 1);
bdata{1}  = sortrows(m2bdata, 1);

monofreqs = unique(adata{1}(:,1));
triplefreqs = adata{1}(:,1);

clear p2adata p1adata z0adata m1adata m2adata p2bdata p1bdata z0bdata m1bdata m2bdata;

% find deltay px --> microns --> Force
% posorforce = 'force';
posorforce = 'pos';

switch posorforce
    case 'force'
    tof = 17.15e-3;
    thold = .5e-3;
    grav = 9.81; 
    case 'pos'
    tof = 1;
    thold = 1;
    grav = 1 ; 
end
gravity = 9.81;
dperpx = 7.4e-6;


% get baselines 
load('DECbaselines');
% usedataabove = [    [1,1,1,0,1];
%                     [1,1,1,0,1];
%                     [1,1,0,0,1] ];
                
usedataabove = [    [1,1,1,0,1];
                    [1,1,1,0,1];
                    [1,1,0,0,1] ];
usedatabelow = [    [1,1,0,0,1];
                    [1,1,0,0,1];
                    [1,1,0,0,1] ];
% usedataabove = [    [0,0,0,0,0];
%                     [0,0,0,0,0];
%                     [0,0,0,0,0] ];
% usedatabelow = [    [0,0,0,0,0];
%                     [0,0,0,0,0];
%                     [0,0,0,0,0] ];    
                
% baseline = baseabove, basebelow
% [ 3x5 ]  =>  -2,-1,0,1,2 l2r,  set123 t2b
%   1x, 2xe, 3y, 4yerror in each cell

meanbaseabovey = zeros(3,5);
errorbaseabovey = zeros(3,5);
meanbasebelowy = zeros(3,5);
errorbasebelowy = zeros(3,5);

meanbaseabovex = zeros(3,5);
errorbaseabovex = zeros(3,5);
meanbasebelowx = zeros(3,5);
errorbasebelowx = zeros(3,5);

xfit = [9.903007321237040e-02, 2.911822570865862e+02];
xfitsigma = [3.203519934231527e-02 ,    5.251255080145256e-02];
yfit = [-1.186367240531855e+00, 2.402143932572746e+02];
yfitsigma = [2.555677105300469e-02 ,    4.189301973437409e-02];


for n=1:1:size(usedataabove,1)
    for m=1:1:size(usedataabove,2)
        clear baselineweighty baselineweightyx;
        if(usedataabove(n,m))

            baselineweighty = 1./baselineabove{n,m}(:,4).^2;
            baselineweightx = 1./baselineabove{n,m}(:,2).^2;

            meanbaseabovey(n,m) = sum(baselineweighty.*baselineabove{n,m}(:,3))./sum(baselineweighty);
            meanbaseabovex(n,m) = sum(baselineweightx.*baselineabove{n,m}(:,1))./sum(baselineweightx);

            errorbaseabovey(n,m) = 1./sqrt(sum(baselineweighty.^2));
            errorbaseabovex(n,m) = 1./sqrt(sum(baselineweightx.^2));
        else    % use fit numbers
            meanbaseabovey(n,m) = yfit(1)*(m-3)+yfit(2);
            meanbaseabovex(n,m) = xfit(1)*(m-3)+xfit(2);

            errorbaseabovey(n,m) = yfitsigma(1)*abs(m-3)+yfitsigma(2);
            errorbaseabovex(n,m) = xfitsigma(1)*abs(m-3)+xfitsigma(2);
        end
        
        if(usedatabelow(n,m))
            baselineweighty = 1./baselinebelow{n,m}(:,4).^2;
            baselineweightx = 1./baselinebelow{n,m}(:,2).^2;

            meanbasebelowy(n,m) = sum(baselineweighty.*baselinebelow{n,m}(:,3))./sum(baselineweighty);
            meanbasebelowx(n,m) = sum(baselineweightx.*baselinebelow{n,m}(:,1))./sum(baselineweightx);

            errorbasebelowy(n,m) = 1./sqrt(sum(baselineweighty.^2));
            errorbasebelowx(n,m) = 1./sqrt(sum(baselineweightx.^2));
        else    % use fit numbers
            meanbasebelowy(n,m) = yfit(1)*(m-3)+yfit(2);
            meanbasebelowx(n,m) = xfit(1)*(m-3)+xfit(2);

            errorbasebelowy(n,m) = yfitsigma(1)*abs(m-3)+yfitsigma(2);
            errorbasebelowx(n,m) = xfitsigma(1)*abs(m-3)+xfitsigma(2);
        end
    end
end

% calculate offsets
deltayabove = cell(5,1);
deltaxabove = cell(5,1);
deltaybelow = cell(5,1);
deltaxbelow = cell(5,1);
forceabove = cell(5,1);
forcebelow = cell(5,1);

% delta pos
for n=1:1:5
    for m=1:1:size(adata{1},1)
            % subtract baseline for pos
            %change in pos      pixel fit       baseline - one per 26frame set
        deltayabove{n}(m,1) = adata{n}(m,4) - meanbaseabovey(ceil(mod(n/26,26)),n);
        deltaxabove{n}(m,1) = adata{n}(m,2) - meanbaseabovex(ceil(mod(n/26,26)),n);
        deltaybelow{n}(m,1) = bdata{n}(m,4) - meanbasebelowy(ceil(mod(n/26,26)),n);
        deltaxbelow{n}(m,1) = bdata{n}(m,2) - meanbasebelowx(ceil(mod(n/26,26)),n);
            % send over error as well
%         deltayabove{n}(m,2) = sqrt(adata{n}(m,5).^2 + errorbaseabovey(ceil(mod(n/26,26)),n)^2);
%         deltaxabove{n}(m,2) = sqrt(adata{n}(m,3).^2 + errorbaseabovex(ceil(mod(n/26,26)),n)^2);
%         deltaybelow{n}(m,2) = sqrt(bdata{n}(m,5).^2 + errorbasebelowy(ceil(mod(n/26,26)),n)^2);
%         deltaxbelow{n}(m,2) = sqrt(bdata{n}(m,3).^2 + errorbasebelowx(ceil(mod(n/26,26)),n)^2);
        
        deltayabove{n}(m,2) = (adata{n}(m,5) + errorbaseabovey(ceil(mod(n/26,26)),n));
        deltaxabove{n}(m,2) = (adata{n}(m,3) + errorbaseabovex(ceil(mod(n/26,26)),n));
        deltaybelow{n}(m,2) = (bdata{n}(m,5) + errorbasebelowy(ceil(mod(n/26,26)),n));
        deltaxbelow{n}(m,2) = (bdata{n}(m,3) + errorbasebelowx(ceil(mod(n/26,26)),n));
    %     deltay{n}(:,2) =data{n}(:,3);
            %convert from px to meters
    end
        deltayabove{n} = deltayabove{n} .* dperpx;
        deltaybelow{n} = deltaybelow{n} .* dperpx;
        deltaxabove{n} = deltaxabove{n} .* dperpx;
        deltaxbelow{n} = deltaxbelow{n} .* dperpx;

    switch whichforce
        case 'y'
            forceabove{n} = deltayabove{n}./tof./thold./grav;
            forcebelow{n} = deltaybelow{n}./tof./thold./grav;
        case 'x'
            forceabove{n} = deltaxabove{n}./tof./thold./grav;
            forcebelow{n} = deltaxbelow{n}./tof./thold./grav;
        case 'xy'
            forceabove{n} = (deltaxabove{n}+deltayabove{n})./tof./thold./grav;
            forcebelow{n} = (deltaxbelow{n}+deltaybelow{n})./tof./thold./grav;
    end

end

% fix bad points
forceabove{4}(46,2) = 1e9;
forceabove{2}(75,2) = 1e9;
forcebelow{4}(1,2) = 1e9;

% average 3 measurements
for n=1:1:5
    forceaboveweight{n} = 1./(forceabove{n}(:,2).^2);
    forcebelowweight{n} = 1./(forcebelow{n}(:,2).^2);
    for m=1:1:numel(monofreqs)
        
        avgforceabove(m,n) = (   ...
            forceaboveweight{n}(3*m-2) * forceabove{n}(3*m-2,1)...
          + forceaboveweight{n}(3*m-1) * forceabove{n}(3*m-1,1) ...
          + forceaboveweight{n}(3*m) * forceabove{n}(3*m,1)   )...
      ./ (  forceaboveweight{n}(3*m-2)...
          + forceaboveweight{n}(3*m-1) ...
          + forceaboveweight{n}(3*m)   );
%         
        avgaboveforceerror(m,n) =  sqrt(...
           ( ((forceabove{n}(3*m-2,1)-avgforceabove(m,n))*(forceabove{n}(3*m-2,2)<10))^2 ...
          +((forceabove{n}(3*m-1,1)-avgforceabove(m,n))*(forceabove{n}(3*m-1,2)<10))^2 ...
          +((forceabove{n}(3*m,1)-avgforceabove(m,n))*(forceabove{n}(3*m,2)<10))^2 )./ 3);
% %       old : new ignores bad pts
%         avgaboveforceerror(m,n) =  sqrt(...
%            ( ((forceabove{n}(3*m-2,1)-avgforceabove(m,n)))^2 ...
%           +((forceabove{n}(3*m-1,1)-avgforceabove(m,n)))^2 ...
%           +((forceabove{n}(3*m,1)-avgforceabove(m,n)))^2 )./ 3);
        
      avgforcebelow(m,n) = (   ...
            forcebelowweight{n}(3*m-2) * forcebelow{n}(3*m-2,1)...
          + forcebelowweight{n}(3*m-1) * forcebelow{n}(3*m-1,1) ...
          + forcebelowweight{n}(3*m) * forcebelow{n}(3*m,1)   )...
      ./ (  forcebelowweight{n}(3*m-2)...
          + forcebelowweight{n}(3*m-1) ...
          + forcebelowweight{n}(3*m)   );
        
        avgbelowforceerror(m,n) =  sqrt(...
           ( ((forcebelow{n}(3*m-2,1)-avgforcebelow(m,n))*(forcebelow{n}(3*m-2,2)<10))^2 ...
          +((forcebelow{n}(3*m-1,1)-avgforcebelow(m,n))*(forcebelow{n}(3*m-1,2)<10))^2 ...
          +((forcebelow{n}(3*m,1)-avgforcebelow(m,n))*(forcebelow{n}(3*m,2)<10))^2 )./ 3);
    end
end

%% Theory curves

% fullfreq  = [5.2e6:1e4:12.2e6];
% fullfreq  = [6e6:1e4:11e6];
fullfreq  = [6e6:5e3:10.75e6];

%estimates
[forestA, INTforestA ] = integratedForce_func_50( fullfreq, Rabi, +5, avgGradRabi, avgdR);
[forestB, INTforestB ] = integratedForce_func_50( fullfreq, Rabi, -5, avgGradRabi, avgdR);

%low bound
[~, INTforceAlow ] = integratedForce_func_50( fullfreq, Rabi-errorRabi, +5, (avgGradRabi-avgGradRabiError), avgdR);
[~, INTforceBlow ] = integratedForce_func_50( fullfreq, Rabi-errorRabi, -5, (avgGradRabi-avgGradRabiError), avgdR);

%high bound
[~, INTforceAhigh ] = integratedForce_func_50( fullfreq, Rabi+errorRabi, +5, (avgGradRabi+avgGradRabiError), avgdR);
[~, INTforceBhigh ] = integratedForce_func_50( fullfreq, Rabi+errorRabi, -5, (avgGradRabi+avgGradRabiError), avgdR);

%estimates
% [forestA, INTforestA ] = intForce_func_move( fullfreq, Rabi, +5, avgGradRabi*gradfact, avgdR);
% [forestB, INTforestB ] = intForce_func_move( fullfreq, Rabi, -5, avgGradRabi*gradfact, avgdR);
% 
% %low bound
% [forceAlow, INTforceAlow ] = intForce_func_move( fullfreq, Rabi-errorRabi, +5, (avgGradRabi-avgGradRabiError)*gradfact, avgdR);
% [forceBlow, INTforceBlow ] = intForce_func_move( fullfreq, Rabi-errorRabi, -5, (avgGradRabi-avgGradRabiError)*gradfact, avgdR);
% 
% %high bound
% [forceAhigh, INTforceAhigh ] = intForce_func_move( fullfreq, Rabi+errorRabi, +5, (avgGradRabi+avgGradRabiError)*gradfact, avgdR);
% [forceBhigh, INTforceBhigh ] = integratedForce_func_50( fullfreq, Rabi+errorRabi, -5, (avgGradRabi+avgGradRabiError)*gradfact, avgdR);


xshade = [fullfreq, fliplr(fullfreq)]./1e6;

for n=1:1:5
    yshadeA{n} = [INTforceAhigh(:,n)', fliplr(INTforceAlow(:,n)')].*1e6;    
    yshadeB{n} = [INTforceBhigh(:,n)', fliplr(INTforceBlow(:,n)')].*1e6;    
end

%% plotting
forcemax = 3.5;
   %    Force/mg  thold  tof      grav  micro
distmax = forcemax*.5e-3*17.15e-3*9.81/1e-6;
% colors = ['mgrkb'];
colors = [[1,0,0];...%r
          [1,0,1];...%m
          [0,.875,0];...%g
          [0,.75,1];...%c
          [0,0,1]];    %b


% colors = ['rrrrr'];
afig=figure('position', [ 70    149   800   800]);
set(afig, 'Renderer', 'Painters');
set(afig,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
hold on;
for n=1:1:5
    %plot all data points
%     errorbar(triplefreqs, forceabove{n}(:,1), forceabove{n}(:,2), '*', 'linestyle', 'none', 'color', colors(n)) ;
    %plot averaged points
    yyaxis left;
    errorbar(monofreqs,avgforceabove(:,n)*1e6, avgaboveforceerror(:,n)*1e6, 'x', 'linestyle', 'none', 'color', colors(6-n,:) ) ;
    plot(fullfreq./1e6, INTforestA(:,n)*1e6, '-', 'color', colors(n,:));
    ylim([-distmax,distmax]);
    ylabel('\Delta y (\mu{}m)')

    ha{n} = fill(xshade,yshadeA{n}, colors(n,:),'EdgeColor','none');
    set(ha{n},'facealpha',.1);
    
    yyaxis right;
    plot(fullfreq./1e6, forestA(:,n), '-.', 'color', colors(n,:));
    ylim([-forcemax,forcemax]);
    ylabel('Force / mg')
    
%     fillobj = fill( xshade, yshade{n}, 'k', 'EdgeColor', 'none');
%     set(fillobj, 'facealpha', .2);
    
end

%dressed state labels
text(8.6,3,'|-->', 'fontsize', 24, 'fontname','cambria math', 'color', colors(1,:));
text(8.2,.8,'|->', 'fontsize', 24, 'fontname','cambria math', 'color', colors(2,:));
text(8.2,.35,'|0>', 'fontsize', 24, 'fontname','cambria math', 'color', colors(3,:));
text(8.2,-.6,'|+>', 'fontsize', 24, 'fontname','cambria math', 'color', colors(4,:));
text(8.6,-3,'|++>', 'fontsize', 24, 'fontname','cambria math', 'color', colors(5,:));
% bare state labels
text(9.7,1,'|+2,N-4>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');
text(9.9,-1,'|-2,N>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');
text(6.1,1,'|-2,N>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');
text(6.1,-1,'|+2,N-4>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');





% title('\delta_0 > 0');
title('From Above');
set(gca,'fontsize',20)
% ylabel('Force / mg')
xlabel('Applied RF frequency');
xlim([fullfreq(1),fullfreq(end)]./1e6);

% grid minor;
grid on;
% axis tight;



bfig = figure('position', [ 900    149   800   800]);
set(bfig, 'Renderer', 'Painters');
set(bfig,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
hold on;
for n=1:1:5
    %plot all data points
%     errorbar(triplefreqs, forceabove{n}(:,1), forceabove{n}(:,2), '*', 'linestyle', 'none', 'color', colors(n)) ;
    yyaxis left  ;  
    errorbar(monofreqs,avgforcebelow(:,n)*1e6, avgbelowforceerror(:,n)*1e6, 'x', 'linestyle', 'none', 'color', colors(n,:)) ;
    plot(fullfreq./1e6, INTforestB(:,n)*1e6, '-', 'color', colors(n,:));
    ylim([-distmax,distmax]);
    ylabel('\Delta y (\mu{}m)')
    
    hb{n} = fill(xshade,yshadeB{n}, colors(n,:),'EdgeColor','none');
    set(hb{n},'facealpha',.1);
        
    yyaxis right;
    plot(fullfreq./1e6, forestB(:,n), '-.', 'color', colors(n,:));
    ylim([-forcemax,forcemax]);  
    ylabel('Force / mg')    
    
%     fillobj = fill( xshade, yshade{n}, 'k', 'EdgeColor', 'none');
%     set(fillobj, 'facealpha', .2);
    
end

%dressed state labels 
text(8.6, 3,'|-->', 'fontsize', 24, 'fontname','cambria math', 'color', colors(1,:));
text(8.2, .8,'|->', 'fontsize', 24, 'fontname','cambria math', 'color', colors(2,:));
text(8.2,.35,'|0>', 'fontsize', 24, 'fontname','cambria math', 'color', colors(3,:));
text(8.2,-.6,'|+>', 'fontsize', 24, 'fontname','cambria math', 'color', colors(4,:));
text(8.6,-3,'|++>', 'fontsize', 24, 'fontname','cambria math', 'color', colors(5,:));
%bare labels
text(9.7,1,'|+2,N-4>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');
text(9.9,-1,'|-2,N>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');
text(6.1,1,'|-2,N>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');
text(6.1,-1,'|+2,N-4>', 'fontsize', 20, 'fontname','cambria math', 'color', 'k');


% title('\delta_0 < 0');
title('From Below');
set(gca,'fontsize',20)

xlabel('Applied RF frequency');
xlim([fullfreq(1),fullfreq(end)]./1e6);

% grid minor;
grid on;
% axis tight;


%% charlie's fill
yyaxis left


% end