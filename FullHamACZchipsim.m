function varargout = FullHamACZchipsim(varargin)
% this is a GUI to visualize AC Zeeman energy amplitude near chip wires.
% Microstrips are simulated with image / reverse currents
% AC skin available

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FullHamACZchipsim_OpeningFcn, ...
                   'gui_OutputFcn',  @FullHamACZchipsim_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function FullHamACZchipsim_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
linkaxes([handles.bplusaxes, handles.bminusaxes, handles.bxaxes, handles.byaxes, handles.tempaxes, handles.btotalaxes]);

%     axes(handles.bplusaxes);cla;
%     axes(handles.bminusaxes);cla;
%     axes(handles.btotalaxes);cla;
%     axes(handles.bxaxes);cla;
%     axes(handles.byaxes);cla;
%     axes(handles.tempaxes);cla;
%     axes(handles.fitx);cla;
%     axes(handles.fity);cla;
%     set(handles.tempreadout, 'string', '0');
%     set(handles.minx, 'string', '0');
%     set(handles.miny, 'string', '0');
%     set(handles.trapfreqtext, 'string', '0');
%     set(handles.timedisp, 'String', 'time');
%     set(handles.pmdisp, 'String', 'time');
%     set(handles.xydisp, 'String', 'time');
%     set(handles.tempdisp, 'String', 'time');
%     set(handles.freqdisp, 'String', 'time');
%     set(handles.totaldisp, 'String', 'time');

go_Callback(handles.go, eventdata, handles);
function varargout = FullHamACZchipsim_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

%% the Action
function go_Callback(~, ~, handles)

set(handles.timedisp, 'String', 'time');
set(handles.pmdisp, 'String', 'time');
set(handles.xydisp, 'String', 'time');
set(handles.tempdisp, 'String', 'time');
set(handles.freqdisp, 'String', 'time');
% set(handles.totaldisp, 'String', 'time');
set(handles.sumdisp, 'String', 'time');

tic; set(handles.working, 'visible', 'on');
drawnow;

whichtrans = get(handles.transtype, 'Value');
    
%% constants
gauss_tesla = 1e-4;         %conversion, * gauss to Tesla (divide for T->G)
% BSfactor = 4*pi*1e-7;       %Biot-Savart factor
BSfactor = 2*1e-7;       %Biot-Savart factor
hbar = 1.054571596e-34 ;     % h bar   J·s
h = hbar * 2 * pi;
massRb = 1.44316060e-25;  %mass of rubidium 87
mu_B = 9.274009994E-24;     %Bohr magneton J*T
muB  = h*1.399624604e6;      %J / G
g = -9.80665;               %m/s^2
k_b = 1.38064852e-23;	    %Boltzmann konstant  J/K
% k_bthreehalf = (3/2)*1.38064852e-23;	    %Boltzmann konstant  J/K
Ahfs = h*3.417341305452145e9;% 5s1/2
gs =  2.0023193043622;
gi = -0.0009951414;

%% gets
phase_left = deg2rad(str2double(get(handles.phaseL, 'string'))); % in degrees, out radians
phase_right = deg2rad(str2double(get(handles.phaseR, 'string'))); % in degrees, out radians
phase_middle = deg2rad(str2double(get(handles.phaseM, 'string'))); % in degrees, out radians
%currents written in milli amp, *1e-3  for Amps
I_left = str2double(get(handles.IL, 'string'))*1e-3*exp(1i*phase_left);    %amps
I_right = str2double(get(handles.IR, 'string'))*1e-3*exp(1i*phase_right);    %amps
I_middle = str2double(get(handles.IM, 'string'))*1e-3*exp(1i*phase_middle);     %amps
% Bz = str2double(get(handles.Bz, 'string'))*gauss_tesla;     %gauss -> tesla
thick = str2double(get(handles.thick, 'string'))*1e-6;  %micrometers
midsep= str2double(get(handles.midsep, 'string'))*1e-6;  %micrometers
wirewidth = str2double(get(handles.wirewidth, 'string'))*1e-6;  %micrometers
halfwidth = wirewidth/2;
% traceH = str2double(get(handles.tracethickness, 'string'))*1e-6;
quantdir = get(handles.quantdir, 'value');
Bdc = str2double(get(handles.Bdc, 'string'));
Bangle = deg2rad(str2double(get(handles.Bangle, 'string')));
wiretype = get(handles.wiretype, 'value');
frf = str2double(get(handles.detunetext, 'string'))*1e6;    % from MHz to Hz
% detuning = str2double(get(handles.detunetext, 'string'));
% detuning = detuning*1e6*2*pi;       % MHz -> Hz (Radians?)
cutfrac = str2double(get(handles.cutoff, 'string'));

booldoplot = get(handles.doplot, 'value');
boolaxes = get(handles.axescheck, 'value');
boolcolorbar = get(handles.colorbarcheck, 'value');
boolscale = get(handles.scalecheck, 'value');
boollog = get(handles.logcheck, 'value');
boollines = get(handles.linecheck, 'value');
boolconlabel = get(handles.conlabelcheck, 'value');
boolmin = get(handles.mincheck, 'value');
boolgravity = get(handles.gravitycheck, 'value');
boolfreqs = get(handles.freqcheck, 'value');
colorcap = str2double(get(handles.colorcap, 'string'));
booloutput = get(handles.outputcheck, 'value');
conspaceno = get(handles.conspace, 'value');
switch conspaceno
    case 1 
        conspace = 1;
    case 2
        conspace = 2;
    case 3 
        conspace = 5;
    case 4 
        conspace = 10;
    case 5 
        conspace = 20;
    case 6 
        conspace = 30;
    case 7 
        conspace = 50;
    case 8 
        conspace = 75;
    case 9
        conspace = 100;
    case 10 
        conspace = 150;
    case 11
        conspace = 200;
    case 12
        conspace = 300;
    case 13
        conspace = 500;
    case 14 
        conspace = 1000;
    case 15
        conspace = 0;
end

%% sizes
xnum = 1+str2double(get(handles.xnum, 'string'));
ynum = 1+str2double(get(handles.ynum, 'string'));
xmin = str2double(get(handles.xmin, 'string'));
ymin = str2double(get(handles.ymin, 'string'));
xmax = str2double(get(handles.xmax, 'string'));
ymax = str2double(get(handles.ymax, 'string'));

x = linspace(xmin*midsep, xmax*midsep, xnum);
y = linspace(ymin*thick, ymax*thick, ynum);
szx = size(x,2);
szy = size(y,2);

%% initialize
clear bx by;
Bx = zeros(szy*szx,1);
By = zeros(szy*szx,1);
Bplus = zeros(szy,szx);
Bminus = zeros(szy,szx);

xplot = x*1e6;      %unit on axis is micrometers
yplot = (y-thick)*1e6;

Bxplot = zeros(szy,szx);
Byplot = zeros(szy,szx);
Bplus = zeros(szy,szx);
Bminus = zeros(szy,szx);
Bplusplot = zeros(szy,szx);
Bminusplot = zeros(szy,szx);
Btot = zeros(szy,szx);
Btotplot = zeros(szy,szx);
Temp = zeros(szy,szx);
Etot = zeros(szy,szx);
gravity = repmat(massRb*g*y',1,szx);

posmat = zeros(szy*szx,2); % [ V xvector V, V yvector V ]
posmat(:,1) = sort(repmat(x', szy,1));
posmat(:,2) = repmat(y', szx,1);

%% B_AC Calculation
switch wiretype         % these all differ by only geometry: how we get B+- from currents
case 1      %3-wire

    Bx = (I_left)*((posmat(:,2)-thick) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)-thick).^2))...
            +(I_right)*((posmat(:,2)-thick) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)-thick).^2))...
           +(I_middle)*(+(posmat(:,2)-thick) ./ (((posmat(:,1)).^2+(posmat(:,2)-thick).^2)));
    Bx = (Bx*BSfactor);

    By = (I_left)*(-(posmat(:,1)+midsep) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)-thick).^2))...
            +(I_right)*(-(posmat(:,1)-midsep) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)-thick).^2))...
           +(I_middle)*(-(posmat(:,1))        ./ ((posmat(:,1)).^2+(posmat(:,2)-thick).^2));
    By = (By*BSfactor);

    Bx = reshape(Bx, szy, szx);
    By = reshape(By, szy, szx);

    Bxplot = log10( abs(Bx));
    Byplot = log10( abs(By)) ;
    Btot   = sqrt( abs( Bx ).^2+ abs( By ).^2 );
    Btotplot = log10(Btot);

    Bplus        =  ( Bx+1i.*By );
    Bminus       =  ( Bx-1i.*By );
    Bplusplot    =  log10( abs(Bplus));
    Bminusplot   =  log10( abs(Bminus));

case 2      %3micro

    Bx = (I_left)*(   ...
        -(posmat(:,2)+thick) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)+thick).^2)...
        +(posmat(:,2)-thick) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)-thick).^2)...
        )+...
    (I_right)*(...
        -(posmat(:,2)+thick) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)+thick).^2)...
        +(posmat(:,2)-thick) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)-thick).^2)...
        )+...
     (I_middle)*(...
        +(posmat(:,2)-thick) ./ ((posmat(:,1)).^2+(posmat(:,2)-thick).^2)...
        -(posmat(:,2)+thick) ./ ((posmat(:,1)).^2+(posmat(:,2)+thick).^2)...
        );
    Bx = Bx*BSfactor;

    By = (I_left)* (  ...
        +(posmat(:,1)+midsep) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)+thick).^2)...
        -(posmat(:,1)+midsep) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)-thick).^2)...
        )+...
        (I_right)*(...
        +(posmat(:,1)-midsep) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)+thick).^2)...
        -(posmat(:,1)-midsep) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)-thick).^2)...
        )+...
        (I_middle)*(...
        -(posmat(:,1)) ./ ((posmat(:,1)).^2+(posmat(:,2)-thick).^2)...
        +(posmat(:,1)) ./ ((posmat(:,1)).^2+(posmat(:,2)+thick).^2)...  
        );
    By = By*BSfactor;

    Bx = reshape(Bx, szy, szx);
    By = reshape(By, szy, szx);

    Bxplot = log10( abs(Bx));
    Byplot = log10( abs(By)) ;
    Btot   = sqrt( abs( Bx ).^2+ abs( By ).^2 );
%         Btot   = sqrt( abs( Bx )^2+ abs( By )^2 + abs(Bz)^2 );
    Btotplot = log10(Btot);

    Bplus        =  ( Bx+1i.*By );
    Bminus       =  ( Bx-1i.*By );
    Bplusplot    =  log10( abs(Bplus));
    Bminusplot   =  log10( abs(Bminus));

case 3      %AC Skin effect
acskingo_Callback(handles.acskingo, [], handles);

tracepositions = get(handles.verticalN, 'UserData');
ACskinI = get(handles.acskingo, 'UserData');
ACskVec = reshape(flipud(ACskinI), numel(ACskinI),1);
set(handles.elpertrace, 'String', num2str(numel(ACskVec))); drawnow;
trpos = zeros(numel(tracepositions{1})*numel(tracepositions{2}),2); % [ V xvector V, V yvector V ]
trpos(:,1) = sort(repmat(tracepositions{1}', numel(tracepositions{2}),1));
trpos(:,2) = repmat(tracepositions{2}, numel(tracepositions{1}),1);    

for n=1:1:numel(ACskVec)
    %leading minus sign  and plus thick-> image current
    Bx =Bx+...
    (I_left*ACskVec(n))*(   ...
        +(posmat(:,2)-thick-trpos(n,2)) ./ ((posmat(:,1)+midsep+trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        )+...
    (I_right*ACskVec(n))*(...
        +(posmat(:,2)-thick-trpos(n,2)) ./ ((posmat(:,1)-midsep-trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        )+...
    (I_middle*ACskVec(n))*(...
        +(posmat(:,2)-thick-trpos(n,2)) ./ ((posmat(:,1)-trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        );
    %leading plus sign  and plus thick-> image current
    By = By+...
    (I_left*ACskVec(n))* (  ...
        -(posmat(:,1)+midsep+trpos(n,1)) ./ ((posmat(:,1)+midsep+trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        )+...
    (I_right*ACskVec(n))*(...
        -(posmat(:,1)-midsep-trpos(n,1)) ./ ((posmat(:,1)-midsep-trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        )+...
    (I_middle*ACskVec(n))*(...
        -(posmat(:,1)+trpos(n,1)) ./ ((posmat(:,1)+trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        );
end

    Bx = Bx*BSfactor;
    By = By*BSfactor;      


    Bx = reshape(Bx, szy, szx);
    By = reshape(By, szy, szx);

    Bxplot = log10( abs(Bx));
    Byplot = log10( abs(By)) ;
    Btot   = sqrt( abs( Bx ).^2+ abs( By ).^2 );
    Btotplot = log10(Btot);

    Bplus        =  ( Bx+1i.*By );
    Bminus       =  ( Bx-1i.*By );
    Bplusplot    =  log10( abs(Bplus));
    Bminusplot   =  log10( abs(Bminus));

case 4      %AC Skin effect
acskingo_Callback(handles.acskingo, [], handles);

tracepositions = get(handles.verticalN, 'UserData');
ACskinI = get(handles.acskingo, 'UserData');
ACskVec = reshape(flipud(ACskinI), numel(ACskinI),1);
set(handles.elpertrace, 'String', num2str(numel(ACskVec))); drawnow;
trpos = zeros(numel(tracepositions{1})*numel(tracepositions{2}),2); % [ V xvector V, V yvector V ]
trpos(:,1) = sort(repmat(tracepositions{1}', numel(tracepositions{2}),1));
trpos(:,2) = repmat(tracepositions{2}, numel(tracepositions{1}),1);    

for n=1:1:numel(ACskVec)
    %leading minus sign  and plus thick-> image current
    Bx =Bx+...
    (I_left*ACskVec(n))*(   ...
        +(posmat(:,2)-thick-trpos(n,2)) ./ ((posmat(:,1)+midsep+trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        -(posmat(:,2)+thick+trpos(n,2)) ./ ((posmat(:,1)+midsep+trpos(n,1)).^2+(posmat(:,2)+thick+trpos(n,2)).^2)...
        )+...
    (I_right*ACskVec(n))*(...
        +(posmat(:,2)-thick-trpos(n,2)) ./ ((posmat(:,1)-midsep-trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        -(posmat(:,2)+thick+trpos(n,2)) ./ ((posmat(:,1)-midsep-trpos(n,1)).^2+(posmat(:,2)+thick+trpos(n,2)).^2)...
        )+...
    (I_middle*ACskVec(n))*(...
        +(posmat(:,2)-thick-trpos(n,2)) ./ ((posmat(:,1)-trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        -(posmat(:,2)+thick+trpos(n,2)) ./ ((posmat(:,1)-trpos(n,1)).^2+(posmat(:,2)+thick+trpos(n,2)).^2)...
        );
    %leading plus sign  and plus thick-> image current
    By = By+...
    (I_left*ACskVec(n))* (  ...
        +(posmat(:,1)+midsep+trpos(n,1)) ./ ((posmat(:,1)+midsep+trpos(n,1)).^2+(posmat(:,2)+thick+trpos(n,2)).^2)...
        -(posmat(:,1)+midsep+trpos(n,1)) ./ ((posmat(:,1)+midsep+trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        )+...
    (I_right*ACskVec(n))*(...
        +(posmat(:,1)-midsep-trpos(n,1)) ./ ((posmat(:,1)-midsep-trpos(n,1)).^2+(posmat(:,2)+thick+trpos(n,2)).^2)...
        -(posmat(:,1)-midsep-trpos(n,1)) ./ ((posmat(:,1)-midsep-trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        )+...
    (I_middle*ACskVec(n))*(...
        +(posmat(:,1)+trpos(n,1)) ./ ((posmat(:,1)+trpos(n,1)).^2+(posmat(:,2)+thick+trpos(n,2)).^2)...  
        -(posmat(:,1)+trpos(n,1)) ./ ((posmat(:,1)+trpos(n,1)).^2+(posmat(:,2)-thick-trpos(n,2)).^2)...
        );
end

    Bx = Bx*BSfactor;
    By = By*BSfactor;      


    Bx = reshape(Bx, szy, szx);
    By = reshape(By, szy, szx);

    Bxplot = log10( abs(Bx));
    Byplot = log10( abs(By)) ;
    Btot   = sqrt( abs( Bx ).^2+ abs( By ).^2 );
    Btotplot = log10(Btot);

    Bplus        =  ( Bx+1i.*By );
    Bminus       =  ( Bx-1i.*By );
    Bplusplot    =  log10( abs(Bplus));
    Bminusplot   =  log10( abs(Bminus));
    
end

%% Breit-Rabi / resonance / C.G. 

qnums = get(handles.transtype, 'UserData');
pluspiminus = qnums(whichtrans, 7);
    % == eight_E_cg_mFham_func(B).m
    % here so it's a single function, it had been calling 2
    % states
    magI = 3/2; vecI = magI:-1:-magI; densq = (2*magI+1);
    magS = 1/2; vecS = magS:-1:-magS;
    magF(2) = magI+magS;
    magF(1) = magI-magS;
    % Hyperfine 
    IdS(2) = .5*(magF(2)*(magF(2)+1) -magI*(magI+1) -magS*(magS+1));
    IdS(1) = .5*(magF(1)*(magF(1)+1) -magI*(magI+1) -magS*(magS+1));
    Hhfs(2) = Ahfs * IdS(2);
    Hhfs(1) = Ahfs * IdS(1);
    % (Hhfs2-Hhfs1)/h;      % 6.834...GHz
    Hhfs = diag([Hhfs(2),Hhfs(2),Hhfs(2),Hhfs(2),Hhfs(2), Hhfs(1),Hhfs(1),Hhfs(1)]);
    % Hhfs = [Hhfs2,Hhfs2,Hhfs2,Hhfs2,Hhfs2, Hhfs1,Hhfs1,Hhfs1]';
    % magnetic
    % construct B in the S, I basis
    % col 1 = i, col 2 = s 
    SIb = [[vecI,vecI];[magS,magS,magS,magS],-[magS,magS,magS,magS]]';
    Hisd = zeros(size(SIb,1),1);
    for s=1:1:size(SIb,1);    Hisd(s) = gi*SIb(s,1) + gs*SIb(s,2); end
    % Spin matrices
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

    % H setup
    B = Bdc;
    % HB = (muB/hbar)*(gs*ms + gi*mi)*B
    HisB = zeros(8,1);
    Hfdiag = zeros(8,8);
    Hoffup = zeros(8,8);
    Hoffdn = zeros(8,8);
    Hftot = zeros(8,8); %#ok<PREALL>
    H_GHz = zeros(8,1);
    E_Hz = zeros(8,1);
    E_Joules = zeros(8,1);
    V= zeros(8,8);
    Vr= zeros(8,8);
    D= zeros(8,8);
    % W= zeros(8,8);
    CGs = zeros(8,8);
    CGi = zeros(8,8);
    HisB = muB*Hisd.*B;

    % self-B terms from I-S, in "F" basis
    Hfdiag(1,1) = (2+2)/densq*HisB(1) + 0;
    Hfdiag(2,2) = (2+1)/densq*HisB(2) + (2-1)/densq*HisB(5);
    Hfdiag(3,3) = (2+0)/densq*HisB(3) + (2-0)/densq*HisB(6);
    Hfdiag(4,4) = (2-1)/densq*HisB(4) + (2+1)/densq*HisB(7);
    Hfdiag(5,5) = 0                   + (2+2)/densq*HisB(8);
    Hfdiag(6,6) = (2-1)/densq*HisB(2) + (2+1)/densq*HisB(5);
    Hfdiag(7,7) = (2+0)/densq*HisB(3) + (2+0)/densq*HisB(6);
    Hfdiag(8,8) = (2+1)/densq*HisB(4) + (2-1)/densq*HisB(7);

    %off- diagonal terms mixing F=1/F=2 from spin up : only +1/0/-1
    Hoffup(2,6)=-sqrt((2+1)*(2-1))/densq*HisB(2);
    Hoffup(3,7)=-sqrt((2+0)*(2-0))/densq*HisB(3);
    Hoffup(4,8)=-sqrt((2-1)*(2+1))/densq*HisB(4);
    Hoffup(6,2)=-sqrt((2+1)*(2-1))/densq*HisB(2);
    Hoffup(7,3)=-sqrt((2+0)*(2-0))/densq*HisB(3);
    Hoffup(8,4)=-sqrt((2-1)*(2+1))/densq*HisB(4);

    %off- diagonal terms mixing F=1/F=2 from spin down : only +1/0/-1
    Hoffdn(2,6)=sqrt((2-1)*(2+1))/densq*HisB(5);
    Hoffdn(3,7)=sqrt((2+0)*(2-0))/densq*HisB(6);
    Hoffdn(4,8)=sqrt((2+1)*(2-1))/densq*HisB(7);
    Hoffdn(6,2)=sqrt((2-1)*(2+1))/densq*HisB(5);
    Hoffdn(7,3)=sqrt((2+0)*(2-0))/densq*HisB(6);
    Hoffdn(8,4)=sqrt((2+1)*(2-1))/densq*HisB(7);

    Hftot = Hfdiag + Hoffup + Hoffdn + Hhfs;

    % B- eigenvalue ( former loop, one Bdc value here, not a Breit Rabi
    % like this was written for
    [V(:,:),D(:,:)] = eig(Hftot(:,:));        % eig can only work on 2-D matrices
    % Vr(:,:) = [V(:,8),V(:,7),V(:,6),V(:,5),V(:,4),V(:,1),V(:,2),V(:,3)];
    %render convenient signs
    % re-order into +2+10-1-2/+10-1 order (from low-to-high Energies) 
    Vr(:,:)=[   V(:,8)* sign(V(1,8)),...
                V(:,7)* sign(V(2,7)),...
                V(:,6)* sign(V(3,6)),...
                V(:,5)* sign(V(4,5)),...
                V(:,4)* sign(V(5,4)),...
                V(:,1)*-sign(V(6,1)),...
                V(:,2)*-sign(V(7,2)),...
                V(:,3)*-sign(V(8,3))];
    % Vr2(:,:) = Vr(:,:).^2;
    % Vr(Vr<1e-6) = 0
    H_GHz(:) = diag(D(:,:))/h/1e9;
    E_Hz(:) = diag(D(:,:))/h;
    E_Joules(:) = diag(D(:,:));
    % I think for these vector-calculation, should use the 2-d mats in the loop 
    % Not obvious how to do this operation across and indep. of B dimension
    % sigma plus RF
    CGs(1,2) = (Vr(:,1)'*SpM*Vr(:,2)); 
    CGs(2,3) = (Vr(:,2)'*SpM*Vr(:,3)); 
    CGs(3,4) = (Vr(:,3)'*SpM*Vr(:,4)); 
    CGs(4,5) = (Vr(:,4)'*SpM*Vr(:,5)); 
    % sigma plus microwave
    CGs(1,6) = (Vr(:,1)'*SpM*Vr(:,6));
    CGs(2,7) = (Vr(:,2)'*SpM*Vr(:,7));
    CGs(3,8) = (Vr(:,3)'*SpM*Vr(:,8));
    % sigma minus RF
    CGs(7,8) = (Vr(:,7)'*SmM*Vr(:,8));
    CGs(6,7) = (Vr(:,6)'*SmM*Vr(:,7));
    % sigma minus Microwave
    CGs(3,6) = (Vr(:,3)'*SmM*Vr(:,6));
    CGs(4,7) = (Vr(:,4)'*SmM*Vr(:,7));
    CGs(5,8) = (Vr(:,5)'*SmM*Vr(:,8));
    % pi trans Microwave
    CGs(2,6) = (Vr(:,2)'*SzM*Vr(:,6));
    CGs(3,7) = (Vr(:,3)'*SzM*Vr(:,7));
    CGs(4,8) = (Vr(:,4)'*SzM*Vr(:,8));
    % copy to the transpose
    CGs(:,:) = CGs(:,:) +CGs(:,:)';

    % sigma plus RF
    CGi(1,2) = (Vr(:,1)'*IpM*Vr(:,2));
    CGi(2,3) = (Vr(:,2)'*IpM*Vr(:,3));
    CGi(3,4) = (Vr(:,3)'*IpM*Vr(:,4));
    CGi(4,5) = (Vr(:,4)'*IpM*Vr(:,5));
    % sigma plus microwave
    CGi(1,6) = (Vr(:,1)'*IpM*Vr(:,6));
    CGi(2,7) = (Vr(:,2)'*IpM*Vr(:,7));
    CGi(3,8) = (Vr(:,3)'*IpM*Vr(:,8));
    % sigma minus RF
    CGi(7,8) = (Vr(:,7)'*ImM*Vr(:,8));
    CGi(6,7) = (Vr(:,6)'*ImM*Vr(:,7));
    % sigma minus Microwave
    CGi(3,6) = (Vr(:,3)'*ImM*Vr(:,6));
    CGi(4,7) = (Vr(:,4)'*ImM*Vr(:,7));
    CGi(5,8) = (Vr(:,5)'*ImM*Vr(:,8));
    % pi trans Microwave
    CGi(2,6) = (Vr(:,2)'*IzM*Vr(:,6));
    CGi(3,7) = (Vr(:,3)'*IzM*Vr(:,7));
    CGi(4,8) = (Vr(:,4)'*IzM*Vr(:,8));
    % copy to the transpose
    CGi(:,:) = CGi(:,:) +CGi(:,:)';
%     end

    CGsg = CGs.*gs;
    CGig = CGi.*gi;
    % CGg = CGig+CGsg; 
    CGg = abs(CGig+CGsg);
%     end
    % assign CG coefs of 21 trans and Energies of 8 states 
    % whichtrans order
    % 1|2,2> to |1,1>  ?+
    % 2|2,1> to |1,0> ?+
    % 3|2,0> to |1,-1> ?+
    % 4|2,1> to |1,1> ?
    % 5|2,0> to |1,0> ?
    % 6|2,-1> to |1,-1> ?
    % 7|2,0> to |1,1> ?-
    % 8|2,-1> to |1,0> ?-
    % 9|2,-2> to |1,-1> ?-
    %10RF= F=2
    %11|++>
    %12|+>
    %13|0>
    %14|->
    %15|- ->
    %16RF= F=1
    %17|+>
    %18|0>
    %19|->

    % but energy is ordered lowest to highest: 1+1 -> 2+2
    % CG coefs should match location in commented array (^~218), symmetric
    % index by whichtrans value, point to CG coef for transin reverse order
    % only needed for microwave quick 2-level calcs
    reCG = [CGg(6,1),CGg(7,2),CGg(8,3),...
            CGg(6,2),CGg(7,3),CGg(8,4),...
            CGg(6,3),CGg(7,4),CGg(8,5)];
    % RF comes in the hamiltonian > 2x2

    % need diff for resonance in microwaves
    reE=[E_Hz(8)-E_Hz(1), E_Hz(7)-E_Hz(2), E_Hz(6)-E_Hz(3), ...
         E_Hz(7)-E_Hz(1), E_Hz(6)-E_Hz(2), E_Hz(5)-E_Hz(3), ...
         E_Hz(6)-E_Hz(1), E_Hz(5)-E_Hz(2), E_Hz(4)-E_Hz(3)  ];

    detuning =0; 

%% Calculate Temp from Bxyz+- fields
if( whichtrans < 10 )       % 2-level micrwave transitions
    detuning = (frf - reE(whichtrans))*2*pi;
    cgcoeff = reCG(whichtrans);
    set(handles.cgcoeff, 'string',num2str(cgcoeff))
    Bminus = abs(Bminus);
    Bplus = abs(Bplus);
    if pluspiminus == 1 
        if(quantdir ==3)
            if(boolgravity) 
                Temp = (gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2)/k_b;
            else
                Temp =  (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2/k_b;
            end
        else
            Temp = 0;
        end
    elseif pluspiminus==0
        if(quantdir ==1)
            if(boolgravity)
                Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2)/k_b;
            else
                Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2/k_b;
            end
        elseif(quantdir ==2)
            if(boolgravity)
                Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2)/k_b;
            else
                Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2/k_b;
            end
        elseif(quantdir ==3)
            if(boolgravity)
                Temp = ( gravity );
            else
                Temp = zeros(szy,szx);
            end
        elseif(quantdir==4)
            if(boolgravity)
                Temp = abs(( gravity+((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar/2).^2))) *hbar)/k_b);
            else
                Temp = abs((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar).^2))) *hbar/2/k_b;
            end
        else
            Temp = 0;
        end
    elseif pluspiminus == -1
        if(quantdir ==3)
             if(boolgravity)
                Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2)/k_b;
             else
                Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2/k_b;
             end
        else 
            Temp = 0;
        end
    else
        Temp = 0;
        end
elseif( 10 < whichtrans && whichtrans < 16 )    % 5-level, F=2
%     Eacz = RFaczEnergy( Bdc, 2, Bplus, Bminus, detuning);

    set(handles.cgcoeff, 'string', 'RF');

    w2     = E_Hz(5)-E_Hz(4);       %
    delta0 = E_Hz(6)-E_Hz(5)-w2;    %
    delta1 = E_Hz(7)-E_Hz(5)-2*w2;  %
    delta2 = E_Hz(8)-E_Hz(5)-3*w2;  %
    
    detuning = frf - w2; 
    
    Rabi = -mu_B*Bminus;
    for n = 1:1:size(Bminus,1)
        for m = 1:1:size(Bminus,2)
            Rabiab = Rabi(n,m)*CGg(2,1)/h;          %Hz, +2 / +1
            Rabibc = Rabi(n,m)*CGg(3,2)/h;          %Hz  +1 /  0
            Rabicd = Rabi(n,m)*CGg(4,3)/h;          %Hz   0 / -1
            Rabide = Rabi(n,m)*CGg(5,4)/h;          %Hz  -1 / -2

           H = h*[2*detuning      Rabide/2        0               0                  0;...
                  conj(Rabide)/2  detuning        Rabicd/2        0                  0;...
                  0               conj(Rabicd)/2  delta0          Rabibc/2            0 ;...
                  0               0               conj(Rabibc)/2  -detuning+delta1   Rabiab/2;...
                  0               0               0                conj(Rabiab)/2     -2*detuning+delta2     ];
%             H = H*h;
            D = eig(H);                      %Calculate eigenvalue and eigenstate

            Eacz{2,5}(n,m) = D(1);      % lowest energy - high field seeking --
            Eacz{2,4}(n,m) = D(2);
            Eacz{2,3}(n,m) = D(3);
            Eacz{2,2}(n,m) = D(4);
            Eacz{2,1}(n,m) = D(5);      % highest energy - low field seeking ++
        end
    end

    if(boolgravity)
        Temp = (gravity + Eacz{2,whichtrans-10})/k_b;
    else
        Temp = Eacz{2,whichtrans-10}/k_b;
    end
elseif( 16 < whichtrans && whichtrans < 20 )    % 5-level, F=1
%     Eacz = RFaczEnergy( Bdc, 1, Bplus, Bminus, detuning);
    set(handles.cgcoeff, 'string', 'RF');
    
    w1    = E_Hz(2)-E_Hz(1);         
    delta = E_Hz(3)-E_Hz(2)-w1;       
    detuning = frf - w1;
    Rabi = -mu_B*Bplus;
    for n = 1:1:size(Bplus,1)
        for m = 1:1:size(Bplus,2)
            Rabifg = Rabi(n,m)*CGg(6,7)/h;          %Hz    +1 to 0
            Rabigh = Rabi(n,m)*CGg(7,8)/h;          %Hz    0 to -1

            H= h*[2*detuning        Rabifg/2            0 ;...
                  conj(Rabifg)/2    detuning           Rabigh/2;...
                  0                 conj(Rabigh)/2     delta     ];
            D = eig(H);                      %Calculate eigenvalue and eigenstate

            Eacz{1,3}(n,m) = D(1);      % lowest energy - high field seeking --
            Eacz{1,2}(n,m) = D(2);
            Eacz{1,1}(n,m) = D(3);      % highest energy - low field seeking ++
        end
    end

    if(boolgravity)
        Temp = (gravity + Eacz{1,whichtrans-16})/k_b;
    else
        Temp = Eacz{1,whichtrans-16}/k_b;
    end
end

set(handles.timedisp, 'String', num2str(round(toc,3)));tic;

%% find the local minimum in Temp
Temp = Temp-min(min(Temp));
Temp = Temp.*1e6;
minisfound = 0;
maxisfound = 0;
rowcounter = 1;
rowmin = zeros(1,szy);
rowmax = zeros(1,szy);
while(~minisfound)
        rowmin(rowcounter) = Temp(rowcounter,round((szx+1)/2));
    if rowcounter == szy
        mincol=1;
        break;
    elseif rowcounter==1 || rowmin(rowcounter)<rowmin(rowcounter-1)
        rowcounter = rowcounter+1;
    elseif rowmin(rowcounter)>rowmin(rowcounter-1)
        minval = rowmin(rowcounter-1);
        minrow = rowcounter-1;
        minisfound = 1;
    else
        mincol=1;
        break;
    end
end

while(~maxisfound)
    rowmax(rowcounter) = Temp(rowcounter,round((szx+1)/2));
    if rowcounter == szy
        break;
    elseif rowcounter==1 || rowmax(rowcounter)>rowmax(rowcounter-1)
        rowcounter = rowcounter+1;
    elseif rowmax(rowcounter)<rowmax(rowcounter-1)
        maxval = rowmax(rowcounter-1);
        maxrow = rowcounter-1;
        maxisfound = 1;
    else rowcounter == szy
        break;
    end
end

if minisfound&&maxisfound
    set(handles.tempreadout, 'String', num2str((maxval-minval)));
else
    set(handles.tempreadout, 'String', 'N/A');
end

[valxrow, locyrow] = min(Temp);
[valy, locx] = min(min(Temp));

valx = valxrow(locx);
locy = locyrow(locx);
locxy = [locx, locy];

set(handles.miny, 'String', num2str(yplot(locy)));
set(handles.minx, 'String', num2str(xplot(locx)));

%% trapping freqs
if(boolfreqs && maxisfound)
    
    warning('off','all');       %sorry
    
    axes(handles.bxaxes);cla;
    axes(handles.byaxes);cla;
    set(handles.bxaxes, 'visible', 'off');
    set(handles.byaxes, 'visible', 'off');
    set(handles.fitx, 'visible', 'on');
    set(handles.fity, 'visible', 'on');
%     cutfrac = .1;

    % setup from inputs
    Etot = Temp.*1e-6.*k_b;     % uK
    % take cross sections around min. 
    yslice = Etot(:,locx);
    xslice = Etot(locy,:)';

    % find inflection point in y
    inflecindexy = find(diff(sign(diff(diff(yslice)))))+3;  %+3 on index since index trims 1 each call
        %this is the wrong thing. I'd rather use the max after the dip
    Ecutoff = yslice(maxrow)*cutfrac;

    xposunder = x(xslice<Ecutoff)';
    yposunder = y(yslice<Ecutoff)';

    xsliceunder = xslice(xslice<Ecutoff);
    ysliceunder = yslice(yslice<Ecutoff);

    horiz_fit = fit(xposunder,xsliceunder,'poly2');
    vert_fit = fit(yposunder,ysliceunder,'poly2');

    horiz_param = coeffvalues(horiz_fit);
    vert_param = coeffvalues(vert_fit);

    w_x = sqrt(2*horiz_param(1)/massRb);
    w_y = sqrt(2*vert_param(1)/massRb);
    f_x = w_x/2/pi;
    f_y = w_y/2/pi;

    f_rb_avg = (f_x + f_y)/2;

    trapfreqround = round(f_rb_avg);
    
    axes(handles.fitx);cla;hold on; legend off;
    plot(xposunder,xsliceunder);    plot(horiz_fit);
    legend off; ylabel('joules'); xlabel('position, um');
    
    axes(handles.fity);cla;hold on; 
    plot(yposunder,ysliceunder);    plot(vert_fit);
    legend off; ylabel('joules'); xlabel('position, um');
    
    set(handles.trapfreqtext, 'string', num2str(trapfreqround));
        
else
    axes(handles.fitx);cla;
    axes(handles.fity);cla;
    set(handles.fitx, 'visible', 'off');
    set(handles.fity, 'visible', 'off');
    set(handles.bxaxes, 'visible', 'on');
    set(handles.byaxes, 'visible', 'on');
    set(handles.trapfreqtext, 'string', 'N/A');
end
set(handles.freqdisp, 'String', num2str(round(toc,3)));tic;

%% plotting 
if(boollines);style='-';else;style='none';end;

if(booldoplot)
colormap jet;

axes(handles.bplusaxes);  
if(boollog);contourf(xplot(:),yplot(:),Bplusplot(:,:),30, 'k','Linestyle', style);title('log10 B +');
else; contourf(xplot(:),yplot(:),abs(Bplus(:,:)),30, 'k', 'Linestyle', style);title('B +'); end;
axis equal;if(boolaxes); axis on; else; axis off; end;if(boolcolorbar); colorbar; end;

axes(handles.bminusaxes); 
if(boollog);contourf(xplot(:),yplot(:),Bminusplot(:,:),30, 'k','Linestyle', style);title('log10 B -');
else; contourf(xplot(:),yplot(:),abs(Bminus(:,:)),30, 'k','Linestyle', style);title('B -'); end;
axis equal;if(boolaxes); axis on; else; axis off; end;if(boolcolorbar); colorbar; end;

set(handles.pmdisp, 'String', num2str(round(toc,3)));tic;

axes(handles.bxaxes); 
if(~boolfreqs)
if(boollog);contourf(xplot(:),yplot(:),Bxplot(:,:),30, 'k','Linestyle', style);title('log10 Bx');
else; contourf(xplot(:),yplot(:),abs(Bx(:,:)),30, 'k','Linestyle', style);title('abs Bx'); end;
axis equal;if(boolaxes); axis on; else; axis off; end;if(boolcolorbar); colorbar; end;
end

axes(handles.byaxes); 
if(~boolfreqs)
if(boollog);contourf(xplot(:),yplot(:),Byplot(:,:),30, 'k','Linestyle', style);title('log10 By');
else; contourf(xplot(:),yplot(:),abs(By(:,:)),30, 'k','Linestyle', style);title('abs By'); end;
axis equal;if(boolaxes); axis on; else; axis off; end;if(boolcolorbar); colorbar; end;
end

set(handles.xydisp, 'String', num2str(round(toc,3)));tic;

axes(handles.btotalaxes); 
if(boollog);contourf(xplot(:),yplot(:),Btotplot(:,:),30, 'k','Linestyle', style);title('log10 B total');
else; contourf(xplot(:),yplot(:),Btot(:,:),30, 'k','Linestyle', style);title('B total'); end;
axis equal;if(boolaxes); axis on; else; axis off; end;if(boolcolorbar); colorbar; end;

else
    cla(handles.btotalaxes);     cla(handles.bxaxes);     cla(handles.byaxes);
    cla(handles.bminusaxes);     cla(handles.bplusaxes);     
end
% plot Temp - always do on plot? maybe not when gifing ? 
if conspace~=0
    if boolscale==0
        conspacing = [0:conspace:floor(max(max(Temp)))*conspace];
    else
        conspacing = [0:conspace:floor(min(colorcap, max(max(Temp)))/conspace)*conspace];
    end
elseif conspace==0
    if boolscale==0
        conspacing=linspace(0, max(max(Temp)), 30);
    else
        conspacing = linspace(0,min(colorcap, max(max(Temp))),30);
    end
end
axes(handles.tempaxes); 
% contourf(xplot,yplot,Temp, 'k', 'LevelList', linspace(0, round(min(colorcap, max(max(Temp)))), 30));title('Temp ( \muK)');
[C,v]=contourf(xplot(:),yplot(:),Temp(:,:), 'k','Linestyle', style, 'LevelList', conspacing);title('Temp ( \muK)');
axis equal;
if(boolaxes); axis on; else; axis off; end;
if(boolcolorbar); colorbar; end;
if(boolconlabel);clabel(C,v);end;
if(boolmin);text(xplot(locx),yplot(locy), 'X', 'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'color','r');end;

set(handles.tempdisp, 'String', num2str(round(toc,3)));
% drawnow;

%% output variables
if(booloutput)
%     assignin('base','xplot',xplot.*1e-6);
%     assignin('base','yplot',yplot.*1e-6);
    assignin('base','xplot',xplot);
    assignin('base','yplot',yplot);
    assignin('base','x',x);
    assignin('base','y',y);
    assignin('base','Bx',Bx);
    assignin('base','By',By);
    assignin('base','Bdc',Bdc);
    assignin('base','Btot',Btot);
    assignin('base','Temp',Temp);   
    assignin('base','E_Hz',E_Hz);   
    assignin('base','CGg',CGg);   
    assignin('base','Ener',Temp.*k_b);   
    assignin('base','Bplus',Bplus);
    assignin('base','Bminus',Bminus);
    assignin('base','logBx',Bxplot);
    assignin('base','logBy',Byplot);
    assignin('base','logBtot',Btotplot);
    assignin('base','logBplus',Bplusplot);
    assignin('base','logBminus',Bminusplot);
    assignin('base','locxy',locxy);
    if( 10 < whichtrans && whichtrans < 16 )    % 5-level, F=2
        assignin('base','Rabi', (mu_B * Bminus)/h);
        assignin('base','EACZ',Eacz{2,whichtrans-10}/h);  
    elseif( 16 < whichtrans && whichtrans < 20 )    % 5-level, F=1
        assignin('base','Rabi', ( mu_B * Bplus)/h );
        assignin('base','EACZ',Eacz{1,whichtrans-16}/h);
    end
end

%% toc time
timer(1) = max([0,str2double(get(handles.timedisp, 'string'))]);
timer(2) = max([0,str2double(get(handles.freqdisp, 'string'))]);
timer(3) = max([0,str2double(get(handles.pmdisp, 'string'))]);
timer(4) = max([0,str2double(get(handles.xydisp, 'string'))]);
% timer(5) = max([0,str2double(get(handles.totaldisp, 'string'))]);
timer(6) = max([0,str2double(get(handles.tempdisp, 'string'))]);
set(handles.sumdisp, 'String', num2str(round(sum(timer),3)));

%% reset button
set(handles.working, 'Visible', 'off');



%% parameters
function dphaseMdown_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phaseM, 'string', num2str(str2double(get(handles.phaseM, 'string'))-change));
if str2double(get(handles.phaseM, 'string')) > 360
    set(handles.phaseM, 'string', num2str(str2double(get(handles.phaseM, 'string'))-360));
elseif str2double(get(handles.phaseM, 'string')) < 0
    set(handles.phaseM, 'string', num2str(str2double(get(handles.phaseM, 'string'))+360));
end
go_Callback(handles.go, eventdata, handles);
function dphaseMup_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phaseM, 'string', num2str(str2double(get(handles.phaseM, 'string'))+change));
if str2double(get(handles.phaseM, 'string')) > 360
    set(handles.phaseM, 'string', num2str(str2double(get(handles.phaseM, 'string'))-360));
elseif str2double(get(handles.phaseM, 'string')) < 0
    set(handles.phaseM, 'string', num2str(str2double(get(handles.phaseM, 'string'))+360));
end
go_Callback(handles.go, eventdata, handles);
function phaseL_Callback(hObject, eventdata, handles)
function phaseL_CreateFcn(hObject, eventdata, handles)
function phaseR_Callback(hObject, eventdata, handles)
function phaseR_CreateFcn(hObject, eventdata, handles)
function dphaseRdown_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phaseR, 'string', num2str(str2double(get(handles.phaseR, 'string'))-change));
if str2double(get(handles.phaseR, 'string')) > 360
    set(handles.phaseR, 'string', num2str(str2double(get(handles.phaseR, 'string'))-360));
elseif str2double(get(handles.phaseR, 'string')) < 0
    set(handles.phaseR, 'string', num2str(str2double(get(handles.phaseR, 'string'))+360));
end
go_Callback(handles.go, eventdata, handles);
function dphaseLup_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phaseL, 'string', num2str(str2double(get(handles.phaseL, 'string'))+change));
if str2double(get(handles.phaseL, 'string')) > 360
    set(handles.phaseL, 'string', num2str(str2double(get(handles.phaseL, 'string'))-360));
elseif str2double(get(handles.phaseL, 'string')) < 0
    set(handles.phaseL, 'string', num2str(str2double(get(handles.phaseL, 'string'))+360));
end
go_Callback(handles.go, eventdata, handles);
function dphaseLdown_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phaseL, 'string', num2str(str2double(get(handles.phaseL, 'string'))-change));
if str2double(get(handles.phaseL, 'string')) > 360
    set(handles.phaseL, 'string', num2str(str2double(get(handles.phaseL, 'string'))-360));
elseif str2double(get(handles.phaseL, 'string')) < 0
    set(handles.phaseL, 'string', num2str(str2double(get(handles.phaseL, 'string'))+360));
end
go_Callback(handles.go, eventdata, handles);
function dphaseRup_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phaseR, 'string', num2str(str2double(get(handles.phaseR, 'string'))+change));
if str2double(get(handles.phaseR, 'string')) > 360
    set(handles.phaseR, 'string', num2str(str2double(get(handles.phaseR, 'string'))-360));
elseif str2double(get(handles.phaseR, 'string')) < 0
    set(handles.phaseR, 'string', num2str(str2double(get(handles.phaseR, 'string'))+360));
end
go_Callback(handles.go, eventdata, handles);
function phaseM_Callback(~, eventdata, handles)
go_Callback(handles.go, eventdata, handles);
function phaseM_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dphase_Callback(~, ~, ~)
function dphase_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function angledown_Callback(~, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.Bangle, 'string', num2str(str2double(get(handles.Bangle, 'string'))-change));
function angleup_Callback(~, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.Bangle, 'string', num2str(str2double(get(handles.Bangle, 'string'))+change));
function verticalN_Callback(hObject, eventdata, handles)
function verticalN_CreateFcn(hObject, eventdata, handles)
function dILdown_Callback(~, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IL, 'string', num2str(str2double(get(handles.IL, 'string'))-change));
% go_Callback(handles.go, eventdata, handles);
function dILup_Callback(~, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IL, 'string', num2str(str2double(get(handles.IL, 'string'))+change));
% go_Callback(handles.go, eventdata, handles);
function IL_Callback(~, eventdata, handles)
% go_Callback(handles.go, eventdata, handles);
function IL_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dIMdown_Callback(~, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IM, 'string', num2str(str2double(get(handles.IM, 'string'))-change));
go_Callback(handles.go, eventdata, handles);
function dIMup_Callback(~, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IM, 'string', num2str(str2double(get(handles.IM, 'string'))+change));
go_Callback(handles.go, eventdata, handles);
function IM_Callback(~, eventdata, handles)
% go_Callback(handles.go, eventdata, handles);
function IM_CreateFcn(hObject, ~, ~)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
function dIRdown_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IR, 'string', num2str(str2double(get(handles.IR, 'string'))-change));
% go_Callback(handles.go, eventdata, handles);
function dIRup_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IR, 'string', num2str(str2double(get(handles.IR, 'string'))+change));
% go_Callback(handles.go, eventdata, handles);
function IR_Callback(hObject, eventdata, handles)
% go_Callback(handles.go, eventdata, handles);
function IR_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Bdc_Callback(hObject, eventdata, handles)
function Bdc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dI_Callback(~, ~, ~)
function dI_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Bz_Callback(hObject, eventdata, handles)
function Bz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function transtype_Callback(hObject, eventdata, handles)
function transtype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function detunetext_Callback(hObject, eventdata, handles)
function detunetext_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function quantdir_Callback(hObject, eventdata, handles)
function quantdir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Bangle_Callback(hObject, eventdata, handles)
function Bangle_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function wiretype_Callback(hObject, eventdata, handles)
function wiretype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% geometry
function xmin_Callback(~, ~, ~)
function xmin_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xmax_Callback(~, ~, ~)
function xmax_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xnum_Callback(~, ~, ~)
function xnum_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ymin_Callback(~, ~, ~)
function ymin_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ymax_Callback(~, ~, ~)
function ymax_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ynum_Callback(~, ~, ~)
function ynum_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function thick_Callback(~, ~, ~)
function thick_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function midsep_Callback(~, ~, ~)
function midsep_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function wirewidth_Callback(~, ~, ~)
function wirewidth_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% AC skin Action
function acskingo_Callback(hObject, eventdata, handles)
%keep secret for now? 
wiretype = get(handles.wiretype, 'value');
freq = str2double(get(handles.detunetext, 'string'))*1e6;    % from MHz to Hz
% freq = str2double(get(handles.frequency, 'string')) ;     
width = str2double(get(handles.wirewidth, 'string'))./1e6 ;            
height = str2double(get(handles.tracethickness, 'string')) ./1e6;           
Nrow = str2double(get(handles.verticalN, 'string')) ;     %<-----set high for precision 
%warning you have to calculate Bfield for every current element. long time if Nrow  is large

cen2cen = 2*str2double(get(handles.thick, 'string'))./1e6 + height;

% constants 
Dmm = .4470491559036622;        %geometric mean of a square from itself
g = str2double(get(handles.conductivity, 'string')) ;       %copper
perm = 4*pi*1e-7;           
E0 = 1;             % voltage across, Volts. 

% creation
dy = height/Nrow;
dx = dy;                    %must be squares
Ncol = floor(width/2/dy);     %if too high, cut width or Nrow
N = Nrow*Ncol;

halfwidth = width/2;
halfheight = height/2;

Z = zeros(N);
X = zeros(N);
D = zeros(N);
ynmat = zeros(N);
ymmat = zeros(N);
Dother = zeros(N);
ydist = zeros(Nrow,1);
xdist = zeros(1,Ncol*2);
% curnt = cell(1,1);
% normcurrent= cell(1,1);
% secondnorm = zeros(1);
% totalcurrent = zeros(1,1);
zmax = 0;
% matrix population D,X
    const = freq*perm*g*dy*dx;

    for n=1:1:N             % n is main, m is others
            yn = ( 1/2 + floor((n-1)/Ncol) )*dy;
            xn = -(1/2 + mod(n-1,Ncol))*dx;
        for m=1:1:N
              ym = ( 1/2 + floor((m-1)/Ncol) )*dy;
              xm = -(1/2 + mod(m-1,Ncol))*dx;
            if(n==m)
                
                switch wiretype
                    case 3
                    X(n,n) = const*log(Dmm*dx* 2*xn );
                    case 4
                    X(n,n) = const*log(Dmm*dx* 2*xn * hypot(2*xn,cen2cen)* cen2cen );
                end
                ynmat(n,m) = yn;
                xnmat(n,m) = xn;
            else                %same quadrant (left)   opposite(right)     image partner       image+opposite
                switch wiretype
                    case 3
                    X(m,n) = const*log(hypot(xn-xm,yn-ym)*hypot(xn+xm,yn-ym));
                    case 4
                    X(m,n) = const*log(hypot(xn-xm,yn-ym)*hypot(xn+xm,yn-ym)*hypot(xn-xm,cen2cen)*hypot(xn+xm,yn-ym+cen2cen));
                end
            end
        end
    end
 
    
    
%%  Calculation
    Z = eye(N)-1i.*X;
    volts = g.*ones(N,1);
    thiscurrent = Z\volts;

    curnt = [fliplr(reshape(thiscurrent, Ncol,Nrow)'),reshape(thiscurrent, Ncol,Nrow)'];
    totalcurrent = sum(sum(abs(curnt)));
    normcurrent = curnt./totalcurrent;
   
    set(handles.acskingo, 'UserData', normcurrent);
    
    ynmat = unique(diag(ynmat));
    xnmat = [unique(diag(xnmat)') -fliplr(unique(diag(xnmat)'))];
    set(handles.verticalN, 'UserData', {xnmat, ynmat});
   
%% axes, plotting options 
function highres_Callback(hObject, eventdata, handles)
set( handles.xnum, 'string', '500');
set( handles.ynum, 'string', '500');
% go_Callback(handles.go, eventdata, handles);
function lowres_Callback(hObject, eventdata, handles)
set( handles.xnum, 'string', '50');
set( handles.ynum, 'string', '50');
% go_Callback(handles.go, eventdata, handles);

function center_Callback(hObject, eventdata, handles)
xmin = str2double(get(handles.xmin, 'string'));
xmax = str2double(get(handles.xmax, 'string'));

set(handles.xmin, 'string', num2str(-max(abs(xmin), abs(xmax))));
set(handles.xmax, 'string', num2str(max(abs(xmin), abs(xmax))));
% go_Callback(handles.go, eventdata, handles);
function tempaxes_CreateFcn(hObject, eventdata, handles)
function btotalaxes_CreateFcn(hObject, eventdata, handles)
function setROI_Callback(hObject, eventdata, handles)

midsep = str2double(get(handles.midsep, 'string'));
thick = str2double(get(handles.thick, 'string'));

Xlimits = get(handles.tempaxes, 'Xlim');
Ylimits = get(handles.tempaxes, 'Ylim');

set(handles.xmin, 'string', num2str(Xlimits(1)/midsep));
set(handles.xmax, 'string', num2str(Xlimits(2)/midsep));
set(handles.ymin, 'string', num2str(Ylimits(1)/thick));
set(handles.ymax, 'string', num2str(Ylimits(2)/thick));
% go_Callback(handles.go, eventdata, handles);
function colorcap_Callback(hObject, eventdata, handles)
function colorcap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function conspace_Callback(hObject, eventdata, handles)
function conspace_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% checks
function gravitycheck_Callback(hObject, eventdata, handles)
function conlabelcheck_Callback(hObject, eventdata, handles)
function linecheck_Callback(hObject, eventdata, handles)
function mincheck_Callback(hObject, eventdata, handles)
function axescheck_Callback(hObject, eventdata, handles)
function colorbarcheck_Callback(hObject, eventdata, handles)
function scalecheck_Callback(hObject, eventdata, handles)
function linkcheck_Callback(hObject, eventdata, handles)
function logcheck_Callback(hObject, eventdata, handles)
function outputcheck_Callback(hObject, eventdata, handles)
function freqcheck_Callback(hObject, eventdata, handles)
function cutoff_Callback(hObject, eventdata, handles)
function cutoff_CreateFcn(hObject, eventdata, handles)
function tracethickness_Callback(hObject, eventdata, handles)
function tracethickness_CreateFcn(hObject, eventdata, handles)
function conductivity_Callback(hObject, eventdata, handles)
function conductivity_CreateFcn(hObject, eventdata, handles)

%% Info

function cucondno_Callback(hObject, eventdata, handles)
function cucondno_CreateFcn(hObject, eventdata, handles)
function aucondno_Callback(hObject, eventdata, handles)
function aucondno_CreateFcn(hObject, eventdata, handles)
function agcondno_Callback(hObject, eventdata, handles)
function agcondno_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in doplot.
function doplot_Callback(hObject, eventdata, handles)
