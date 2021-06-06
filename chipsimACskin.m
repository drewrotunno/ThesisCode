function varargout = chipsimACskin(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chipsimACskin_OpeningFcn, ...
                   'gui_OutputFcn',  @chipsimACskin_OutputFcn, ...
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

function chipsimACskin_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
linkaxes([handles.bplusaxes, handles.bminusaxes, handles.bxaxes, handles.byaxes, handles.tempaxes, handles.tempaxes]);

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

function varargout = chipsimACskin_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

%% parameters
function dphasedown_Callback(~, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phase, 'string', num2str(str2double(get(handles.phase, 'string'))-change));
if str2double(get(handles.phase, 'string')) > 360
    set(handles.phase, 'string', num2str(str2double(get(handles.phase, 'string'))-360));
elseif str2double(get(handles.phase, 'string')) < 0
    set(handles.phase, 'string', num2str(str2double(get(handles.phase, 'string'))+360));
end
% go_Callback(handles.go, eventdata, handles);
function dphaseup_Callback(~, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.phase, 'string', num2str(str2double(get(handles.phase, 'string'))+change));
if str2double(get(handles.phase, 'string')) > 360
    set(handles.phase, 'string', num2str(str2double(get(handles.phase, 'string'))-360));
elseif str2double(get(handles.phase, 'string')) < 0
    set(handles.phase, 'string', num2str(str2double(get(handles.phase, 'string'))+360));
end
% go_Callback(handles.go, eventdata, handles);
function phase_Callback(~, ~, ~)
function phase_CreateFcn(hObject, ~, ~)
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
% go_Callback(handles.go, eventdata, handles);
function angleup_Callback(~, eventdata, handles)
change = str2double(get(handles.dphase, 'string'));
set(handles.Bangle, 'string', num2str(str2double(get(handles.Bangle, 'string'))+change));
% go_Callback(handles.go, eventdata, handles);

function verticalN_Callback(hObject, eventdata, handles)
function verticalN_CreateFcn(hObject, eventdata, handles)
function frequency_Callback(hObject, eventdata, handles)
function frequency_CreateFcn(hObject, eventdata, handles)


function dILdown_Callback(~, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IL, 'string', num2str(str2double(get(handles.IL, 'string'))-change));
go_Callback(handles.go, eventdata, handles);
function dILup_Callback(~, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IL, 'string', num2str(str2double(get(handles.IL, 'string'))+change));
go_Callback(handles.go, eventdata, handles);
function IL_Callback(~, ~, ~)
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
function IM_Callback(~, ~, ~)
function IM_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dIRdown_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IR, 'string', num2str(str2double(get(handles.IR, 'string'))-change));
go_Callback(handles.go, eventdata, handles);
function dIRup_Callback(hObject, eventdata, handles)
change = str2double(get(handles.dI, 'string'));
set(handles.IR, 'string', num2str(str2double(get(handles.IR, 'string'))+change));
go_Callback(handles.go, eventdata, handles);
function IR_Callback(hObject, eventdata, handles)
function IR_CreateFcn(hObject, eventdata, handles)
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
function detuneslide_Callback(hObject, eventdata, handles)
detunevalue = get(handles.detuneslide, 'value');
set(handles.detunetext, 'string', num2str(detunevalue));
function detuneslide_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function detunetext_Callback(hObject, eventdata, handles)
set(handles.detuneslide, 'value', str2double(get(handles.detunetext, 'string')));
% go_Callback(handles.go, eventdata, handles);
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

%% the action
function acskingo_Callback(hObject, eventdata, handles)
freq = str2double(get(handles.frequency, 'string')) ;     

width = str2double(get(handles.wirewidth, 'string'))./1e6 ;            
height = str2double(get(handles.tracethickness, 'string')) ./1e6;           
Nrow = str2double(get(handles.verticalN, 'string')) ;     %<-----set for precision

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
ydist = zeros(Nrow,numel(freq));
xdist = zeros(numel(freq),Ncol*2);
curnt = cell(numel(freq),1);
normcurrent= cell(numel(freq),1);
secondnorm = zeros(numel(freq));
totalcurrent = zeros(numel(freq),1);
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
                X(n,n) = const*log(Dmm*dx* 2*xn * hypot(2*xn,cen2cen)* cen2cen );
%                 D(n,n) = Dmm*dx;
%                 Dother(n,n) = 2*xn;
                ynmat(n,m) = yn;
%                 ymmat(n,m) = ym;
                xnmat(n,m) = xn;
%                 xmmat(n,m) = xm;
            else                %same quadrant (left)   opposite(right)  image partner       image+opposite
                X(m,n) = const*log(hypot(xn-xm,yn-ym)*hypot(xn+xm,yn-ym)*hypot(xn-xm,cen2cen)*hypot(xn+xm,yn-ym+cen2cen));
%                 D(m,n) = hypot(xn-xm,yn-ym);
%                 Dother(m,n) = hypot(xn+xm,yn-ym);
%                 ynmat(n,m) = yn;
%                 ymmat(n,m) = ym;
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

    
    
function go_Callback(~, ~, handles)

set(handles.timedisp, 'String', 'time');
set(handles.pmdisp, 'String', 'time');
set(handles.xydisp, 'String', 'time');
set(handles.tempdisp, 'String', 'time');
set(handles.freqdisp, 'String', 'time');
set(handles.totaldisp, 'String', 'time');

tic; set(handles.working, 'visible', 'on');
drawnow;

%% get C.G. numbers
whichtrans = get(handles.transtype, 'Value');
qnums = get(handles.transtype, 'UserData');

Fstart = qnums(whichtrans, 1);
mstart = qnums(whichtrans, 2);
Fprime = qnums(whichtrans, 3); %only for RF, F+ = 2 in others
mprime = qnums(whichtrans, 4);
Inuc = qnums(whichtrans, 5); % denom always 2I+1 = 4, here for completeness
eta = qnums(whichtrans, 6);
pluspiminus = qnums(whichtrans, 7);
cgcoeff = 100; % catch all, replaced in ifs

%% C.G. coeffeicients
if 1<= whichtrans&&whichtrans <=3
    cgcoeff = sqrt((2 + mstart)*(2 + mprime)) / 4;
elseif 4<= whichtrans&&whichtrans <=6
    cgcoeff = 2*sqrt((2 + mstart)*(2 - mprime)) / 4;
elseif 7<= whichtrans&&whichtrans <=9
    cgcoeff = sqrt((2 - mstart)*(2 - mprime)) / 4;
elseif whichtrans==10
    cgcoeff = 0;
elseif 11<= whichtrans &&whichtrans<=14
    cgcoeff = eta*sqrt((Fstart + mstart)*(Fstart - mprime)) / 4;
elseif 15<= whichtrans &&whichtrans<=19
    cgcoeff = eta*mstart / 4;
elseif 20<= whichtrans &&whichtrans<=23
    cgcoeff = eta*sqrt((Fstart - mstart)*(Fstart + mprime)) / 4;
elseif 24<= whichtrans &&whichtrans<=25
    cgcoeff = eta*sqrt((Fstart + mstart)*(Fstart - mprime)) / 4;
elseif 26<= whichtrans &&whichtrans<=28
    cgcoeff = eta*mstart / 4;
elseif 29<= whichtrans &&whichtrans<=31
    cgcoeff = eta*sqrt((Fstart - mstart)*(Fstart + mprime)) / 4;
else
    cgcoeff = 0;    %catchall, shouldn't ever be this
end

set(handles.cgcoeff, 'string', num2str(cgcoeff));
    
%% constants
gauss_tesla = 1e-4;         %conversion, * gauss to Tesla (divide for T->G)
% BSfactor = 4*pi*1e-7;       %Biot-Savart factor
BSfactor = 2*1e-7;       %Biot-Savart factor
massRb = 1.44316060e-25;  %mass of rubidium 87
mu_B = 9.274009994E-24;     %Bohr magneton J*T
g = -9.80665;               %m/s^2
k_b = 1.38064852e-23;	    %Boltzmann konstant  J/K
k_bthreehalf = (3/2)*1.38064852e-23;	    %Boltzmann konstant  J/K
hbar = 1.054571596e-34 ;     % h bar   J·s

%% gets
phase = deg2rad(str2double(get(handles.phase, 'string'))); % in degrees, out radians
I_left = str2double(get(handles.IL, 'string')) ;    %amps
I_right = str2double(get(handles.IR, 'string')) ;    %amps
I_middle = str2double(get(handles.IM, 'string'))*exp(i*phase);     %amps
% Bz = str2double(get(handles.Bz, 'string'))*gauss_tesla;     %gauss -> tesla
thick = str2double(get(handles.thick, 'string'))*1e-6;  %micrometers
midsep= str2double(get(handles.midsep, 'string'))*1e-6;  %micrometers
wirewidth = str2double(get(handles.wirewidth, 'string'))*1e-6;  %micrometers
halfwidth = wirewidth/2;
traceH = str2double(get(handles.tracethickness, 'string'))*1e-6;
quantdir = get(handles.quantdir, 'value');
Bangle = deg2rad(str2double(get(handles.Bangle, 'string')));
wiretype = get(handles.wiretype, 'value');
detuning = str2double(get(handles.detunetext, 'string'));
detuning = detuning*1e6*2*pi;       % MHz -> Hz
cutfrac = str2double(get(handles.cutoff, 'string'));

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
Tempplot = zeros(szy,szx);
Etot = zeros(szy,szx);
gravity = repmat(massRb*g*y',1,szx);

posmat = zeros(szy*szx,2); % [ V xvector V, V yvector V ]
posmat(:,1) = sort(repmat(x', szy,1));
posmat(:,2) = repmat(y', szx,1);


%% B Calculation

switch wiretype
case 1      %3-wire
    
        Bx = (I_left)*((posmat(:,2)-thick) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)-thick).^2))...
                +(I_right)*((posmat(:,2)-thick) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)-thick).^2))...
               +(I_middle)*(+(posmat(:,2)-thick) ./ (((posmat(:,1)).^2+(posmat(:,2)-thick).^2)));
        Bx = abs(Bx*BSfactor);

        By = (I_left)*(-(posmat(:,1)+midsep) ./ ((posmat(:,1)+midsep).^2+(posmat(:,2)-thick).^2))...
                +(I_right)*(-(posmat(:,1)-midsep) ./ ((posmat(:,1)-midsep).^2+(posmat(:,2)-thick).^2))...
               +(I_middle)*(-(posmat(:,1))        ./ ((posmat(:,1)).^2+(posmat(:,2)-thick).^2));
        By = abs(By*BSfactor);

        Bx = reshape(Bx, szy, szx);
        By = reshape(By, szy, szx);
    
        Bxplot = log10( abs(Bx));
        Byplot = log10( abs(By)) ;
        Btot   = sqrt( abs( Bx ).^2+ abs( By ).^2 );
        Btotplot = log10(Btot);

        Bplus        =  abs( Bx+1i.*By )./sqrt(2);
        Bminus       =  abs( Bx-1i.*By )./sqrt(2);
        Bplusplot    =  log10( Bplus);
        Bminusplot   =  log10( Bminus);

        if pluspiminus == 1 
            if(quantdir ==3)
                if(boolgravity)
                    Temp = (gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2)/k_b;
                else
                    Temp =  (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2/k_b;
                end
                    Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus==0
            if(quantdir ==1)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir ==2)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir==4)
                if(boolgravity)
                    Temp = abs(( gravity+((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar/2).^2))) *hbar)/k_b);
                else
                    Temp = abs((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar).^2))) *hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus == -1
            if(quantdir ==3)
                 if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2)/k_b;
                 else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2/k_b;
                 end
                Tempplot = log10(Temp);
            else 
                Temp = 0;
                Tempplot = log10(Temp);
            end
        else
            Temp = 0;
        end

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

        Bplus  =  abs( Bx+1i.*By )./sqrt(2);
        Bminus =  abs( Bx-1i.*By )./sqrt(2);
        Bplusplot  =  log10( Bplus);
        Bminusplot =  log10( Bminus);


        if pluspiminus == 1 
            if(quantdir ==3)
                if(boolgravity)
                    Temp = (gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2)/k_b;
                else
                    Temp =  (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2/k_b;
                end
                    Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus==0
            if(quantdir ==1)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir ==2)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir==4)
                if(boolgravity)
                    Temp = abs(( gravity+((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar/2).^2))) *hbar)/k_b);
                else
                    Temp = abs((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar).^2))) *hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus == -1
            if(quantdir ==3)
                 if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2)/k_b;
                 else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2/k_b;
                 end
                Tempplot = log10(Temp);
            else 
                Temp = 0;
                Tempplot = log10(Temp);
            end
        else
            Temp = 0;
        end

case 3      %split 3 wire

        Bx = (I_left/2)*((posmat(:,2)-thick) ./ ((posmat(:,1)+midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                           +(posmat(:,2)-thick) ./ ((posmat(:,1)+midsep-halfwidth).^2+(posmat(:,2)-thick).^2))...
                +(I_right/2)*((posmat(:,2)-thick) ./ ((posmat(:,1)-midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                           +(posmat(:,2)-thick) ./ ((posmat(:,1)-midsep-halfwidth).^2+(posmat(:,2)-thick).^2))...
               +(I_middle/2)*((posmat(:,2)-thick) ./ (((posmat(:,1)+halfwidth).^2+(posmat(:,2)-thick).^2))...
                           +(posmat(:,2)-thick) ./ (((posmat(:,1)-halfwidth).^2+(posmat(:,2)-thick).^2)));
        Bx = abs(Bx*BSfactor);

        By = (I_left/2)*(-(posmat(:,1)+midsep+halfwidth) ./ ((posmat(:,1)+midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)+midsep-halfwidth) ./ ((posmat(:,1)+midsep-halfwidth).^2+(posmat(:,2)-thick).^2))...
                +(I_right/2)*(-(posmat(:,1)-midsep+halfwidth) ./ ((posmat(:,1)-midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)-midsep-halfwidth) ./ ((posmat(:,1)-midsep-halfwidth).^2+(posmat(:,2)-thick).^2))...
               +(I_middle/2)*(-(posmat(:,1)+halfwidth)        ./ ((posmat(:,1)+halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)-halfwidth)        ./ ((posmat(:,1)-halfwidth).^2+(posmat(:,2)-thick).^2));
        By = abs(By*BSfactor);

        Bx = reshape(Bx, szy, szx);
        By = reshape(By, szy, szx);
        
        Bxplot = log10( abs(Bx));
        Byplot = log10( abs(By)) ;
        Btot   = sqrt( abs( Bx ).^2+ abs( By ).^2 );
%         Btot   = sqrt( abs( Bx )^2+ abs( By )^2 + abs(Bz)^2 );
        Btotplot = log10(Btot);

        Bplus  =  abs( Bx+1i.*By )./sqrt(2);
        Bminus =  abs( Bx-1i.*By )./sqrt(2);
        Bplusplot  =  log10( Bplus);
        Bminusplot =  log10( Bminus);

        if pluspiminus == 1 
            if(quantdir ==3)
                if(boolgravity)
                    Temp = (gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2)/k_b;
                else
                    Temp =  (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2/k_b;
                end
                    Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus==0
            if(quantdir ==1)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir ==2)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir==4)
                if(boolgravity)
                    Temp = abs(( gravity+((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar/2).^2))) *hbar)/k_b);
                else
                    Temp = abs((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar).^2))) *hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus == -1
            if(quantdir ==3)
                 if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2)/k_b;
                 else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2/k_b;
                 end
                Tempplot = log10(Temp);
            else 
                Temp = 0;
                Tempplot = log10(Temp);
            end
        else
            Temp = 0;
        end
    
case 4      %split 3-micro

        Bx = (I_left/2)*((posmat(:,2)-thick) ./ ((posmat(:,1)+midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                           +(posmat(:,2)-thick) ./ ((posmat(:,1)+midsep-halfwidth).^2+(posmat(:,2)-thick).^2)...
                           +(posmat(:,2)+thick) ./ ((posmat(:,1)+midsep+halfwidth).^2+(posmat(:,2)+thick).^2)...
                           +(posmat(:,2)+thick) ./ ((posmat(:,1)+midsep-halfwidth).^2+(posmat(:,2)+thick).^2))...
                +(I_right/2)*((posmat(:,2)-thick) ./ ((posmat(:,1)-midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                           +(posmat(:,2)-thick) ./ ((posmat(:,1)-midsep-halfwidth).^2+(posmat(:,2)-thick).^2)...
                           +(posmat(:,2)+thick) ./ ((posmat(:,1)-midsep+halfwidth).^2+(posmat(:,2)+thick).^2)...
                           +(posmat(:,2)+thick) ./ ((posmat(:,1)-midsep-halfwidth).^2+(posmat(:,2)+thick).^2))...
               +(I_middle/2)*((posmat(:,2)-thick) ./ (((posmat(:,1)+halfwidth).^2+(posmat(:,2)-thick).^2))...
                           +(posmat(:,2)-thick) ./ (((posmat(:,1)-halfwidth).^2+(posmat(:,2)-thick).^2))...
                           +(posmat(:,2)+thick) ./ (((posmat(:,1)+halfwidth).^2+(posmat(:,2)+thick).^2))...
                           +(posmat(:,2)+thick) ./ (((posmat(:,1)-halfwidth).^2+(posmat(:,2)+thick).^2)));
        Bx = abs(Bx*BSfactor);

        By = (I_left/2)*(-(posmat(:,1)+midsep+halfwidth) ./ ((posmat(:,1)+midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)+midsep-halfwidth) ./ ((posmat(:,1)+midsep-halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)+midsep+halfwidth) ./ ((posmat(:,1)+midsep+halfwidth).^2+(posmat(:,2)+thick).^2)...
                            -(posmat(:,1)+midsep-halfwidth) ./ ((posmat(:,1)+midsep-halfwidth).^2+(posmat(:,2)+thick).^2))...
                +(I_right/2)*(-(posmat(:,1)-midsep+halfwidth) ./ ((posmat(:,1)-midsep+halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)-midsep-halfwidth) ./ ((posmat(:,1)-midsep-halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)-midsep+halfwidth) ./ ((posmat(:,1)-midsep+halfwidth).^2+(posmat(:,2)+thick).^2)...
                            -(posmat(:,1)-midsep-halfwidth) ./ ((posmat(:,1)-midsep-halfwidth).^2+(posmat(:,2)+thick).^2))...
               +(I_middle/2)*(-(posmat(:,1)+halfwidth)        ./ ((posmat(:,1)+halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)-halfwidth)        ./ ((posmat(:,1)-halfwidth).^2+(posmat(:,2)-thick).^2)...
                            -(posmat(:,1)+halfwidth)        ./ ((posmat(:,1)+halfwidth).^2+(posmat(:,2)+thick).^2)...
                            -(posmat(:,1)-halfwidth)        ./ ((posmat(:,1)-halfwidth).^2+(posmat(:,2)+thick).^2));
        By = abs(By*BSfactor);

        Bx = reshape(Bx, szy, szx);
        By = reshape(By, szy, szx);
        
        Bxplot = log10( abs(Bx));
        Byplot = log10( abs(By)) ;
        Btot   = sqrt( abs( Bx ).^2+ abs( By ).^2 );
%         Btot   = sqrt( abs( Bx )^2+ abs( By )^2 + abs(Bz)^2 );
        Btotplot = log10(Btot);

        Bplus  =  abs( Bx+1i.*By )./sqrt(2);
        Bminus =  abs( Bx-1i.*By )./sqrt(2);
        Bplusplot  =  log10( Bplus);
        Bminusplot =  log10( Bminus);

        if pluspiminus == 1 
            if(quantdir ==3)
                if(boolgravity)
                    Temp = (gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2)/k_b;
                else
                    Temp =  (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2/k_b;
                end
                    Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus==0
            if(quantdir ==1)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir ==2)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir==4)
                if(boolgravity)
                    Temp = abs(( gravity+((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar/2).^2))) *hbar)/k_b);
                else
                    Temp = abs((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar).^2))) *hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus == -1
            if(quantdir ==3)
                 if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2)/k_b;
                 else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2/k_b;
                 end
                Tempplot = log10(Temp);
            else 
                Temp = 0;
                Tempplot = log10(Temp);
            end
        else
            Temp = 0;
        end

case 5      %AC Skin effect
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

        Bplus        =  abs( Bx+1i.*By )./sqrt(2);
        Bminus       =  abs( Bx-1i.*By )./sqrt(2);
        Bplusplot    =  log10( Bplus);
        Bminusplot   =  log10( Bminus);

        if pluspiminus == 1 
            if(quantdir ==3)
                if(boolgravity)
                    Temp = (gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2)/k_b;
                else
                    Temp =  (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2/k_b;
                end
                    Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus==0
            if(quantdir ==1)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir ==2)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir==4)
                if(boolgravity)
                    Temp = abs(( gravity+((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar/2).^2))) *hbar)/k_b);
                else
                    Temp = abs((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar).^2))) *hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus == -1
            if(quantdir ==3)
                 if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2)/k_b;
                 else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2/k_b;
                 end
                Tempplot = log10(Temp);
            else 
                Temp = 0;
                Tempplot = log10(Temp);
            end
        else
            Temp = 0;
        end

        
case 6      %AC Skin effect
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

        Bplus        =  abs( Bx+1i.*By )./sqrt(2);
        Bminus       =  abs( Bx-1i.*By )./sqrt(2);
        Bplusplot    =  log10( Bplus);
        Bminusplot   =  log10( Bminus);

        if pluspiminus == 1 
            if(quantdir ==3)
                if(boolgravity)
                    Temp = (gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2)/k_b;
                else
                    Temp =  (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bplus/hbar).^2))*hbar/2/k_b;
                end
                    Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus==0
            if(quantdir ==1)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bx/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir ==2)
                if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2)/k_b;
                else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*By/hbar).^2))*hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            elseif(quantdir==4)
                if(boolgravity)
                    Temp = abs(( gravity+((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar/2).^2))) *hbar)/k_b);
                else
                    Temp = abs((-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*cos(Bangle)*Bx/hbar/2).^2)) + (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*sin(Bangle)*By/hbar).^2))) *hbar/2/k_b;
                end
                Tempplot = log10(Temp);
            else
                Temp = 0;
                Tempplot = log10(Temp);
            end
        elseif pluspiminus == -1
            if(quantdir ==3)
                 if(boolgravity)
                    Temp = ( gravity+(-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2)/k_b;
                 else
                    Temp = (-abs(detuning) + sqrt(detuning^2 + (mu_B*cgcoeff*Bminus/hbar).^2))*hbar/2/k_b;
                 end
                Tempplot = log10(Temp);
            else 
                Temp = 0;
                Tempplot = log10(Temp);
            end
        else
            Temp = 0;
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
    elseif rowcounter==1 | rowmin(rowcounter)<rowmin(rowcounter-1)
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
if(boolfreqs)
    
    warning('off','all');       %sorry
    
    axes(handles.bxaxes);cla;
    axes(handles.byaxes);cla;
    set(handles.bxaxes, 'visible', 'off');
    set(handles.byaxes, 'visible', 'off');
    set(handles.fitx, 'visible', 'on');
    set(handles.fity, 'visible', 'on');
%     cutfrac = .1;

    % setup from inputs
    Etot = Temp.*1e-6.*k_b;
    % take cross sections around min. 
    yslice = Etot(:,locx);
    xslice = Etot(locy,:)';

    % find inflection point in y
    inflecindexy = find(diff(sign(diff(diff(yslice)))))+3;  %+3 on index since index trims 1 each call
    Ecutoff = yslice(inflecindexy)*cutfrac;

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
end
set(handles.freqdisp, 'String', num2str(round(toc,3)));tic;

%% plotting 
colormap jet;
% colorbar;
if(boollines);style='-';else;style='none';end;

axes(handles.bplusaxes);  
if(boollog);contourf(xplot(:),yplot(:),Bplusplot(:,:),30, 'k','Linestyle', style);title('log10 B +');
else; contourf(xplot(:),yplot(:),Bplus(:,:),30, 'k', 'Linestyle', style);title('B +'); end;
axis equal;if(boolaxes); axis on; else; axis off; end;if(boolcolorbar); colorbar; end;

axes(handles.bminusaxes); 
if(boollog);contourf(xplot(:),yplot(:),Bminusplot(:,:),30, 'k','Linestyle', style);title('log10 B -');
else; contourf(xplot(:),yplot(:),Bminus(:,:),30, 'k','Linestyle', style);title('B -'); end;
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

drawnow;

%% output variables
if(booloutput)
    assignin('base','xplot',xplot.*1e-6);
    assignin('base','yplot',yplot.*1e-6);
    assignin('base','Bx',Bx);
    assignin('base','By',By);
    assignin('base','Btot',Btot);
    assignin('base','Temp',Temp);   
    assignin('base','Ener',Temp.*k_b);   
    assignin('base','Bplus',Bplus);
    assignin('base','Bminus',Bminus);
    assignin('base','logBx',Bxplot);
    assignin('base','logBy',Byplot);
    assignin('base','logBtot',Btotplot);
    assignin('base','logTemp',Tempplot);
    assignin('base','logBplus',Bplusplot);
    assignin('base','logBminus',Bminusplot);
    assignin('base','locxy',locxy);
end

%% reset button
set(handles.working, 'Visible', 'off');
set(handles.totaldisp, 'String', num2str(toc));


%% axes, plotting options 
function highres_Callback(hObject, eventdata, handles)
set( handles.xnum, 'string', '500');
set( handles.ynum, 'string', '500');
% go_Callback(handles.go, eventdata, handles);
function lowres_Callback(hObject, eventdata, handles)
set( handles.xnum, 'string', '50');
set( handles.ynum, 'string', '50');
% go_Callback(handles.go, eventdata, handles);
function resetroi_Callback(hObject, eventdata, handles)

set(handles.xmin, 'string', num2str(-2));
set(handles.xmax, 'string', num2str(2));
set(handles.ymin, 'string', num2str(.5));
set(handles.ymax, 'string', num2str(6));
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



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)

function edit35_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
