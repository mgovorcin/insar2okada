function varargout = insar2okada(varargin)
% INSAR2OKADA MATLAB code for insar2okada.fig
%      INSAR2OKADA, by itself, creates a new INSAR2OKADA or raises the existing
%      singleton*.
%
%      H = INSAR2OKADA returns the handle to a new INSAR2OKADA or the handle to
%      the existing singleton*.
%
%      INSAR2OKADA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSAR2OKADA.M with the given input arguments.
%
%      INSAR2OKADA('Property','Value',...) creates a new INSAR2OKADA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before insar2okada_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to insar2okada_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help insar2okada

% Last Modified by GUIDE v2.5 12-Jul-2019 09:18:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @insar2okada_OpeningFcn, ...
                   'gui_OutputFcn',  @insar2okada_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before insar2okada is made visible.
function insar2okada_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to insar2okada (see VARARGIN)

% Choose default command line output for insar2okada
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes insar2okada wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = insar2okada_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in insardata.
function insardata_Callback(hObject, eventdata, handles)
% hObject    handle to insardata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat');
insardata = load([pathname filename]);
handles.insardata = insardata;
guidata(hObject, handles)
fprintf(['INSAR Data:' [pathname filename] '\n'])

function lonMin_Callback(hObject, eventdata, handles)
% hObject    handle to lonMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lonMin as text
%        str2double(get(hObject,'String')) returns contents of lonMin as a double
handles.LonMin = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function lonMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lonMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lonMax_Callback(hObject, eventdata, handles)
% hObject    handle to lonMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lonMax as text
%        str2double(get(hObject,'String')) returns contents of lonMax as a double
handles.LonMax = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function lonMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lonMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LatMin_Callback(hObject, eventdata, handles)
% hObject    handle to LatMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LatMin as text
%        str2double(get(hObject,'String')) returns contents of LatMin as a double
handles.LatMin = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function LatMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LatMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function latMax_Callback(hObject, eventdata, handles)
% hObject    handle to latMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latMax as text
%        str2double(get(hObject,'String')) returns contents of latMax as a double
handles.LatMax = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function latMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in insarcrop.
function insarcrop_Callback(hObject, eventdata, handles)
% hObject    handle to insarcrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

insardata = handles.insardata;
if isnan(handles.LonMin) || isnan(handles.LonMax) || isnan(handles.LatMin) || isnan(handles.LatMax)
    fprintf('ERROR!! Define geo boundary box parameters \n')
else
LonMin = handles.LonMin;
LonMax = handles.LonMax;
LatMin = handles.LatMin;
LatMax = handles.LatMax;

insardata1 = insardata;

% Apply bounding box and remove data points outside the AOI
    iOutBox = find(insardata.Lon > LonMax | insardata.Lon < LonMin | insardata.Lat < LatMin | insardata.Lat > LatMax);
    if sum(iOutBox)>0
        insardata1.Phase(iOutBox) = [];
        insardata1.Lat(iOutBox) = [];
        insardata1.Lon(iOutBox) = [];
        insardata1.Heading(iOutBox) = [];
        insardata1.Inc(iOutBox) = [];
    end
handles.insarcrop = insardata1;
end
guidata(hObject, handles);
fprintf(['Finished cropping InSAR Data\n'])


function lambdaInput_Callback(hObject, eventdata, handles)
% hObject    handle to lambdaInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambdaInput as text
%        str2double(get(hObject,'String')) returns contents of lambdaInput as a double
handles.lambda = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function lambdaInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambdaInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function georefLon_Callback(hObject, eventdata, handles)
% hObject    handle to georefLon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of georefLon as text
%        str2double(get(hObject,'String')) returns contents of georefLon as a double
handles.georefLon = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function georefLon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to georefLon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function georefLat_Callback(hObject, eventdata, handles)
% hObject    handle to georefLat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of georefLat as text
%        str2double(get(hObject,'String')) returns contents of georefLat as a double
handles.georefLat = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function georefLat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to georefLat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in insarConvert.
function insarConvert_Callback(hObject, eventdata, handles)
% hObject    handle to insarConvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

insardata = handles.insarcrop;
%insardata = handles.insardata;

if isnan(handles.georefLon) || isnan(handles.georefLat);
  georeferencePoint = [ (min(insardata.Lon)+max(insardata.Lon))/2  (min(insardata.Lat)+max(insardata.Lat))/2 ]
  fprintf('Center of InSAR data is taken as geoReferencePoint for converstion \n')
else
georeferencePoint = [handles.georefLon handles.georefLat];
end

ll = [single(insardata.Lon) single(insardata.Lat)];   % Create Longitude and Latitude 2Xn matrix
    xy = llh2local(ll', georeferencePoint);    % Transform from geografic to local coordinates
    
    nPointsThis = size(ll, 1);   % Calculate length of current InSAR data vector
    xy = double([(1:nPointsThis)', xy'*1000]);   % Add ID number column to xy matrix with local coordinates
handles.xy = xy;
guidata(hObject, handles)
fprintf(['InSAR Data converted to local Cartesian coordinate system\n'])

% --- Executes on button press in insarPlot.
function insarPlot_Callback(hObject, eventdata, handles)
% hObject    handle to insarPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wavelength = handles.lambda;
xy = handles.xy;
insardata = handles.insarcrop;
phase = insardata.Phase;
convertedPhase = (phase / (4*pi) )  * wavelength;    % Convert phase from radians to m
los = single(-convertedPhase);  % Convert to Line-of-sigth displacement in m
handles.dlos = los;
guidata(hObject, handles)

cmap.seismo = colormap_cpt('GMT_seis.cpt', 100);    % GMT 'Seismo' colormap for wrapped data
cmap.redToBlue = colormap_cpt('polar.cpt', 100);    % Red to Blue colormap for unwrapped data
flag = get(handles.wrapped_flag,'Value');

axes(handles.insarplot);
cla;
 plotInsarW(xy, los, wavelength, cmap,flag);
 
    function plotInsarW(xy, los, wavelength, cmap,flag)
        edge = round(min(abs(diff(xy(:,3)))))+2; % Size of patch set to minumum distance between points
    if edge < 50
        edge = 50;
    end
    xs = [xy(:,2)'; xy(:,2)'+edge; xy(:,2)'+edge; xy(:,2)']; % Coordinates of four vertex of patch
    ys = [xy(:,3)'; xy(:,3)'; xy(:,3)'+edge; xy(:,3)'+edge];
     
   
    % Display wrapped interferogram at 5.6 cm wavelength 
    h1 = patch(xs, ys, 'r');
    
    if flag == 1
    colormap(cmap.seismo);
    set(h1, 'facevertexcdata', mod(los,wavelength/2), 'facecolor', 'flat', 'edgecolor', 'none')
    caxis([0 0.025])
    t = title('Wrapped InSAR Data','FontSize', 14);
    else
    colormap(cmap.redToBlue);
    set(h1, 'facevertexcdata', los, 'facecolor', 'flat', 'edgecolor', 'none')
    c = max(abs([min(los), max(los)])); % Calculate maximu value for symmetric colormap
    caxis([-c c])
    t = title('Unwrapped InSAR Data','FontSize', 14);
    end
    axis equal; axis tight;
    ax = gca;
    grid on
    ax.Layer = 'top';
    ax.Box = 'on';
    ax.LineWidth = 1.0;
    ax.GridLineStyle = '--';
    cbar = colorbar; ylabel(cbar,'Line-of-sight displacement m','FontSize', 10); 
    xlabel('X distance from local origin (m)','FontSize', 10)
    ylabel('Y distance from local origin (m)','FontSize', 10)
    
    %set(t,'Position',get(t,'Position')+[0 500 0]);
    drawnow



function edit_length_Callback(hObject, eventdata, handles)
% hObject    handle to edit_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_length as text
%        str2double(get(hObject,'String')) returns contents of edit_length as a double
handles.okadaLength = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_width_Callback(hObject, eventdata, handles)
% hObject    handle to edit_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_width as text
%        str2double(get(hObject,'String')) returns contents of edit_width as a double
handles.okadaWidth = str2double(get(hObject,'String'));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_depth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_depth as text
%        str2double(get(hObject,'String')) returns contents of edit_depth as a double
handles.okadaDepth = str2double(get(hObject,'String'));
set(handles.Gmidpointx,'String',num2str(0));
set(handles.Gmidpointy,'String',num2str(0));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_strike_Callback(hObject, eventdata, handles)
% hObject    handle to edit_strike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_strike as text
%        str2double(get(hObject,'String')) returns contents of edit_strike as a double
handles.okadaStrike = str2double(get(hObject,'String'));
set(handles.Gmidpointx,'String',num2str(0));
set(handles.Gmidpointy,'String',num2str(0));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_strike_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_strike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dip_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dip as text
%        str2double(get(hObject,'String')) returns contents of edit_dip as a double
handles.okadaDip = str2double(get(hObject,'String'));
set(handles.Gmidpointx,'String',num2str(0));
set(handles.Gmidpointy,'String',num2str(0));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_strikeslip_Callback(hObject, eventdata, handles)
% hObject    handle to edit_strikeslip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_strikeslip as text
%        str2double(get(hObject,'String')) returns contents of edit_strikeslip as a double
handles.okadaStrikeslip = str2double(get(hObject,'String'));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_strikeslip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_strikeslip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dipslip_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dipslip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dipslip as text
%        str2double(get(hObject,'String')) returns contents of edit_dipslip as a double
handles.okadaDipslip = str2double(get(hObject,'String'));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_dipslip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dipslip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in okadainversion.
function okadainversion_Callback(hObject, eventdata, handles)
% hObject    handle to okadainversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
L = handles.okadaLength;
W = handles.okadaWidth;
D = handles.okadaDepth;
Strike = handles.okadaStrike;
Dip = handles.okadaDip;
SS = handles.okadaStrikeslip;
DS = handles.okadaDipslip;
Smx = str2double(get(handles.Smidpointx,'String'));
Smy = str2double(get(handles.Smidpointy,'String'));
Gmx = str2double(get(handles.Gmidpointx,'String'));
Gmy = str2double(get(handles.Gmidpointy,'String'));
xy = handles.xy;
los = handles.dlos;
wavelength = handles.lambda;
cmap.seismo = colormap_cpt('GMT_seis.cpt', 100);    % GMT 'Seismo' colormap for wrapped data
cmap.redToBlue = colormap_cpt('polar.cpt', 100);    % Red to Blue colormap for unwrapped data

insardata = handles.insarcrop;
[X Y] = createGrid(xy);
data = [(1:1:size(X(:),1))' X(:) Y(:)];
if Gmx == 0 &&  Gmy == 0
[Gmx Gmy] = findGmidpoint(L,W,D,Strike,Dip,SS,DS,Smx,Smy);
set(handles.Gmidpointx,'String',num2str(Gmx));
set(handles.Gmidpointy,'String',num2str(Gmy));
else
[Smx Smy] = findSmidpoint(L,W,D,Strike,Dip,SS,DS,Gmx,Gmy);
set(handles.Smidpointx,'String',num2str(Smx));
set(handles.Smidpointy,'String',num2str(Smy));
end
Green = [L W D Dip Strike Gmx Gmy SS DS 0]';
U = okada(Green,data(:,2:3), 0.25);

[m M] = faultg(L,W,D,Strike,Dip,Gmx,Gmy,Smx,Smy);
[Mo Mw] = gedeticmoment(L,W,SS,DS);
handles.green = Green;
handles.faultg = {m M Smx Smy Gmx Gmy};
guidata(hObject, handles)

Heading = mean(insardata.Heading);
Inc = mean(insardata.Inc);

ULos = U2los(U,Heading, Inc);

MaxMD = max(ULos);
MaxD =max(los);
MinMD = min(ULos);
MinD =min(los);
flag = get(handles.wrapped_flag,'Value');

axes(handles.okadamodel);
cla;
plotInsarW(data, ULos', wavelength, cmap, flag);
hold on
plot(Smx,Smy,'k^','MarkerFaceColor','black','LineWidth',2)
plot(Gmx,Gmy,'kd','MarkerFaceColor','black','LineWidth',2)
plot(m(1,:),m(2,:),'k--','LineWidth',2)
plot(M(1,:),M(2,:),'k-','LineWidth',2)
hold off

axis tight, axis equal; 
colorbar
cbar = colorbar; ylabel(cbar,'Line-of-sight displacement m','FontSize', 10);
t = title('Okada Model','FontSize', 14);


axes(handles.insarplot)
insarPlot_Callback(hObject, eventdata, handles)
hold on
midsurface = plot(Smx,Smy,'r^','MarkerFaceColor','red');
midground = plot(Gmx,Gmy,'b*','MarkerFaceColor','blue');
faultgground = plot(m(1,:),m(2,:),'k--','LineWidth',2);
faultgsurface = plot(M(1,:),M(2,:),'k-','LineWidth',2);
hold off
refreshdata

Magtxt = sprintf('Mo = %.2e Nm   Mw = %.2f \nMaxDisp = %.2fm   MaxModelDisp = %.2fm\nMinDisp = %.2fm   MinModelDisp = %.2fm',Mo,Mw, MaxD, MaxMD, MinD, MinMD);
set(handles.text24, 'String', '');
set(handles.text24,'String',Magtxt);


axes(handles.fault)
cla
plot(0,0,'r*')
hold on
if Dip<0
plot((D/tand(abs(Dip))),-D,'b*')
plot([0,(D/tand(abs(Dip)))],[0,-D],'k--')
plot((D/tand(abs(Dip)))+W*cosd(abs(Dip)),(-D)-W*sind(abs(Dip)),'b*')
plot([(D/tand(abs(Dip))),(D/tand(abs(Dip)))+W*cosd(abs(Dip))],[-D,-D-W*sind(abs(Dip))],'r-','LineWidth',2)
else
plot((D/tand(abs(Dip))),-D,'b*')
plot([0,(D/tand(abs(Dip)))],[0,-D],'k--')
plot((D/tand(abs(Dip)))-W*cosd(Dip),(-D)+W*sind(Dip),'b*')
plot([(D/tand(abs(Dip))),(D/tand(abs(Dip)))-W*cosd(abs(Dip))],[-D,(-D)+W*sind(Dip)],'r-','LineWidth',2)   
end
axis tight; axis equal;
ylabel('Depth [m]')
xlabel('Distance from the fault on surface [m]')
hold off


    function [X Y] = createGrid(xy)
        min_x = min(xy(:,2));
        max_x = max(xy(:,2));
        min_y = min(xy(:,3));
        max_y = max(xy(:,3));
        
        [X,Y]=meshgrid(min_x:200:max_x,min_y:200:max_y);
   
   function [x1 y1] = findGmidpoint(L,W,D,Strike,Dip,SS,DS,MX,MY)
               Strike = 90-Strike; %rotate towards north
               if Strike > 180
                 Strike = Strike -180;
               else
                  Strike = Strike;
               end
       
       
               if Dip > 0 
                    x1 = MX + D/tand(abs(Dip))*cosd(Strike-90);
                    y1 = MY + D/tand(abs(Dip))*sind(Strike-90);
 
                else
                    x1 = MX + D/tand(abs(Dip))*cosd(Strike+90);
                    y1 = MY + D/tand(abs(Dip))*sind(Strike+90);
               end
                
               
               function [x1 y1] = findSmidpoint(L,W,D,Strike,Dip,SS,DS,MX,MY)
               Strike = 90-Strike; %rotate towards north
               if Strike > 180
                 Strike = Strike -180;
               else
                  Strike = Strike;
               end
       
       
               if Dip > 0 
                    x1 = MX - D/tand(abs(Dip))*cosd(Strike-90);
                    y1 = MY - D/tand(abs(Dip))*sind(Strike-90);
 
                else
                    x1 = MX - D/tand(abs(Dip))*cosd(Strike+90);
                    y1 = MY - D/tand(abs(Dip))*sind(Strike+90);
                end
               
                      
             
   function U = okada(Green,data, nu)    
            U = disloc(Green,data',nu);
           
       function ULos = U2los(U,Heading, Inc)
            UEast = -cosd(Heading).* sind(Inc); % East unit vector
            UNorth = sind(Heading).* sind(Inc); % North unit vector
            UVert = cosd(Inc); % Vertical unit vector
            
            ULos = UEast.* U(1,:) + ...
                UNorth.* U(2,:) + ...             % Convert to line of sight displacement
                UVert.* U(3,:);
            
           function [m M] = faultg(L,W,D,Strike,Dip,x1,y1,MX,MY)
               Strike = 90-Strike;
               if Strike >180
                       Strike = Strike -180;
                   else
                       Strike = Strike;
                   end
               
               
                m(1,:)=[x1-L/2*cosd(Strike);x1+L/2*cosd(Strike)]';
                m(2,:)=[y1-L/2*sind(Strike);y1+L/2*sind(Strike)]';
                  
                                
                   if Dip>0
                    m(1,3) = m(1,2) + W*cosd(abs(Dip))*cosd(Strike+90);
                    m(1,4) = m(1,1) + W*cosd(abs(Dip))*cosd(Strike+90);
                    m(1,5) = m(1,1);
                    
                    m(2,3) = m(2,2) + W*cosd(abs(Dip))*sind(Strike+90);
                    m(2,4) = m(2,1) + W*cosd(abs(Dip))*sind(Strike+90);
                    m(2,5) = m(2,1);
                    
                   else
                    m(1,3) = m(1,2) + W*cosd(abs(Dip))*cosd(Strike+90);
                    m(1,4) = m(1,1) + W*cosd(abs(Dip))*cosd(Strike+90);
                    m(1,5) = m(1,1);
                    
                    m(2,3) = m(2,2) + W*cosd(abs(Dip))*sind(Strike+90);
                    m(2,4) = m(2,1) + W*cosd(abs(Dip))*sind(Strike+90);
                    m(2,5) = m(2,1);
                   end
                
                
                M(1,:)=[MX-L/2*cosd(Strike);MX+L/2*cosd(Strike)];
                M(2,:)=[MY-L/2*sind(Strike);MY+L/2*sind(Strike)];
                
                
                
               function [Mo Mw] = gedeticmoment(L,W,SS,DS)
                   Mo = 30e9 * (abs(SS)+abs(DS))*L*W;
                   Mw = (2/3)*(log10(Mo)) - 6.03; 
                


% --- Executes on button press in modelresiduals.
function modelresiduals_Callback(hObject, eventdata, handles)
% hObject    handle to modelresiduals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
green = handles.green;
insardata = handles.insarcrop;
Heading = mean(insardata.Heading);
Inc = mean(insardata.Inc);
wavelength = handles.lambda;
xy = handles.xy;
insardata = handles.insarcrop;
los = handles.dlos;
MX = handles.faultg{3};
MY = handles.faultg{4};
x1 = handles.faultg{5};
y1 = handles.faultg{6};
m =  handles.faultg{1};
M =  handles.faultg{2};

U = okada(green,xy(:,2:3), 0.25);
ULos = U2los(U,Heading, Inc);
handles.ulos = ULos;
guidata(hObject, handles)

Ures = los - ULos';

cmap.seismo = colormap_cpt('GMT_seis.cpt', 100);    % GMT 'Seismo' colormap for wrapped data
cmap.redToBlue = colormap_cpt('polar.cpt', 100);    % Red to Blue colormap for unwrapped data
flag = get(handles.wrapped_flag,'Value');
axes(handles.residual);
cla;

plotInsarW(xy, Ures, wavelength, cmap,flag);
 axis tight; axis equal;
hold on
midsurface = plot(MX,MY,'r^','MarkerFaceColor','red');
midground = plot(x1,y1,'b*','MarkerFaceColor','blue');
faultgground = plot(m(1,:),m(2,:),'k--','LineWidth',2);
faultgsurface = plot(M(1,:),M(2,:),'k-','LineWidth',2);
hold off

if ~isfield(handles,'statbox')
msgbox('Draw closed polygon using mouse');
display('Draw closed polygon using mouse');
axes(handles.residual);
polyMask = impoly;

pos=getPosition(polyMask);
handles.statbox = pos;
guidata(hObject, handles)
fprintf('Stat box is defined\n')
else
pos= handles.statbox;
hold on
fill(pos(:,1),pos(:,2),'black','FaceColor','none')
hold off
end

in=inpolygon(xy(:,2),xy(:,3),pos(:,1),pos(:,2));
ixsubset = find(in==1);

stat_res = mad(Ures(ixsubset,:));
Stat_txt=sprintf('MAD residual %.4g m',stat_res);
set(handles.text31, 'String', '');
set(handles.text31,'String',Stat_txt);

 


% --- Executes on button press in loadfaults.
function loadfaults_Callback(hObject, eventdata, handles)
% hObject    handle to loadfaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.shp');
faultdata = shaperead([pathname filename],'UseGeoCoords',true);
n = size(faultdata,1);

for i = 1:n
    faults{i}(:,1:2) = [faultdata(i).Lon' faultdata(i).Lat'];
    faults{i}(isnan(faults{i}(:,1)),:) = [];
end
faultdata = faults;
handles.faultdata = faultdata;
guidata(hObject, handles)
fprintf(['Fault Data:' [pathname filename] '\n'])


% --- Executes on button press in plotfaultdata.
function plotfaultdata_Callback(hObject, eventdata, handles)
% hObject    handle to plotfaultdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
faultdata = handles.faultdata;
georeferencePoint = [handles.georefLon handles.georefLat];
n = size(faultdata,2);
for i = 1:n
ll{i} = [single(faultdata{i}(:,1)) faultdata{i}(:,2)];   % Create Longitude and Latitude 2Xn matrix
    xy{i} = (llh2local(ll{i}', georeferencePoint))'*1000;    % Transform from geografic to local coordinates
end

axes(handles.insarplot)
hold on
for ii = 1: size(xy,2)
    ff = plot(xy{ii}(:,1),xy{ii}(:,2),'k:','LineWidth',1.5);
end
hold off
refreshdata

axes(handles.okadamodel)
hold on
for ii = 1: size(xy,2)
    ff = plot(xy{ii}(:,1),xy{ii}(:,2),'k:','LineWidth',1.5);
end
hold off
refreshdata
guidata(hObject, handles)

axes(handles.residual)
hold on
for ii = 1: size(xy,2)
    ff = plot(xy{ii}(:,1),xy{ii}(:,2),'k:','LineWidth',1.5);
end
hold off
refreshdata
guidata(hObject, handles)


% --- Executes on button press in export_kml.
function export_kml_Callback(hObject, eventdata, handles)
% hObject    handle to export_kml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ULos = handles.ulos;
xy = handles.xy;
georeferencePoint = [handles.georefLon handles.georefLat];
MX = handles.faultg{3};
MY = handles.faultg{4};
x1 = handles.faultg{5};
y1 = handles.faultg{6};
m =  handles.faultg{1};
M =  handles.faultg{2};

ll_points = (local2llh(xy(:,2:3)'/1000, georeferencePoint))';
ll_mpsurf = (local2llh([MX MY]'/1000,georeferencePoint))';
ll_mpground = (local2llh([x1 y1]'/1000,georeferencePoint))';
ll_fsurf = (local2llh(M/1000,georeferencePoint))';
ll_fground = (local2llh(m/1000,georeferencePoint))';

model = [ll_points, ULos'];
model(abs(model(:,3))<0.001,:) = [];
model = model((1:10:length(model(:,3))),:);

FPS = geoshape(ll_fsurf(:,2),ll_fsurf(:,1));
FPG = geoshape(ll_fground(:,2),ll_fground(:,1));
MPS = geopoint(ll_mpsurf(:,2),ll_mpsurf(:,1));
MPG = geopoint(ll_mpground(:,2),ll_mpground(:,1));
kmlwrite('fault_FPS.kml',FPS,'Color','blue')
kmlwrite('fault_MPS.kml',MPS,'Color','blue')
kmlwrite('fault_FPG.kml',FPG,'Color','red')
kmlwrite('fault_MPG.kml',MPG,'Color','red')
gescatter('model.kml',model(:,1),model(:,2),model(:,3));


% --- Executes on button press in Smidpoint_manual.
function Smidpoint_manual_Callback(hObject, eventdata, handles)
% hObject    handle to Smidpoint_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.insarplot);
[mx,my] = getpts;
set(handles.Smidpointx,'String',num2str(mx));
set(handles.Smidpointy,'String',num2str(my));
set(handles.Gmidpointx,'String',num2str(0));
set(handles.Gmidpointy,'String',num2str(0));



function Smidpointy_Callback(hObject, eventdata, handles)
% hObject    handle to Smidpointy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Smidpointy as text
%        str2double(get(hObject,'String')) returns contents of Smidpointy as a double
Smidpointy = str2double(get(hObject,'String'));
set(handles.Smidpointy,'String',Smidpointy);
set(handles.Gmidpointx,'String',num2str(0));
set(handles.Gmidpointy,'String',num2str(0));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function Smidpointy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Smidpointy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Smidpointx_Callback(hObject, eventdata, handles)
% hObject    handle to Smidpointx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Smidpointx as text
%        str2double(get(hObject,'String')) returns contents of Smidpointx as a double
Smidpointx = str2double(get(hObject,'String'));
set(handles.Smidpointx,'String',Smidpointx);
set(handles.Gmidpointx,'String',num2str(0));
set(handles.Gmidpointy,'String',num2str(0));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function Smidpointx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Smidpointx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gmidpointx_Callback(hObject, eventdata, handles)
% hObject    handle to Gmidpointx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gmidpointx as text
%        str2double(get(hObject,'String')) returns contents of Gmidpointx as a double
Gmidpointx = str2double(get(hObject,'String'));
set(handles.Gmidpointx,'String',Gmidpointx);
set(handles.Smidpointx,'String',num2str(0));
set(handles.Smidpointy,'String',num2str(0));
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function Gmidpointx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gmidpointx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gmidpointy_Callback(hObject, eventdata, handles)
% hObject    handle to Gmidpointy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gmidpointy as text
%        str2double(get(hObject,'String')) returns contents of Gmidpointy as a double
Gmidpointy = str2double(get(hObject,'String'));
set(handles.Gmidpointy,'String',Gmidpointy);
set(handles.Smidpointx,'String',num2str(0));
set(handles.Smidpointy,'String',num2str(0));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function Gmidpointy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gmidpointy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Gmidpoint_manual.
function Gmidpoint_manual_Callback(hObject, eventdata, handles)
% hObject    handle to Gmidpoint_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.insarplot);
[mx,my] = getpts;
set(handles.Gmidpointx,'String',num2str(mx));
set(handles.Gmidpointy,'String',num2str(my));
set(handles.Smidpointx,'String',num2str(0));
set(handles.Smidpointy,'String',num2str(0));


% --------------------------------------------------------------------
function uitoggletool4_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm = datacursormode( handles.figure1 );
if isequal(eventdata.Source.State,'on')
dcm.enable = 'on'; 
set(dcm,'UpdateFcn',@updateDataCursor);
elseif isequal(eventdata.Source.State,'off')
dcm.enable = 'off';
end


function txt = updateDataCursor(empt,event_obj)
    pos = get(event_obj, 'Position' );
    target =get(event_obj,'Target');
    los = target.FaceVertexCData;
    txt = { ['X: ', num2str( round( pos( 1 ), 3 ) )],...
          ['Y: ', num2str( round( pos( 2 ), 3 ) )],...
    ['Value: ', num2str( round( los( 3 ), 3 ) )] };


% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat');
if 'filename'~=0
load_handles = load([pathname filename]);
if isstruct(load_handles)
fig = handles.figure1;
handles = load_handles;
handles.figure1 = fig;
guidata(hObject, handles)
fprintf(['Load Parameter Data:' [pathname filename] '\n'])
else
    fprintf(['Error: Wrong input file\n'])
end
end


% --------------------------------------------------------------------
function uipushtool4_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 pathname=uigetdir;
 if pathname ~= 0
 saveParameter = handles;
 %saveParameter.figure1 = [];
 dir = [pathname,'/insar2okada_parmdata.mat'];
 save(dir,'saveParameter')
 end
