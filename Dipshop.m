function varargout = Dipshop(varargin)
% DIPSHOP MATLAB code for Dipshop.fig
%      DIPSHOP, by itself, creates a new DIPSHOP or raises the existing
%      singleton*.
%
%      H = DIPSHOP returns the handle to a new DIPSHOP or the handle to
%      the existing singleton*.
%
%      DIPSHOP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIPSHOP.M with the given input arguments.
%
%      DIPSHOP('Property','Value',...) creates a new DIPSHOP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Dipshop_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Dipshop_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Dipshop

% Last Modified by GUIDE v2.5 05-Jun-2021 18:58:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Dipshop_OpeningFcn, ...
                   'gui_OutputFcn',  @Dipshop_OutputFcn, ...
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


% --- Executes just before Dipshop is made visible.
function Dipshop_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Dipshop (see VARARGIN)

% Choose default command line output for Dipshop
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Dipshop wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Dipshop_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Start_Callback(hObject, eventdata, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function openFile_Callback(hObject, eventdata, handles)
% hObject    handle to openFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[imgfilename imgpathname]=uigetfile({'*.jpg;*.png;*.bmp'},'Select a RGB image');

if imgfilename
    imgdata=imread([imgpathname '/' imgfilename]); 
%   Windows?????????????????? 
%   imgdata=imread([imgpathname '/' imgfilename]);
    axes(handles.axes1) 
    imshow(imgdata)
    handles.imgfilename=imgfilename;
    handles.imgdata=imgdata;    
end
guidata(hObject,handles)


% --------------------------------------------------------------------
function saveGray_Callback(hObject, eventdata, handles)
% hObject    handle to saveGray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile({'*.jpg','JPEG(*.jpg)';
                                 '*.bmp','Bitmap(*.bmp)';
                                 '*.gif','GIF(*.gif)';
                                 '*.*',  'All Files (*.*)'},'Save Picture','Untitled');
if FileName==0
    return;
else
    imwrite(handles.imggray,[PathName,FileName]);
end;


% --------------------------------------------------------------------
function saveOutcome_Callback(hObject, eventdata, handles)
% hObject    handle to saveOutcome (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile({'*.jpg','JPEG(*.jpg)';
                                 '*.bmp','Bitmap(*.bmp)';
                                 '*.gif','GIF(*.gif)';
                                 '*.*',  'All Files (*.*)'},'Save Picture','Untitled');

if FileName==0
    return;
else
    h=getframe(handles.axes3);
    imwrite(h.cdata,[PathName,FileName]);
end;


% --------------------------------------------------------------------
function editImage_Callback(hObject, eventdata, handles)
% hObject    handle to editImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function enlargeImage_Callback(hObject, eventdata, handles)
% hObject    handle to enlargeImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????1???'};
name='??????'; 
numlines=1;
defaultanswer={'1'};
anss=inputdlg(prompt,name,numlines,defaultanswer);

i=str2num(anss{1});
%??????????????????????????????
g=imresize(handles.imgdata,i);
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function shrinkImage_Callback(hObject, eventdata, handles)
% hObject    handle to shrinkImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????1???'};
name='??????'; 
numlines=1;
defaultanswer={'1'};
anss=inputdlg(prompt,name,numlines,defaultanswer);

i=str2num(anss{1});
%??????????????????????????????
g=imresize(handles.imgdata,1/i);
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function rotateImage_Callback(hObject, eventdata, handles)
% hObject    handle to rotateImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cwRotate_Callback(hObject, eventdata, handles)
% hObject    handle to cwRotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'????????????????????????????????????'};
name='??????'; 
numlines=1;
defaultanswer={'90'};
anss=inputdlg(prompt,name,numlines,defaultanswer);

i=str2num(anss{1});
%??????????????????
g=imrotate(handles.imgdata,360-i)
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function ccwRotate_Callback(hObject, eventdata, handles)
% hObject    handle to ccwRotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'????????????????????????????????????'};
name='??????'; 
numlines=1;
defaultanswer={'90'};
anss=inputdlg(prompt,name,numlines,defaultanswer);

i=str2num(anss{1});
%??????????????????
g=imrotate(handles.imgdata,i)
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function flipImage_Callback(hObject, eventdata, handles)
% hObject    handle to flipImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filpUD_Callback(hObject, eventdata, handles)
% hObject    handle to filpUD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%????????????????????????
g=flipud(handles.imgdata);
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function flipLR_Callback(hObject, eventdata, handles)
% hObject    handle to flipLR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%????????????????????????
g=fliplr(handles.imgdata);
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function cropImage_Callback(hObject, eventdata, handles)
% hObject    handle to cropImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%??????????????????
g = imcrop(handles.imgdata);
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function grayTrans_Callback(hObject, eventdata, handles)
% hObject    handle to grayTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function showGrayImage_Callback(hObject, eventdata, handles)
% hObject    handle to showGrayImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%???????????????
imgoutput=rgb2gray(handles.imgdata);
axes(handles.axes2);
imshow(imgoutput)
colormap(handles.axes2,gray(256))
handles.imggray=imgoutput;
guidata(hObject,handles)


% --------------------------------------------------------------------
function imageNeg_Callback(hObject, eventdata, handles)
% hObject    handle to imageNeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


f=rgb2gray(handles.imgdata);
%?????????????????????
g=255-f;
axes(handles.axes3);
imshow(g);

guidata(hObject, handles);


% --------------------------------------------------------------------
function binarize_Callback(hObject, eventdata, handles)
% hObject    handle to binarize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function setThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'???????????????'};
name='?????????'; 
numlines=1;
defaultanswer={'0.5'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
th=str2num(anss{1});


handles.binary_thresh = th;
guidata(hObject,handles)
a=im2double(handles.imggray);
%?????????????????????
g=im2bw(a,handles.binary_thresh);
axes(handles.axes3)
imshow(g)


% --------------------------------------------------------------------
function logTrans_Callback(hObject, eventdata, handles)
% hObject    handle to logTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'?????????i,0<i<2.35'};
name='????????????'; 
numlines=1;
defaultanswer={'1'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
i=str2num(anss{1});

imggrey=rgb2gray(handles.imgdata);
f=mat2gray(imggrey);
%????????????
g=log(f+i);
axes(handles.axes3);
imshow(g);

guidata(hObject, handles);


% --------------------------------------------------------------------
function powerlawTrans_Callback(hObject, eventdata, handles)
% hObject    handle to powerlawTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'???????????????'};
name='????????????'; 
numlines=1;
defaultanswer={'1'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
i=str2num(anss{1});

imggrey=rgb2gray(handles.imgdata);
f=mat2gray(imggrey);
%????????????
g=f.^i;
axes(handles.axes3);
imshow(g);

guidata(hObject, handles);

% --------------------------------------------------------------------
function histoEqualiz_Callback(hObject, eventdata, handles)
% hObject    handle to histoEqualiz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function showHisto_Callback(hObject, eventdata, handles)
% hObject    handle to showHisto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes4)
%???????????????????????????
[nk,rk]=imhist(handles.imggray);
bar(rk,nk,0.1);


% --------------------------------------------------------------------
function equalize_Callback(hObject, eventdata, handles)
% hObject    handle to equalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%??????????????????????????????
g=histeq(handles.imggray); 

axes(handles.axes3);
imshow(g);
handles.after_img = g;
guidata(hObject,handles);

%??????????????????????????????
axes(handles.axes4)
[new0,x0]=imhist(handles.after_img); 
bar(x0,new0,0.3);%0.3??????????????????


% --------------------------------------------------------------------
function spatialFilter_Callback(hObject, eventdata, handles)
% hObject    handle to spatialFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function smoSpatialFilter_Callback(hObject, eventdata, handles)
% hObject    handle to smoSpatialFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function averagingFilter_Callback(hObject, eventdata, handles)
% hObject    handle to averagingFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????????????????'};
name='?????????????????????'; 
numlines=1;
defaultanswer={'5'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
x1=str2num(anss{1});

%????????????????????????
m=fspecial('average',[x1,x1]);
%??????????????????????????????????????????
g=imfilter(handles.imggray,m);  
axes(handles.axes3);
imshow(g)

% --------------------------------------------------------------------
function medianFilter_Callback(hObject, eventdata, handles)
% hObject    handle to medianFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????????????????'};
name='?????????????????????'; 
numlines=1;
defaultanswer={'5'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
x1=str2num(anss{1});

%????????????????????????
g=medfilt2(handles.imggray,[x1,x1]);  
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function shaSpatialFilter_Callback(hObject, eventdata, handles)
% hObject    handle to shaSpatialFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function laplacianFilter_Callback(hObject, eventdata, handles)
% hObject    handle to laplacianFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'????????????????????????????????????'};
name='????????????(0-1)';
numlines=1;
defaultanswer={'0.2'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
x1=str2num(anss{1});

%??????????????????????????????
m=fspecial('laplacian',x1);
%????????????????????????????????????????????????
g=imfilter(handles.imggray,m);
axes(handles.axes3);imshow(g);


% --------------------------------------------------------------------
function gradientFilter_Callback(hObject, eventdata, handles)
% hObject    handle to gradientFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%????????????
[dx,dy]=gradient(double(handles.imggray));
g=sqrt(dx.^2+dy.^2);
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function sobelMethod_Callback(hObject, eventdata, handles)
% hObject    handle to sobelMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%???Sobel???????????????????????????????????????
g = edge(handles.imggray,'Sobel'); 
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function robertMethod_Callback(hObject, eventdata, handles)
% hObject    handle to robertMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%???Roberts???????????????????????????????????????
g = edge(handles.imggray,'Roberts'); 
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function frequencyFilter_Callback(hObject, eventdata, handles)
% hObject    handle to frequencyFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fourierTrans_Callback(hObject, eventdata, handles)
% hObject    handle to fourierTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function forwardFourierTrans_Callback(hObject, eventdata, handles)
% hObject    handle to forwardFourierTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


f=handles.imggray;
[M,N]=size(f);
y=fft2(f);
g=log(abs(fftshift(y))+1);
axes(handles.axes4);
imshow(g,[8,12]);

handles.ft=y;
handles.ftimg=log(abs(fftshift(y))+1);
guidata(hObject, handles);

% --------------------------------------------------------------------
function inverseFourierTrans_Callback(hObject, eventdata, handles)
% hObject    handle to inverseFourierTrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=handles.ft;
[M,N]=size(f);
y=ifft2(f);
axes(handles.axes3);
imshow(y,[min(min(abs(y))),max(max(abs(y)))]);


% --------------------------------------------------------------------
function lowPassFilter_Callback(hObject, eventdata, handles)
% hObject    handle to lowPassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function highPassFilter_Callback(hObject, eventdata, handles)
% hObject    handle to highPassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function idealLPFilter_Callback(hObject, eventdata, handles)
% hObject    handle to idealLPFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'???????????????????????????????????????' };
name='??????????????????????????????';
numlines=1;
defaultanswer={'300'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});
 
%??????????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=double(D<=D0);
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=fftshift(G);
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);


% --------------------------------------------------------------------
function butterworthLPFilter_Callback(hObject, eventdata, handles)
% hObject    handle to butterworthLPFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'?????????????????????????????????????????????????????????' '?????????????????????????????????????????????'};
name='?????????????????????????????????????????????????????????';
numlines=1;
defaultanswer={'300' '1'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});
n=str2num(anss{2});

%????????????????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=1./(1+(D./D0).^(2*n));
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=fftshift(G);
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);
handles.imdata=g;
guidata(hObject, handles);

% --------------------------------------------------------------------
function gaussianLPFilter_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianLPFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'?????????????????????????????????????????????????????????' };
name='????????????????????????????????????????????????';
numlines=1;
defaultanswer={'300'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});

%??????????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=exp(-(D.^2)./(2*(D0^2)));
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=fftshift(G);
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);
handles.imdata=g;
guidata(hObject, handles);


% --------------------------------------------------------------------
function idealLHFilter_Callback(hObject, eventdata, handles)
% hObject    handle to idealLHFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'???????????????????????????????????????' };
name='??????????????????????????????';
numlines=1;
defaultanswer={'10'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});

%??????????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=double(D>=D0);
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=fftshift(G);
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);
handles.imdata=g;
guidata(hObject, handles);


% --------------------------------------------------------------------
function butterworthHPFilter_Callback(hObject, eventdata, handles)
% hObject    handle to butterworthHPFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'?????????????????????????????????????????????????????????' '?????????????????????????????????????????????'};
name='?????????????????????????????????????????????????????????';
numlines=1;
defaultanswer={'10' '1'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});
n=str2num(anss{2});

%????????????????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=1-1./(1+(D./D0).^(2*n));
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=real(fftshift(G));
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);
handles.imdata=g;
guidata(hObject, handles);

% --------------------------------------------------------------------
function gaussianHPFilter_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianHPFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'?????????????????????????????????????????????????????????' };
name='????????????????????????????????????????????????';
numlines=1;
defaultanswer={'10'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});

%??????????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=1-exp(-(D.^2)./(2*(D0^2)));
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=real(fftshift(G));
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);
handles.imdata=g;
guidata(hObject, handles);


% --------------------------------------------------------------------
function selectiveFilter_Callback(hObject, eventdata, handles)
% hObject    handle to selectiveFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function bandpassFilter_Callback(hObject, eventdata, handles)
% hObject    handle to bandpassFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????????????????' '??????????????????????????????????????????'};
name='????????????????????????????????????';
numlines=1;
defaultanswer={'5' '15'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});
D1=str2num(anss{2});


%????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=double(D>=D0&D<=D1);
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=fftshift(G);
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);
handles.imdata=g;
guidata(hObject, handles);

% --------------------------------------------------------------------
function bandrejectFilter_Callback(hObject, eventdata, handles)
% hObject    handle to bandrejectFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????????????????' '??????????????????????????????????????????'};
name='????????????????????????????????????';
numlines=1;
defaultanswer={'5' '15'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
D0=str2num(anss{1});
D1=str2num(anss{2});

%????????????
r=handles.imggray;
[M,N]=size(r);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>M/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
H=1-double(D>=D0&D<=D1);
F=fft2(r,size(H,1),size(H,2));
G=H.*F;

%???????????????
f=fftshift(G);
axes(handles.axes4);
imshow(f);

%??????????????????
g=real(ifft2(G));
g=g(1:size(r,1),1:size(r,2));
g=uint8(g);
axes(handles.axes3);
imshow(g);
handles.imdata=g;
guidata(hObject, handles);


% --------------------------------------------------------------------
function imageDegradation_Callback(hObject, eventdata, handles)
% hObject    handle to imageDegradation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function probaDensityFunction_Callback(hObject, eventdata, handles)
% hObject    handle to probaDensityFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gaussianNoise_Callback(hObject, eventdata, handles)
% hObject    handle to gaussianNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????' '??????????????????????????????' };
name='??????????????????';
numlines=1;
defaultanswer={'0' '0.01'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
x1=str2num(anss{1});
x2=str2num(anss{2});

%????????????
g=imnoise(handles.imggray,'gaussian',x1,x2);
axes(handles.axes3)
imshow(g)

% --------------------------------------------------------------------
function saltpepperNoise_Callback(hObject, eventdata, handles)
% hObject    handle to saltpepperNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%?????????????????????
prompt={'??????????????????????????????'};
name='??????????????????';
numlines=1;
defaultanswer={'0.05'};
anss=inputdlg(prompt,name,numlines,defaultanswer);
x1=str2num(anss{1});

%??????????????????
g=imnoise(handles.imggray,'salt & pepper',x1);
axes(handles.axes3)
imshow(g)


% --------------------------------------------------------------------
function poissonNoise_Callback(hObject, eventdata, handles)
% hObject    handle to poissonNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%??????????????????
g=imnoise(handles.imggray,'poisson');
axes(handles.axes3)
imshow(g)
