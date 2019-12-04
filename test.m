%% By : Simiao Salvador da Gama
%% Source : githb.com/nukhugama
function varargout = test(varargin) 
% TEST MATLAB code for test.fig
%      TEST, by itself, creates a new TEST or raises the existing
%      singleton*. 
%
%      H = TEST returns the handle to a new TEST or the handle to 
%      the existing singleton*. 
%
%      TEST('CALLBACK',hObject,eventData,handles,...) call    s the local
%      function named CALLBACK in TEST.M with the given input arguments. 
%
%      TEST('Property','Value',...) creates a new TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_OpeningFcn via varargin. 
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test

% Last Modified by GUIDE v2.5 13-May-2019 01:49:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_OpeningFcn, ...
                   'gui_OutputFcn',  @test_OutputFcn, ...
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


% --- Executes just before test is made visible.
function test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test (see VARARGIN)

% Choose default command line output for test
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in importImage.
function importImage_Callback(hObject, eventdata, handles)
% hObject    handle to importImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of importImage
        


% --- Executes on button press in browseImage.
function browseImage_Callback(hObject, eventdata, handles)
% hObject    handle to browseImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % mencari gambar yang berupa jpg, png dan bmp;
       [filename pathname] = uigetfile({'*.*'}, 'Pilih Gambar');
       gambar = strcat(pathname, filename);
       
        % Display gambar
        axes(handles.showImage);
        imshow(gambar,'InitialMagnification', 'fit');
     
        
        
        infoFilename = imfinfo(gambar);  
        
        SIZE = infoFilename.FileSize;
        Size = SIZE/1024;
        set(handles.ukuranGasli,'String',[num2str(Size),' kb']);
        
        handles.gambar = gambar;
        handles.infoFilename = infoFilename;
        guidata(hObject, handles); 
        title('Gambar Asli');



function edit1_Callback(hObject, eventdata, handles) 
% hObject    handle to showImage (see GCBO) 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of showImage as text
%        str2double(get(hObject,'String')) returns contents of showImage as a double


% --- Executes during object creation, after setting all properties.
function showImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in deleteImg.
function deleteImg_Callback(hObject, eventdata, handles)
% hObject    handle to deleteImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Dialbox Untuk memilih menghapus gambar atau tidak
    choice = questdlg('Are sure to delete the image ? ...','Yes','Cancel' )
    
    %Proses penghapusan gambar tau tidak 
    if (strcmp('Yes',choice))
       cla (handles.showImage, 'reset');
       cla (handles.mKecil, 'reset');
       cla (handles.mBesar, 'reset');
       cla (handles.showCrop, 'reset');
       cla (handles.showHisteq,'reset');
       cla (handles.showHistogram,'reset');
       cla (handles.showConvolution,'reset');
       cla (handles.ukuranGasli,'reset');
       cla (handles.editLoss,'reset');
       
    elseif (strcmp('Cancel',choice))
        return;
    end;
    
    


% --- Executes on button press in btnExit.
function btnExit_Callback(hObject, eventdata, handles)
% hObject    handle to btnExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Dialbox Untuk memilih keluar program atau tidak
    choice = questdlg('Are sure to close the app ? ...','Yes','Cancel' )
    
    % Statement untuk keluar program atau tidak
    if (strcmp('Yes',choice))
        close;
        
    elseif (strcmp('Cancel',choice))
        return;
    end;
    
  


% --- Executes on button press in btnGrayscale.
function btnGrayscale_Callback(hObject, eventdata, handles)
% hObject    handle to btnGrayscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        
        % Membaca gambar yang sudah dipilih ke dalam variabel dptGambar
        dptGambar = imread(handles.gambar); 
        
        % Mendefenisi RGB pada matriks gambar
        red = dptGambar(:,:,1);
        green = dptGambar(:,:,2); 
        blue = dptGambar(:,:,3);
        
        % Membuat Opersai supaya mengkonversi gambar ke grayscale
        imgGray= 0.2*red + 0.2*green + 0.6*blue; 
        
        % Output gambar yang sudah diubah ke grayscale;
        
        % Display gambar yang sudah dioperasikan
        axes(handles.showImage);
        imshow(imgGray);
        title('Grayscale Image');
        
%         axes(handles.mKecil);
%         axes(handles.mBesar);   
        
      
%        

% --- Executes on button press in btnTambah.
function btnTambah_Callback(hObject, eventdata, handles)
% hObject    handle to btnTambah (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        % Membaca gambar yang sudah dipilih ke dalam variabel dptGambar
        dptgambar = imread(handles.gambar);
        %Proses Penambahan Gambar
        
        gambarTambah =dptgambar + 100 ;
        axes(handles.showImage);
        imshow(gambarTambah);
        title('Penambahan Image');
        %dptgambar = gambarTambah;
       
       % Display gambar yang sudah dioperasikan
        
        
%         if (handles.mKecil==1)
%              axes(handles.mKecil) 
%              imshow(gambarTambah);
%         end;
%         
%         if (handles.mBesar==1)
%              axes(handles.mBesar);
%              imshow(gambarTambah);
%         end;
       

% --- Executes on button press in btnKurang.
function btnKurang_Callback(hObject, eventdata, handles)
% hObject    handle to btnKurang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        % Membaca gambar yang sudah dipilih ke dalam variabel dptGambar
        dptGambar = imread(handles.gambar);
        % Proses Pengurangan
        gambarKurang = dptGambar - 100 ;
        
        % Display gambar yang sudah dioperasikan 
        axes(handles.showImage);
        imshow(gambarKurang);
        title('Pengurangan Image');

        
%         if (handles.mKecil==1)
%              axes(handles.mKecil)  
%              imshow(gambarKurang);
% 
%         end;
%         
%         if (handles.mBesar==1)
%              axes(handles.mBesar);
%              imshow(gambarKurang);
% 
%         end;
       

% --- Executes on button press in btnBagi.
function btnBagi_Callback(hObject, eventdata, handles)
% hObject    handle to btnBagi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        % Membaca gambar yang sudah dipilih ke dalam variabel dptGambar
        dptGambar = imread(handles.gambar);
        % Proses Pembagian
        gambarBagi = dptGambar/5 ;
        
        % Display gambar yang sudah dioperasikan
        axes(handles.showImage);
        imshow(gambarBagi);
        title('Pembagian Image');
       
%         if (handles.mKecil==1)
%              axes(handles.mKecil)     
%         end;
%         
%         if (handles.mBesar==1)
%              axes(handles.mBesar);     
%         end;

% --- Executes on button press in btnKali.
function btnKali_Callback(hObject, eventdata, handles)
% hObject    handle to btnKali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        % Membaca gambar yang sudah dipilih ke dalam variabel dptGambar
        dptGambar = imread(handles.gambar);
        %Proses Perkalian
        gambarKali = dptGambar * 5 ;
        
         % Display gambar yang sudah dioperasikan
        axes(handles.showImage);
        imshow(gambarKali);
        title('Perkalian Image');
        
        
%         if (handles.mKecil==1)
%              axes(handles.mKecil)    
%         end;
%         
%         if (handles.mBesar==1)
%              axes(handles.mBesar);     
%         end;
%       
        


        
%    Fungsi ntuk Memperkecil atau memperbesar pada gambar
function [out] = memperBesar_Kecil(im, out_dims)

    %// Membuat beberapa variabel
    in_rows = size(im,1);
    in_cols = size(im,2);
    out_rows = out_dims(1);
    out_cols = out_dims(2);

    %// Variabel S_R = R / R'        
    S_R = in_rows / out_rows;
    
    %// Let S_C = C / C'
    S_C = in_cols / out_cols;

    %// Mendefinisi grid dari koordinasi pada gambar
    %// Generate (x,y) pada tiap titik di dalam gambar
    [cf, rf] = meshgrid(1 : out_cols, 1 : out_rows);

    %// Variabel r_f = r'*S_R for r = 1,...,R'
    %// Variabel c_f = c'*S_C for c = 1,...,C'
    rf = rf * S_R;
    cf = cf * S_C;

    %// variabel  r diisi floor(rf) dan c diisi floor(cf)
    r = floor(rf);
    c = floor(cf);

    %// Any values out of range, cap
    r(r < 1) = 1;
    c(c < 1) = 1;
    r(r > in_rows - 1) = in_rows - 1;
    c(c > in_cols - 1) = in_cols - 1;

    %// Variabel delta_R = rf - r dan delta_C = cf - c
    delta_R = rf - r;
    delta_C = cf - c;

    
    %// Ambil column yang indeks besar pada tiap titik yang kita mau
    %//  access
    in1_ind = sub2ind([in_rows, in_cols], r, c);
    in2_ind = sub2ind([in_rows, in_cols], r+1,c);
    in3_ind = sub2ind([in_rows, in_cols], r, c+1);
    in4_ind = sub2ind([in_rows, in_cols], r+1, c+1);       

    %// Sekarang menginterpolarisasi
    %// Menuju ke tiap channel untuk mendeteksi kasus warna
    %// Membuat output image dimana terdapat class input yang sama
    out = zeros(out_rows, out_cols, size(im, 3));
    out = cast(out, class(im));

    for idx = 1 : size(im, 3)
        chan = double(im(:,:,idx)); %// Get i'th channel
        %// Interpolarisasi ke channel
        tmp = chan(in1_ind).*(1 - delta_R).*(1 - delta_C) + ...
                       chan(in2_ind).*(delta_R).*(1 - delta_C) + ...
                       chan(in3_ind).*(1 - delta_R).*(delta_C) + ...
                       chan(in4_ind).*(delta_R).*(delta_C);
        out(:,:,idx) = cast(tmp, class(im));
    end

        
% --- Executes on button press in btnZoomIn.
function btnZoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to btnZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
      
     % Membaca gambar yang sudah dipilih ke dalam variabel dptGambar
      dptGambar = imread(handles.gambar);
      
      % Menyimpan ukuran Gambar tersebut ke dalam matriks(baris dan kolom)
      [baris, kolom] = size(dptGambar); 
      
%       disp(baris)
%       disp('=====');
%       disp(kolom);
%       disp('=====');

      % Memperbesar 2 kali dari baris dan kolom ukuran gambar tersebut
      a = (baris*2);
      b = (kolom*2); 
      
     % Dalam statement ini saya membatasi size baris dan kolom tidak lebih
     % 2000
     if (a>=2000)
        a=2000;
     end
     
     if (b>=2000)
        b = 2000;
     end;
     
      % Proses perbesaran Gambar 
      perbesar = memperBesar_Kecil(dptGambar,[a b]);
      % Display gambar yang sudah diperbesar
      axes(handles.mBesar);
      imshow(perbesar,'InitialMagnification', 'fit');
      title('Memperbesar Image');
           
 


% --- Executes on button press in btnZoomOut.
function btnZoomOut_Callback(hObject, eventdata, handles) 
% hObject    handle to btnZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%        
      % membaca gambar yang sudah dipilih
      dptGambar = imread(handles.gambar);
      
      % Menyimpan ukuran Gambar tersebut ke dalam matriks(baris dan kolom)
      [baris, kolom] = size(dptGambar);
      

    % Dalam statement ini saya membatasi size baris dan kolom tidak lebih
    % dari 1000
     if (baris>=1000)
        baris=1000;
     end
     
     if (kolom>=1000)
        kolom = 1000;
     end
     
     % Proses perkecilan baris dan kolom
     a = round(baris/2);
     b = round(kolom/2);
     

     % Memperkecil gambar
      perkecil = memperBesar_Kecil(dptGambar,[a b]);
      
      %Display gambar yang sudah  diperkecil
      axes(handles.mKecil);
      imshow(perkecil,'InitialMagnification', 'fit');
      title('Memperkecil Image');
     
      

% --- Executes during object creation, after setting all properties.
function mBesar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mBesar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mBesar


% --- Executes during object creation, after setting all properties.
function mKecil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mKecil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mKecil


% --- Executes during object creation, after setting all properties.
function showCrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate showCrop


% --- Executes on button press in btnCrop.
function btnCrop_Callback(hObject, eventdata, handles)
% hObject    handle to btnCrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % membaca gambar yang sudah dipilih
    I = imread(handles.gambar);
   

    % Pesan, apakah mau melakukan crop pada gambar tersebut
    message = sprintf('Pilih titik yang yg maun dicrop..!');
    uiwait(msgbox(message));
    
    % Proses Crop baru dimulai ketika buttonnya dipress
    k = waitforbuttonpress;
    
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up terdeteksi
    point1 = point1(1,1:2);              % extract ke x dan y
    point2 = point2(1,1:2);
    
    p1 = min(point1,point2);             % kalkulasi lokasi
    offset = abs(point1-point2);         % dan dimensi 
    x = round([p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)]) 
    y = round([p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)]);
   
    % proses crop 
    im = I(y(1):y(3), x(1):x(2));

    imageCrop= uint8(im);
    
    % output into axes 
    axes(handles.showCrop); 
    imshow(imageCrop,[]);
    title('Crop Image');


% --- Executes on button press in Histo.
function Histo_Callback(hObject, eventdata, handles)
% hObject    handle to Histo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca gambar 
    I = imread(handles.gambar);
    
    %Ubah rgb ke grayscale
    Im=rgb2gray(I);
    
   
    [baris,kolom]=size(Im);  %Menentukan ukuran baris dan kolom pada gambar
    z=zeros(1,256);        % data untuk matriks z di assign nol  
    for i=1:baris
        for j=1:kolom
            b=Im(i,j);
            z(b+1)=z(b+1)+1;
        end
    end
    
    N=sum(z); % Nilai N di assign jumlah total dari z
    p=zeros(1,256); % data untuk matriks p di assign nol  
    s=zeros(1,256); % data untuk matriks s di assign nol
    kolom=zeros(1,256); % data untuk array kolom di assign nol
    baris=zeros(1,256); % data untuk array baris di assign nol
    for k=1:256
        p(k)=z(k)/N;
        if k==1
            kolom(k)=p(k);
            s(k)=kolom(k)*255;
            baris(k)=floor(s(k));
        else
            kolom(k)=kolom(k-1)+p(k); 
            s(k)=kolom(k)*255;
            baris(k)=floor(s(k)); 
        end
    end
   
    axes(handles.showHistogram); 
    %Plotting baris dan z 
    plot(baris,z); 
    title('Histogram');
   



% --- Executes on button press in Histeq.
function Histeq_Callback(hObject, eventdata, handles)
% hObject    handle to Histeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca gambar 
    Im = imread(handles.gambar);
    % Menentukan ukuran baris dan kolom pada gambar;
%     [baris, kolom] = size(Im);
%     Menentukan banyaknya pixel
    numofpixels=size(Im,1)*size(Im,2);
    
    imgHistogramed=uint8(zeros(size(Im,1),size(Im,2)));
% frequency pdf dan cdf serta beberapa variables
    freq=zeros(256,1);
    pdf=zeros(256,1);
    cdf=zeros(256,1);
    cum=zeros(256,1);
    output=zeros(256,1);

    % Proses Pembuatan Histeq
    for i=1:size(Im,1)
        for j=1:size(Im,2)
            value=Im(i,j); %nilai RGB pada tiap posisi diisi di variable value
            freq(value+1)=freq(value+1)+1; % Menentukan nilai frekuensi  
            pdf(value+1)=freq(value+1)/numofpixels; % Menentukan nilai pdf(probability density function)
        end
    end

    sum=0;

    for i=1:size(pdf) 
       sum=sum+freq(i); 
       cum(i)=sum; % Kumulatif di assign dengan nilai sum
       cdf(i)=cum(i)/numofpixels; % Mencari cdf (fungsi kumulatif)
       output(i)=round(cdf(i)*255);
    end

    for i=1:size(Im,1)
        for j=1:size(Im,2)
                imgHistogramed(i,j)=output(Im(i,j)+1);
        end
    end
    
    
    % Proses melakukan histogramnya yang sudah diperbesar
    [baris,kolom]=size(imgHistogramed);
    % Membuat nilai matriks zero untuk matriks z,
    z=zeros(1,256);  
    for u=1:baris
        for v=1:kolom
            b=imgHistogramed(u,v);
            z(b+1)=z(b+1)+1;
            zTot = z(b+1);
        end
    end
    
    N = zTot; % Variabel N di assign dengan nilai zTot
    p=zeros(1,256); % data untuk matriks p di assign nol  
    s=zeros(1,256); % data untuk matriks s di assign nol
    kolom=zeros(1,256); % data untuk array kolom di assign nol
    baris=zeros(1,256); % data untuk array baris di assign nol
    for k=1:256
        p(k)=z(k)/N; 
        if k==1
            kolom(k)=p(k);
            s(k)=kolom(k)*255;
            baris(k)=floor(s(k)); % ambil floornya 
        else
            kolom(k)=kolom(k-1)+p(k); 
            s(k)=kolom(k)*255;
            baris(k)=floor(s(k)); % ambil floornya
        end
    end; 
   
   % Display Histeq
    axes(handles.showHisteq);
    plot(baris,z);
    title('Histeq');
    
%     bb1 = Im(1,1);
%     ba1 = Im(1,1);
%     for i=1:size(Im,1)
%         for j=1:size(Im,2)
%             if (Im(i,j)>=ba1)
%                 ba1 = Im(i,j);   
%             end
%             if (Im(i,j)<=bb1)
%                 bb1 = Im(i,j);
%             end;
%             
%         end;    
%     end;
    
%     Input manual batas atas dan batas bawah baru

%     bb2 = bb1;
%     ba2 = ba1 + 150;
%    
%     for a=1:256
%         
%         for b=1:size(Im,2)
%             % x nilai lama
%             x = Im(a,b);
%             %x dengan nilai baru
%             x_baru(a,b) = bb2 + ((x - bb1)* ((ba2-bb2)/(ba1-bb1)))
%         end;
%     end;
%     axes(handles.showHisteq);
%     plot(x_baru);
   
    
    


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text6.
function text6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function showHisteq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showHisteq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate showHisteq 


% --- Executes during object creation, after setting all properties.
function showHistogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called 

% Hint: place code in OpeningFcn to populate showHistogram  


% --- Executes on button press in btnBlur.
function btnBlur_Callback(hObject, eventdata, handles)
% hObject    handle to btnBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca gambar 
    rgbImage = imread(handles.gambar);
    % Menentukan ukuran baris, kolom dan banyaknya warna 
    [rows, columns, numberOfColorBands] = size(rgbImage);

    % Mengekstraksi saluran warna merah, hijau dan biru 
    redChannel = rgbImage(:, :, 1);
    greenChannel = rgbImage(:, :, 2);
    blueChannel = rgbImage(:, :, 3);
    
    % Membuat Gaussian Kernel 
    m=15; 
    n=15;
    sigma=17;
    [h1 h2]=meshgrid(-(m-1)/2:(m-1)/2, -(n-1)/2:(n-1)/2);
    hg= exp(-(h1.^2+h2.^2)/(2*sigma^2));            %Fungsi Gaussian 
    h=hg/sum(hg(:))


    % Membelitkan tiga warna RGB di setiap channel masing masing 
    redBlurred = conv2(redChannel, h);
    greenBlurred = conv2(greenChannel, h);
    blueBlurred = conv2(blueChannel, h);

    % Kombinasi kembali warna yang terpisah ke dalam satu chanel(saluran)
    % atau kembali ke gambar RGB
    rgbImage2 = cat(3, uint8(redBlurred), uint8(greenBlurred), uint8(blueBlurred));
    
    % Display gambar blur.
    axes(handles.showConvolution);
    imshow(rgbImage2); 
    title('Blur Image');
    


% --- Executes on button press in btnEdge.
function btnEdge_Callback(hObject, eventdata, handles)
% hObject    handle to btnEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% EDGE DETECTION menggunakan Sobel Method
    % Membaca Gambar
    Im = imread(handles.gambar);
    Im=double(rgb2gray(Im));
    

    for i=1:size(Im,1)-2
        for j=1:size(Im,2)-2
        %Sobel mask untuk arah x
            mx=((2*Im(i+2,j+1)+Im(i+2,j)+Im(i+2,j+2))-(2*Im(i,j+1)+Im(i,j)+Im(i,j+2)));
            %Sobel mask untuk arah y
            my=((2*Im(i+1,j+2)+Im(i,j+2)+Im(i+2,j+2))-(2*Im(i+1,j)+Im(i,j)+Im(i+2,j)));

            Im(i,j)=sqrt(mx.^2+my.^2);
        end
    end
    


%     Thresh=100;
%     B=max(Im,Thresh);
%     B(B==round(Thresh))=0;
    B=uint8(Im);
    % Display Hasil Edge 
    axes(handles.showConvolution)
    imshow(B);
    title('Edge Detection');
    
    
    



% --- Executes on button press in btnSharp.
function btnSharp_Callback(hObject, eventdata, handles)
% hObject    handle to btnSharp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca Gambar
    Im = imread(handles.gambar); 
    % Diubah kedalam bentuk double
    Im = im2double(Im);

    lap = [0 -1 0; -1 5 -1; 0 -1 0]; %// Isi matriks untuk sharp
    resp = imfilter(Im, lap, 'conv'); %// Memfilter gambar dengan lap 

    %// Normalisasi gambar 
    minR = min(resp(:));
    maxR = max(resp(:));
    resp = (resp - minR) / (maxR - minR);

    %// Menambah ke gambar yang asli;
    sharpened = Im + resp;

    %// Normalisasi hasil Sharppened
    minA = min(sharpened(:));
    maxA = max(sharpened(:));
    sharpened = (sharpened - minA) / (maxA - minA);

    %//  Ubah performasi peningkatan linear constrant
    sharpened = imadjust(sharpened, [60/255 200/255], [0 1]);
    
    % Display Gambar Sharp
    axes(handles.showConvolution) 
    imshow([sharpened]); 
    title('Sharp Image'); 
    
    


% --- Executes on button press in btnBaseTreshold.
function btnBaseTreshold_Callback(hObject, eventdata, handles) 
% hObject    handle to btnBaseTreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca Gambar
    Im = imread(handles.gambar);
    % Ubah ke dalam tipe bentukan double
    Im = double(Im);
    [baris,kolom,dim] = size(Im);
    region = inputdlg({'Masukan Nilai R : ', 'Masukan Nilai G :', 'Masukan Nilai B :'})
     
    % Ubah nilai R, G dan B dari string ke number 
    R=str2num(region{1});
    G=str2num(region{2})
    B=str2num(region{3})
    
    for i=1:baris
        for j=1:kolom
            if (Im(i,j,1)>R) && (Im(i,j,2)>G) && (Im(i,j,3)<=B) % Isi dengan warna putih 
                Im(i,j,1) = 255;
                Im(i,j,2) = 255;
                Im(i,j,3) = 255;
            else  % isi dengan warna hitam
                Im(i,j,1) = 0;
                Im(i,j,2) = 0;
                Im(i,j,3) = 0;
            end;
        end;
    end;
    
    baseTresh = cat(3, uint8(Im(:,:,1)),uint8(Im(:,:,2)), uint8(Im(:,:,3)))
    
    % Display Gambar Seed Region Growth 
    axes(handles.showImage)
    imshow(baseTresh); 
    title('Base Treshold'); 

    
    
% Fungsi Seed Region Growth 
function [J] = seedRegionGrowth(I,x,y,reg_maxdist)
   %Jika nilai reg_maxdist=0 maka reg_maxdist=0.2
    if(exist('reg_maxdist','var')==0) 
        reg_maxdist=0.2; 
    end
    %Jika y=0 
    if(exist('y','var')==0) 
        imshow(I,[]); 
        [y,x]=getpts; 
        y=round(y(1)); 
        x=round(x(1)); 
    end
    J = zeros(size(I)); % Output 
    Isizes = size(I); % Dimensi input image
    reg_mean = I(x,y); % Rata rata segmented region
    reg_size = 1; % Jumlah pixel di dalam region
    % Free memory 
    neg_free = 10000;
    neg_pos=0;
    neg_list = zeros(neg_free,3); 
    pixdist=0; % Jarak antara regio pixel yg terbaru ke rata rata regio 
    % Neighbor locations (footprint)
    neigb=[-1 0; 1 0; 0 -1;0 1];
    % Mulai regiongrowing sampai jarak antara region memungkinkan pixel
    % baru lebih tinggi dibanding treshold yg sebenarnya
    while(pixdist<reg_maxdist && reg_size<numel(I))
        % Tambahkan pixel tetangga yang baru (new neighbors pixels)
        for j=1:4,
            % Kalkulasi koordinasi ketetanggan
            xn = x +neigb(j,1); yn = y +neigb(j,2);
            % Mengecek Jika ada tetangga yang ada di dalam atau luar gambar
            ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));

            % Tambahkan neighbor(tetangga)jika dalamnya sudah tidak
            % termasuk bagian dari luas segmen (segment area)
            if(ins&&(J(xn,yn)==0)) 
                    neg_pos = neg_pos+1;
                    neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
            end
        end
        % Tambahkan Block baru dari memori yang masih kosong
        if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end

        % Tambahkan pixel dengan intensitas yang paling mendekati rata-rata
        % dari suatu region ke region lainnya
        dist = abs(neg_list(1:neg_pos,3)-reg_mean); 
        [pixdist, index] = min(dist); 
        J(x,y)=2; reg_size=reg_size+1; 

        % Menghitung rata-rata dari wilayah(region) 
        reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1); 

        % Mengimpan koordinasi pixel x dan y (tambahkan proses untuk neighbour) 
        x = neg_list(index,1); y = neg_list(index,2);

        % Menghapus pixel dari neighbour(check) list
        neg_list(index,:)=neg_list(neg_pos,:); 
        neg_pos=neg_pos-1;
    end
    % Mengembalikan area segmen menjadi logikal matriks 
    J=J>1;
    
    
    


% --- Executes on button press in btnSRegionGrowth.
function btnSRegionGrowth_Callback(hObject, eventdata, handles)
% hObject    handle to btnSRegionGrowth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca gambar
    Im = imread(handles.gambar);
    % Ubah gambar dalam bentuk double
    Im = im2double(Im);
    % Pilih nilai x dan y , (dimana gambar itu sendiri mulai akan disegmentasi) melalui dialogbox 
    region = inputdlg({'Masukan Nilai X : ', 'Masukan Nilai Y :' })
     
    % Ubah nilai x dan y dari string ke number 
    x=str2num(region{1});
    y=str2num(region{2}); 
%     
%     J = seedRegionGrowth(Im); 
    J = seedRegionGrowth (Im,x,y); 
    
    % Display Gambar Seed Region Growth 
    axes(handles.showImage)
    imshow(Im+J); 
    title('Seed Region Growth'); 

 function G = dilasi(F, H, hotx, hoty)
% DILASI Berguna untuk melaksanakan operasi dilasi.
%     Masukan:
%        F = citra yang akan dikenai dilasi
%        H = elemen pentruksur
%        (hy, hx) koordinat pusat piksel
 
    [th, lh]=size(H);
    [tf, lf]=size(F);

    if nargin < 3
        hotx = round(lh/2);
        hoty = round(th/2);
    end

    Xh = [];
    Yh = [];
    jum_anggota = 0;

    % Menentukan koordinat piksel bernilai 1 pada H
    for baris = 1 : th
        for kolom = 1 : lh
            if H(baris, kolom) == 1
                jum_anggota = jum_anggota + 1;
                Xh(jum_anggota) = -hotx + kolom;
                Yh(jum_anggota) = -hoty + baris;
            end
        end
    end

    G = zeros(tf, lf); % Nolkan semua pada hasil dilasi

    % Memproses dilasi
    for baris = 1 : tf
        for kolom = 1 : lf
            for indeks = 1 : jum_anggota
                if F(baris, kolom) == 1
                    xpos = kolom + Xh(indeks);
                    ypos = baris + Yh(indeks);
                    if (xpos >= 1) && (xpos <= lf) && ...
                       (ypos >= 1) && (ypos <= tf)
                        G(ypos, xpos) = 1;
                    end
                end
            end
        end
    end


% --- Executes on button press in btnDilasi.
function btnDilasi_Callback(hObject, eventdata, handles)
% hObject    handle to btnDilasi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca Gambar
    Im = imread(handles.gambar);
    ImBW = im2bw(Im, 0.5);
    
    % Elemen pentruksur
    H = [1 0
         2 4];
%     H = ones(7)
    
    outDil = dilasi(ImBW, H); 
    
    % Display Gambar Dilasi 
    axes(handles.showImage)
    imshow(outDil); 
    title('Dilasi');
    
    
function G = erosi(F, H, hotx, hoty)
% DILASI Berguna untuk melaksanakan operasi erosi.
%     Masukan:
%        F = citra yang akan dikenai erosi
%        H = elemen pentruksur
%        (hy, hx) koordinat pusat piksel
 
    [th, lh]=size(H);
    [tf, lf]=size(F);

    if nargin < 3
        hotx = round(lh/2);
        hoty = round(th/2);
    end

    Xh = [];
    Yh = [];
    jum_anggota = 0; 

    % Menentukan koordinat piksel bernilai 1 pada H
    for baris = 1 : th
        for kolom = 1 : lh
            if H(baris, kolom) == 1
                jum_anggota = jum_anggota + 1;
                Xh(jum_anggota) = -hotx + kolom;
                Yh(jum_anggota) = -hoty + baris;
            end
        end
    end

    G = zeros(tf, lf); % Nolkan semua pada hasil erosi
    % Memproses erosi
    for baris = 1 : tf
        for kolom = 1 : lf
            cocok = true;
            for indeks = 1 : jum_anggota
                xpos = kolom + Xh(indeks);
                ypos = baris + Yh(indeks);
                if (xpos >= 1) && (xpos <= lf) && (ypos >= 1) && (ypos <= tf)
                    if F(ypos, xpos) ~= 1
                        cocok = false;
                        break;
                    end;
                else
                     cocok = false;
                end
            end
            if cocok
                G(baris, kolom) = 1;
            end
        end
    end

% --- Executes on button press in btnErosi.
function btnErosi_Callback(hObject, eventdata, handles)
% hObject    handle to btnErosi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % Membaca Gambar
    Im = imread(handles.gambar);
    ImBW = im2bw(Im, 0.1);
    
    H = ones(4);
    G = erosi(ImBW, H); 
    
    
    % Display Gambar Erosi
    axes(handles.showImage)
    imshow(G, [0 1]); 
    title('Erosi'); 
    
    
   
  

% --- Executes on button press in btnLossless.
function btnLossless_Callback(hObject, eventdata, handles)
% hObject    handle to btnLossless (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % memanggil variabel Img & fileinfo yang ada pada lokasi handles
    Img = imread(handles.gambar);
    fileinfo =  handles.infoFilename;

    % melakukan lossless kompresi terhadap citra asli
    Img_comp = zeros(size(Img));
    for n = 1:3
        I = im2double(Img(:,:,n));
        T = dctmtx(8);
        B = blkproc(I,[8 8],'P1*x*P2',T,T');
        mask = [1   1   1   1   0   0   0   0
            1   1   1   0   0   0   0   0
            1   1   0   0   0   0   0   0
            1   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0
            0   0   0   0   0   0   0   0];
        B2 = blkproc(B,[8 8],'P1.*x',mask);
        Img_comp(:,:,n) = blkproc(B2,[8 8],'P1*x*P2',T',T);
    end

    
    % menyimpan file citra hasil kompresi
    file_comp = ['Hasil Kompresi.',fileinfo.Format];
    imwrite(Img_comp,file_comp);
    
    % melihat ukuran file citra hasil kompresi dalam satuan kb
    fileinfo_comp = imfinfo(file_comp);
    SIZE = fileinfo_comp.FileSize;
    Size = SIZE/1024;
    set(handles.editLoss,'String',[num2str(Size),' kb']);
    
    % menampilkan file citra hasil kompresi ke dalam axes
    axes(handles.showComp);
    imshow( Img_comp);
    title('Lossless Image');



function ukuranGasli_Callback(hObject, eventdata, handles)
% hObject    handle to ukuranGasli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ukuranGasli as text
%        str2double(get(hObject,'String')) returns contents of ukuranGasli as a double


% --- Executes during object creation, after setting all properties.
function ukuranGasli_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ukuranGasli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editLoss_Callback(hObject, eventdata, handles)
% hObject    handle to editLoss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLoss as text
%        str2double(get(hObject,'String')) returns contents of editLoss as a double


% --- Executes during object creation, after setting all properties.
function editLoss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLoss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
