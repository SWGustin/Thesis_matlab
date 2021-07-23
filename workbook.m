clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
apparatus = 'PEL'; %['ARPEL'] % this defines which geometry we are working with

speed_isoplane_threshold = .5; %the speed at which the isoplane is drawn

visualize = false; %should we create graphs this time

%define the reduction_step_size value to reduce data size 
%(data will be 1/reduction_step_size in each dimension)
reduction_reduction_step_size_size = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("-----------------------------------------------------------------------")
%store if we're working with a pel or arpel
if strcmp(apparatus,'PEL')
    fprintf('running for pel\n')
    %percentage of interest for PEL
    %valid percentag65
    percent = 100;
    recip_prcnt = 100-percent;
    %locations of files with wildcard
    dataDir = sprintf('D:\\Sam\\Documents\\school\\Grad\\THESIS\\Data\\SAM_GUSTIN_PEL\\PEL_%01d-%01d_x_*mm_77V_16K\\PIV_MPd(4x16x16_75%%ov_ImgCorr)\\Avg Vel vector field',percent,recip_prcnt);
    %store locations of the numberOfSlicesices
    Y_numberOfSlicesice_loc = [0 5 10 15 20 25];
    
elseif strcmp(apparatus, 'ARPEL')
    fprintf('running for arpel\n'); 
    %percentage of interest for ARPEL
    prcnt1 = 50;
    prcnt2 = 0;
    prcnt3 = 0;
    prcnt4 = 50;
    dataDir = sprintf('D:\\Sam\\Documents\\school\\THESIS\\Data\\Grad\\SAM_GUSTIN\\SAM_GUSTIN_ARPEL\\ARPEL_%01d_%01d_%01d_%01d_x_*mm_77V_16K\\PIV_MPd(4x16x16_75%%ov_ImgCorr)\\Avg Vel vector field', prcnt1, prcnt2, prcnt3, prcnt4);  
    %store locations of the numberOfSlicesices
    Y_numberOfSlicesice_loc = [0 10 22 30];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         LOAD RAW VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this loads the vectors into a 2d array by file and number
primaryVelocityArray = loadarrayvec(dataDir,'B00001_AvgV.vc7');

%fix loading issue caused by a naming convention causing...
%the above to load two numberOfSlicesices (numberOfSlicesice 2 and 6) out of order
if strcmp(apparatus, 'PEL')
    temp = primaryVelocityArray(6);
    temp2 = primaryVelocityArray(2);
    primaryVelocityArray(6) = temp;
    for i = 6:-1:3
        primaryVelocityArray(i) = primaryVelocityArray(i-1);
    end
    primaryVelocityArray(2) = temp;
        %primaryVelocityArray(6) = temp2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CONSISTENCY CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%verify the right number of numberOfSlicesices loaded
if not(isequal(size(primaryVelocityArray),size(Y_numberOfSlicesice_loc')))
    fprintf('the wrong number of files/folder was read.  the program found')
    size(primaryVelocityArray)
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PREALLOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get dimensions of the velocity array
[primaryArrayRowNum, primaryArrayColNum] = size(primaryVelocityArray(1,1).vx);

%preallocate memory
numberOfSlices = length(Y_numberOfSlicesice_loc);


DataTensor=zeros(6,primaryArrayRowNum,primaryArrayColNum,numberOfSlices);  %688x550x6
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MANUAL MESHGRID EQUIVALENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b, c] = meshgrid(1:2,1:2,4:7);


vect = primaryVelocityArray(1,1);
%this only needs to be done once
%reformat x & z locations
x = repmat(vect.x,[1,primaryArrayColNum]);
x = reshape(x,[primaryArrayRowNum,primaryArrayColNum]);
z = repelem(vect.y',primaryArrayRowNum)';
z = reshape(z, [primaryArrayRowNum, primaryArrayColNum]);
y = zeros(primaryArrayRowNum,primaryArrayColNum,numberOfSlices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                BUILD 3D SPEED, POSITION, AND VELOCITY FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(primaryVelocityArray)
    for j = 1:size(primaryVelocityArray(1))
        vect = primaryVelocityArray(i,j);
        
        speedTensor(:,:,i)= (vect.vx.^2+vect.vy.^2).^(0.5)/1;% 1 is the normalization factor for making the measurement dimensionless
        
        xPositionTensor(:,:,i) = x;
        zPositionTensor(:,:,i) = z;
        vxPositionTensor(:,:,i) = vect.vx;
        vzPositionTensor(:,:,i) = vect.vy;
        
        y(:,:,i) = ones(size(speedTensor(:,:,i))).*Y_numberOfSlicesice_loc(i);

        transparency(:,:,i) = speedTensor(:,:,i)>speed_isoplane_threshold;

   end
end

if not(visualize)
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        REDUCE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%y_reduced = y(1:reduction_step_size:end,1:reduction_step_size:end,:);
%spd = speedTensor(1:reduction_step_size:end,1:reduction_step_size:end,:);
%x_reduced = xPositionTensor(1:reduction_step_size:end,1:reduction_step_size:end,:);
%z_reduced = zPositionTensor(1:reduction_step_size:end,1:reduction_step_size:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     CALC TRANSPARENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trans = transparency(1:reduction_step_size:end,1:reduction_step_size:end,size(primaryVelocityArray));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        INTERPOLATE IN Y AXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interp_reduction_step_size = 1;

refV = primaryVelocityArray(1,1);
if strcmp(apparatus,'PEL')
   desired_y_vect = Y_numberOfSlicesice_loc(1):interp_reduction_step_size:Y_numberOfSlicesice_loc(6);
elseif strcmp(apparatus, 'ARPEL')
    desired_y_vect = Y_numberOfSlicesice_loc(1):interp_reduction_step_size:Y_numberOfSlicesice_loc(4);
end

[Xq,Yq,Zq] = meshgrid(refV.x, desired_y_vect, refV.y);
[Xi,Yi,Zi] = meshgrid(refV.x, refV.y, Y_numberOfSlicesice_loc);

Xi = unique(Xi);
Yi = unique(Yi);
Zi = unique(Zi);

%Xi = refV.x';
%Yi = refV.y';
%Zi = Y_numberOfSlicesice_loc;

speedTensor = permute(speedTensor, [2 1 3]);

interp_speed = interp3(Xi, Yi, Zi, speedTensor, Xq, Yq, Zq, 'spline');
[Xi,Yi,Zi] = size(interp_speed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FORMAT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [0 0 1200 1000]);
hold on
tset = 'test';
title(tset);
%ylim([min(Y_numberOfSlicesice_loc),max(Y_numberOfSlicesice_loc)]);
%zlim([-5,105]);
%xlim([30, inf]);
colorbar
caxis([0,1])
xlabel('x position [mm]');
ylabel('y position [mm]');
zlabel('z position [mm]');

%make a raw isosurface
speedTensor = permute(speedTensor,[2 1 3]);
interp_speed = permute(interp_speed, [2 3 1]);
[Xq,Yq,Zq] = meshgrid(refV.y, refV.x, desired_y_vect);

[f1,v1] = isosurface(Xq, Yq, Zq, interp_speed,speed_threshold);
[f2,v2,e2] = isocaps(Xq,Yq,Zq,interp_speed,speed_threshold);

%[f1,v1] = isosurface(xPositionTensor, y, zPositionTensor, speedTensor,speed_threshold);
%[f2,v2,e2] = isocaps(xPositionTensor,y,zPositionTensor,speedTensor,speed_threshold);
p1 = patch('Faces',f1,'Vertices',v1);
p1.EdgeColor = 'none';
p1.FaceColor = 'blue';
%p1.FaceAlpha = 0.5;
p2 = patch('Faces',f2,'Vertices',v2,'FaceVertexCData',e2);
p2.FaceColor = 'interp';
p2.EdgeColor = 'none';
camlight(135,135);


%pbaspect([2,1,1]);
view(-45,45);
[Sx, Sy, Sz] = meshgrid(10:10:100, 0:5:200, 0:10:200);

xPositionTensor = unique(xPositionTensor);
zPositionTensor = unique(zPositionTensor);
y = unique(y);
z=y;
yPositionTensor = zPositionTensor;

vxPositionTensor = permute(vxPositionTensor, [2 3 1]);
vyPositionTensor = permute(vyPositionTensor, [2 3 1]);
vzPositionTensor = permute(vzPositionTensor, [2 3 1]);
disp("-------")
size(xPositionTensor)
size(yPositionTensor)
size(z)
size(vxPositionTensor)
size(vyPositionTensor)
size(vzPositionTensor)
%size(y)

%h2 = streamline(xPositionTensor,y,zPositionTensor,vyPositionTensor, vzPositionTensor, vxPositionTensor, Sx, Sy,Sz);

%h2 = streamline(zPositionTensor,xPositionTensor,y,vxPositionTensor, ones(size(vxPositionTensor))*0, vzPositionTensor, Sx, Sy,Sz);
h2 = streamline(xPositionTensor, z, yPositionTensor, vyPositionTensor, vzPositionTensor, vxPositionTensor, Sx, Sy,Sz);

%set the streamline color
set(h2, 'color', [1 0 0]);

axis tight equal
load('variables.mat') %load variables

%get dimensions of the velocity array
[nrows, ncols] = size(first_v.vx);

%store locations of the slices
Y_slice_loc = [0 11 15 20 25];

%build coherent vectorized dataset
for i = 1:size(v)
    for j = 1:size(v(1))
        %assign the ariablename v_working to v(i,j)
        v_working = v(i,j);
        %reformat x and z locations
        x_working = repmat(v_working.x,[1,ncols]);1
        z_working = repelem(v_working.y',nrows)';
        
        %reduce the vectors
        x_working = x_working(1:step:end);
        z_working = z_working(1:step:end);
        
        %build y slice
        y_working = ones(size(x_working)).*Y_slice_loc(i);
        
        %reformat the velocity matrices 
        vx_working = reshape(v_working.vx, [1,numel(v_working.vx)]);
        vz_working = reshape(v_working.vy, [1,numel(v_working.vy)]);
        
        %reduce the vecotrs
        vx_working = vx_working(1:step:end);
        vz_working = vz_working(1:step:end);

        %accumulate the results
        X = [X x_working];
        Z = [Z z_working];
        Y = [Y y_working];
        vx = [vx vx_working];
        vz = [vz vz_working];
        
    end
end
vy = zeros(size(vz));

figure('Renderer', 'painters', 'Position', [0 0 1200 1000]);
quiver3(X,Y,Z,vx,vy,vz,5);