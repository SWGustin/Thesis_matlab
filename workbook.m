clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
apparatus = 'PEL'; %['ARPEL'] % this defines which geometry we are working with

speedIsoplaneThreshold = .25; %the speed at which the isoplane is drawn

visualize = true; %should we create graphs this time

%define the reductionStepSize value to reduce data size 
%(data will be 1/reductionStepSize in each dimension)
reductionStepSize = 8;

%step between y layers in the interpolation (units of mm)
interpStep = 1;
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
    ySliceLocs = [0 5 10 15 20 25];
    
elseif strcmp(apparatus, 'ARPEL')
    fprintf('running for arpel\n'); 
    %percentage of interest for ARPEL
    prcnt1 = 50;
    prcnt2 = 0;
    prcnt3 = 0;
    prcnt4 = 50;
    dataDir = sprintf('D:\\Sam\\Documents\\school\\THESIS\\Data\\Grad\\SAM_GUSTIN\\SAM_GUSTIN_ARPEL\\ARPEL_%01d_%01d_%01d_%01d_x_*mm_77V_16K\\PIV_MPd(4x16x16_75%%ov_ImgCorr)\\Avg Vel vector field', prcnt1, prcnt2, prcnt3, prcnt4);  
    %store locations of the numberOfSlicesices
    ySliceLocs = [0 10 22 30];
end

%Characteristic length (Used in "Build 3d Arrays" section)
L = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         LOAD RAW VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this loads the vectors into a 2d array by file and number
unmodifiedFlowField = loadarrayvec(dataDir,'B00001_AvgV.vc7');

%fix loading issue caused by a naming convention causing...
%the above to load two numberOfSlicesices (numberOfSlicesice 2 and 6) out of order
if strcmp(apparatus, 'PEL')
    temp = unmodifiedFlowField(6);
    temp2 = unmodifiedFlowField(2);
    unmodifiedFlowField(6) = temp;
    for i = 6:-1:3
        unmodifiedFlowField(i) = unmodifiedFlowField(i-1);
    end
    unmodifiedFlowField(2) = temp;
        %unmodifiedFlowField(1,1)(6) = temp2;
end

clearvars temp temp2
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CONSISTENCY CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%verify the right number of numberOfSlicesices loaded
if not(isequal(size(unmodifiedFlowField),size(ySliceLocs')))
    fprintf('the wrong number of files/folder was read.  the program found')
    size(unmodifiedFlowField)
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              RENAME AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = fieldnames(unmodifiedFlowField(1,1));
f{strmatch('y',f,'exact')} = 'z';
f{strmatch('vy',f,'exact')} = 'vz';
f{strmatch('unitvy',f,'exact')} = 'unitvz';
f{strmatch('namey',f,'exact')} = 'namez';
f{strmatch('namevy',f,'exact')} = 'namevz';
f{strmatch('ysign',f,'exact')} = 'zsign';

for i = 1:length(unmodifiedFlowField)
    v = struct2cell(unmodifiedFlowField(i));
    flowSlice = cell2struct(v,f);
    flowField(i) = flowSlice;
    flowField(i).zsign = 'Z axis upward';
    flowField(i).namevz = 'u_z';
    flowField(i).namez = 'z';
end
clearvars f v i unmodifiedFlowField flowSlice;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        REDUCE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%these are always the same so we can do them just once
flowField(1).x = reduceData(reductionStepSize,flowField(1).x);
flowField(1).y = ySliceLocs;
flowField(1).z = reduceData(reductionStepSize,flowField(1).z);

%these differ
% for i = 1:length(flowField)
%     flowField(i).vx = reduceData(reductionStepSize,flowField(i).vx);
%     flowField(i).vy = zeros(size(flowField(i).vx));
%     flowField(i).vz = reduceData(reductionStepSize,flowField(i).vz);
% end

for i = 1:length(flowField)
    flowField(i).vx = blockRMS(reductionStepSize,flowField(i).vx);
    flowField(i).vy = zeros(size(flowField(i).vx));
    flowField(i).vz = blockRMS(reductionStepSize,flowField(i).vz);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PREALLOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get dimensions of the velocity array
[primaryArrayRowNum, primaryArrayColNum] = size(flowField(1,1).vx);

%preallocate memory
numberOfSlices = length(ySliceLocs);

speedTensor = zeros(primaryArrayRowNum,primaryArrayColNum,numberOfSlices);        %688x550x6   
%transparencyTensor = zeros(primaryArrayRowNum,primaryArrayColNum,numberOfSlices);      %688x550x6
xVelocityTensor = zeros(primaryArrayRowNum,primaryArrayColNum,numberOfSlices);        %688x550x6
yVelocityTensor = zeros(primaryArrayRowNum,primaryArrayColNum,numberOfSlices);        %688x550x6
zVelocityTensor = zeros(primaryArrayRowNum,primaryArrayColNum,numberOfSlices);        %688x550x6

clearvars ySliceLocs primaryArrayRowNum primaryArrayColNum numberOfSlices;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                BUILD 3D SPEED, POSITION, AND VELOCITY FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get axis vectors
xAxis = flowField(1).x;
yAxis = flowField(1).y;
zAxis = flowField(1).z;

for i = 1:length(flowField)         
% L is the normalization factor for making the measurement dimensionless
    speedTensor(:,:,i)= (flowField(i).vx.^2+flowField(i).vz.^2).^(0.5)./L;% 
    xVelocityTensor(:,:,i) = flowField(i).vx./L;
    zVelocityTensor(:,:,i) = flowField(i).vy./L;
%	transparencyTensor(:,:,i) = speedTensor(:,:,i)>speedIsoplaneThreshold;
end

%permute to switch z and y
permOrder = [1 3 2]; %switch z and y
speedTensor = permute(speedTensor, permOrder);
xVelocityTensor = permute(xVelocityTensor, permOrder);
yVelocityTensor = permute(yVelocityTensor, permOrder);
zVelocityTensor = permute(zVelocityTensor, permOrder);
%transparencyTensor = permute(transparencyTensor, permOrder);
%dispose of the uneccessary overhead
clearvars flowSlice flowField permOrder;

if not(visualize)
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        INTERPOLATE IN Y AXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xref, yref, zref] = meshgrid(yAxis, xAxis, zAxis);

yAxis = [yAxis(1):interpStep:yAxis(end)];
[vxInterp, vyInterp, vzInterp] ...
        = meshgrid(xAxis,yAxis,zAxis);

speedTensor = interp3(xref,yref,zref, speedTensor, vyInterp, vxInterp, vzInterp, 'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FORMAT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f1,v1] = isosurface(xAxis, yAxis, zAxis, speedTensor, speedIsoplaneThreshold);
[f2,v2,e2] = isocaps(xAxis, yAxis, zAxis, speedTensor, speedIsoplaneThreshold);

figure('Renderer', 'painters', 'Position', [0 0 1200 1000]);
title('test');
ylim([yAxis(1)-1,yAxis(end)+1]);
zlim([-5,105]);
xlim([20, inf]);
colorbar
%caxis([0,1])
xlabel('x position [mm]');
ylabel('y position [mm]');
zlabel('z position [mm]');
%make a raw isosurface

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

clearvars *
return;
%h2 = streamline(xPositionTensor,y,zPositionTensor,yVelocityTensor, zVelocityTensor, xVelocityTensor, Sx, Sy,Sz);

%h2 = streamline(zPositionTensor,xPositionTensor,y,xVelocityTensor, ones(size(xVelocityTensor))*0, zVelocityTensor, Sx, Sy,Sz);
%h2 = streamline(xPositionTensor, z, yPositionTensor, yVelocityTensor, zVelocityTensor, xVelocityTensor, Sx, Sy,Sz);

%set the streamline color
set(h2, 'color', [1 0 0]);

axis tight equal
load('variables.mat') %load variables

%get dimensions of the velocity array
[primaryArrayRowNum, primaryArrayColNum] = size(first_v.vx);

%store locations of the slices
Y_slice_loc = [0 11 15 20 25];

%build coherent vectorized dataset
for i = 1:size(v)
    for j = 1:size(v(1))
        %assign the ariablename v_working to v(i,j)
        v_working = v(i,j);
        %reformat x and z locations
        x_working = repmat(v_working.x,[1,primaryArrayColNum]);1
        z_working = repelem(v_working.y',primaryArrayRowNum)';
        
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
load('variables.mat') %load variables

%get dimensions of the velocity array
[primaryArrayRowNum, primaryArrayColNum] = size(first_v.vx);

%store locations of the slices
Y_slice_loc = [0 11 15 20 25];

%build coherent vectorized dataset
for i = 1:size(v)
    for j = 1:size(v(1))
        %assign the ariablename v_working to v(i,j)
        v_working = v(i,j);
        %reformat x and z locations
        x_working = repmat(v_working.x,[1,primaryArrayColNum]);
        z_working = repelem(v_working.y',primaryArrayRowNum)';
        
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

