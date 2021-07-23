clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%specify which apparatus you're using
%apparatus = 'ARPEL';
apparatus = 'PEL';

speed_threshold = .5;

%define the step value to reduce data size 
%(data will be 1/step in each dimension)
step = 15;

if strcmp(apparatus,'PEL')
    fprintf('running for pel\n')
    
    %percentage of interest for PEL
    %valid percentag65
    percent = 100;
    recip_prcnt = 100-percent;
    
    %locations of files with wildcard
    dirs = sprintf('D:\\Sam\\Documents\\school\\Grad\\THESIS\\Data\\SAM_GUSTIN_PEL\\PEL_%01d-%01d_x_*mm_77V_16K\\PIV_MPd(4x16x16_75%%ov_ImgCorr)\\Avg Vel vector field',percent,recip_prcnt);

    %store locations of the slices
    Y_slice_loc = [0 5 10 15 20 25];
    
elseif strcmp(apparatus, 'ARPEL')
    fprintf('running for arpel\n'); 
    %percentage of interest for ARPEL
    prcnt1 = 50;
    prcnt2 = 0;
    prcnt3 = 0;
    prcnt4 = 50;

    dirs = sprintf('D:\\Sam\\Documents\\school\\THESIS\\Data\\Grad\\SAM_GUSTIN\\SAM_GUSTIN_ARPEL\\ARPEL_%01d_%01d_%01d_%01d_x_*mm_77V_16K\\PIV_MPd(4x16x16_75%%ov_ImgCorr)\\Avg Vel vector field', prcnt1, prcnt2, prcnt3, prcnt4);
    
    %store locations of the slices
    Y_slice_loc = [0 10 22 30];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         LOAD RAW VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this loads the vectors into a 2d array by file and number
v_primary = loadarrayvec(dirs,'B00001_AvgV.vc7','verbose');

%fix loading issue
if strcmp(apparatus, 'PEL')
    temp = v_primary(6);
    temp2 = v_primary(2);
    v_primary(6) = temp;
    for i = 6:-1:3
        v_primary(i) = v_primary(i-1);
    end
    v_primary(2) = temp;
        %v_primary(6) = temp2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CONSISTENCY CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%verify the right number of slices loaded
if not(isequal(size(v_primary),size(Y_slice_loc')))
    fprintf('the wrong number of files/folder was read.  the program found')
    size(v_primary)
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PREALLOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get dimensions of the velocity array
[nrows, ncols] = size(v_primary(1,1).vx);

%preallocate memory
sL = size(Y_slice_loc);
speedAmalg = zeros(nrows,ncols,sL(2));
transparency = zeros(nrows,ncols,sL(2));
x_amalg = zeros(nrows,ncols,sL(2));
z_amalg = zeros(nrows,ncols,sL(2));
vx_amalg = zeros(nrows,ncols,sL(2));
vz_amalg = zeros(nrows,ncols,sL(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       MANUAL MESHGRID EQUIVALENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b, c] = meshgrid(1:2,1:2,4:7);


vect = v_primary(1,1);
%this only needs to be done once
%reformat x & z locations
x = repmat(vect.x,[1,ncols]);
x = reshape(x,[nrows,ncols]);
z = repelem(vect.y',nrows)';
z = reshape(z, [nrows, ncols]);
y = zeros(nrows,ncols,sL(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                BUILD 3D SPEED, POSITION, AND VELOCITY FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(v_primary)
    for j = 1:size(v_primary(1))
        vect = v_primary(i,j);
        
        speedAmalg(:,:,i)= (vect.vx.^2+vect.vy.^2).^(0.5)/1;% 1 is the normalization factor for making the measurement dimensionless
        
        x_amalg(:,:,i) = x;
        z_amalg(:,:,i) = z;
        vx_amalg(:,:,i) = vect.vx;
        vz_amalg(:,:,i) = vect.vy;
        
        y(:,:,i) = ones(size(speedAmalg(:,:,i))).*Y_slice_loc(i);

        transparency(:,:,i) = speedAmalg(:,:,i)>speed_threshold;

   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        REDUCE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%y_reduced = y(1:step:end,1:step:end,:);
%spd = speedAmalg(1:step:end,1:step:end,:);
%x_reduced = x_amalg(1:step:end,1:step:end,:);
%z_reduced = z_amalg(1:step:end,1:step:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     CALC TRANSPARENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trans = transparency(1:step:end,1:step:end,size(v_primary));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        INTERPOLATE IN Y AXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interp_step = 1;

refV = v_primary(1,1);
if strcmp(apparatus,'PEL')
   desired_y_vect = Y_slice_loc(1):interp_step:Y_slice_loc(6);
elseif strcmp(apparatus, 'ARPEL')
    desired_y_vect = Y_slice_loc(1):interp_step:Y_slice_loc(4);
end

[Xq,Yq,Zq] = meshgrid(refV.x, desired_y_vect, refV.y);
[Xi,Yi,Zi] = meshgrid(refV.x, refV.y, Y_slice_loc);

Xi = unique(Xi);
Yi = unique(Yi);
Zi = unique(Zi);

%Xi = refV.x';
%Yi = refV.y';
%Zi = Y_slice_loc;

speedAmalg = permute(speedAmalg, [2 1 3]);

interp_speed = interp3(Xi, Yi, Zi, speedAmalg, Xq, Yq, Zq, 'spline');
[Xi,Yi,Zi] = size(interp_speed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FORMAT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [0 0 1200 1000]);
hold on
tset = 'test';
title(tset);
%ylim([min(Y_slice_loc),max(Y_slice_loc)]);
%zlim([-5,105]);
%xlim([30, inf]);
colorbar
caxis([0,1])
xlabel('x position [mm]');
ylabel('y position [mm]');
zlabel('z position [mm]');

%make a raw isosurface
speedAmalg = permute(speedAmalg,[2 1 3]);
interp_speed = permute(interp_speed, [2 3 1]);
[Xq,Yq,Zq] = meshgrid(refV.y, refV.x, desired_y_vect);

[f1,v1] = isosurface(Xq, Yq, Zq, interp_speed,speed_threshold);
[f2,v2,e2] = isocaps(Xq,Yq,Zq,interp_speed,speed_threshold);

%[f1,v1] = isosurface(x_amalg, y, z_amalg, speedAmalg,speed_threshold);
%[f2,v2,e2] = isocaps(x_amalg,y,z_amalg,speedAmalg,speed_threshold);
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

x_amalg = unique(x_amalg);
z_amalg = unique(z_amalg);
y = unique(y);
z=y;
y_amalg = z_amalg;

vx_amalg = permute(vx_amalg, [2 3 1]);
vy_amalg = permute(vy_amalg, [2 3 1]);
vz_amalg = permute(vz_amalg, [2 3 1]);
disp("-------")
size(x_amalg)
size(y_amalg)
size(z)
size(vx_amalg)
size(vy_amalg)
size(vz_amalg)
%size(y)

%h2 = streamline(x_amalg,y,z_amalg,vy_amalg, vz_amalg, vx_amalg, Sx, Sy,Sz);

%h2 = streamline(z_amalg,x_amalg,y,vx_amalg, ones(size(vx_amalg))*0, vz_amalg, Sx, Sy,Sz);
h2 = streamline(x_amalg, z, y_amalg, vy_amalg, vz_amalg, vx_amalg, Sx, Sy,Sz);

%set the streamline color
set(h2, 'color', [1 0 0]);

axis tight equal