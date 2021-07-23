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
        x_working = repmat(v_working.x,[1,ncols]);
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
