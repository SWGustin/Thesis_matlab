function result = blockRMS(step,data)

%look to implement rms rather than just picking a value
if ndims(data) == 2
    if size(data(1)) == 1 %1D
        result = zeros(size(data(1:step:end)));
        j = 1;
        for i = 1:step:length(data)
            noOfCells = min(length(data)-i-step,step);
            result(j) = (sum(data(i:i+noOfCells)).^2/noOfCells)^0.5;
            j = j + 1;
        end
    else
        disp 2d
        
    end
elseif ndims(data) == 3
   result = zeros(1:step:size(data,1),1:step:size(data,2),1:step:size(data,3));
   
end



