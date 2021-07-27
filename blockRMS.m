function result = blockRMS(step,data)

    [xLimIn, yLimIn] =  size(data);
    result = zeros(ceil(xLimIn/step),ceil(yLimIn/step));
    [xLimOut, yLimOut] =  size(result);
    for x = 1:xLimOut
        for y = 1:yLimOut
                xmin = (x-1)*step+1;
                xmax = min(x*step,xLimIn);
                ymin = (y-1)*step+1;
                ymax = min(y*step,yLimIn);
                result(x,y) = (sum(sum(data(xmin:xmax,ymin:ymax).^2))...
                /numel(data(xmin:xmax,ymin:ymax)))^0.5;
        end
    end
end