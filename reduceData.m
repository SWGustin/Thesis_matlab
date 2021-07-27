function reduced = reduceData(step,data)
if ndims(data) == 3
    reduced = data(1:step:end,:,1:step:end);
elseif ndims(data) == 2
    reduced = data(1:step:end,1:step:end);
elseif ndims(data) == 1
    reduced = data(1:step:end);
end
            