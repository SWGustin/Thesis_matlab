%Constants

step = 20; %step size for reducing data size 

%location of data
datafilePath = "D:\Sam\Documents\school\Grad\SAM_GUSTIN\SAM_GUSTIN_THESIS_PILOT\top down\regular\100-0_PEL*\PIV_MPd(2x16x16_75%ov_ImgCorr)*\Avg_StdDev','B00001_AvgV.vc7";

%define our position arrays
X = [];
Y = [];
Z = [];

%define velocity arrays
vx = [];
vy = [];
vz = [];




%-------------------------------
%save the whole thing
save D:\Sam\Documents\school\Grad\THESIS\variables.mat