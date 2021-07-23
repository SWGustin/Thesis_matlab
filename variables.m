%Constants

step = 100;


%this loads the vectors into a 2d array by file and number
filePath = "D:\Sam\Documents\school\Grad\SAM_GUSTIN\SAM_GUSTIN_THESIS_PILOT\top down\regular\100-0_PEL*\PIV_MPd(2x16x16_75%ov_ImgCorr)*\Avg_StdDev','B00001_AvgV.vc7";

%define our arrays
X = [];
Y = [];
Z = [];

vx = [];
%vy = [];
vz = [];

save variables.mat