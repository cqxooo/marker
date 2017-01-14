function K = calib_conic(imdata, K0)
cps = [];
for imLoop=1:length(imdata)-1
    cp  = getCpCenter(imdata(imLoop),imdata(imLoop+1),imLoop);
    cps = [cps cp];
end
[error, IAC, K] = calib_4dof(cps, K0);