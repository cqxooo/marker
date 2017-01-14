function [E Proj] = getE_Proj(data, K)
for imLoop = 1:length(data)-1
    [lcor rcor] = findcorres(data(imLoop), data(imLoop+1));
    E(:,:,imLoop) = getE(lcor, rcor, K);
%    [E(:,:,imLoop), lcor, rcor] = ransacE(lcor, rcor, K);
   if imLoop == 1
       [Proj(:,:,imLoop) Proj(:,:,imLoop+1)] = getProj(lcor, rcor, E(:,:,imLoop), K);
   else
       [~, Proj(:,:,imLoop+1)] = getProj(lcor, rcor, E(:,:,imLoop), K, Proj(:,:,imLoop));
   end
end
