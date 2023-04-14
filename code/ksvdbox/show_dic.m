function [] = show_dic(opt)
dic = opt.dic;
faceW = opt.blocksize; 
faceH = size(dic,1)/faceW; 
numPerLine = 17; 
ShowLine = 17; 

Y = zeros(faceH*ShowLine,faceW*numPerLine); 
idx = 1;
for i=0:ShowLine-1 
    for j=0:numPerLine-1 
         Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(dic(:,idx),[faceH,faceW])'; 
         idx = idx + 1;
    end 
end 

imagesc(Y);
