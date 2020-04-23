function PatchVec = Patching(im,PatchM)
%[Nim,Mim] = size(im);
%imPlusN = PatchM - mod(Nim,PatchM);
%imPlusM = PatchM - mod(Mim,PatchM);
%im = [zeros(Nim,floor(imPlusN/2)) im zeros(Nim,ceil(imPlusN/2))];
%im = [zeros(floor(imPlusM/2),Mim+imPlusM); im ;zeros(ceil(imPlusM/2),Mim+imPlusM)];

imSz = size(im);
patchSz = [PatchM PatchM];
xIdxs = [1:patchSz(2):imSz(2) imSz(2)+1];
yIdxs = [1:patchSz(1):imSz(1) imSz(1)+1];
patches = cell(length(yIdxs)-1,length(xIdxs)-1);
for i = 1:length(yIdxs)-1
    Isub = im(yIdxs(i):yIdxs(i+1)-1,:);
    for j = 1:length(xIdxs)-1
        patches{i,j} = Isub(:,xIdxs(j):xIdxs(j+1)-1);
    end
end
PatchVec = zeros(patchSz(1)*patchSz(2),size(patches,1)*size(patches,2));
for i = 1:size(PatchVec,2)
    PatchVec(:,i) = reshape(patches{i},size(PatchVec,1),1);
end

