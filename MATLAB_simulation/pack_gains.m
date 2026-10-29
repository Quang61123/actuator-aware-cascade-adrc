function G=pack_gains(cells)
classes=[1 1 2 3 3 4];G=zeros(37,1);
for j=1:6,G(6*(j-1)+(1:6))=cells{classes(j)};end
G(37)=1;
end
