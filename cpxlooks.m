function mlookigram=cpxlooks(igram,nlooks)

nr=size(igram,2);
nr_ml=floor(nr/nlooks);

for i=1:nr_ml
    k=(i-1)*nlooks+1:i*nlooks;
    temp=igram(:,k);
    mlookigram(:,i)=sum(temp,2);
end
mlookigram=mlookigram/nlooks;
return