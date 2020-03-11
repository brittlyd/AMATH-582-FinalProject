%Rename the fields from U and V to u and v
names=fieldnames(data);
for n=1:length(names)
    names2=lower(fieldnames(data.(names{n})));
    %replace the fieldnames
    data.(names{n})=cell2struct(struct2cell(data.(names{n})),names2);
    for k=length(names2)-2:length(names2)
        names3=lower(fieldnames(data.(names{n}).(names2{k})));
        %replace the fieldnames
        data.(names{n}).(names2{k})=...
            cell2struct(struct2cell(data.(names{n}).(names2{k})),names3);
    end
        
end
save('up1_1 Crop','data')