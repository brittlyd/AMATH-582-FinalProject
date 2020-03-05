%Rename the fields from U and V to u and v
names=fieldnames(up1_1FULL);
for n=1:length(names)
    names2=lower(fieldnames(up1_1FULL.(names{n})));
    %replace the fieldnames
    up1_1FULL.(names{n})=cell2struct(struct2cell(up1_1FULL.(names{n})),names2);
end
