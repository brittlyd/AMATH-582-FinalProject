%Rename the fields from U and V to u and v
names=fieldnames(down1_1FULL);
for n=1:length(names)
    names2=lower(fieldnames(down1_1FULL.(names{n})));
    %replace the fieldnames
    down1_1FULL.(names{n})=cell2struct(struct2cell(down1_1FULL.(names{n})),names2);
end
save('down1_1FULL','down1_1FULL')