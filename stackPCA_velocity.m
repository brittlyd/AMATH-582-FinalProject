function [X,ind_nan_total,ind_nan_no_mask,count] = stackPCA_velocity(U,V,mask_ind)
%%% Subtract mean
mean_U = repmat(mean(U,1,'omitnan'),size(U,1),1);
mean_V = repmat(mean(V,1,'omitnan'),size(V,1),1);
U = U-mean_U;
V = V-mean_V;

%%% Reshape into a vector
U_vec = reshape(U,[],1);
V_vec = reshape(V,[],1);

%%% Find location of all non-NaN values in the array (for later)
ind_nan_total = isnan(U_vec);

%%% Remove indices for mask
% U_vec(mask_ind) = [];
% V_vec(mask_ind) = [];

count = 0;
for j = 1:length(mask_ind)
   if isnan(U_vec(mask_ind(j)))
       U_vec(mask_ind(j)) = [];
       count = count+1;
   end
end

%%% Replace NaNs with outliers
ind_nan_no_mask = isnan(U_vec);
U_vec(ind_nan_no_mask) = randi([5 10],1,1);
V_vec(ind_nan_no_mask) = randi([5 10],1,1);

%%% Concatenate U and V
X = [U_vec; V_vec];
end