function [U_new,V_new] = unstackPCA_velocity(X,U,V,mask_ind,ind_nan_total,ind_nan_no_mask)
n = length(X)/2;

U_vec = X(1:n);
V_vec = X(n+1:end);

U_new_vec = zeros(size(U,1)*size(U,2),1);
V_new_vec = U_new_vec;

U_new_vec(ind_nan_total) = NaN;
V_new_vec(ind_nan_total) = NaN;

U_new_vec(~isnan(U_new_vec)) = U_vec(~ind_nan_no_mask);
V_new_vec(~isnan(U_new_vec)) = V_vec(~ind_nan_no_mask);

U_new = reshape(U_new_vec,size(U));
V_new = reshape(V_new_vec,size(V));

end