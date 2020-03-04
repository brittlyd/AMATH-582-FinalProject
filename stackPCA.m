%%% This function takes in one or two arrays (u,v) as well as a 
%%% vector of linear indices for the mask on the frame (mask_ind) and
%%% outputs a column vector, X, that contains the reshaped data from the u
%%% (and v) array(s) where the indices containing the mask are removed and
%%% the remaining NaNs are replaced with outliers (very large, positive
%%% values). 

function [X] = stackPCA(u,v,mask_ind)
%%% Subtract mean of the frame from each element
mean_u = repmat(mean(u,'all','omitnan'),size(u));
u = u-mean_u;

%%% Reshape into a column vector
u_vec = reshape(u,[],1);

%%% Remove indices for mask
u_vec(mask_ind) = [];

%%% Replace remaining NaNs with outliers
ind_nan_no_mask_u = isnan(u_vec);
u_vec(ind_nan_no_mask_u) = randi([5 10],sum(ind_nan_no_mask_u),1);

if isempty(v) %%% if there is no second array
    X = u_vec; %%% return the vector 
else 
    %%% Do the same thing for v if it exists
    mean_v = repmat(mean(v,'all','omitnan'),size(v));
    v = v-mean_v;
    
    v_vec = reshape(v,[],1);
    v_vec(mask_ind) = [];
    
    ind_nan_no_mask_v = isnan(v_vec);
    v_vec(ind_nan_no_mask_v) = randi([5 10],sum(ind_nan_no_mask_v),1);
    
    %%% Concatenate U and V
    X = [u_vec; v_vec];
end