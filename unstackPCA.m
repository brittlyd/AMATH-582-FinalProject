%%% This function takes in a column vector X, original array(s), u (and v),
%%% and a vector of indices indicating the location of the mask and returns
%%% new array(s), u_new (and v_new), of the same size of u (and v).

function [u_new,v_new] = unstackPCA(X,nx,ny,mask_log, vExist)
if vExist %%% if there are two arrays in X
    %%% Determine length of each array
    n = length(X)/2;

    %%% Extract data
    u_vec = X(1:n);
    v_vec = X(n+1:end);

    %%% Create new vectors of full size
    u_new_vec = zeros(nx*ny,1);
    v_new_vec = zeros(nx*ny,1);

    %%% Reshape mask logicals into a column vector
    mask_log=reshape(mask_log,nx*ny,1);
    
    %%% Replace mask indices with NaNs
    u_new_vec(mask_log) = NaN;
    v_new_vec(mask_log) = NaN;

    %%% Fill in the rest of the vectors with data
    u_new_vec(~isnan(u_new_vec)) = u_vec;
    v_new_vec(~isnan(u_new_vec)) = v_vec;

    %%% Reshape back into arrays of same size as u and v
    u_new = reshape(u_new_vec,[nx ny]);
    v_new = reshape(v_new_vec,[nx ny]);

else %%% if there is only one array in X
    %%% Do the same thing but for one array
    u_new_vec = zeros(nx*ny,1);
    u_new_vec(mask_log) = NaN;
    u_new_vec(~isnan(u_new_vec)) = X;
    u_new = reshape(u_new_vec,[nx ny]);
    v_new = [];
end