%%% This function takes in one or two arrays (u,v) as well as a 
%%% vector of linear indices for the mask on the frame (mask_ind) and
%%% outputs a column vector, X, that contains the reshaped data from the u
%%% (and v) array(s) where the indices containing the mask are removed and
%%% the remaining NaNs are replaced with outliers (very large, positive
%%% values). 

function [X,mask_log] = stackPCA2(u,v,nx,ny,mask_ind,fl)
%%% Turns mask_ind into logical. Need to reshpape the mask_ind vector as ny
%%% by nx. Also need to transpose. This works because of the way
%%% matlab indexes matricies and how I wrote everything in my PIV
%%% processing code. Commented out, but can plot after to show the mask (yellow area). The result
%%% mask_log is a matrix of logicals that is nx by ny which is the same size as the
%%% velocity fields


mask_log=zeros(nx*ny,1);
mask_log(mask_ind)=1;
mask_log=logical(reshape(mask_log,[ny,nx])');
if fl
    mask_log=fliplr(mask_log);
end
%figure
% pcolor(mask_log')
% axis equal
% figure
% pcolor(u')
% axis equal

%%% Remove indices for mask
u(mask_log) = [];

%%% Reshape into a column vector
u_vec = reshape(u,[],1);

%%% Replace remaining NaNs with outliers
ind_nan_no_mask_u = isnan(u_vec);
p = -1 + 2*rand(sum(ind_nan_no_mask_u),1);
u_vec(ind_nan_no_mask_u) = p;

if isempty(v) %%% if there is no second array
    X = u_vec; %%% return the vector 
else 
    %%% Do the same thing for v if it exists
    
    %%% Remove indices for mask
    v(mask_log) = [];

    v_vec = reshape(v,[],1);
    
    ind_nan_no_mask_v = isnan(v_vec);
    q = -1 + 2*rand(sum(ind_nan_no_mask_v),1);
    v_vec(ind_nan_no_mask_v) = q;
    
    %%% Concatenate U and V
    X = [u_vec; v_vec];
end