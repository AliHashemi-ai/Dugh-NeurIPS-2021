function Bnorm = blobnorm(B,varargin)
% Bnorm = blobnorm(B,groups,options)
% Normalize the spatial basis functions used in SFLEX, ARD, etc...
%  Input:
%       B -> NdxNk (number of spatial basis functions). Matrix containing
%               the spatial basis function in each of the columns.
%       groups -> 1xNk. Vector indexing the groups of the basis functions.
%               For example: if the first 2 columns are basis with variance
%               0.1 the next 3 columns with variance 0.2 and the last one
%               with variance 0.5. then 'groups' would be: [1 1 2 2 2 3]
%       options.norm -> string. The normalization will be with respect to this norm
%               The naming is the same as the second argument of matlab's "norm"
%               function.
%       options.norm_group -> Bool. Normalize with respect to groups (true)
%               or normalize each column separately (false).
%  Output:
%       Bnorm -> NdxNk. Normalized basis functions.
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 15 Aug 2013


p = inputParser;

def_norm = 1;
def_ngroup = false;
def_groups = [];
addParamValue(p,'norm',def_norm);
addParamValue(p,'ngroup',def_ngroup);
addParamValue(p,'groups',def_groups);

parse(p,varargin{:})
options = p.Results;


if options.ngroup
    c = unique(options.groups);
    Bnorm = [];
    for i=1:length(c);
        idx_cur = find(options.groups == c(i));
        Bact = B(:,idx_cur);
        Btemp = blobnorm(Bact,'norm',options.norm,'ngroup',false);
        Bnorm = [Bnorm Btemp];
    end        
else 
    if sum(options.norm == 1) == 1 % Norm l1
        norm_term = sum(abs(B));
        norm_term = repmat(norm_term,size(B,1),1);
    elseif sum(options.norm == 2) == 1 % Norm l2
        norm_term = sqrt(sum(B.^2));
        norm_term = repmat(norm_term,size(B,1),1);
    else
        try
            for i=1:size(B,2)
                norm_term(1,i) = norm(B(:,i),options.norm);
            end
            norm_term = repmat(norm_term,size(B,1),1);
        catch
            error('Invalid norm selected')
        end
    end
    Bnorm = B./norm_term;
end