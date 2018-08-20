function [VIM, alphas] = dynGENIE3(TS_data,time_points,alpha,SS_data,input_idx,tree_method,K,ntrees)
%Computation of tree-based weights for all putative edges.
%
%
%[VIM, alphas] = dynGENIE3(TS_data,time_points) learns p tree models from
%the time series expression data contained in TS_data, where p is the
%number of columns (genes) in the matrices of TS_data. A weight is then
%assigned to each edge directed from any gene in TS_data to any other gene.
%   * TS_data is a 1-by-n cell array where each cell is a matrix containing
%   expression values related to one time series experiment. Each line of a
%   matrix corresponds to a time point and each column corresponds to a
%   gene. All matrices in TS_data must have the same number of
%   columns/genes.
%   * time_points is a 1-by-n cell array where the k-th cell is a vector
%   containing the time points corresponding to the k-th experiment/cell of
%   TS_data. 
%The function dynGENIE3() returns:
%   * VIM: a matrix of size p x p. VIM(i,j) is the weight of the edge
%   directed from the i-th gene of TS_data to the j-th gene. VIM(i,i) is
%   set to zero for all i.
%   * alphas: a vector in which the i-th element is the degradation rate
%   of the i-th gene (as ordered in TS_data).
%
%[VIM, alphas] = dynGENIE3(TS_data,time_points,alpha) learns the p tree
%models using the gene degradation rates specified by the input argument
%alpha. alpha can be either:
%   - 'from_data' (default value). In that case, the degradation rate of
%   each gene is estimated from the data, by assuming an exponential decay
%   between the highest and lowest observed expression values.
%   - a vector of positive numbers, where the i-th element is the
%   degradation rate of the i-th gene (as ordered in TS_data).
%   - a single positive number. In that case, all the genes are assumed to
%   have the same degradation rate alpha.
%
%[VIM, alphas] = dynGENIE3(TS_data,time_points,alpha,SS_data) learns the p
%tree models jointly from TS_data and SS_data. SS_data is a matrix
%containing gene expression values at steady-state level. Each line
%corresponds to an experiment and each column corresponds to a gene. The
%number of genes/columns in SS_data must be equal to the number of genes in
%TS_data. By default, SS_data = [] (empty matrix), i.e. no steady-state
%data are used.
%
%[VIM, alphas] = dynGENIE3(TS_data,time_points,alpha,SS_data,input_idx)
%only uses as candidate regulators the genes whose index (as ordered in
%TS_data) is in input_idx. input_idx is a vector of length <= p. VIM(i,:)
%such that i is not in input_idx is set to zero. The default vector
%contains the indices of all the genes in TS_data.
%
%[VIM, alphas] =
%dynGENIE3(TS_data,time_points,alpha,SS_data,input_idx,tree_method)
%specifies which tree precedure is used. Available methods are:
%   * 'RF' - Random Forests (Default method) 
%   * 'ET' - Extra Trees
%
%[VIM, alphas] =
%dynGENIE3(TS_data,time_points,alpha,SS_data,input_idx,tree_method,K)
%specifies the number K of randomly selected attributes at each node of one
%tree. Possible values of K:
%   * 'sqrt': K = square root of the number of candidate regulators
%   (Default value)
%   * 'all': K = number of candidate regulators
%   * any numerical value
%
%[VIM, alphas] =
%dynGENIE3(TS_data,time_points,alpha,SS_data,input_idx,tree_method,K,ntrees
%) specifies the number of trees grown in each ensemble. 
%Default value: 1000.

%%
tic;

%% Check input arguments
error(nargchk(2,8,nargin));

if ~iscell(TS_data) || (size(TS_data,1) ~=1 && size(TS_data,2) ~=1)
    error('Input argument TS_data must be a 1-by-n cell array where each element is a two-dimensional matrix corresponding to a time series experiment. The element(i,j) of each of these matrices is the expression of the j-th gene at the i-th time point.')
end

for k=1:length(TS_data)
    if length(size(TS_data{k})) ~= 2
        error('Input argument TS_data must be a 1-by-n cell array where each element is a two-dimensional matrix corresponding to a time series experiment. The element(i,j) of each of these matrices is the expression of the j-th gene at the i-th time point.')
    end
end

ngenes = size(TS_data{1},2);

for k=2:length(TS_data)
    
    if ngenes ~= size(TS_data{k},2)
        error('The number of columns/genes must be the same in every matrix of TS_data.')
    end
end

if ~iscell(time_points) || (size(time_points,1) ~=1 && size(time_points,2) ~=1)
    error('Input argument time_points must be a 1-by-n cell array, where n is the number of time series experiments in TS_data.')
end

if length(time_points) ~= length(TS_data)
    error('Input argument time_points must be a 1-by-n cell array, where n is the number of time series experiments in TS_data.')
end

for k=1:length(time_points)
    if ~isvector(time_points{k})
        error('Each array in time_points must be a vector.')
    end
end

for k=1:length(time_points)
    if length(time_points{k}) ~= size(TS_data{k},1)
        error('The length of the k-th vector of time_points must be equal to the number of rows in the i-th matrix of TS_data.')
    end
end

if nargin > 2

    if ~isscalar(alpha) && ~isvector(alpha) && ~isa(alpha,'char')
        error('Input argument alpha must be either "from_data", a vector of positive numbers or a single positive number.')
    end

    if isa(alpha,'char') && ~strcmp(alpha,'from_data')
        error('Input argument alpha must be either "from_data", a vector of positive numbers or a single positive number.')
    end

    if isscalar(alpha) && alpha < 0
        error('Input argument alpha must be either "from_data", a vector of positive numbers or a single positive number.')
    end

    if isvector(alpha) && isnumeric(alpha) && length(alpha) > 1
        if length(alpha) ~= ngenes
            error('When input argument alpha is a vector, this must be a vector of length p, where p is the number of genes.')
        end
        for i=1:length(alpha)
            if alpha(i) < 0
                error('Input argument alpha must be either "from_data", a vector of positive numbers or a single positive number.')
            end
        end
    end
end

if nargin > 3
    if ~isempty(SS_data)
        if length(size(SS_data)) ~= 2 || size(SS_data,1) < 2 || size(SS_data,2) < 2
            error('Input argument SS_data must be either an empty matrix or a two-dimensional matrix, where SS_data(i,j) is the steady-state expression of the j-th gene in the i-th condition.')
        end

        if size(SS_data,2) ~= ngenes
            error('The number of genes/colums of input argument SS_data must be equal to the number of genes/colums in the matrices of TS_data.')
        end
    end
end

if nargin > 4 && sum(ismember(input_idx,1:ngenes)) ~= length(input_idx)
    error('Input argument input_idx must be a vector containing integers between 1 and p, where p is the number of genes in TS_data.')
end

if nargin > 5 && sum(strcmp(tree_method,{'RF' 'ET'})) == 0
    error('Input argument tree_method must be ''RF'' (Random Forests) or ''ET'' (Extra-Trees).')
end

if nargin > 6 && ((isa(K,'char') && sum(strcmp(K,{'sqrt' 'all'})) == 0) || (isa(K,'numeric') && K <=0))
    error('Input argument K must be ''sqrt'', ''all'' or a strictly positive integer.')
end

if nargin > 7 && (~isa(ntrees,'numeric') || (isa(ntrees,'numeric') && ntrees <=0))
    error('Input argument ntrees must be a strictly positive integer.')
end

%% Print parameters
if nargin < 6
    print_method = 'RF';
else
    print_method = tree_method;
end
    
if nargin < 7
    print_K = 'sqrt';
else
    print_K = K;
    if isa(K,'numeric')
        print_K = num2str(K);
    end
end

if nargin < 8
    print_ntrees = 1000;
else
    print_ntrees = ntrees;
end

fprintf('Tree method: %s\n',print_method)
fprintf('K: %s\n',print_K)
fprintf('Number of trees: %d\n',print_ntrees)


%% Degradation rates
if nargin < 3
    alphas = estimate_decay_rates(TS_data,time_points);
elseif isscalar(alpha)
    alphas = zeros(1,ngenes) + alpha;
elseif isvector(alpha) && isnumeric(alpha)
    alphas = alpha;
else
    alphas = estimate_decay_rates(TS_data,time_points);
end

fprintf('alpha min: %.3e\n',min(alphas))
fprintf('alpha max: %.3e\n\n',max(alphas))

%% Learn an ensemble of trees for each gene and compute importance


VIM = zeros(ngenes,ngenes);

for i=1:ngenes
    fprintf('Gene %d/%d...\n',i,ngenes)
    
    if nargin < 4
        VIM(i,:) = dynGENIE3_single(TS_data,time_points,alphas(i),i);
    elseif nargin == 4
        VIM(i,:) = dynGENIE3_single(TS_data,time_points,alphas(i),i,SS_data);
    elseif nargin == 5
        VIM(i,:) = dynGENIE3_single(TS_data,time_points,alphas(i),i,SS_data,input_idx);
    elseif nargin == 6
        VIM(i,:) = dynGENIE3_single(TS_data,time_points,alphas(i),i,SS_data,input_idx,tree_method);
    elseif nargin == 7
        VIM(i,:) = dynGENIE3_single(TS_data,time_points,alphas(i),i,SS_data,input_idx,tree_method,K);
    else
        VIM(i,:) = dynGENIE3_single(TS_data,time_points,alphas(i),i,SS_data,input_idx,tree_method,K,ntrees);
    end
end

VIM = VIM';

%% 
toc;