function vi = dynGENIE3_single(TS_data,time_points,alpha,output_idx,SS_data,input_idx,tree_method,K,ntrees)
%Computation of tree-based weights for putative edges directed towards a 
%specified target gene.
%
%vi = dynGENIE3_single(TS_data,time_points,alpha,output_idx) learns a tree
%model from the time series expression data contained in TS_data, and
%assigns a weight to each edge directed from a candidate regulator to the
%target gene. 
%   * TS_data is a 1-by-n cell array where each cell is a matrix containing
%   expression values related to one time series experiment. Each line of a
%   matrix corresponds to a time point and each column corresponds to a
%   gene. All matrices in TS_data must have the same number of
%   columns/genes.
%   * time_points is a 1-by-n cell array where the k-th cell is a vector
%   containing the time points corresponding to the k-th experiment/cell of
%   TS_data. 
%   * alpha is a positive number specifying the degradation rate of the
%   target gene
%   * output_idx is the (column) index of the target gene in
%   the matrices of TS_data.
%vi is a vector of length p, where p is the number of columns/genes in the
%matrices of TS_data. vi(i) is the weight of the edge directed from
%the i-th gene to the target gene. vi(output_idx) is set to zero.
%
%vi = dynGENIE3_single(TS_data,time_points,alpha,output_idx,SS_data) learns
%the tree model jointly from TS_data and SS_data. SS_data is a matrix
%containing gene expression values at steady-state level. Each line
%corresponds to an experiment and each column corresponds to a gene. The
%number of genes/columns in SS_data must be equal to the number of genes in
%TS_data. 
%By default, SS_data = [] (empty matrix), i.e. no steady-state data are
%used.
%
%vi =
%dynGENIE3_single(TS_data,time_points,alpha,output_idx,SS_data,input_idx)
%only uses as candidate regulators the genes whose index (as ordered in
%TS_data) is in input_idx. input_idx is a vector of length <= p. vi(i) such
%that i is not in input_idx is set to zero. The default vector contains the
%indices of all the genes in TS_data.
%
%vi =
%dynGENIE3_single(TS_data,time_points,alpha,output_idx,SS_data,input_idx,tr
%ee_method) specifies which tree precedure is used. Available methods are:
%   * 'RF' - Random Forests (Default method) 
%   * 'ET' - Extra Trees
%
%vi =
%dynGENIE3_single(TS_data,time_points,alpha,output_idx,SS_data,input_idx,tr
%ee_method,K) specifies the number K of randomly selected attributes at
%each node of one tree. Possible values of K:
%   * 'sqrt' - K = square root of the number of candidate regulators
%   (Default value)
%   * 'all' - K = number of candidate regulators
%   * any numerical value
%
%vi =
%dynGENIE3_single(TS_data,time_points,alpha,output_idx,SS_data,input_idx,tr
%ee_method,K,ntrees) specifies the number of trees grown in the ensemble.
%Default value: 1000.


%% Check input arguments
error(nargchk(4,9,nargin));

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

if ~isscalar(alpha)
    error('Input argument alpha must be a single number.')
end

if alpha < 0
    error('Input argument alpha must be a positive number.')
end

if ~isscalar(output_idx)
    error('Input argument output_idx must be one integer.')
end

if ~ismember(output_idx,1:ngenes)
    error('Input argument output_idx must be an integer between 1 and p, where p is the number of genes in TS_data.')
end

if nargin > 4
    if ~isempty(SS_data)
        if length(size(SS_data)) ~= 2 || size(SS_data,1) < 2 || size(SS_data,2) < 2
            error('Input argument SS_data must be either an empty matrix or a two-dimensional matrix, where SS_data(i,j) is the steady-state expression of the j-th gene in the i-th condition.')
        end

        if size(SS_data,2) ~= ngenes
            error('The number of genes/colums of input argument SS_data must be equal to the number of genes/colums in the matrices of TS_data.')
        end
    end
end

if nargin > 5 && sum(ismember(input_idx,1:ngenes)) ~= length(input_idx)
    error('Input argument input_idx must be a vector containing integers between 1 and p, where p is the number of genes in TS_data.')
end

if nargin > 6 && sum(strcmp(tree_method,{'RF' 'ET'})) == 0
    error('Input argument tree_method must be ''RF'' (Random Forests) or ''ET'' (Extra-Trees).')
end

if nargin > 7 && ((isa(K,'char') && sum(strcmp(K,{'sqrt' 'all'})) == 0) || (isa(K,'numeric') && K <=0))
    error('Input argument K must be ''sqrt'', ''all'' or a strictly positive integer.')
end

if nargin > 8 && (~isa(ntrees,'numeric') || (isa(ntrees,'numeric') && ntrees <=0))
    error('Input argument ntrees must be a strictly positive integer.')
end


%% Re-order time points in increasing order
for k=1:length(time_points)
    tp = time_points{k};
    [tp_sort,order] = sort(tp);
    time_points{k} = tp_sort;
    expr_data = TS_data{k};
    TS_data{k} = expr_data(order,:);
end


%% Indices of input genes
if nargin >= 6
    input_idx = unique(input_idx);
else
    % Default: all genes are candidate regulators
    input_idx = 1:ngenes;
end

ninputs = length(input_idx);

%% Construct learning sample

h = 1;

% Time series data
nexp = length(TS_data);

nsamples_time = 0;
for k=1:nexp
    nsamples_time = nsamples_time + size(TS_data{k},1);
end

input_matrix_time = zeros(nsamples_time-h*nexp,ninputs);
output_vect_time = zeros(nsamples_time-h*nexp,1);

nsamples_count = 0;
for k=1:nexp
    current_time_points = time_points{k};
    current_timeseries = TS_data{k};
    time_diff_current = current_time_points(h+1:end) - current_time_points(1:end-h);
    current_timeseries_input = current_timeseries(1:end-h,input_idx);
    current_timeseries_output = (current_timeseries(h+1:end,output_idx) - current_timeseries(1:end-h,output_idx)) ./ time_diff_current' + alpha*current_timeseries(1:end-h,output_idx);
    nsamples_current = size(current_timeseries_input,1);
    input_matrix_time(nsamples_count+1:nsamples_count+nsamples_current,:) = current_timeseries_input;
    output_vect_time(nsamples_count+1:nsamples_count+nsamples_current) = current_timeseries_output;
    nsamples_count = nsamples_count + nsamples_current;
end


% Add steady-state data (if any)
if nargin >= 5 && ~isempty(SS_data)
    input_all = [SS_data(:,input_idx);input_matrix_time];
    output_all = [SS_data(:,output_idx)*alpha;output_vect_time];
else
    input_all = input_matrix_time;
    output_all = output_vect_time;
end

clear input_matrix_time output_vect_time

%% Tree parameters

% Default parameters: Random Forests, K=sqrt(number of input genes),
% 1000 trees in the ensemble
if nargin < 7 || (nargin >= 7 && strcmp(tree_method,'RF'))
    ok3ensparam = init_rf();
    if nargin >= 8
        if strcmp(K,'all')
            ok3ensparam=init_rf(ninputs);
        elseif isa(K,'numeric')
            ok3ensparam=init_rf(round(K));
        end
    end
elseif nargin >= 7 && strcmp(tree_method,'ET')
    ok3ensparam = init_extra_trees();
    if nargin >= 8
        if strcmp(K,'all')
            ok3ensparam=init_extra_trees(ninputs);
        elseif isa(K,'numeric')
            ok3ensparam=init_extra_trees(round(K));
        end
    end
end

% Number of trees in the ensemble
if nargin < 9
    ok3ensparam.nbterms = 1000;
else
    ok3ensparam.nbterms = ntrees;
end


%% Learning of tree model
[tree,varimp]=rtenslearn_c(single(input_all),single(output_all),int32(1:length(output_all)),[],ok3ensparam,0);
% Some variable importance might be slightly negative due to some rounding
% error
varimp(varimp<0) = 0;
vi = zeros(1,ngenes);
vi(input_idx) = varimp';
vi(output_idx) = 0;

% Normalize variable importances
vi = vi / sum(vi);
