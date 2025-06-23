from sklearn.tree import BaseDecisionTree
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
import numpy as np
import time
from operator import itemgetter
from multiprocessing import Pool
from itertools import combinations

def compute_feature_importances(estimator):
    
    """Computes variable importances from a trained tree-based model.
    """
    
    if isinstance(estimator, BaseDecisionTree):
        return estimator.tree_.compute_feature_importances(normalize=False)
    else:
        importances = [e.tree_.compute_feature_importances(normalize=False)
                       for e in estimator.estimators_]
        importances = np.array(importances)
        return np.sum(importances,axis=0) / len(estimator)



def get_link_list(VIM,gene_names=None,regulators='all',maxcount='all',file_name=None):
    
    """Gets the ranked list of (directed) regulatory links.
    
    Parameters
    ----------
    
    VIM: numpy array
        Array as returned by the function dynGENIE3(), in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene. 
        
    gene_names: list of strings, optional
        List of length p, where p is the number of rows/columns in VIM, containing the names of the genes. The i-th item of gene_names must correspond to the i-th row/column of VIM. When the gene names are not provided, the i-th gene is named Gi.
        default: None
        
    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names), and the returned list contains only edges directed from the candidate regulators. When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
        
    maxcount: 'all' or positive integer, optional
        Writes only the first maxcount regulatory links of the ranked list. When maxcount is set to 'all', all the regulatory links are written.
        default: 'all'
        
    file_name: string, optional
        Writes the ranked list of regulatory links in the file file_name.
        default: None
        
        
    
    Returns
    -------
    
    The list of regulatory links, ordered according to the edge score. Auto-regulations do not appear in the list. Regulatory links with a score equal to zero are randomly permuted. In the ranked list of edges, each line has format:
        
        regulator   target gene     score of edge
    """
    
    # Check input arguments      
    if not isinstance(VIM,np.ndarray):
        raise ValueError('VIM must be a square array')
    elif VIM.shape[0] != VIM.shape[1]:
        raise ValueError('VIM must be a square array')
        
    ngenes = VIM.shape[0]
        
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
        
    if regulators != 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')
        
    if maxcount != 'all' and not isinstance(maxcount,int):
        raise ValueError('input argument maxcount must be "all" or a positive integer')
        
    if file_name is not None and not isinstance(file_name,str):
        raise ValueError('input argument file_name must be a string')
    
    

    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = list(range(ngenes))
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]
               
    nTFs = len(input_idx)
    
    # Get the non-ranked list of regulatory links
    vInter = [(i,j,score) for (i,j),score in np.ndenumerate(VIM) if i in input_idx and i!=j]
    
    # Rank the list according to the weights of the edges        
    vInter_sort = sorted(vInter,key=itemgetter(2),reverse=True)
    nInter = len(vInter_sort)
    
    # Random permutation of edges with score equal to 0
    flag = 1
    i = 0
    while flag and i < nInter:
        (TF_idx,target_idx,score) = vInter_sort[i]
        if score == 0:
            flag = 0
        else:
            i += 1
            
    if not flag:
        items_perm = vInter_sort[i:]
        items_perm = np.random.permutation(items_perm)
        vInter_sort[i:] = items_perm
        
    # Write the ranked list of edges
    nToWrite = nInter
    if isinstance(maxcount,int) and maxcount >= 0 and maxcount < nInter:
        nToWrite = maxcount
        
    if file_name:
    
        outfile = open(file_name,'w')
    
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('%s\t%s\t%.6f\n' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('G%d\tG%d\t%.6f\n' % (TF_idx+1,target_idx+1,score))
            
        
        outfile.close()
        
    else:
        
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('%s\t%s\t%.6f' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('G%d\tG%d\t%.6f' % (TF_idx+1,target_idx+1,score))
                
                
def estimate_degradation_rates(TS_data,time_points):
    
    """
    For each gene, the degradation rate is estimated by assuming that the gene expression x(t) follows:
    x(t) =  A exp(-alpha * t) + C_min,
    between the highest and lowest expression values.
    C_min is set to the minimum expression value over all genes and all samples.
    """
    
    ngenes = TS_data[0].shape[1]
    nexp = len(TS_data)
    
    C_min = TS_data[0].min()
    if nexp > 1:
        for current_timeseries in TS_data[1:]:
            C_min = min(C_min,current_timeseries.min())
    
    alphas = np.zeros((nexp,ngenes))
    
    for (i,current_timeseries) in enumerate(TS_data):
        current_time_points = time_points[i]
        
        for j in range(ngenes):
            
            idx_min = np.argmin(current_timeseries[:,j])
            idx_max = np.argmax(current_timeseries[:,j])
            
            xmin = current_timeseries[idx_min,j]
            xmax = current_timeseries[idx_max,j]

            if xmin != xmax:
            
                tmin = current_time_points[idx_min]
                tmax = current_time_points[idx_max]

                xmin = max(xmin-C_min,1e-6)
                xmax = max(xmax-C_min,1e-6)

                xmin = np.log(xmin)
                xmax = np.log(xmax)

                alphas[i,j] = (xmax - xmin) / abs(tmin - tmax)
                
    alphas = alphas.max(axis=0)
 
    return alphas
        

    
         


def dynGENIE3(TS_data,time_points,alpha='from_data',SS_data=None,gene_names=None,regulators='all',tree_method='RF',K='sqrt',ntrees=1000,compute_quality_scores=False,save_models=False,nthreads=1):

    '''Computation of tree-based scores for all putative regulatory links.

    Parameters
    ----------

    TS_data: list of numpy arrays
        List of arrays, where each array contains the gene expression values of one time series experiment. Each row of an array corresponds to a time point and each column corresponds to a gene. The i-th column of each array must correspond to the same gene.
    
    time_points: list of one-dimensional numpy arrays
        List of n vectors, where n is the number of time series (i.e. the number of arrays in TS_data), containing the time points of the different time series. The i-th vector specifies the time points of the i-th time series of TS_data.
    
    alpha: either 'from_data', a positive number or a vector of positive numbers
        Specifies the degradation rate of the different gene expressions. 
        When alpha = 'from_data', the degradation rate of each gene is estimated from the data, by assuming an exponential decay between the highest and lowest observed expression values.
        When alpha is a vector of positive numbers, the i-th element of the vector must specify the degradation rate of the i-th gene.
        When alpha is a positive number, all the genes are assumed to have the same degradation rate alpha.
        default: 'from_data'
    
    SS_data: numpy array, optional
        Array containing steady-state gene expression values. Each row corresponds to a steady-state condition and each column corresponds to a gene. The i-th column/gene must correspond to the i-th column/gene of each array of TS_data.
        default: None

    gene_names: list of strings, optional
        List of length p containing the names of the genes, where p is the number of columns/genes in each array of TS_data. The i-th item of gene_names must correspond to the i-th column of each array of TS_data (and the i-th column of SS_data when SS_data is not None).
        default: None
    
    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names). When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'
    
    tree-method: 'RF' or 'ET', optional
        Specifies which tree-based procedure is used: either Random Forest ('RF') or Extra-Trees ('ET')
        default: 'RF'
    
    K: 'sqrt', 'all' or a positive integer, optional
        Specifies the number of selected attributes at each node of one tree: either the square root of the number of candidate regulators ('sqrt'), the number of candidate regulators ('all'), or any positive integer.
        default: 'sqrt'
     
    ntrees: positive integer, optional
        Specifies the number of trees grown in an ensemble.
        default: 1000

    compute_quality_scores: boolean, optional
        Indicates if the scores assessing the edge ranking quality must be computed or not. These scores are:
        - the score of prediction of out-of-bag samples, i.e. the Pearson correlation between the predicted and true output values. To be able to compute this score, Random Forests must be used (i.e. parameter tree_method must be set to 'RF').
        - the stability score, measuring the average stability among the top-5 candidate regulators returned by each tree of a forest.
        default: False

    save_models: boolean, optional
        Indicates if the tree models (one for each gene) must be saved or not.

    nthreads: positive integer, optional
        Number of threads used for parallel computing
        default: 1
    
    
    Returns
    -------

    A tuple (VIM, alphas, prediction_score, stability_score, treeEstimators).

    VIM: array in which the element (i,j) is the score of the edge directed from the i-th gene to the j-th gene. All diagonal elements are set to zero (auto-regulations are not considered). When a list of candidate regulators is provided, all the edges directed from a gene that is not a candidate regulator are set to zero.
 
    alphas: vector in which the i-th element is the degradation rate of the i-th gene.
    
    prediction_score: prediction score on out-of-bag samples (averaged over all genes and all trees). prediction_score is an empty list if compute_quality_scores is set to False or if tree_method is not set to 'RF'.
    
    stability_score: stability score (averaged over all genes). stability_score is an empty list if compute_quality_scores is set to False.
    
    treeEstimators: list of tree models, where the i-th model is the model predicting the expression of the i-th gene. treeEstimators is an empty list if save_models is set to False.

    '''

    time_start = time.time()

    # Check input arguments
    if not isinstance(TS_data,(list,tuple)):
        raise ValueError('TS_data must be a list of arrays, where each row of an array corresponds to a time point/sample and each column corresponds to a gene')
    
    for expr_data in TS_data:
        if not isinstance(expr_data,np.ndarray):
            raise ValueError('TS_data must be a list of arrays, where each row of an array corresponds to a time point/sample and each column corresponds to a gene')
    
    ngenes = TS_data[0].shape[1]

    if len(TS_data) > 1:
        for expr_data in TS_data[1:]:
            if expr_data.shape[1] != ngenes:
                raise ValueError('The number of columns/genes must be the same in every array of TS_data.')
                
                
    if not isinstance(time_points,(list,tuple)):
        raise ValueError('time_points must be a list of n one-dimensional arrays, where n is the number of time series experiments in TS_data')
    
    if len(time_points) != len(TS_data):
        raise ValueError('time_points must be a list of n one-dimensional arrays, where n is the number of time series experiments in TS_data')
    
    for tp in time_points:
        if (not isinstance(tp,(list,tuple,np.ndarray))) or (isinstance(tp,np.ndarray) and tp.ndim > 1):
            raise ValueError('time_points must be a list of n one-dimensional arrays, where n is the number of time series in TS_data')
        
    for (i,expr_data) in enumerate(TS_data):
        if len(time_points[i]) != expr_data.shape[0]:
            raise ValueError('The length of the i-th vector of time_points must be equal to the number of rows in the i-th array of TS_data')

    if alpha != 'from_data':
        if not isinstance(alpha,(list,tuple,np.ndarray,int,float)):
            raise ValueError("input argument alpha must be either 'from_data', a positive number or a vector of positive numbers")
        
        if isinstance(alpha,(int,float)) and alpha < 0:
            raise ValueError("the degradation rate(s) specified in input argument alpha must be positive")
        
        if isinstance(alpha,(list,tuple,np.ndarray)):
            if isinstance(alpha,np.ndarray) and alpha.ndim > 1:
                raise ValueError("input argument alpha must be either 'from_data', a positive number or a vector of positive numbers")
            if len(alpha) != ngenes:
                raise ValueError('when input argument alpha is a vector, this must be a vector of length p, where p is the number of genes')
            for a in alpha:
                if a < 0:
                    raise ValueError("the degradation rate(s) specified in input argument alpha must be positive")

    if SS_data is not None:
        if not isinstance(SS_data,np.ndarray):
            raise ValueError('SS_data must be an array in which each row corresponds to a steady-state condition/sample and each column corresponds to a gene')
        
        if SS_data.ndim != 2:
            raise ValueError('SS_data must be an array in which each row corresponds to a steady-state condition/sample and each column corresponds to a gene')
    
        if SS_data.shape[1] != ngenes:
            raise ValueError('The number of columns/genes in SS_data must by the same as the number of columns/genes in every array of TS_data.')
        

    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
    
    if regulators != 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')        
    
    if tree_method != 'RF' and tree_method != 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')
    
    if K != 'sqrt' and K != 'all' and not isinstance(K,int): 
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
    
    if isinstance(K,int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')

    if not isinstance(ntrees,int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')
    
    if not isinstance(compute_quality_scores, bool):
        raise ValueError('input argument compute_quality_scores must be a boolean (True or False)')
    
    if not isinstance(save_models, bool):
        raise ValueError('input argument save_models must be a boolean (True or False)')
        
    if not isinstance(nthreads,int):
        raise ValueError('input argument nthreads must be a stricly positive integer')
    elif nthreads <= 0:
        raise ValueError('input argument nthreads must be a stricly positive integer')
    
    

    
    # Re-order time points in increasing order
    for (i,tp) in enumerate(time_points):
        tp = np.array(tp, dtype=np.float32)
        indices = np.argsort(tp)
        time_points[i] = tp[indices]
        expr_data = TS_data[i]
        TS_data[i] = expr_data[indices,:]
        
    if alpha == 'from_data':
        alphas = estimate_degradation_rates(TS_data,time_points)
    elif isinstance(alpha,(int,float)):
        alphas = np.zeros(ngenes) + float(alpha)    
    else:
        alphas = [float(a) for a in alpha]

                
    print('Tree method: ' + str(tree_method))
    print('K: ' + str(K))
    print('Number of trees: ' + str(ntrees))
    print('alpha min: ' + str(min(alphas)))
    print('alpha max: ' + str(max(alphas)))
    print('\n')
                

    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = list(range(ngenes))
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]


    # Learn an ensemble of trees for each target gene and compute scores for candidate regulators
    VIM = np.zeros((ngenes,ngenes))

    if compute_quality_scores:
        if tree_method == 'RF':
            prediction_score = np.zeros(ngenes)
        else:
            prediction_score = []
        stability_score = np.zeros(ngenes)
    else:
        prediction_score = []
        stability_score = []
    
    if save_models:
        treeEstimators = [0] * ngenes
    else:
        treeEstimators = []

    if nthreads > 1:
        print('running jobs on %d threads' % nthreads)

        input_data = [[TS_data,time_points,SS_data,i,alphas[i],input_idx,tree_method,K,ntrees,compute_quality_scores,save_models] for i in range(ngenes)]

        pool = Pool(nthreads)
        alloutput = pool.map(wr_dynGENIE3_single, input_data)

        for out in alloutput:
            i = out[0]
        
            (vi,prediction_score_i,stability_score_i,treeEstimator) = out[1]
            VIM[i,:] = vi
        
            if compute_quality_scores:
                if tree_method == 'RF':
                    prediction_score[i] = prediction_score_i
                stability_score[i] = stability_score_i
        
            if save_models:
                treeEstimators[i] = treeEstimator

    else:
        print('running single threaded jobs')
        for i in range(ngenes):
            print('Gene %d/%d...' % (i+1,ngenes))
        
            (vi,prediction_score_i,stability_score_i,treeEstimator) = dynGENIE3_single(TS_data,time_points,SS_data,i,alphas[i],input_idx,tree_method,K,ntrees,compute_quality_scores,save_models)
            VIM[i,:] = vi
        
            if compute_quality_scores:
                if tree_method == 'RF':
                    prediction_score[i] = prediction_score_i
                stability_score[i] = stability_score_i
            
            if save_models:
                treeEstimators[i] = treeEstimator

    VIM = np.transpose(VIM)
    if compute_quality_scores:
        if tree_method == 'RF':
            prediction_score = np.mean(prediction_score)
        stability_score = np.mean(stability_score)

    time_end = time.time()
    print("Elapsed time: %.2f seconds" % (time_end - time_start))

    return VIM, alphas, prediction_score, stability_score, treeEstimators



def wr_dynGENIE3_single(args):
    return([args[3], dynGENIE3_single(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10])])
    


def dynGENIE3_single(TS_data,time_points,SS_data,output_idx,alpha,input_idx,tree_method,K,ntrees,compute_quality_scores,save_models):

    h = 1 # lag (in number of time points) used for the finite approximation of the derivative of the target gene expression
    ntop = 5 # number of top-ranked candidate regulators over which to compute the stability score

    ngenes = TS_data[0].shape[1]
    nexp = len(TS_data)
    nsamples_time = sum([expr_data.shape[0] for expr_data in TS_data]) 
    ninputs = len(input_idx)

    # Construct learning sample 

    # Time-series data
    input_matrix_time = np.zeros((nsamples_time-h*nexp,ninputs))
    output_vect_time = np.zeros(nsamples_time-h*nexp)
    
    # Data for the computation of the prediction score on out-of-bag samples
    output_vect_time_present = np.zeros(nsamples_time-h*nexp)
    output_vect_time_future = np.zeros(nsamples_time-h*nexp)
    time_diff = np.zeros(nsamples_time-h*nexp)

    nsamples_count = 0

    for (i,current_timeseries) in enumerate(TS_data):
        current_time_points = time_points[i]
        npoints = current_timeseries.shape[0]
        time_diff_current = current_time_points[h:] - current_time_points[:npoints-h]
        current_timeseries_input = current_timeseries[:npoints-h,input_idx]
        current_timeseries_output = (current_timeseries[h:,output_idx] - current_timeseries[:npoints-h,output_idx]) / time_diff_current + alpha*current_timeseries[:npoints-h,output_idx]
        nsamples_current = current_timeseries_input.shape[0]
        input_matrix_time[nsamples_count:nsamples_count+nsamples_current,:] = current_timeseries_input
        output_vect_time[nsamples_count:nsamples_count+nsamples_current] = current_timeseries_output
        output_vect_time_present[nsamples_count:nsamples_count+nsamples_current] = current_timeseries[:npoints-h,output_idx]
        output_vect_time_future[nsamples_count:nsamples_count+nsamples_current] = current_timeseries[h:,output_idx]
        time_diff[nsamples_count:nsamples_count+nsamples_current] = time_diff_current
        nsamples_count += nsamples_current
    
    # Steady-state data (if any)
    if SS_data is not None:

        input_matrix_steady = SS_data[:,input_idx]
        output_vect_steady = SS_data[:,output_idx] * alpha
    
        # Concatenation
        input_all = np.vstack([input_matrix_steady,input_matrix_time])
        output_all = np.concatenate((output_vect_steady,output_vect_time))
        
        del input_matrix_time
        del output_vect_time
        del input_matrix_steady
        del output_vect_steady
    
    else:
  
        input_all = input_matrix_time
        output_all = output_vect_time
        
        del input_matrix_time
        del output_vect_time


    # Parameters of the tree-based method

    # Whether or not to compute the prediction score of out-of-bag samples
    if compute_quality_scores and tree_method =='RF':
        oob_score = True
    else:
        oob_score = False
    
    # Parameter K of the tree-based method
    if (K == 'all') or (isinstance(K,int) and K >= len(input_idx)):
        max_features = "auto"
    else:
        max_features = K

    if tree_method == 'RF':
        treeEstimator = RandomForestRegressor(n_estimators=ntrees,max_features=max_features,oob_score=oob_score)
    elif tree_method == 'ET':
        treeEstimator = ExtraTreesRegressor(n_estimators=ntrees,max_features=max_features,oob_score=oob_score)

    # Learn ensemble of trees
    treeEstimator.fit(input_all,output_all)

    # Compute importance scores
    feature_importances = compute_feature_importances(treeEstimator)
    vi = np.zeros(ngenes)
    vi[input_idx] = feature_importances
    vi[output_idx] = 0

    # Normalize importance scores
    vi_sum = np.sum(vi)
    if vi_sum > 0:
        vi = vi / vi_sum
        
        
    # Ranking quality scores
    prediction_score_oob = []
    stability_score = []
        
    if compute_quality_scores:
        
        if tree_method == 'RF':
            
            # Prediction of out-of-bag samples
        
            if SS_data is not None:
            
                nsamples_SS = SS_data.shape[0]
            
                # Samples coming from the steady-state data
                oob_prediction_SS = treeEstimator.oob_prediction_[:nsamples_SS]
                output_pred_SS = oob_prediction_SS / alpha
            
                # Samples coming from the time series data
                oob_prediction_TS = treeEstimator.oob_prediction_[nsamples_SS:]
                output_pred_TS = (oob_prediction_TS - alpha*output_vect_time_present) * time_diff + output_vect_time_present
            
                output_pred = np.concatenate((output_pred_SS,output_pred_TS))
                output_true = np.concatenate((SS_data[:,output_idx],output_vect_time_future))
            
                (prediction_score_oob,tmp) = np.pearsonr(output_pred,output_true)
            
            else:
                oob_prediction_TS = treeEstimator.oob_prediction_
                output_pred_TS = (oob_prediction_TS - alpha*output_vect_time_present) * time_diff + output_vect_time_present
   
                (prediction_score_oob,tmp) = np.pearsonr(output_pred_TS,output_vect_time_future)
            
            
        # Stability score
   
        # Importances returned by each tree
        importances_by_tree = np.array([e.tree_.compute_feature_importances(normalize=False) for e in treeEstimator.estimators_])
        if output_idx in input_idx:
            idx = input_idx.index(output_idx)
            # Remove importances of target gene
            importances_by_tree = np.delete(importances_by_tree,idx,1)
            
        # Add some jitter to avoir numerical errors
        importances_by_tree = importances_by_tree + np.random.uniform(low=1e-12,high=1e-11,size=importances_by_tree.shape)
            
        if np.sum(importances_by_tree) > 0:
        
            # Ranking of candidate regulators
            ranking_by_tree = [importances_by_tree[i,:].argsort()[::-1] for i in range(ntrees)]
            top_by_tree = [set(r[:ntop]) for r in ranking_by_tree]
    
            # Stability score computed over the top-ranked candidate regulators
            stability_score = np.mean([len(top_by_tree[i].intersection(top_by_tree[j])) for (i,j) in combinations(range(ntrees),2)]) / float(ntop)
            
                
        # Variance of output is too small --> no forest was built and all the importances are zero    
        else:
            stability_score = 0.0
            
    if save_models: 
        return vi, prediction_score_oob, stability_score, treeEstimator
    else:
        return vi, prediction_score_oob, stability_score, []
    
    
    
    
    
    
    
    
def dynGENIE3_predict_doubleKO(expr_WT,treeEstimators,alpha,gene_names,regulators,KO1_gene,KO2_gene,nTimePoints,deltaT):
    
    '''Prediction of gene expressions in a double knockout experiment.

    Parameters
    ----------

    expr_WT: vector containing the gene expressions in the wild-type.
    
    treeEstimators: list of tree models, as returned by the function dynGENIE3(), where the i-th model is the model predicting the expression of the i-th gene. 
        The i-th model must correspond to the i-th gene in expr_WT.
    
    alpha: a positive number or a vector of positive numbers
        Specifies the degradation rate of the different gene expressions. 
        When alpha is a vector of positives, the i-th element of the vector must specify the degradation rate of the i-th gene.
        When alpha is a positive number, all the genes are assumed to have the same degradation rate.
    
    gene_names: list of strings
        List containing the names of the genes. The i-th item of gene_names must correspond to the i-th gene in expr_WT.
    
    regulators: list of strings
        List containing the names of the candidate regulators. When regulators is set to 'all', any gene can be a candidate regulator.
        Note that the candidate regulators must be the same as the ones used when calling the function dynGENIE3().
    
    KO1_gene: name of the first knocked-out gene.
    
    KO2_gene: name of the second knocked-out gene.
    
    nTimePoints: integer
        Specifies the number of time points for which to make a prediction.
    
    deltaT: a positive number
        Specifies the (constant) time interval between two predictions.
    
    
    
    Returns
    -------

    An array in which the element (t,i) is the predicted expression of the i-th gene at the t-th time point.
    The first row of the array contains the initial gene expressions (i.e. the expressions in expr_WT), where the expressions of the two knocked-out genes are set to 0.

    '''
    
    
    time_start = time.time()
    
    # Check input arguments
    if not isinstance(expr_WT,np.ndarray) or expr_WT.ndim > 1:
        raise ValueError("input argument expr_WT must be a vector of numbers")
        
    ngenes = len(expr_WT)
    
    if len(treeEstimators) != ngenes:
        raise ValueError("input argument treeEstimators must contain p tree models, where p is the number of genes in expr_WT")
    
    if not isinstance(alpha,(list,tuple,np.ndarray,int,float)):
        raise ValueError("input argument alpha must be a positive number or a vector of positive numbers")
        
    if isinstance(alpha,(int,float)) and alpha < 0:
        raise ValueError("the degradation rate(s) specified in input argument alpha must be positive")
        
    if isinstance(alpha,(list,tuple,np.ndarray)):
        if isinstance(alpha,np.ndarray) and alpha.ndim > 1:
            raise ValueError("input argument alpha must be a positive number or a vector of positive numbers")
        if len(alpha) != ngenes:
            raise ValueError('when input argument alpha is a vector, this must be a vector of length p, where p is the number of genes')
        for a in alpha:
            if a < 0:
                raise ValueError("the degradation rate(s) specified in input argument alpha must be positive")
     
    if not isinstance(gene_names,(list,tuple)):
        raise ValueError('input argument gene_names must be a list of gene names')
    elif len(gene_names) != ngenes:
        raise ValueError('input argument gene_names must be a list of length p, where p is the number of genes in expr_WT')

    if regulators != 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        sIntersection = set(gene_names).intersection(set(regulators))
        if not sIntersection:
            raise ValueError('The genes must contain at least one candidate regulator')
            
    if not (KO1_gene in gene_names):
        raise ValueError('input argument KO1_gene was not found in gene_names')
        
    if not (KO2_gene in gene_names):
        raise ValueError('input argument KO2_gene was not found in gene_names')

    if not isinstance(nTimePoints,int) or nTimePoints < 1:
        raise ValueError("input argument nTimePoints must be a strictly positive integer")
                   
    if not isinstance(deltaT,(int,float)) or deltaT < 0:
        raise ValueError("input argument deltaT must be a positive number")
        
    KO1_idx = gene_names.index(KO1_gene)
    KO2_idx = gene_names.index(KO2_gene)

    geneidx = list(range(ngenes))
    geneidx.remove(KO1_idx)
    geneidx.remove(KO2_idx)
    
    if isinstance(alpha,(int,float)):
        alphas = np.zeros(ngenes) + float(alpha)    
    else:
        alphas = [float(a) for a in alpha]
    
    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = list(range(ngenes))
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]
                
        
    # Predict time series
    
    print('Predicting time series...')
       
    TS_predict = np.zeros((nTimePoints+1,ngenes))
    TS_predict[0,:] = expr_WT
    TS_predict[0,KO1_idx] = 0
    TS_predict[0,KO2_idx] = 0
    
    for t in range(1,nTimePoints+1):
        new_expr = [(treeEstimators[i].predict(TS_predict[t-1,input_idx].reshape(1,-1)) - alphas[i]*TS_predict[t-1,i]) * deltaT + TS_predict[t-1,i] for i in geneidx] 
        TS_predict[t,geneidx] = np.array(new_expr,dtype=np.float32).flatten()

    time_end = time.time()
    print("Elapsed time: %.2f seconds" % (time_end - time_start))

    return TS_predict
