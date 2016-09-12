#include "mex.h"
#include "PriorityVector.hpp"


#define IND2D(i,j,ni,nj)    ((i) + (j)*(ni))
#define ROW(k,ni,nj)        ((k) % (ni))
#define COL(k,ni,nj)        ((k) / (ni))

#define ABS(x) (((x)>=0) ? (x) : -(x))


/* Prototype */

mxArray* GMP(
    const mxArray* P,       /* dxK, patches             */
    const mxArray* D,       /* dxN, dictionary          */
    const mxArray* L0,      /*  K,  column sparsity     */
    const mxArray* L0t,     /*  N , row sparsity        */
    size_t K,               /* number of patches        */
    size_t N,               /* number of atoms          */
    size_t niter            /* number of iterations     */
);


/* Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Check argument counts and types */
    if(nrhs != 5)
        mexErrMsgTxt("Error: mexGMP requires 5 arguments: P, D, L0, L0t, niter");
    if(nlhs != 1)
        mexErrMsgTxt("Error: mexGMP returns 1 value");
    for(int i=0; i<nrhs; ++i)
        if(mxIsSparse(prhs[i]) || mxIsComplex(prhs[i]) || !mxIsDouble(prhs[i]))
            mexErrMsgTxt("Error: mexGMP arguments must be double real non-sparse");

    /* Names and dimensions*/
    const mxArray* P  = prhs[0];
    const mxArray* D  = prhs[1];
    const mxArray* L0 = prhs[2];
    const mxArray* L0t  = prhs[3];
    const mxArray* iterMat = prhs[4];
    const size_t K = mxGetN(P);
    const size_t N = mxGetN(D);

    /* Check arguments sizes */
    if(mxGetNumberOfElements(iterMat) != 1)
        mexErrMsgTxt("Error in mexGMP: argument 'niter' must be scalar");
    if(mxGetNumberOfElements(L0) != K)
        mexErrMsgTxt("Error in mexGMP: invalid dimension for argument L0");
    if(mxGetNumberOfElements(L0t) != N)
        mexErrMsgTxt("Error in mexGMP: invalid dimension for argument L0t");
    if(mxGetM(P) != mxGetM(D))
        mexErrMsgTxt("Error in mexGMP: incompatible dimension for arguments P and D");

    /* Call the processing function */
    const size_t niter = mxGetScalar(prhs[4]);
    plhs[0] = GMP(P, D, L0, L0t, K, N, niter);
}


/* Compute the inner product A'*B between the columns of A and B */

mxArray* columns_inner_product(const mxArray* A, const mxArray* B) {
    mxArray* A_ = const_cast<mxArray*>(A);
    mxArray* B_ = const_cast<mxArray*>(B);

    /* Transpose */
    mxArray* At = NULL;
    mexCallMATLAB(1, &At, 1, &A_, "transpose");

    /* Multiply */
    mxArray* prod = NULL;
    mxArray* tab[] = {At, B_};
    mexCallMATLAB(1, &prod, 2, tab, "mtimes");

    /* Return */
    return prod;
}


/* Perform the matching pursuit */

mxArray* GMP(
    const mxArray* P, const mxArray* D,
    const mxArray* L0_, const mxArray* L0t_,
    size_t K, size_t N, size_t niter)
{
    /* Duplicate / create matrices */
    mxArray* W = mxCreateDoubleMatrix(N, K, mxREAL);
    mxArray* DP = columns_inner_product(D, P);
    mxArray* DD = columns_inner_product(D, D);
    mxArray* L0 = mxDuplicateArray(L0_);
    mxArray* L0t = mxDuplicateArray(L0t_);

    /* Matrix to vectors */
    double* W_vec   = mxGetPr(W);
    double* L0_vec  = mxGetPr(L0);
    double* L0t_vec = mxGetPr(L0t);
    double* DP_vec  = mxGetPr(DP);
    const double* DD_vec  = mxGetPr(DD);

    /* Build the heap */
    PriorityVector abs_DP(K*N);
    for (size_t n = 0; n < N; n++)
        if(L0t_vec[n] > 0)
            for (size_t k = 0; k < K; k++)
                if (L0_vec[k] > 0) {
                    const size_t q = IND2D(n, k, N, K);
                    abs_DP.setPriority(q, ABS(DP_vec[q]));
                }

    /* Iterate */
    for(size_t step = 0; step < niter && abs_DP.count(); step++) {

        /* Max. inner product index */
        const size_t q_min = abs_DP.getFirstIndex();
        const size_t n_min = ROW(q_min, N,K);
        const size_t k_min = COL(q_min, N,K);
        const double dW = DP_vec[q_min];

        /* Update L0 and L0t */
        if(W_vec[q_min] == 0) {
            L0_vec[k_min]--;
            if(L0_vec[k_min] <= 0)
                for(size_t n = 0; n < N; n++)
                    if(abs_DP.hasPriority(IND2D(n, k_min, N, K)))
                        abs_DP.removePriority(IND2D(n, k_min, N, K));
            L0t_vec[n_min]--;
            if(L0t_vec[n_min] <= 0)
                for(size_t k = 0; k < K; k++)
                    if(abs_DP.hasPriority(IND2D(n_min, k, N, K)))
                        abs_DP.removePriority(IND2D(n_min, k, N, K));
        }

        /* Update the element and inner products*/
        W_vec[q_min] += dW;
        for(size_t n = 0; n < N; n++) {
            const size_t q = IND2D(n, k_min, N, K);
            DP_vec[q] -= dW * DD_vec[IND2D(n, n_min, N, N)];
            if(abs_DP.hasPriority(q))
                abs_DP.setPriority(q, ABS(DP_vec[q]));
        }
    }

    /* Finished */
    return W;
}

