#include <mex.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	mwIndex  D = mxGetM(prhs[0]);
	mwIndex M = mxGetN(prhs[0]);
	double temp, accum, sigma2, sigma3, summ, *sr;
    mwIndex oldnzmax, *irs, *jcs;
	  double *X = mxGetPr(prhs[0]);
    double sigma = *mxGetPr(prhs[1]);
    mwIndex nzmax = M;
    plhs[0] = mxCreateSparse(M, M, nzmax, mxREAL);
    sr  = mxGetPr(plhs[0]);
    irs = mxGetIr(plhs[0]);
    jcs = mxGetJc(plhs[0]);
    
    sigma2 = -1 / (2*sigma*sigma);
    sigma3 = 9*sigma*sigma;
    summ = 0;
    int k = 0;
	for (int m=0; m<M; m++)
	{
        jcs[m] = k;
		for (int n=0; n<=m; n++)
		{
			accum = 0.0;
			for (int d=0; d<D; d++)
			{
				temp = X[d + m * D] - X[d + n * D];
				accum += temp * temp;
			}
            if (accum <= sigma3)
            {   
                if (k >= nzmax){
                    mwIndex oldnzmax = nzmax;
                    nzmax = (int) ceil((double) nzmax * 1.25);
                    if (oldnzmax == nzmax)
                    {
                        nzmax = (int) ceil((double) nzmax * 1.1);
                    }
                    mxSetNzmax(plhs[0], nzmax); 
                    mxSetPr(plhs[0], (double*) mxRealloc(sr, nzmax*sizeof(double)));
                    mxSetIr(plhs[0], (mwIndex*) mxRealloc(irs, nzmax*sizeof(mwIndex)));
                    sr  = mxGetPr(plhs[0]);
                    irs = mxGetIr(plhs[0]);
                }
                if (n == m) {
                    temp = 0.5;
                }
                else
                {
                    temp = exp(accum * sigma2);
                }
                sr[k] = temp;
                irs[k] = n;
                k++;
                summ += temp;
            }
		}   
	}
    jcs[M] = k;
    
    summ *= 2.0 / ((double) M);
    
    for (int i=0; i<k; i++)
    {
        sr[i] /= summ;
    }
}
