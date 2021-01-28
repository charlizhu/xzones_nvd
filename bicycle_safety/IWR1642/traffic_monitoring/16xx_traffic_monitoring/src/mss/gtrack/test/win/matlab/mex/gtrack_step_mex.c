#include "mex.h"
#include "gtrack.h"

void mexFunction ( int         nlhs,
                   mxArray     *plhs[],
                   int         nrhs,
                   const       mxArray *prhs[]);

void mexFunction ( int         nlhs,
                   mxArray     *plhs[],
                   int         nrhs,
                   const       mxArray *prhs[])                   
{
    GTRACK_measurementPoint *point;
    GTRACK_measurementVariance *var;
    uint32_t *vid;
    GTRACK_targetDesc tDescr[100];
	uint8_t mIndex[1000];

    double *pDouble;
    float *pSingle;

	uint64_t *pU64;
    uint16_t *pU16T;
	uint8_t *pU8T;

    uint16_t tNum;
    void *handle;
    uint16_t mNum;
    
    float *dest;
    float *src;
    
    size_t mxM, mxN;
    uint32_t n, k;
    uint32_t benchmark[32];
            
    
    /* check proper input and output */
    if(nrhs!=4)
        mexErrMsgIdAndTxt( "MATLAB:phonebook:invalidNumInputs", "Four inputs required");
    else if(nlhs != 5)
        mexErrMsgIdAndTxt( "MATLAB:phonebook:maxlhs", "Five outputs required");
	
	/* First Argument */
    if (mxGetM(prhs[0]) == 0 || mxGetN(prhs[0]) == 0)
        mexErrMsgIdAndTxt( "MATLAB:phonebook:invalidNumInputs", "First Input can not be empty");
    pU64 = (uint64_t *)mxGetPr(prhs[0]);
    handle = (void *)*pU64;

	/* Second Argument */
    mxM = mxGetM(prhs[1]);
    mxN = mxGetN(prhs[1]);

    if ((mxM == 0) || (mxN == 0))
		point = NULL;
	else {
	    pSingle = mxGetPr(prhs[1]);
		point = (GTRACK_measurementPoint *)pSingle;
	}

	/* Third Argument */
    if (mxGetM(prhs[2]) == 0 || mxGetN(prhs[2]) == 0)
		var = NULL;
	else {	
		pSingle = mxGetPr(prhs[2]);
		var = (GTRACK_measurementVariance *)pSingle;
	}

	/* Forth Argument */
    if (mxGetM(prhs[3]) == 0 || mxGetN(prhs[3]) == 0)
        mexErrMsgIdAndTxt( "MATLAB:phonebook:invalidNumInputs", "Forth Input can not be empty");
    pSingle = mxGetPr(prhs[3]);
    mNum = (uint16_t)*pSingle;
  
    gtrack_step(handle, point, var, mNum, tDescr, &tNum, &mIndex, benchmark);

	if(tNum == 0) {    
        plhs[0] = mxCreateNumericMatrix(0, 0, mxUINT16_CLASS, mxREAL);
        plhs[1] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);
        plhs[3] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);
    }
    else {
        /* First is an array of UIDs */  
        plhs[0] = mxCreateNumericMatrix(1, tNum, mxUINT16_CLASS, mxREAL);
        pU16T = (uint16_t *)mxGetPr(plhs[0]);
    
        /* Copy data into the mxArray uid[tnum]*/
        for (n = 0; n < tNum; n++ ) {
            *pU16T++ = tDescr[n].uid;
        }
        mxSetM(plhs[0], 1);
        mxSetN(plhs[0], tNum);

        /* Second is an array of State inforamtion: S[6,tnum]*/
        plhs[1] = mxCreateNumericMatrix(6, tNum, mxSINGLE_CLASS, mxREAL);
        dest = (float *)mxGetPr(plhs[1]);

        /* Copy data into the mxArray */
        for (n = 0; n < tNum; n++ ) {
            src = tDescr[n].S;
            for (k = 0; k < 6; k++ ) {
                *dest++ = *src++;
            }
        }
        mxSetM(plhs[1], 6);
        mxSetN(plhs[1], tNum);    
    
        /* Third is an array of iC_inv[3,3] = Inverse of Group Innovation Covariance Matrices */
        plhs[2] = mxCreateNumericMatrix(9, tNum, mxSINGLE_CLASS, mxREAL);
        dest = (float *)mxGetPr(plhs[2]);
        /* Copy data into the mxArray */
        for (n = 0; n < tNum; n++ ) {
            src = tDescr[n].EC;
		    for (k = 0; k < 9; k++ ) {
                *dest++ = *src++;
            }
        }
        mxSetM(plhs[2], 9);
        mxSetN(plhs[2], tNum);  

        /* Forth is an array of Gain values used in association process */
        plhs[3] = mxCreateNumericMatrix(1, tNum, mxSINGLE_CLASS, mxREAL);
        dest = (float *)mxGetPr(plhs[3]);
        /* Copy data into the mxArray */
        for (n = 0; n < tNum; n++ ) {
		    *dest++ = tDescr[n].G;
        }
        mxSetM(plhs[3], 1);
        mxSetN(plhs[3], tNum);
    }
    

    /* Fifth is an array of association indices */
    if(mNum == 0) {
        plhs[4] = mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL);
    }
    else {
        plhs[4] = mxCreateNumericMatrix(1, mNum, mxUINT8_CLASS, mxREAL);
        pU8T = (uint8_t *)mxGetPr(plhs[4]);

        /* Copy data into the mxArray uid[tnum]*/
        for (n = 0; n < mNum; n++ ) {
            *pU8T++ = mIndex[n];
        }
        mxSetM(plhs[4], 1);
        mxSetN(plhs[4], mNum);
    }
}
