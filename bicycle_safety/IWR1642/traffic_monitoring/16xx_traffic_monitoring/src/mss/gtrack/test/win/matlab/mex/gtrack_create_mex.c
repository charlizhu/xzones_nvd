#include "memory.h"
#include "mex.h"
#include "gtrack.h"

#define ROWS 2
#define COLUMNS 2
#define ELEMENTS 4

void mx2c_config(GTRACK_moduleConfig *config, const mxArray *mx);   

void mexFunction ( int         nlhs,
                   mxArray     *plhs[],
                   int         nrhs,
                   const       mxArray *prhs[]);

void mexFunction ( int         nlhs,
                   mxArray     *plhs[],
                   int         nrhs,
                   const       mxArray *prhs[])                   
{

    GTRACK_moduleConfig config;
	GTRACK_advancedParameters advParams;

	GTRACK_gatingParams gatingParams;
	GTRACK_stateParams stateParams;
	GTRACK_varParams variationParams;
	GTRACK_allocationParams allocationParams;
	GTRACK_sceneryParams sceneryParams;

    void        *handle;
    int32_t     errCode;
    UINT64_T    *dynamicData;                 /* pointer to dynamic data */

        
    /* check proper input and output */
    if(nrhs!=1)
        mexErrMsgIdAndTxt( "MATLAB:phonebook:invalidNumInputs",
                "One input required.");
    else if(nlhs > 1)
        mexErrMsgIdAndTxt( "MATLAB:phonebook:maxlhs",
                "Too many output arguments.");
    else if(!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:phonebook:inputNotStruct",
                "Input must be a structure.");

	memset(&config, 0, sizeof(GTRACK_moduleConfig));
	memset(&advParams, 0, sizeof(GTRACK_advancedParameters));
	config.advParams = &advParams;

	memset(&gatingParams, 0, sizeof(GTRACK_gatingParams));
	memset(&stateParams, 0, sizeof(GTRACK_stateParams));
	memset(&variationParams, 0, sizeof(GTRACK_varParams));
	memset(&allocationParams, 0, sizeof(GTRACK_allocationParams));
	memset(&sceneryParams, 0, sizeof(GTRACK_sceneryParams));

	config.advParams->gatingParams = &gatingParams;
	config.advParams->stateParams = &stateParams;
	config.advParams->variationParams = &variationParams;
	config.advParams->allocationParams = &allocationParams;
	config.advParams->sceneryParams = &sceneryParams;

	mx2c_config(&config, prhs[0]);
    
    handle = gtrack_create(&config, &errCode);
    

    /* Create a local array and load data */
    dynamicData = (UINT64_T *)mxCalloc(1, sizeof(UINT64_T));
	*dynamicData = handle;

    /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
    plhs[0] = mxCreateNumericMatrix(0, 0, mxUINT64_CLASS, mxREAL);

    /* Point mxArray to dynamicData */
    mxSetData(plhs[0], dynamicData);
    mxSetM(plhs[0], 1);
    mxSetN(plhs[0], 1);
#if 0
    /* Create a 0-by-0 mxArray; the memory is allocated in C world, here we need to pass a pointer */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    /* Allocate local variable to hold a void pointer */
    mxHandle = mxCalloc(1, sizeof(UINT64_T));
	mxHandle = handle;
    /* Point mxArray to dynamicData */
    mxSetData(plhs[0], mxHandle);
	/* Do not call mxFree(mxHandle) because plhs[0] points to it */
#endif
}

/* Parse configuration from matlab to C */
void mx2c_config(GTRACK_moduleConfig *config, const mxArray *mx) 
{ 
    mxArray *field;
	mxArray *params;
	char *pdata;
	double *srcData;
	double value;

	size_t m,n;

    if( !mxIsStruct(mx) || mxGetNumberOfElements(mx) == 0 ) { 
        mexErrMsgTxt("Input must be a non-empty struct"); 
    } 
    
    params = mxGetField(mx, 0, "advParams"); 
    if(params == NULL ) { 
		config->advParams = 0;
    } 
	else if(!mxIsStruct(params) || mxGetNumberOfElements(params) == 0 ) { 
        mexErrMsgTxt("advParams must be a non-empty struct"); 
    } 
	else {
		
		field = mxGetField(params, 0, "gating"); 
		if(field != NULL ) {
			pdata =	mxGetData(field);
			if (mxGetM(field) != 1 && mxGetN(field) != 4)
				mexErrMsgTxt( "Gating Parameters shall have 4 values");

			srcData = (double *)pdata;
			config->advParams->gatingParams->volume = (float) *srcData++;
			config->advParams->gatingParams->limits.length = (float) *srcData++;
			config->advParams->gatingParams->limits.width = (float) *srcData++;
			config->advParams->gatingParams->limits.vel = (float) *srcData;
		}
		else {
			/* Default gating params shall be used */
			config->advParams->gatingParams = 0;
		}

		field = mxGetField(params, 0, "allocation"); 
		if(field != NULL ) { 
			pdata =	mxGetData(field);
			if (mxGetM(field) != 1 && mxGetN(field) != 6)
				mexErrMsgTxt( "Allocation Parameters shall have 5 values");

			srcData = (double *)pdata;				
			config->advParams->allocationParams->snrThre = (float) *srcData++;
			config->advParams->allocationParams->snrThreObscured = (float) *srcData++;
			config->advParams->allocationParams->velocityThre = (float) *srcData++;
			config->advParams->allocationParams->pointsThre = (UINT16_T) *srcData++;
			config->advParams->allocationParams->maxDistanceThre = (float) *srcData++;
			config->advParams->allocationParams->maxVelThre = (float) *srcData;
		} 
		else {
			/* Default allocation params shall be used */
			config->advParams->allocationParams = 0;
		}

		field = mxGetField(params, 0, "thresholds"); 
		if(field != NULL ) { 
			pdata =	mxGetData(field);
			if (mxGetM(field) != 1 && mxGetN(field) != 5)
				mexErrMsgTxt( "State Parameters shall have 5 values");
			srcData = (double *)pdata;
			config->advParams->stateParams->det2actThre = (UINT16_T) *srcData++;
			config->advParams->stateParams->det2freeThre = (UINT16_T) *srcData++;
			config->advParams->stateParams->active2freeThre = (UINT16_T) *srcData++;
			config->advParams->stateParams->static2freeThre = (UINT16_T) *srcData++;
			config->advParams->stateParams->exit2freeThre = (UINT16_T) *srcData;
		}
		else {
			/* Default thresholds params shall be used */
			config->advParams->stateParams = 0;
		}


		field = mxGetField(params, 0, "variations"); 
		if(field != NULL ) { 
			pdata =	mxGetData(field);
			if (mxGetM(field) != 1 && mxGetN(field) != 3)
				mexErrMsgTxt( "Variation Parameters shall have 3 values");

			srcData = (double *)pdata;				
			config->advParams->variationParams->lengthStd = (float) *srcData++;
			config->advParams->variationParams->widthStd = (float) *srcData++;
			config->advParams->variationParams->dopplerStd = (float) *srcData;
		}
		else {
			/* Default measurements params shall be used */
			config->advParams->variationParams = 0;
		}

		field = mxGetField(params, 0, "scenery"); 
		if(field != NULL ) { 
			pdata =	mxGetData(field);
			if (mxGetM(field) != 1)
				mexErrMsgTxt( "Scenery Parameters shall have single dimensional array");
			srcData = (double *)pdata;

            config->advParams->sceneryParams->numBoundaryBoxes = (uint8_t) *srcData++;
            if(config->advParams->sceneryParams->numBoundaryBoxes > GTRACK_MAX_BOUNDARY_BOXES)
                mexErrMsgTxt( "Scenery Parameters has too many boundary boxes");

            if (mxGetN(field) <= 1+4*config->advParams->sceneryParams->numBoundaryBoxes+1)
				mexErrMsgTxt( "Scenery Parameters not enough parameters for boundary box");

            for(n=0; n<config->advParams->sceneryParams->numBoundaryBoxes; n++)  {
    			config->advParams->sceneryParams->boundaryBox[n].left = (float) *srcData++;
    			config->advParams->sceneryParams->boundaryBox[n].right = (float) *srcData++;
    			config->advParams->sceneryParams->boundaryBox[n].bottom = (float) *srcData++;
    			config->advParams->sceneryParams->boundaryBox[n].top = (float) *srcData++;
    		}
            config->advParams->sceneryParams->numStaticBoxes = (uint8_t) *srcData++;

            if(config->advParams->sceneryParams->numStaticBoxes > GTRACK_MAX_STATIC_BOXES)
                mexErrMsgTxt( "Error: Scenery Parameters has too many static boxes");

            if (mxGetN(field) != 1+4*config->advParams->sceneryParams->numBoundaryBoxes+1+4*config->advParams->sceneryParams->numStaticBoxes)
				mexErrMsgTxt( "Scenery Parameters wrong unumber of inputs");

            for(n=0; n<config->advParams->sceneryParams->numStaticBoxes; n++)  {
    			config->advParams->sceneryParams->staticBox[n].left = (float) *srcData++;
    			config->advParams->sceneryParams->staticBox[n].right = (float) *srcData++;
    			config->advParams->sceneryParams->staticBox[n].bottom = (float) *srcData++;
    			config->advParams->sceneryParams->staticBox[n].top = (float) *srcData++;
    		}
        }
		else {
			/* Default measurements params shall be used */
			config->advParams->sceneryParams = 0;
		}
	}

    field = mxGetField(mx, 0, "stateVectorType"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'stateVectorType'"); 
    } 
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->stateVectorType = (GTRACK_STATE_VECTOR_TYPE)value;

    field = mxGetField( mx, 0, "verbose"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'verbose'");  
    } 
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->verbose = (GTRACK_VERBOSE_TYPE)value;    
     
    field = mxGetField( mx, 0, "maxNumPoints"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'maxNumPoints'"); 
    } 
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->maxNumPoints = (uint16_t)value;    
    
    field = mxGetField( mx, 0, "maxNumTracks"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'maxNumTracks'"); 
    } 
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->maxNumTracks = (uint16_t)value;  
    
    field = mxGetField( mx, 0, "maxRadialVelocity"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'maxRadialVelocity'"); 
    } 
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->maxRadialVelocity = (float)value;
    
    field = mxGetField( mx, 0, "radialVelocityResolution"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'radialVelocityResolution'"); 
    } 
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->radialVelocityResolution = (float)value;

    field = mxGetField( mx, 0, "maxAcceleration"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'maxAcceleration'"); 
    }
    if (mxGetM(field) != 1 && mxGetN(field) != 2) {
        mexErrMsgTxt( "max Acceleration shall have 2 values");
    }

    pdata =	mxGetData(field);
    srcData = (double *)pdata;				
    config->maxAccelerationX = (float) *srcData++;
    config->maxAccelerationY = (float) *srcData++;

        
    field = mxGetField( mx, 0, "initialRadialVelocity"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'initialRadialVelocity'"); 
    }
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->initialRadialVelocity = (float)value;

	
	field = mxGetField( mx, 0, "deltaT"); 
    if( field == NULL ) { 
        mexErrMsgTxt("Input struct must have field 'deltaT'"); 
    }
	pdata =	mxGetData(field);
	value = *(double *)pdata;
    config->deltaT = (float)value;
}
