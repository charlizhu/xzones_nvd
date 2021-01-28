#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gtrack.h"

/* This test application (traffic monitoring), wants to modify default parameters */
GTRACK_sceneryParams appSceneryParams = {
	0,{{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f}},0,{{0.f,0.f,0.f,0.f},{0.f,0.f,0.f,0.f}} 	/* no boundary boxes, static boxes */
};
GTRACK_gatingParams appGatingParams = {
     16.f, {12.f, 8.f, 0.f}    /* TM: 16 gating volume, Limits are set to 8m in length, 2m in width, 0 no limit in doppler */
};
GTRACK_stateParams appStateParams = {
     3U, 3U, 5U, 5U, 5U              /* TM: 3 frames to activate, 3 to forget, 5 to delete */
};
GTRACK_allocationParams appAllocationParams = {
     -1.f, -1.f, 1.f, 3U, 4.f, 2.f           /* TM: any SNRs, 1m/s minimal velocity, 3 points, 4m in distance, 2m/c in velocity */
};
/* Using standard deviation of uniformly distributed variable in the range [a b]: 1/sqrt(12)*(b-a) */
GTRACK_varParams appVariationParams = {
     /* Standard deviation of uniformly distributed number in range [a b]: sqrt(1/12)*(b-a) */
     1.f/3.46f, 1.f/3.46f, 1.f     /* TM: 1m height, 1m in width, 2 m/s for doppler */
};


extern unsigned int memoryBytesUsed;

int main(int argc, char **argv)
{
	GTRACK_moduleConfig config;    
	GTRACK_advancedParameters advParams;

	GTRACK_measurementPoint points[1000];
	GTRACK_measurementVariance variances[1000];
	GTRACK_targetDesc targetDescr[100];
	uint32_t vid[1000];


	uint32_t gtrackStartTime, benchmarkCycles;
	uint32_t benchmarks[GTRACK_BENCHMARK_SIZE];

	int32_t errCode;


	void *hTrackModule;
	unsigned int tick, gtick;	
	uint16_t loop;	
	uint16_t mNum;	
	uint16_t tNum;	
	uint16_t vNum = 0;	
	

	unsigned int n, k;
	FILE *fCloud;
	char fileName[120];

	memset((void *)&config, 0, sizeof(GTRACK_moduleConfig));
	config.stateVectorType = GTRACK_STATE_VECTORS_2DA; // Track two dimensions with acceleration 
	config.verbose = GTRACK_VERBOSE_NONE;
	config.deltaT = 0.04f; // 40ms ticks
	config.maxRadialVelocity = 12.5f; // Radial velocity from sensor is limited to +/- maxURV (in m/s)
	config.radialVelocityResolution = 0.4f; // Radial velocity resolution (in m/s)
	config.maxAccelerationX = 20.f; // Target acceleration is not excedding 12m/s2 in lateral direction
	config.maxAccelerationY = 20.f; // Target acceleration is not excedding 12m/s2 in longitudinal direction
	config.maxNumPoints = 1000;
	config.maxNumTracks = 40;
	config.initialRadialVelocity = -20; // Expected target radial velocity at the moment of detection, m/s
	
    /* Here, we want to set allocation, gating, and threshold parameters, leaving the rest to default */
	memset((void *)&advParams, 0, sizeof(GTRACK_advancedParameters));
	advParams.allocationParams = &appAllocationParams;
	advParams.gatingParams = &appGatingParams;
	advParams.stateParams = &appStateParams;
    advParams.sceneryParams = &appSceneryParams;
    advParams.variationParams = &appVariationParams;

	config.advParams = &advParams;

	hTrackModule = gtrack_create(&config, &errCode);

	printf("%d Bytes of data memory used\n", memoryBytesUsed);

	for(loop=0; loop<10; loop++) {

		sprintf(fileName, "../../../test/vectors/usecases/urv20/pointCloud1500_%03u.dat", loop+1);
		fCloud = fopen(fileName, "r");

		for(tick=0; tick<1500; tick++) {
			
			gtick = loop*1500+tick+1;

			fscanf(fCloud, "NUM=%hu\n",&mNum); 
			

			// Limit the points
			if(mNum > config.maxNumPoints)
				mNum = config.maxNumPoints;

			for(n=0; n<mNum; n++) {
				fscanf(fCloud, "%f,%f,%f,%f,%f,%f,%u\n", 
					&points[n].range, &points[n].angle, &points[n].doppler, 
					&variances[n].rangeVar, &variances[n].angleVar, &variances[n].dopplerVar,
					&vid[n]); 
				points[n].snr = 0;
			}
#if 0
				// Radial velocity alliasing
				rv = pointCloud[n].doppler; 
				rv = fmod((pointCloud[n].doppler + config.maxURadialVelocity),(2*config.maxURadialVelocity));
				if(rv < 0)
					rv+= 2*config.maxURadialVelocity;
				rv -= config.maxURadialVelocity;
				pointCloud[n].doppler = rv;
#endif

		    gtrackStartTime = gtrack_getCycleCount();
			gtrack_step(hTrackModule, points, variances, mNum, targetDescr, &tNum, 0, benchmarks);
	        benchmarkCycles = gtrack_getCycleCount() - gtrackStartTime;      
			printf("Frame #%u, %u Target(s), %u Measurements, %u Tracks, %u Ticks =(", gtick, vNum, mNum, tNum, benchmarkCycles);
            for(k=0; k<GTRACK_BENCHMARK_SIZE; k++) {
                printf("%u, ", benchmarks[k] - gtrackStartTime);
                gtrackStartTime = benchmarks[k];
				
			}
            printf("%u)\n", benchmarks[k] - gtrackStartTime);
		}
		fclose(fCloud);
	}
	gtrack_delete(hTrackModule);
	printf("%d Bytes of data memory used\n", memoryBytesUsed);
	
}
