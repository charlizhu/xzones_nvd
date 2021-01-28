/*
 *   @file  gtrackApp.c
 *
 *   @brief
 *      Gtrack Application Unit Test code
 *
 *  \par
 *  NOTE:
 *      (C) Copyright 2016 Texas Instruments, Inc.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *    Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *    Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 *    Neither the name of Texas Instruments Incorporated nor the names of
 *    its contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**************************************************************************
 *************************** Include Files ********************************
 **************************************************************************/

/* Standard Include Files. */
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>

/* BIOS/XDC Include Files. */
#include <xdc/std.h>
#include <xdc/cfg/global.h>
#include <xdc/runtime/IHeap.h>
#include <xdc/runtime/System.h>
#include <xdc/runtime/Error.h>
#include <xdc/runtime/Memory.h>
#include <ti/sysbios/BIOS.h>
#include <ti/sysbios/knl/Task.h>
#include <ti/sysbios/knl/Event.h>
#include <ti/sysbios/knl/Semaphore.h>
#include <ti/sysbios/knl/Clock.h>
#include <ti/sysbios/heaps/HeapBuf.h>
#include <ti/sysbios/heaps/HeapMem.h>

/* mmWave SK Include Files: */
#include <ti/common/sys_common.h>
#include <ti/drivers/soc/soc.h>
#include <ti/drivers/esm/esm.h>
#include <ti/utils/cycleprofiler/cycle_profiler.h>

#include <ti/alg/gtrack/gtrack.h>
#include "math.h"
#include "float.h"

GTRACK_measurementPoint points[1000];
GTRACK_measurementVariance variances[1000];
uint32_t vid[1000];

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
     -1.f, -1.f, 1.f, 3U, 4.f, 2.f           /* TM: any SNR, any SNR even if obscured, 1m/s minimal velocity, 3 points, 4m in distance, 2m/c in velocity */
};
/* Using standard deviation of uniformly distributed variable in the range [a b]: 1/sqrt(12)*(b-a) */
GTRACK_varParams appVariationParams = {
     /* Standard deviation of uniformly distributed number in range [a b]: sqrt(1/12)*(b-a) */
     1.f/3.46, 1.f/3.46, 1.f     /* TM: 1m height, 1m in width, 2 m/s for doppler */
};

#if defined(__TI_COMPILER_VERSION__)
#pragma DATA_SECTION(points, ".data:pointCloudArray")
#pragma DATA_SECTION(variances, ".data:pointCloudArray")
#endif

extern uint32_t memoryBytesUsed;

/**
 *  @b Description
 *  @n
 *      System Initialization Task which initializes the various
 *      components in the system.
 *
 *  @retval
 *      Not Applicable.
 */
void Test_initTask(UArg arg0, UArg arg1)
{
	GTRACK_moduleConfig config;
	GTRACK_targetDesc targetDescr[40];

	uint32_t vunique[20];
    Memory_Stats            startMemoryStats;
    Memory_Stats            peakMemoryStats;
    Memory_Stats            endMemoryStats;

    GTRACK_advancedParameters advParams;

	void *hTrackModule;
	uint32_t tick, gtick;	
	uint16_t loop;	
	uint16_t mNum;	
	uint16_t tNum;	
	uint16_t vNum = 0;	

	uint32_t n, k;
	uint16_t vidFound;

	FILE *fCloud;
	char fileName[120];
    
	uint32_t gtrackStartTime, benchmarkCycles;
	uint32_t benchmarks[GTRACK_BENCHMARK_SIZE];

	int32_t errCode;
	
	
	/* Get the heap statistics at the beginning of the tests */
    HeapMem_getStats (heap0, &startMemoryStats);
    System_printf ("Total Heap Size = %d bytes\n", startMemoryStats.totalSize);
    System_printf ("HeapMem: Used = %d bytes, Free = %d bytes\n", startMemoryStats.totalSize - startMemoryStats.totalFreeSize, startMemoryStats.totalFreeSize);
    
	memset((void *)&config, 0, sizeof(GTRACK_moduleConfig));

	config.stateVectorType = GTRACK_STATE_VECTORS_2DA; // Track two dimensions with acceleration 
	config.verbose = GTRACK_VERBOSE_NONE;
	config.deltaT = 0.04f; // 40ms ticks
	config.maxRadialVelocity = 12.5f; // Radial velocity from sensor is limited to +/- maxURV (in m/s)
	config.radialVelocityResolution = 0.4f; // Radial velocity resolution (in m/s)
	config.maxAccelerationX = 5.f; // Target acceleration is not excedding 12m/s2 in lateral direction
	config.maxAccelerationY = 5.f; // Target acceleration is not excedding 12m/s2 in longitudinal direction	config.maxNumPoints = 150;
	config.maxNumTracks = 20;
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
	
    /* Get the heap statistics at the beginning of the tests */
    HeapMem_getStats (heap0, &peakMemoryStats);
    System_printf ("HeapMem: Used = %d bytes, Free = %d bytes\n", peakMemoryStats.totalSize - peakMemoryStats.totalFreeSize, peakMemoryStats.totalFreeSize);
	System_printf("Gtrack is using %d bytes of data memory\n", memoryBytesUsed);
	
	Cycleprofiler_init();
	gtrackStartTime = gtrack_getCycleCount();
	for(loop=0; loop<1000; loop++) {
        /* one milisecond */
        SOC_microDelay(1000);
	}
    benchmarkCycles = gtrack_getCycleCount() - gtrackStartTime;
    System_printf("Calibration: 1 seconds is %u Ticks\n", benchmarkCycles);

    gtrackStartTime = gtrack_getCycleCount();
    benchmarkCycles = gtrack_getCycleCount() - gtrackStartTime;
    System_printf("Calibration: CCNT read is %u Ticks\n", benchmarkCycles);

	for(loop=0; loop<10; loop++) {

		sprintf(fileName, "../../../test/vectors/usecases/urv20/pointCloud1500_%03u.dat", loop+1);

        fCloud = fopen(fileName, "r");
		if (fCloud==NULL)
		{
			System_printf("Can not open file %s\n",fileName);
			
		}

		for(tick=0; tick<1500; tick++) {
			
			gtick = loop*1500+tick+1;

			fscanf(fCloud, "NUM=%hu\n",&mNum); 

			if(mNum > config.maxNumPoints)
				mNum = config.maxNumPoints;
			    
			for(n=0; n<mNum; n++) {
				fscanf(fCloud, "%f,%f,%f,%f,%f,%f,%u\n", 
					&points[n].range, &points[n].angle, &points[n].doppler, 
					&variances[n].rangeVar, &variances[n].angleVar, &variances[n].dopplerVar,
					&vid[n]);
				points[n].snr = 0;
			}
			
			for(n=0, vNum = 0; n<mNum; n++) {
				vidFound = 0;
				for(k=0; k<vNum; k++) {
					if(vunique[k] == vid[n]) {
						vidFound = 1;
						break;
					}
				}
				if(vidFound == 0)
					vunique[vNum++] = vid[n];
			}

		    gtrackStartTime = gtrack_getCycleCount();
			gtrack_step(hTrackModule, points, variances, mNum, targetDescr, &tNum, 0, benchmarks);
		    benchmarkCycles = gtrack_getCycleCount() - gtrackStartTime;
			System_printf("Frame #%u, %u Target(s), %u Measurements, %u Tracks, %u Ticks =(", gtick, vNum, mNum, tNum, benchmarkCycles);
            for(k=0; k<GTRACK_BENCHMARK_SIZE; k++) {
                System_printf("%u, ", benchmarks[k] - gtrackStartTime);
                gtrackStartTime = benchmarks[k];
            }
            System_printf("%u)\n", benchmarks[k] - gtrackStartTime);
		}
		fclose(fCloud);
	}
	gtrack_delete(hTrackModule);
	System_printf("Gtrack is using %d bytes of data memory used\n", memoryBytesUsed);
	
    /* Get the heap statistics at the beginning of the tests */
    HeapMem_getStats (heap0, &endMemoryStats);
    System_printf ("HeapMem: Used = %d bytes, Free = %d bytes\n", endMemoryStats.totalSize - peakMemoryStats.totalFreeSize, endMemoryStats.totalFreeSize);
	return;
}
