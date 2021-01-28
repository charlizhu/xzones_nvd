/**
 *   @file  task_mbox.c
 *
 *   @brief
 *     MSS main implementation of the millimeter wave Demo
 *
 *  Copyright (C) 2017 Texas Instruments Incorporated - http://www.ti.com/
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
#include <math.h>


/* BIOS/XDC Include Files. */
#include <xdc/runtime/System.h>

/* mmWave SDK Include Files: */
#include <ti/common/sys_common.h>
#include <ti/drivers/uart/UART.h>
#include <ti/utils/cli/cli.h>
#include <ti/alg/gtrack/gtrack.h>
#include <ti/drivers/osal/MemoryP.h>
#include "float.h"

/* Demo Include Files */
#include <mmw_mss.h>


/**************************************************************************
 *************************** Local Definitions ****************************
 **************************************************************************/

/**************************************************************************
 *************************** Global Definitions ***************************
 **************************************************************************/

typedef enum {
    TRACKING_DEFAULT_PARAM_SET = 0,
    TRACKING_TRAFFIC_MONITORING_PARAM_SET,
    TRACKING_PEOPLE_COUNTING_PARAM_SET,
    TRACKING_OUTDOOR_PARAM_SET,	
    TRACKING_CEILING_MOUNT_PARAM_SET	
} TRACKING_ADVANCED_PARAM_SET;

typedef enum {
    TRACKING_PARAM_SET_TM = 0,
    TRACKING_PARAM_SET_PC,
    TRACKING_PARAM_SET_OUTDOOR,
    TRACKING_PARAM_SET_CEILING_MOUNT
} TRACKING_ADVANCED_PARAM_SET_TABLE;

/**
 * @brief
 *  Global Variable for tracking information required by the mmw Demo
 */
extern MmwDemo_MSS_MCB  gMmwMssMCB;

unsigned int gGtrackMemoryUsed = 0;

/* @TODO: These functions need to be abstracted to the DPC */
void *gtrack_alloc(unsigned int numElements, unsigned int sizeInBytes)
{
	gGtrackMemoryUsed += numElements*sizeInBytes;
    return MemoryP_ctrlAlloc(numElements*sizeInBytes, 0);
}
void gtrack_free(void *pFree, unsigned int sizeInBytes)
{
	gGtrackMemoryUsed -= sizeInBytes;
	MemoryP_ctrlFree(pFree,sizeInBytes);
}
void gtrack_log(GTRACK_VERBOSE_TYPE level, const char *format, ...)
{
#if 0
	va_list args;
    va_start(args, format);
	vprintf(format, args);
	va_end(args);
#endif
}

/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for tracking configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
int32_t MmwDemo_CLITrackingCfg (int32_t argc, char* argv[])
{
    GTRACK_moduleConfig         config;


    TRACKING_ADVANCED_PARAM_SET trackingParamSet;

    //Memory_Stats            startMemoryStats;
    //Memory_Stats            endMemoryStats;

    if (argc >= 1) {
        gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.trackerEnabled = (uint16_t) atoi (argv[1]);
    }

    if(gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.trackerEnabled != 1) {
        return 0;
    }

    /* Sanity Check: Minimum argument check */
    if (argc != 9)
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }
#if 0
    System_printf("Debug: Heap before creating a tracker\n");
    MmwDemo_printHeapStats();
#endif
    /* Initialize CLI configuration: */
    memset ((void *)&config, 0, sizeof(GTRACK_moduleConfig));

    trackingParamSet = (TRACKING_ADVANCED_PARAM_SET) atoi (argv[2]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.trackingParamSet = trackingParamSet;

    switch(trackingParamSet)
    {
        case TRACKING_DEFAULT_PARAM_SET:
            // Do not configure advanced parameters, use library default parameters
            config.advParams = 0;
            break;

        case TRACKING_TRAFFIC_MONITORING_PARAM_SET:
            /* Initialize CLI configuration: */
            config.initialRadialVelocity = -8.0f; 	// for TM, detected targets are approaching
            config.maxAcceleration[0] = 2.0f;     	// for TM, maximum acceleration in lateral direction is set to 2m/s2
            config.maxAcceleration[1] = 20.0f;     	// for TM, maximum acceleration is longitudinal direction set to 20m/s2
            config.maxAcceleration[2] = 0.f;		// ignored
            break;

        case TRACKING_PEOPLE_COUNTING_PARAM_SET:
            /* Initialize CLI configuration: */
            config.initialRadialVelocity = 0;	//For PC, detected target velocity is unknown
            config.maxAcceleration[0] = 1.0f;
            config.maxAcceleration[1] = 1.0f;
            config.maxAcceleration[2] = 1.0f;
            break;

        case TRACKING_OUTDOOR_PARAM_SET:
            /* Initialize CLI configuration: */
            config.initialRadialVelocity = 0; 	// for OUTDOOR, detected targets velocity is unknown 
            config.maxAcceleration[0] = 1.0f;
            config.maxAcceleration[1] = 1.0f;
            config.maxAcceleration[2] = 1.0f;
            break;			

        case TRACKING_CEILING_MOUNT_PARAM_SET:
            /* Initialize CLI configuration: */
            config.initialRadialVelocity = 0;
            config.maxAcceleration[0] = 0.1f;
            config.maxAcceleration[1] = 0.1f;
            config.maxAcceleration[2] = 0.1f;
            break;


        default:
            CLI_write ("Error: Invalid usage of the CLI command\n");
            return -1;
    }
#ifdef GTRACK_3D
    config.stateVectorType      	= GTRACK_STATE_VECTORS_3DA; // Track three dimensions with acceleration
#else
    config.stateVectorType      	= GTRACK_STATE_VECTORS_2DA; // Track two dimensions with acceleration
#endif
    config.verbose              	= GTRACK_VERBOSE_NONE;
    config.maxNumPoints         	= (uint16_t) atoi(argv[3]);
    config.maxNumTracks         	= (uint16_t) atoi(argv[4]);
    config.maxRadialVelocity   		= (float) atoi(argv[5]) *0.1f;
#ifndef GTRACK_3D
    config.radialVelocityResolution	= (float) atoi(argv[6]) *0.001f;
#endif
    config.deltaT               	= (float) atoi(argv[7]) *0.001f;

    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sensorAzimuthTilt = (90-atoi(argv[8]))*3.14f/180;
#if 0
    if(gMmwMCB.cfg.trackingCfg.config.maxNumPoints != 0) {
        pointCloudSize = gMmwMCB.cfg.trackingCfg.config.maxNumPoints*sizeof(GTRACK_measurementPoint);
        if(gMmwMCB.pointCloud != NULL) {
            MemoryP_ctrlFree(gMmwMCB.pointCloud, pointCloudSize);
        }
        if(gMmwMCB.pointCloudTlv != NULL) {
            pointCloudTlvSize = sizeof(MmwDemo_output_message_tl) + gMmwMCB.cfg.trackingCfg.config.maxNumPoints*sizeof(GTRACK_measurementPoint);
            MemoryP_ctrlFree(gMmwMCB.pointCloudTlv, pointCloudTlvSize);
        }
    }
    if(gMmwMCB.cfg.trackingCfg.config.maxNumTracks != 0) {
        targetListSize = sizeof(MmwDemo_output_message_tl) + gMmwMCB.cfg.trackingCfg.config.maxNumTracks*sizeof(GTRACK_targetDesc);
        targetIndexSize = sizeof(MmwDemo_output_message_tl) + gMmwMCB.cfg.trackingCfg.config.maxNumPoints*sizeof(uint8_t);
        if(gMmwMCB.targetDescrHandle != NULL) {
            /* Free Target List Arrays */
            if(gMmwMCB.targetDescrHandle->tList[0] != NULL)
                MemoryP_ctrlFree(gMmwMCB.targetDescrHandle->tList[0], targetListSize);
            if(gMmwMCB.targetDescrHandle->tList[1] != NULL)
                MemoryP_ctrlFree(gMmwMCB.targetDescrHandle->tList[1], targetListSize);

            /* Free Target Index Arrays */
            if(gMmwMCB.targetDescrHandle->tIndex[0] != NULL)
                MemoryP_ctrlFree(gMmwMCB.targetDescrHandle->tIndex[0], targetIndexSize);
            if(gMmwMCB.targetDescrHandle->tIndex[1] != NULL)
                MemoryP_ctrlFree(gMmwMCB.targetDescrHandle->tIndex[1], targetIndexSize);
        }
    }
    if(gMmwMCB.targetDescrHandle != NULL) {
        MemoryP_ctrlFree(gMmwMCB.targetDescrHandle, sizeof(MmwDemo_targetDescrHandle));
    }
#endif
    /* Save Configuration to use later */
    memcpy((void *)&gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gtrackModuleConfig, (void *)&config, sizeof(GTRACK_moduleConfig));
    //memcpy((void *)&gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gTrackAdvParams, (void *)&advParams, sizeof(GTRACK_advancedParameters));

    return 0;
}

int32_t MmwDemo_CLIStaticBoundaryBoxCfg (int32_t argc, char* argv[])
{
    /* Sanity Check: Minimum argument check */
#ifdef GTRACK_3D
    if (argc != 7)
#else
    if (argc != 5)
#endif
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }

    /* Initialize the ADC Output configuration: */
    //memset ((void *)&sceneryParams, 0, sizeof(GTRACK_sceneryParams));
    
    /* Populate configuration: */
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.numStaticBoxes = 1;
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.staticBox[0].x1 = (float) atof (argv[1]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.staticBox[0].x2 = (float) atof (argv[2]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.staticBox[0].y1 = (float) atof (argv[3]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.staticBox[0].y2 = (float) atof (argv[4]);
    
#ifdef GTRACK_3D
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.staticBox[0].z1 = (float) atof (argv[5]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.staticBox[0].z2 = (float) atof (argv[6]);
#endif

    return 0;
    
}

int32_t MmwDemo_CLIBoundaryBoxCfg (int32_t argc, char* argv[])
{
    /* Sanity Check: Minimum argument check */
#ifdef GTRACK_3D
    if (argc != 7)
#else
    if (argc != 5)
#endif
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }

    /* Initialize the ADC Output configuration: */
    //memset ((void *)&sceneryParams, 0, sizeof(GTRACK_sceneryParams));
    
    /* Populate configuration: */
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.numBoundaryBoxes = 1;
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.boundaryBox[0].x1 = (float) atof (argv[1]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.boundaryBox[0].x2 = (float) atof (argv[2]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.boundaryBox[0].y1 = (float) atof (argv[3]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.boundaryBox[0].y2 = (float) atof (argv[4]);
#ifdef GTRACK_3D
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.boundaryBox[0].z1 = (float) atof (argv[5]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.sceneryParams.boundaryBox[0].z2 = (float) atof (argv[6]);
#endif
    
    return 0;
    
}

/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for GatingParam configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
int32_t MmwDemo_CLIGatingParamCfg (int32_t argc, char* argv[])
{
    /* Sanity Check: Minimum argument check */
#ifdef GTRACK_3D
    if (argc != 6)
#else
    if (argc != 5)
#endif
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }
    /* Initialize the ADC Output configuration: */
    //memset ((void *)&gatingParams, 0, sizeof(GTRACK_gatingParams));
    
    /* Populate configuration: */
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gatingParams.gain = (float) atof (argv[1]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gatingParams.limits.width = (float) atof (argv[2]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gatingParams.limits.depth = (float) atof (argv[3]);
#ifdef GTRACK_3D
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gatingParams.limits.height = (float) atof (argv[4]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gatingParams.limits.vel = (float) atof (argv[5]);
#else
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.gatingParams.limits.vel = (float) atof (argv[4]);
#endif

    return 0;
    
}


/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for StateParam configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
int32_t MmwDemo_CLIStateParamCfg (int32_t argc, char* argv[])
{
    /* Sanity Check: Minimum argument check */
    if (argc != 6)
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }

    /* Initialize the ADC Output configuration: */
    //memset ((void *)&stateParams, 0, sizeof(GTRACK_stateParams));
    
    /* Populate configuration: */
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.stateParams.det2actThre = (uint16_t) atoi (argv[1]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.stateParams.det2freeThre= (uint16_t) atoi (argv[2]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.stateParams.active2freeThre= (uint16_t) atoi (argv[3]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.stateParams.static2freeThre= (uint16_t) atoi (argv[4]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.stateParams.exit2freeThre= (uint16_t) atoi (argv[5]);

    return 0;
    
}

/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for AllocationParam configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
int32_t MmwDemo_CLIAllocationParamCfg (int32_t argc, char* argv[])
{
    /* Sanity Check: Minimum argument check */
#ifdef GTRACK_3D
    if (argc != 7)
#else
    if (argc != 6)
#endif
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }

    /* Initialize the ADC Output configuration: */
    //memset ((void *)&allocationParams, 0, sizeof(GTRACK_allocationParams));

    /* Populate configuration: */
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.snrThre = (float) atof (argv[1]);
#ifdef GTRACK_3D
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.snrThreObscured = (float) atof (argv[2]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.velocityThre = (float) atof (argv[3]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.pointsThre = (uint16_t) atoi (argv[4]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.maxDistanceThre = (float) atof (argv[5]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.maxVelThre = (float) atof (argv[6]);
#else
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.velocityThre = (float) atof (argv[2]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.pointsThre = (uint16_t) atoi (argv[3]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.maxDistanceThre = (float) atof (argv[4]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.allocationParams.maxVelThre = (float) atof (argv[5]);
#endif

    return 0;
    
}

int32_t MmwDemoCLIMaxAccelerationParamCfg(int32_t argc, char* argv[])
{
#ifdef GTRACK_3D
    if (argc != 4)
#else
    if (argc != 3)
#endif
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }
    
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.accelerationParams[0] = (float) atof (argv[1]);
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.accelerationParams[1] = (float) atof (argv[2]);
#ifdef GTRACK_3D    
    gMmwMssMCB.trackerCfg.trackerDpuCfg.staticCfg.accelerationParams[2] = (float) atof (argv[3]);
#endif
    
    return 0;

}
