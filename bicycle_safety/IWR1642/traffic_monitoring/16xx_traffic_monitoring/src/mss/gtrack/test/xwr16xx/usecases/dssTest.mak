###################################################################################
# GTRACK Usecase Unit Test on DSS Makefile
###################################################################################
.PHONY: usecaseDssTest usecaseDssTestClean

###################################################################################
# Setup the VPATH:
###################################################################################
vpath %.c src
vpath %.c test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases
vpath %.c test/common

###################################################################################
# The GTRACK Unit test requires additional libraries
###################################################################################
GTRACK_USECASE_DSS_TEST_STD_LIBS = $(C674_COMMON_STD_LIB)						\
						  -llibgtrack_$(MMWAVE_SDK_DEVICE_TYPE).$(C674_LIB_EXT)
GTRACK_USECASE_DSS_TEST_LOC_LIBS = $(C674_COMMON_LOC_LIB)						\
						  -ilib

###################################################################################
# Unit Test Files
###################################################################################
GTRACK_USECASE_DSS_TEST_CFG		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/dss.cfg
GTRACK_USECASE_DSS_TEST_CMD		= $(MMWAVE_SDK_INSTALL_PATH)/ti/platform/$(MMWAVE_SDK_DEVICE_TYPE)
GTRACK_USECASE_DSS_TEST_CONFIGPKG	= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/dss_configPkg_$(MMWAVE_SDK_DEVICE_TYPE)
GTRACK_USECASE_DSS_TEST_MAP		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/$(MMWAVE_SDK_DEVICE_TYPE)_gtrack_usecase_dss.map
GTRACK_USECASE_DSS_TEST_OUT		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/$(MMWAVE_SDK_DEVICE_TYPE)_gtrack_usecase_dss.$(C674_EXE_EXT)
GTRACK_USECASE_DSS_TEST_BIN		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/$(MMWAVE_SDK_DEVICE_TYPE)_gtrack_usecase_dss.bin
GTRACK_USECASE_DSS_TEST_APP_CMD	= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/dss_gtrack_linker.cmd
GTRACK_USECASE_DSS_TEST_SOURCES	= main_dss.c			\
							  gtrackApp.c		\
							  gtrackAlloc.c		\
							  gtrackLog.c

GTRACK_USECASE_DSS_TEST_DEPENDS = $(addprefix $(PLATFORM_OBJDIR)/, $(GTRACK_USECASE_DSS_TEST_SOURCES:.c=.$(C674_DEP_EXT)))
GTRACK_USECASE_DSS_TEST_OBJECTS = $(addprefix $(PLATFORM_OBJDIR)/, $(GTRACK_USECASE_DSS_TEST_SOURCES:.c=.$(C674_OBJ_EXT)))

###################################################################################
# RTSC Configuration:
###################################################################################
usecaseDssRTSC:
	@echo 'Configuring RTSC packages...'
	$(XS) --xdcpath="$(XDCPATH)" xdc.tools.configuro $(C674_XSFLAGS) -o $(GTRACK_USECASE_DSS_TEST_CONFIGPKG) $(GTRACK_USECASE_DSS_TEST_CFG)
	@echo 'Finished configuring packages'
	@echo ' '

###################################################################################
# Build Unit Test:
###################################################################################
usecaseDssTest: BUILD_CONFIGPKG=$(GTRACK_USECASE_DSS_TEST_CONFIGPKG)
usecaseDssTest: C674_CFLAGS += --cmd_file=$(BUILD_CONFIGPKG)/compiler.opt
usecaseDssTest: buildDirectories usecaseDssRTSC $(GTRACK_USECASE_DSS_TEST_OBJECTS)
	$(C674_LD) $(C674_LDFLAGS) $(GTRACK_USECASE_DSS_TEST_LOC_LIBS) $(GTRACK_USECASE_DSS_TEST_STD_LIBS) 	\
	-l$(GTRACK_USECASE_DSS_TEST_CONFIGPKG)/linker.cmd --map_file=$(GTRACK_USECASE_DSS_TEST_MAP) 		\
	$(GTRACK_USECASE_DSS_TEST_OBJECTS) $(PLATFORM_C674X_LINK_CMD) $(GTRACK_USECASE_DSS_TEST_APP_CMD) 	\
	$(C674_LD_RTS_FLAGS) -o $(GTRACK_USECASE_DSS_TEST_OUT)
	@echo '******************************************************************************'
	@echo 'Built the GTRACK Usecase DSS Unit Test '
	@echo '******************************************************************************'

###################################################################################
# Cleanup Unit Test:
###################################################################################
usecaseDssTestClean:
	@echo 'Cleaning the GTRACK Usecase DSS Unit Test objects'
	@$(DEL) $(GTRACK_USECASE_DSS_TEST_OBJECTS) $(GTRACK_USECASE_DSS_TEST_OUT) $(GTRACK_USECASE_DSS_TEST_BIN)
	@$(DEL) $(GTRACK_USECASE_DSS_TEST_MAP) $(GTRACK_USECASE_DSS_TEST_DEPENDS)
	@echo 'Cleaning the GTRACK Usecase DSS Unit RTSC package'
	@$(DEL) $(GTRACK_USECASE_DSS_TEST_CONFIGPKG)
	@$(DEL) $(PLATFORM_OBJDIR)

###################################################################################
# Dependency handling
###################################################################################
-include $(GTRACK_USECASE_DSS_TEST_DEPENDS)

