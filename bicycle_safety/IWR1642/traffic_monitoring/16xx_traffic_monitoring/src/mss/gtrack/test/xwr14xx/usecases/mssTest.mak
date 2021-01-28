###################################################################################
# GTRACK Usecase Unit Test on MSS Makefile
###################################################################################
.PHONY: usecaseMssTest usecaseMssTestClean

###################################################################################
# Setup the VPATH:
###################################################################################
vpath %.c src
vpath %.c test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases
vpath %.c test/common

###################################################################################
# The GTRACK Unit test requires additional libraries
###################################################################################
GTRACK_USECASE_MSS_TEST_STD_LIBS = $(R4F_COMMON_STD_LIB)						\
						  -llibgtrack_$(MMWAVE_SDK_DEVICE_TYPE).$(R4F_LIB_EXT)
GTRACK_USECASE_MSS_TEST_LOC_LIBS = $(R4F_COMMON_LOC_LIB)						\
						  -ilib

###################################################################################
# Unit Test Files
###################################################################################
GTRACK_USECASE_MSS_TEST_CFG		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/mss.cfg
GTRACK_USECASE_MSS_TEST_CMD		= $(MMWAVE_SDK_INSTALL_PATH)/ti/platform/$(MMWAVE_SDK_DEVICE_TYPE)
GTRACK_USECASE_MSS_TEST_CONFIGPKG	= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/mss_configPkg_$(MMWAVE_SDK_DEVICE_TYPE)
GTRACK_USECASE_MSS_TEST_MAP		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/$(MMWAVE_SDK_DEVICE_TYPE)_gtrack_usecase_mss.map
GTRACK_USECASE_MSS_TEST_OUT		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/$(MMWAVE_SDK_DEVICE_TYPE)_gtrack_usecase_mss.$(R4F_EXE_EXT)
GTRACK_USECASE_MSS_TEST_BIN		= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/$(MMWAVE_SDK_DEVICE_TYPE)_gtrack_usecase_mss.bin
GTRACK_USECASE_MSS_TEST_APP_CMD	= test/$(MMWAVE_SDK_DEVICE_TYPE)/usecases/mss_gtrack_linker.cmd
GTRACK_USECASE_MSS_TEST_SOURCES	= main.c					\
							  gtrackApp.c				\
							  gtrackAlloc.c				\
							  gtrackLog.c
GTRACK_USECASE_MSS_TEST_DEPENDS 	 = $(addprefix $(PLATFORM_OBJDIR)/, $(GTRACK_USECASE_MSS_TEST_SOURCES:.c=.$(R4F_DEP_EXT)))
GTRACK_USECASE_MSS_TEST_OBJECTS 	 = $(addprefix $(PLATFORM_OBJDIR)/, $(GTRACK_USECASE_MSS_TEST_SOURCES:.c=.$(R4F_OBJ_EXT)))

###################################################################################
# RTSC Configuration:
###################################################################################
usecaseMssRTSC:
	@echo 'Configuring RTSC packages...'
	$(XS) --xdcpath="$(XDCPATH)" xdc.tools.configuro $(R4F_XSFLAGS) -o $(GTRACK_USECASE_MSS_TEST_CONFIGPKG) $(GTRACK_USECASE_MSS_TEST_CFG)
	@echo 'Finished configuring packages'
	@echo ' '

###################################################################################
# Build Unit Test:
###################################################################################
usecaseMssTest: BUILD_CONFIGPKG=$(GTRACK_USECASE_MSS_TEST_CONFIGPKG)
usecaseMssTest: R4F_CFLAGS += --cmd_file=$(BUILD_CONFIGPKG)/compiler.opt
usecaseMssTest: buildDirectories usecaseMssRTSC $(GTRACK_USECASE_MSS_TEST_OBJECTS)
	$(R4F_LD) $(R4F_LDFLAGS) $(GTRACK_USECASE_MSS_TEST_LOC_LIBS) $(GTRACK_USECASE_MSS_TEST_STD_LIBS) 	\
	-l$(GTRACK_USECASE_MSS_TEST_CONFIGPKG)/linker.cmd --map_file=$(GTRACK_USECASE_MSS_TEST_MAP) 		\
	$(GTRACK_USECASE_MSS_TEST_OBJECTS) $(PLATFORM_R4F_LINK_CMD) $(GTRACK_USECASE_MSS_TEST_APP_CMD) 	\
	$(R4F_LD_RTS_FLAGS) -o $(GTRACK_USECASE_MSS_TEST_OUT)
	@echo '******************************************************************************'
	@echo 'Built the GTRACK Usecase MSS Unit Test '
	@echo '******************************************************************************'

###################################################################################
# Cleanup Unit Test:
###################################################################################
usecaseMssTestClean:
	@echo 'Cleaning the GTRACK Usecase MSS Unit Test objects'
	@$(DEL) $(GTRACK_USECASE_MSS_TEST_OBJECTS) $(GTRACK_USECASE_MSS_TEST_OUT) $(GTRACK_USECASE_MSS_TEST_BIN)
	@$(DEL) $(GTRACK_USECASE_MSS_TEST_MAP) $(GTRACK_USECASE_MSS_TEST_DEPENDS)
	@echo 'Cleaning the GTRACK Usecase MSS Unit RTSC package'
	@$(DEL) $(GTRACK_USECASE_MSS_TEST_CONFIGPKG)
	@$(DEL) $(PLATFORM_OBJDIR)

###################################################################################
# Dependency handling
###################################################################################
-include $(GTRACK_USECASE_MSS_TEST_DEPENDS)

