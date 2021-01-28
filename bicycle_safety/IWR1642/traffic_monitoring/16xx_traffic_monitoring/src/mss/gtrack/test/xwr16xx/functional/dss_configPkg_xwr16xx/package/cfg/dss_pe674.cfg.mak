# invoke SourceDir generated makefile for dss.pe674
dss.pe674: .libraries,dss.pe674
.libraries,dss.pe674: package/cfg/dss_pe674.xdl
	$(MAKE) -f package/cfg/dss_pe674.src/makefile.libs

clean::
	$(MAKE) -f package/cfg/dss_pe674.src/makefile.libs clean

