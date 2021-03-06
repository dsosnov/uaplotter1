SRC=../src
OBJ=../lib
INTERFACE=../interface
CMS=../../UATree/UADataFormat
TOT=../../TOTEMdataFormat
CXXFLAGS+=$(shell root-config --cflags) -Wall -I$(SRC) -I$(INTERFACE) -I$(CMS)/interface -I$(TOT)/interface -ggdb3 -fPIC -std=gnu++14 -O2
LDFLAGS+=$(shell root-config --glibs)  -ggdb3 
LINKD=uaplot_LinkDef.h
HEAD=$(filter-out $(INTERFACE)/$(LINKD) $(INTERFACE)/uaplot_dict.h, $(wildcard $(INTERFACE)/*.h))
CODE=$(wildcard $(SRC)/ua*.cc)
OBJS=$(subst $(SRC),$(OBJ),$(CODE:.cc=.o))

all: $(OBJ) $(OBJ)/libuaplotter.so $(CMS)/lib/libUADataFormat.so $(TOT)/lib/libTOTEMdataFormat.so docs
necessary: $(OBJ) $(OBJ)/libuaplotter.so

$(OBJ)/libuaplotter.so: $(OBJS) $(OBJ)/uaplot_dict.o 
	$(CXX) -shared -o $@ $^ $(LDFLAGS) #order important for Ubuntu 
	if [ ! -s $(OBJ)/uaplot_dict_rdict.pcm ] ; then ln -s $(SRC)/uaplot_dict_rdict.pcm $(OBJ)/ ; fi
#	@echo ""

$(OBJ)/%.o: $(SRC)/%.cc $(HEAD)
	$(CXX) $(CXXFLAGS) -c $< -o $@
#	@echo ""

$(SRC)/uaplot_dict.cc: $(HEAD)
	@echo LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}
	cd $(SRC); rootcint -f uaplot_dict.cc -c -p -I$(INTERFACE) -I$(CMS)/interface -I$(TOT)/interface -D_POSIX_C_SOURCE=200809L $(subst $(SRC)/,,$(HEAD)) $(LINKD)
#	@echo ""

# $(LINKD).h:
#       #compatibility with cmssw
#       mv $(LINKD)h  $(LINKD).h
#       @echo ""

$(OBJ):
	mkdir -p $(OBJ)

info:
	@echo HEADERS: $(HEAD)
	@echo SOURCES: $(CODE)
	@echo OBJECTS: $(OBJS)

clean:
	rm -rf $(OBJ) $(SRC)/uaplot_dict*

clean_all: clean
	cd $(CMS)/utilities; make clean
	cd $(TOT)/utilities; make clean

docs:
	@doxygen

$(CMS)/lib/libUADataFormat.so:
	cd $(CMS)/utilities; $(MAKE)

$(TOT)/lib/libTOTEMdataFormat.so:
	cd $(TOT)/utilities; $(MAKE)
