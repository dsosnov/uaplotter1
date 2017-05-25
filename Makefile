SRC=src
OBJ=lib
CMS=/home/kkuzn/pA13/CMSTotem/CMSdataFormat
TOT=/home/kkuzn/pA13/CMSTotem/TOTEMdataFormat
RD=/data/SWnew/FairSoft/tools/root/bin
export LD_LIBRARY_PATH=/data/SWnew/FairSoftInst/lib/root
CXXFLAGS+=$(shell $(RD)/root-config --cflags) -Wall -I$(SRC) -I$(CMS)/src -I$(TOT)/src -ggdb3 -fPIC
LDFLAGS+=$(shell $(RD)/root-config --glibs)  -ggdb3 
LINKD=uaplot_LinkDef.h
HEAD=$(filter-out $(SRC)/$(LINKD) $(SRC)/uaplot_dict.h, $(wildcard $(SRC)/*.h))
CODE=$(wildcard $(SRC)/ua*.cc)
OBJS=$(subst $(SRC),$(OBJ),$(CODE:.cc=.o))



all: $(OBJ) $(OBJ)/libuaplotter.so 
$(OBJ)/libuaplotter.so: $(OBJS) $(OBJ)/uaplot_dict.o 
	$(CXX) -shared -o $@ $^ $(LDFLAGS) #order important for Ubuntu 
	@echo ""

$(OBJ)/%.o: $(SRC)/%.cc $(HEAD)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo ""

$(SRC)/uaplot_dict.cc: $(HEAD)
	@echo LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}
	cd $(SRC);$(RD)/rootcint -f uaplot_dict.cc -c -p  -I$(CMS)/src -I$(TOT)/src $(subst $(SRC)/,,$(HEAD)) $(LINKD)
	@echo ""
# 
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
	rm -rf $(OBJS) $(OBJ)/libuaplotter.so $(SRC)/uaplot_dict.cc $(SRC)/uaplot_dict.h
