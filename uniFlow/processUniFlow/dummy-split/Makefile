#
# Project variables
# =================
#
NAME			:= vp
LIBNAME			:= VP

#
# ROOT parameters for building
# ============================
#
# Use the script `root-config` to extract some variables we need when
# building.  Note, we massage some fo the return values from the
# script to better suite our needs
#
ROOT_CONFIG		:= $(shell which root-config)
ROOT_MAJOR		:= $(basename \
			     $(subst /,,\
			       $(dir $(shell $(ROOT_CONFIG) --version))))
ROOT_CXX		:= $(shell $(ROOT_CONFIG) --cxx)
ROOT_LD			:= $(shell $(ROOT_CONFIG) --ld)
ROOT_CPPFLAGS		:= -I$(shell $(ROOT_CONFIG) --incdir)
ROOT_CXXFLAGS		:= $(filter-out $(ROOT_CPPFLAGS),\
			     $(shell $(ROOT_CONFIG) --cflags))
ROOT_LIBS		:= $(shell $(ROOT_CONFIG) \
			     --noldflags --noauxlibs --libs)
ROOT_LDFLAGS		:= $(filter-out $(ROOT_LIBS),\
			     $(shell $(ROOT_CONFIG) --libs))
ROOT_EXE		:= $(shell root-config --bindir)/root.exe
# What dictionary generator we should use depend on ROOT major version
ifeq ($(ROOT_MAJOR),6)
ROOT_DICTGEN		:= $(shell root-config --bindir)/rootcling
else
ROOT_DICTGEN		:= $(shell root-config --bindir)/rootcint -v4 -f
endif

#
# System specific flags
# =====================
#
# What flag do we need for Position-Independent-Code (PIC).  Most use
# `-fPIC`.
#
PIC_FLAG		:= -fPIC

#
# Final values
# ============
#
CXX			:= $(ROOT_CXX)
LD			:= $(ROOT_LD)
SO			:= $(ROOT_LD)
CPPFLAGS		:= $(ROOT_CPPFLAGS) -I.
CXXFLAGS		:= $(ROOT_CXXFLAGS) $(PIC_FLAG) -c
LDFLAGS			:= $(ROOT_LDFLAGS)
SOFLAGS			:= $(ROOT_LDFLAGS) -shared
LIBS			:= $(ROOT_LIBS)

#
# Our code
# ========
#
HEADERS			:= Base.h		\
			   Other.h		\
			   Child1.h		\
			   Child2.h		\
			   GrandChild1.h	\
			   GrandChild2.h	\
			   VPLinkDef.h		\
			   VPInclude.h
SOURCES			:= Base.cxx		\
			   Other.cxx		\
			   Child2.cxx		\
			   GrandChild1.cxx
SCRIPTS			:= Load.C Run.C
DICT			:= VPDict.cxx
OBJECTS			:= $(DICT:.cxx=.o) $(SOURCES:.cxx=.o)

#
# Generic rules (a.k.a. pattern rules)
# ====================================
#
# Rule to make an object file from a source file
#
%.o:	%.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $<

# Make dictionary source code
%Dict.cxx:%Include.h %LinkDef.h
	$(ROOT_DICTGEN) $@ -c $(CPPFLAGS) $^

# Make dictionary header
%Dict.h:%Dict.cxx
	@if test -f $<; then rm -f $<; $(MAKE) $(MAKEFLAGS) $< ; fi

# Make loadable module
%.so:%Dict.o
	$(SO) $(SOFLAGS) -o $@ $(LIBS) $^

#
# Specific rules and dependencies
# ===============================
#
all:	$(LIBNAME).so

$(LIBNAME).so:	$(OBJECTS)
$(LIBNAME)Include.h: $(filter-out $(LIBNAME)Include.h, $(HEADERS))
$(LIBNAME)Dict.cxx: $(LIBNAME)Include.h $(LIBNAME)LinkDef.h

clean:
	rm -f *~ *ACLiC* *_*.d *Dict.* VP.so *.so *.o $(NAME).tar.gz *.pcm

show:
	@echo "ROOT_CONFIG    = $(ROOT_CONFIG)"
	@echo "ROOT_MAJOR     = $(ROOT_MAJOR)"
	@echo "ROOT_CXX       = $(ROOT_CXX)"
	@echo "ROOT_LD        = $(ROOT_LD)"
	@echo "ROOT_CPPFLAGS  = $(ROOT_CPPFLAGS)"
	@echo "ROOT_CXXFLAGS  = $(ROOT_CXXFLAGS)"
	@echo "ROOT_LIBS      = $(ROOT_LIBS)"
	@echo "ROOT_LDFLAGS   = $(ROOT_LDFLAGS)"
	@echo "ROOT_DICTGEN   = $(ROOT_DICTGEN)"
	@echo "DICT           = $(DICT)"
	@echo "SOURCES	      = $(SOURCES)"
	@echo "OBJECTS	      = $(OBJECTS)"


test_so:
	$(MAKE) clean
	$(MAKE) $(LIBNAME).so
	echo "gSystem->Load(\"$(LIBNAME).so\");gROOT->Macro(\"Run.C\")" > bla
	$(ROOT_EXE) < bla
	rm -f bla


test_aclic:	$(SOURCES) $(HEADERS)
	$(MAKE) clean
	echo "gROOT->Macro(\"Load.C\");gROOT->Macro(\"Run.C\")" > bla
	$(ROOT_EXE) < bla
	rm -f bla

dist:
	mkdir $(NAME)
	cp Makefile $(SOURCES) $(HEADERS) $(SCRIPTS) $(NAME)/
	tar -czvf $(NAME).tar.gz $(NAME)/
	rm -rf $(NAME)

distcheck:	dist
	tar -xzvf $(NAME).tar.gz
	$(MAKE) -C $(NAME) test_so
	$(MAKE) -C $(NAME) test_aclic
	rm -rf $(NAME)

#
# EOF
#
