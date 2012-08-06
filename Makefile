#========================================
# USER RESERVED -- The following are reserved for users to set on the
# command line.  Makefiles should not set these.  These variables are
# for C/C++ compilation, and linking.
## -pg is for use with gprof.  Use on CFLAGS and on LDFLAGS
#CFLAGS		= -Wall -pg
#CFLAGS         = -Winline
#JFLAGS		=
#LDFLAGS	= -pg

# OPTIMIZE with the -O option.  Override from the command line for
# building debug versions.
#
#OPTFLAGS	= -O3 -funroll-loops -DNDEBUG=1

#========================================
### For Seqan support (required):
SEQAN_CFLAGS 	= -I./seqan
SEQAN_LDFLAGS 	=
SEQAN_LIBS	=

#========================================
## For the HMMoC-BFloat-Algebra library (required):
ALGEBRA_CFLAGS 	= -I./HMMoC-BFloat-Algebra
ALGEBRA_LDFLAGS = -L./HMMoC-BFloat-Algebra
ALGEBRA_LIBS	= -lHMMoC-BFloat-Algebra

#========================================
## For Boost support (required):
BOOST_CFLAGS 	= -I./boost-include
BOOST_LDFLAGS 	= -L./boost-lib
BOOST_LIBS	= -lboost_serialization -lboost_graph -lboost_filesystem -lboost_system -lboost_regex -lboost_program_options

###==============================================

#========================================
### For muscle support (optional):
# You must uncomment these for muscle support.
#MUSCLE_CPPSRC = $(sort $(wildcard muscle/*.cpp))
#MUSCLE_CPPOBJ_TMP = $(subst .cpp,.o,$(MUSCLE_CPPSRC))
#MUSCLE_CPPOBJ = $(subst muscle/main.o,,$(MUSCLE_CPPOBJ_TMP)) ;
#MUSCLE_CFLAGS = -D__HAVE_MUSCLE
#MUSCLE_LDFLAGS =
#MUSCLE_LIBS =

#========================================
### For HMMer (and squid) support (optional):
#HMMER_CFLAGS 	= -I./hmmer/src -I./hmmer/squid
#HMMER_LDFLAGS 	= -L./hmmer/src -L./hmmer/squid
#HMMER_LIBS	= -lsquid -lhmmer
HMMER_CFLAGS 	=
HMMER_LDFLAGS 	=
HMMER_LIBS	=

#========================================
### For ALGLIB support (optional):
ALGLIB_CPPSRC = $(sort $(wildcard alglib/*.cpp))
ALGLIB_CPPOBJ_TMP = $(subst .cpp,.o,$(ALGLIB_CPPSRC))
ALGLIB_CPPOBJ = $(subst alglib/main.o,,$(ALGLIB_CPPOBJ_TMP)) ;
ALGLIB_CFLAGS =
ALGLIB_LDFLAGS =
ALGLIB_LIBS =

###==============================================

INCS = Prolific.hpp \
Ambiguous.hpp \
Sequence.hpp \
MultinomialDistribution.hpp \
ProfileHMM.hpp \
Profile.hpp \
ProfileTree.hpp \
Fasta.hpp \
Random.hpp \
DynamicProgramming.hpp \
Parameters.hpp \
AminoAcid20.hpp

PROLIFICTESTS_INCS = $(INCS)

PROLIFICTESTS_OBJS = ProlificTests.o

PROLIFICTESTS_SOURCES = ProlificTests.cpp

default: prolifictests

prolifictests: $(PROLIFICTESTS_SOURCES) $(PROLIFICTESTS_INCS) $(PROLIFICTESTS_OBJS) $(MUSCLE_CPPOBJ)
	     $(CXX_LINK) -o prolifictests $(PROLIFICTESTS_OBJS) $(MUSCLE_CPPOBJ)

all: prolifictests

## Recompile if the includes are modified ...
$(PROLIFICTESTS_OBJS): $(PROLIFICTESTS_SOURCES) $(PROLIFICTESTS_INCS)

.PHONY: clean
clean:
	rm -f prolifictests $(PROLIFICTESTS_OBJS)

#========================================
# FILE EXTENSIONS.  Extensions and prefixes for different types of
# files change from platform to platform.  Hide these in macros so
# that we can more easily cut and paste between makefiles.
o		= .o
EXE_SFX		= 
SCRIPT_SFX 	= 
LIB_PFX		= lib
LIB_SFX		= .a
LIB_SHARED_SFX	= .so
TMPLIB		= libtemp.a

# FILE TOOLS
AR 	= ar qv
CHMOD 	= chmod
CP	= cp
GREP	= grep
MKDIR 	= mkdir
MUNCH 	= stepmunch
MV	= mv
NM 	= nm
RANLIB	= ranlib
RM 	= rm -f
RMDIR 	= rm -rf
STRIP	= strip
UNZIP 	= unzip
ZIP 	= zip


#========================================
# ANSI C Compile and Link
#
CC		= gcc
CC_COMPILE	= $(CC) -c $(OPTFLAGS) $(CFLAGS) $(CC_CFLAGS) $(CC_SYSCFLAGS)
CC_LINK		= $(CC) $(LDFLAGS) $(CC_LDFLAGS) $(CC_SYSLDFLAGS) $(CC_LIBS)
CC_CFLAGS 	= $(ALGEBRA_CFLAGS) $(BOOST_CFLAGS) $(SEQAN_CFLAGS) $(ALGLIB_CFLAGS) $(HMMER_CFLAGS) $(MUSCLE_CFLAGS)
CC_LDFLAGS	= $(ALGEBRA_LDFLAGS) $(BOOST_LDFLAGS) $(SEQAN_LDFLAGS) $(ALGLIB_LDFLAGS) $(HMMER_LDFLAGS) $(MUSCLE_LDFLAGS)
CC_LIBS		= $(ALGEBRA_LIBS) $(BOOST_LIBS) $(SEQAN_CLIBS) $(ALGLIB_CLIBS) $(HMMER_CLIBS) $(MUSCLE_CLIBS)

# Global system things used for compilation, static linking, etc.
CC_SYSCFLAGS 	= -I.
CC_SYSLDFLAGS 	=
CC_SYSLIBS	=

#========================================
# C++ Compile and Link
#
CXX		= g++
CXX_COMPILE	= $(CXX) -c  $(OPTFLAGS) $(CFLAGS) $(CXX_CFLAGS) $(CXX_SYSCFLAGS)
CXX_LINK	= $(CXX) $(LDFLAGS) $(CXX_LDFLAGS) $(CXX_SYSLDFLAGS) $(CXX_LIBS)
CXX_CFLAGS 	= $(ALGEBRA_CFLAGS) $(BOOST_CFLAGS) $(SEQAN_CFLAGS) $(ALGLIB_CFLAGS) $(HMMER_CFLAGS) $(MUSCLE_CFLAGS)
CXX_LDFLAGS	= $(ALGEBRA_LDFLAGS) $(BOOST_LDFLAGS) $(SEQAN_LDFLAGS) $(ALGLIB_LDFLAGS) $(HMMER_LDFLAGS) $(MUSCLE_LDFLAGS)
CXX_LIBS	= $(ALGEBRA_LIBS) $(BOOST_LIBS) $(SEQAN_LIBS) $(ALGLIB_LIBS) $(HMMER_LIBS) $(MUSCLE_LIBS)

# The force flags are used for C/C++ compilers that select the
# language based on the file naming conventions.  Some C++ source
# may be in files with C naming conventions.
CXX_FORCE	= 

# System Flags -- Things for static linking or making sure that the
# compiler understands that a file is a C++ file or whatever.  These
# usually change from platform to platform.
CXX_SYSCFLAGS 	= -I.
CXX_SYSLDFLAGS 	= 
CXX_SYSLIBS	= 

# Compilation Rules -- Repeat the rules for all of the different
# naming conventions.
#
.cxx.o:	; $(CXX_COMPILE) $<
.cpp.o:	; $(CXX_COMPILE) $<
.cc.o:	; $(CXX_COMPILE) $<
.C.o:	; $(CXX_COMPILE) $<

.cxx:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cpp:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cc:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.C:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)

# for legacy reasons also compile .c as c++
.c.o:	; $(CXX_COMPILE) $(CXX_FORCE) $<
.c:	
	$(CXX_COMPILE) $(CXX_FORCE) $<

