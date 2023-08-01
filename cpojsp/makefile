CCC = g++
CPOHOME = /storagespace/store1/CPAIExpe/src
CONCERTDIR = $(CPOHOME)/concert
CPLEXDIR = $(CPOHOME)/cplex
CFLAGS = -DIL_STD -O -DNDEBUG -I$(CPOHOME)/cpoptimizer/include -I$(CONCERTDIR)/include -fPIC -fstrict-aliasing -pedantic -Wall -fexceptions -Wno-long-long -m64 -DILOUSEMT -D_REENTRANT -DILM_REENTRANT
LDFLAGS = -L$(CPOHOME)/cpoptimizer/lib/x86-64_sles10_4.1/static_pic -lcp -L$(CPLEXDIR)/lib/x86-64_sles10_4.1/static_pic -lcplex -L$(CONCERTDIR)/lib/x86-64_sles10_4.1/static_pic -lconcert  -lpthread 
JAVA = $(JAVA_HOME)/bin/java
JAVAC = $(JAVA_HOME)/bin/javac
JAVAFLAGS = -d64
CPCLASSPATH = $(CPOHOME)/cpoptimizer/lib/ILOG.CP.jar
LIBRARYPATH = $(CPOHOME)/cpoptimizer/bin/x86-64_sles10_4.1:$(CPLEXDIR)/bin/x86-64_sles10_4.1
CPPEXDIR = src
JAVAEXDIR = $(CPOHOME)/cpoptimizer/examples/src/java

all: cpp

cpp: cpo_sched 

clean:
	rm -f cpo_sched  

FORCE:

cpo_sched: $(CPPEXDIR)/cpo_sched.cpp $(CPPEXDIR)/osp.hpp $(CPPEXDIR)/jsp.hpp $(CPPEXDIR)/jet.hpp $(CPPEXDIR)/reader.hpp
	$(CCC) -o cpo_sched $(CFLAGS) $(CPPEXDIR)/cpo_sched.cpp $(LDFLAGS)

