// -------------------------------------------------------------- -*- C++ -*-
// File: ./examples/src/cpp/sched_jobshop.cpp
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
//
// 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5725-A06 5725-A29
// Copyright IBM Corporation 1990, 2013. All Rights Reserved.
//
// Note to U.S. Government Users Restricted Rights:
// Use, duplication or disclosure restricted by GSA ADP Schedule
// Contract with IBM Corp.
// --------------------------------------------------------------------------

/* ------------------------------------------------------------

Problem Description
-------------------

In the classical Job-Shop Scheduling problem a finite set of jobs is
processed on a finite set of machines. Each job is characterized by a
fixed order of operations, each of which is to be processed on a
specific machine for a specified duration.  Each machine can process
at most one operation at a time and once an operation initiates
processing on a given machine it must complete processing
uninterrupted.  The objective of the problem is to find a schedule
that minimizes the makespan of the schedule.

------------------------------------------------------------ */

#include <ilcp/cp.h>
#include "reader.hpp"
#include "osp.hpp"
#include "jsp.hpp"
#include "jet.hpp"
#include "sds.hpp"
#include "jtl.hpp"




class FileError: public IloException {
public:
  FileError() : IloException("Cannot open data file") {}
};

class ArgError: public IloException {
public:
  ArgError() : IloException("Wrong arguments") {}
};

int main(int argc, const char* argv[]){
  IloEnv env;
  try {
    if (argc != 6 && argc != 4) {
      env.out() << "usage: " << argv[0] << " <type> <file> <timeLimit> [-seed <seed>]" << std::endl;
      throw ArgError();
    }

    const char* filename = argv[2];
    IloInt timeLimit = atoi(argv[3]);
    IloInt seed = 12345;
    if(argc==6) seed = atoi(argv[5]);

    std::ifstream file(filename);
    if (!file){
      env.out() << "The file \"" << filename << "\" does not exists or cannot be read" << std::endl;
      throw FileError();
    }


    IloIntervalVarArray tasks(env);
    IloIntervalVarArray ends(env);
    IloIntervalSequenceVarArray disjuncts(env);


    Instance data;
    IloModel model(env);

    
    
    if(!strncmp(argv[1], "osp", 3)) {
      data.osp_readData(filename);
      osp_setup(env, model, data, tasks, disjuncts);
    } else if(!strncmp(argv[1], "jsp", 3)) {
      data.jsp_readData(filename);
      jsp_setup(env, model, data, tasks, disjuncts);
    } else if(!strncmp(argv[1], "jla", 3)) {
      data.jla_readData(filename);
      jsp_setup(env, model, data, tasks, disjuncts);
    } else if(!strncmp(argv[1], "jet", 3)) {
      data.jet_readData(filename);
      jet_setup(env, model, data, tasks, ends, disjuncts);
    } else if(!strncmp(argv[1], "dyn", 3)) {
      data.dyn_readData(filename, 100);
      dyn_setup(env, model, data, tasks, ends, disjuncts);
    } else if(!strncmp(argv[1], "sds", 3)) {
      data.sds_readData(filename);
      sds_setup(env, model, data, tasks, disjuncts);
    } else if(!strncmp(argv[1], "jtl", 3)) {
      data.jtl_readData(filename);
      jtl_setup(env, model, data, tasks, disjuncts);
    } else if(!strncmp(argv[1], "now", 3)) {
      data.now_readData(filename);
      now_setup(env, model, data, tasks, disjuncts);
    } else {
      env.out() << "Unknown type \"" << argv[1] << "\"" << std::endl;
      throw ArgError();
    }






    IloCP cp(model);

    cp.setParameter(IloCP::Workers, 1);
    cp.setParameter(IloCP::TimeLimit, timeLimit);
    cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet);
    //cp.setParameter(IloCP::LogVerbosity, IloCP::Terse);
    cp.setParameter(IloCP::RandomSeed, seed);
    cp.setParameter(IloCP::NoOverlapInferenceLevel, IloCP::Extended);


    cp.out()   << std::left << std::setw(20) << "c Instance: " << std::right << std::setw(20) << filename << std::endl;


    bool ok = false;

    if(argv[1][3] == 's') {
      cp.out()   << "c Sequence search" << std::endl;
      IloSearchPhaseArray phaseArray(env);
      phaseArray.add(IloSearchPhase(env, tasks));
      phaseArray.add(IloSearchPhase(env, disjuncts));
      ok = cp.solve(phaseArray);
    } else {
      cp.out()   << "c Default search" << std::endl;
      ok = cp.solve();
    }



    double norm_obj = 0;
    if(ends.getSize() == data.nJobs()) {
      for(int i=0; i<data.nJobs(); ++i) {
	double incr = data.getRealCost(i, cp.getValue(IloEndOf(ends[i])));
	//std::cout << incr << std::endl;
	norm_obj += incr;
      }
      norm_obj /= data.getNormalizer();
      //std::cout << norm_obj << std::endl;
    }
      

    if (ok) {
      cp.out() << std::left << std::setw(20) << "d UPPERBOUND " << std::right << std::setw(20) << cp.getObjValue() << std::endl;
      cp.out() << std::left << std::setw(20) << "d OBJECTIVE " << std::right << std::setw(20) << cp.getObjValue() << std::endl;
      cp.out() << std::left << std::setw(20) << "d NORMOBJECTIVE " << std::right << std::setw(20) << norm_obj << std::endl;
      cp.out() << std::left << std::setw(20) << "d REALTIME " << std::right << std::setw(20) << cp.getInfo(IloCP::TotalTime) << std::endl;
      cp.out() << std::left << std::setw(20) << "d RUNTIME " << std::right << std::setw(20) << cp.getInfo(IloCP::SolveTime) << std::endl;
      cp.out() << std::left << std::setw(20) << "d NODES " << std::right << std::setw(20) << cp.getInfo(IloCP::NumberOfChoicePoints) << std::endl;
      cp.out() << std::left << std::setw(20) << "d BACKTRACKS " << std::right << std::setw(20) << cp.getInfo(IloCP::NumberOfFails) << std::endl;
      cp.out() << std::left << std::setw(20) << "d FAILS " << std::right << std::setw(20) << cp.getInfo(IloCP::NumberOfFails) << std::endl;
      cp.out() << std::left << std::setw(20) << "d NODES/s " << std::right << std::setw(20) << (cp.getInfo(IloCP::SolveTime)>0 ? cp.getInfo(IloCP::NumberOfChoicePoints)/cp.getInfo(IloCP::SolveTime) : 0) << std::endl;
      cp.out() << std::left << std::setw(20) << "d BACKTRACKS/s " << std::right << std::setw(20) << (cp.getInfo(IloCP::SolveTime)>0 ? cp.getInfo(IloCP::NumberOfFails)/cp.getInfo(IloCP::SolveTime) : 0) << std::endl; 
      cp.out() << std::left << std::setw(20) << "d FAILS/s " << std::right << std::setw(20) << (cp.getInfo(IloCP::SolveTime)>0 ? cp.getInfo(IloCP::NumberOfFails)/cp.getInfo(IloCP::SolveTime) : 0) << std::endl; 
      cp.out() << std::left << std::setw(20) << "d SOLUTIONS " << std::right << std::setw(20) << cp.getInfo(IloCP::NumberOfSolutions) << std::endl;
      cp.out() << std::left << std::setw(20) << "d OPTIMAL " << std::right << std::setw(20) << (cp.getInfo(IloCP::FailStatus) == IloCP::SearchHasFailedNormally) << std::endl;
      cp.out() << "s SATISFIABLE\nv 0\n";
    } 


    for(int i=0; i<data.nTasks(); ++i) {
      std::cout << cp.getValue(IloStartOf(tasks[i])) << std::endl;
      
    }



  } catch(IloException& e){
    env.out() << " ERROR: " << e << std::endl;
  }
  env.end();
  return 0;
}
