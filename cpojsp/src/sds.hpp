
#include <ilcp/cp.h>


void sds_setup(IloEnv& env, IloModel& model, Instance& data, 
	       IloIntervalVarArray& tasks,
	       IloIntervalSequenceVarArray& disjuncts) {
    IloIntervalVarArray2 machines(env, data.nMachines());
    for (IloInt j = 0; j < data.nMachines(); j++)
      machines[j] = IloIntervalVarArray(env);
    IloIntExprArray ends(env);
    int t=0;
    for (IloInt i = 0; i < data.nJobs(); i++) {
      IloIntervalVar prec;
      for (IloInt j = 0; j < data.nMachines(); j++) {
        IloInt m = data.getMachine(t), d = data.getDuration(t);
        IloIntervalVar ti(env, d);
        machines[m].add(ti);
        if (0 != prec.getImpl())
          model.add(IloEndBeforeStart(env, prec, ti));
        prec = ti;
	tasks.add(ti);
	++t;
      }
      ends.add(IloEndOf(prec));
    }
    for (IloInt j = 0; j < data.nMachines(); j++) {

      int n_tasks = data.nTasksInMachine(j);
      IloTransitionDistance setupTimes(env, n_tasks);
      IloIntArray           type      (env, n_tasks);
      for(int a = 0; a<n_tasks; ++a) {
	type[a] = a;
      }
      for(int a = 0; a<n_tasks; ++a) {
	for(int b = 0; b<n_tasks; ++b) //if(a!=b) 
	  {	
	    int t1 = data.getMachineTask(j, a);
	    int t2 = data.getMachineTask(j, b);
	    setupTimes.setValue(a,b,data.getSetupTime(j, t1, t2));
	    //std::cout << " " << (data.getSetupTime(j, t1, t2)) ;
	  }
	//std::cout << std::endl;
      }
      //std::cout << std::endl;
      
      IloIntervalSequenceVar s(env, machines[j], type);
      disjuncts.add(s);
      model.add(IloNoOverlap(env, s, setupTimes, IloTrue));
    }


    for (IloInt i = 0; i < data.nTasks(); i++) {
      model.add( IloStartOf(tasks[i]) >= data.getReleaseDate(i) );
    }
    //model.add(IloNoOverlap(env, machines[j]));

    IloObjective objective = IloMinimize(env,IloMax(ends));
    model.add(objective);
}
