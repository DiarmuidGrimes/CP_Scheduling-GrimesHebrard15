
#include <ilcp/cp.h>


void jtl_setup(IloEnv& env, IloModel& model, Instance& data, 
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
        if (0 != prec.getImpl()) {
          model.add(IloEndBeforeStart(env, prec, ti));
	  model.add(IloStartBeforeEnd(env, ti, prec, -data.getMaxLag(i, j)));
	}
        prec = ti;
	tasks.add(ti);
	++t;
      }
      ends.add(IloEndOf(prec));
    }
    for (IloInt j = 0; j < data.nMachines(); j++) {
      IloIntervalSequenceVar s(env, machines[j]);
      disjuncts.add(s);
      model.add(IloNoOverlap(env, s));
    }
    //model.add(IloNoOverlap(env, machines[j]));

    IloObjective objective = IloMinimize(env,IloMax(ends));
    model.add(objective);
}

void now_setup(IloEnv& env, IloModel& model, Instance& data, 
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
          model.add(IloStartAtEnd(env, ti, prec));
        prec = ti;
	tasks.add(ti);
	++t;
      }
      ends.add(IloEndOf(prec));
    }
    for (IloInt j = 0; j < data.nMachines(); j++) {
      IloIntervalSequenceVar s(env, machines[j]);
      disjuncts.add(s);
      model.add(IloNoOverlap(env, s));
    }
    //model.add(IloNoOverlap(env, machines[j]));

    IloObjective objective = IloMinimize(env,IloMax(ends));
    model.add(objective);
}
