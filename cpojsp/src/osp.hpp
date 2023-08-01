
#include <ilcp/cp.h>


void osp_setup(IloEnv& env, IloModel& model, Instance& data, 
	       IloIntervalVarArray& tasks,
	       IloIntervalSequenceVarArray& disjuncts) {
    IloInt i, j;
    IloIntervalVarArray2 machines(env, data.nMachines());
    for (j = 0; j < data.nMachines(); j++)
      machines[j] = IloIntervalVarArray(env);
    IloIntExprArray ends(env);
    for (i = 0; i < data.nTasks(); i++) {
      IloIntervalVar ti(env, data.getDuration(i));
      for (j = 0; j < data.nMachines(i); j++) {
        machines[data.getMachine(i,j)].add(ti);
        ends.add(IloEndOf(ti));
      }
      tasks.add(ti);
    }
    for (j = 0; j < data.nMachines(); j++) {
      IloIntervalSequenceVar s(env, machines[j]);
      disjuncts.add(s);
      model.add(IloNoOverlap(env, s));
    }
    //model.add(IloNoOverlap(env, machines[j]));

    IloObjective objective = IloMinimize(env,IloMax(ends));
    model.add(objective);
}
