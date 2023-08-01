
#include <ilcp/cp.h>
#include <iostream>

IloIntExpr TardinessCost(IloIntervalVar task, IloInt dd, IloInt weight) {
  return weight * IloMax(0, IloEndOf(task)-dd);
}

IloIntExpr EarlinessCost(IloIntervalVar task, IloInt dd, IloInt weight) {
  return weight * IloMax(0, dd-IloEndOf(task));
}

void jet_setup(IloEnv& env, IloModel& model, Instance& data, 
	       IloIntervalVarArray& tasks,
	       //IloIntVarArray& ends,
	       IloIntervalVarArray& ends,
	       IloIntervalSequenceVarArray& disjuncts) {

    IloIntervalVarArray2 machines(env, data.nMachines());
    for (IloInt j = 0; j < data.nMachines(); j++)
      machines[j] = IloIntervalVarArray(env);
    IloIntExprArray penalties(env);
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

      //ends.add(IloEndOf(prec));
      ends.add(prec);

      penalties.add(EarlinessCost(prec, data.getJobDueDate(i), data.getJobEarlyCost(i)));
      penalties.add(TardinessCost(prec, data.getJobDueDate(i), data.getJobLateCost(i)));
    }
    for (IloInt j = 0; j < data.nMachines(); j++) {
      IloIntervalSequenceVar s(env, machines[j]);
      disjuncts.add(s);
      model.add(IloNoOverlap(env, s));
    }
    //model.add(IloNoOverlap(env, machines[j]));


    IloObjective objective = IloMinimize(env,IloSum(penalties));
    model.add(objective);
}


void dyn_setup(IloEnv& env, IloModel& model, Instance& data, 
	       IloIntervalVarArray& tasks,
	       IloIntervalVarArray& ends,
	       IloIntervalSequenceVarArray& disjuncts) {

  IloIntervalVarArray2 machines(env, data.nMachines());
  for (IloInt j = 0; j < data.nMachines(); j++)
    machines[j] = IloIntervalVarArray(env);
  IloIntExprArray penalties(env);

  std::cout << "ntasks = " << data.nTasks() << ", num release dates: " << data.release_date.size() << std::endl;
  
  int t;
  for (IloInt i = 0; i < data.nJobs(); i++) {
    IloIntervalVar prec;
    for (IloInt j = 0; j < data.nTasksInJob(i); j++) {
      t = data.getJobTask(i,j);
      IloInt m = data.getMachine(t), d = data.getDuration(t);
      IloIntervalVar ti(env, d);


      if(j==0) {
	model.add(IloStartOf(ti) >= data.getReleaseDate(t));
	std::cout << data.getReleaseDate(t) << std::endl;
      }
      //std::cout << "add t" << t << " (" << d << ") in machine " << m << std::endl;

      machines[m].add(ti);
      if (0 != prec.getImpl()) {
	model.add(IloEndBeforeStart(env, prec, ti));
	//std::cout << "add prec t" << data.getJobTask(i,j-1) << " < t" << t << std::endl;
      }
      prec = ti;
      tasks.add(ti);
    }

    //std::cout <<  data.getJobEarlyCost(i) << " / " << data.getJobLateCost(i) << std::endl;
    
    //ends.add(IloEndOf(prec));
    ends.add(prec);

    //penalties.add(data.getJobLateCost(i) * IloAbs(IloEndOf(prec)-data.getJobDueDate(i)));
    penalties.add(EarlinessCost(prec, data.getJobDueDate(i), data.getJobEarlyCost(i)));
    penalties.add(TardinessCost(prec, data.getJobDueDate(i), data.getJobLateCost(i)));
  }
  for (IloInt j = 0; j < data.nMachines(); j++) {
    if(data.nTasksInMachine(j) > 0) {
      IloIntervalSequenceVar s(env, machines[j]);
      disjuncts.add(s);
      model.add(IloNoOverlap(env, s));
    }
  }
  // model.add(IloNoOverlap(env, machines[j]));

  
  IloObjective objective = IloMinimize(env,IloSum(penalties));
  model.add(objective);
}
