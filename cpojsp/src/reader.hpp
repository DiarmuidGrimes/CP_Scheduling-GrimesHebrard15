 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//_DEBUGSCHED = 1

#ifdef _DEBUGSCHED
  #define DBG(fmt, args...) printf("dbg - l%d: "fmt,__LINE__,args)
#else
  #define DBG(fmt, args...)
#endif

#define INFTY 0x0fffffff



class pair_ {
public:
    
  int element;
  int rank;

  pair_(const int elt, const int rk) : element(elt), rank(rk) {}
};

struct Term {
  int X;
  int Y;
  int k;
};
class Instance {

public:

  //dtp clauses (atoms are terms like x-y <= k)
  int dtp_nodes;
  std::vector< std::vector < Term > > dtp_clauses;
  std::vector< int > dtp_weights;

  std::vector< int >              release_date;

private:

  std::vector< std::vector<int> > tasks_in_job;
  std::vector< std::vector<int> > jobs_of_task;
  std::vector< std::vector<int> > tasks_in_machine;
  std::vector< std::vector<int> > machines_of_task;
  std::vector< int >              duration;
  //  std::vector< int >              release_date;
  std::vector< int >              due_date;

  std::vector< std::vector< pair_ > > task_rank_in_machine;
  std::vector< std::vector< pair_ > > task_rank_in_job;

  int***  setup_time;
  int**   time_lag[2];
  int*    jsp_duedate;
  int*    jsp_latecost;
  int*    jsp_earlycost;
  double* jsp_floatcost;

  int max_makespan;


  int addTask(const int dur, const int job, const int machine) ;
  void addTaskToJob(const unsigned int index, const unsigned int j) ;
  void addTaskToMachine(const unsigned int index, const unsigned int j) ;



  int get_machine_rank(const int ti, const int mj) const {
    int rk=-1, mc;
    for(unsigned int i=0; i<task_rank_in_machine[ti].size(); ++i) {
      mc = task_rank_in_machine[ti][i].element;
      if(mc == mj) {
	rk = task_rank_in_machine[ti][i].rank;
	break;
      }
    }
    return rk;
  }

  int get_job_rank(const int ti, const int jj) const {
    int rk=-1, jb;
    for(unsigned int i=0; i<task_rank_in_job[ti].size(); ++i) {
      jb = task_rank_in_job[ti][i].element;
      if(jb == jj) {
	rk = task_rank_in_job[ti][i].rank;
	break;
      }
    }
    return rk;
  }


public:

  Instance();
  virtual ~Instance();


  void osp_readData( const char* filename );
  void sds_readData( const char* filename );
  void jtl_readData( const char* filename );
  void now_readData( const char* filename );
  void jla_readData( const char* filename );
  void tsp_readData( const char* filename );
  void fsp_readData( const char* filename );
  void jsp_readData( const char* filename );
  void jet_readData( const char* filename );
  void dtp_readData( const char* filename );
  void dyn_readData( const char* filename, const int p );


  std::ostream& print(std::ostream& os);

  int  nJobs() const {return tasks_in_job.size();}
  int  nJobs(const int i) const {return jobs_of_task[i].size();}

  int  nMachines() const {return tasks_in_machine.size();}
  int  nMachines(const int i) const {return machines_of_task[i].size();}

  int  nTasks() const {return duration.size();}
  int  nTasksInJob(const int j) const {return tasks_in_job[j].size();}
  int  nTasksInMachine(const int j) const {return tasks_in_machine[j].size();}

  int  getJobTask(const int i, const int j) const {return tasks_in_job[i][j];}
  int  getMachineTask(const int i, const int j) const {return tasks_in_machine[i][j];}
  int  getLastTaskofJob(const int i) const {return tasks_in_job[i][nTasksInJob(i)-1];}


  int  getJob(const int i, const int j=0) const {return jobs_of_task[i][j];}
  int  getMachine(const int i, const int j=0) const {return machines_of_task[i][j];}

  int  getDuration(const int i) const {return duration[i];}
  int  getReleaseDate(const int i) const {return release_date[i];}
  void setReleaseDate(const int i, const int d) {release_date[i] = d;}
  int  getDueDate(const int i) const {return due_date[i];}

  bool hasSetupTime() const {return setup_time != NULL;}
  int  getSetupTime(const int k, const int i, const int j) const;

  bool hasTimeLag() const {return time_lag[0] != NULL;}
  int  getMinLag(const int i, const int j) const {return time_lag[0][i][j];}
  int  getMaxLag(const int i, const int j) const {return time_lag[1][i][j];}

  bool hasJobDueDate() const {return jsp_duedate != NULL;}
  int  getJobDueDate(const int i) const {return jsp_duedate[i];}

  bool hasEarlyCost() const {return jsp_earlycost != NULL; }
  bool hasLateCost() const {return jsp_latecost != NULL; }
  int  getJobEarlyCost(const int i) const {return jsp_earlycost[i];}
  int  getJobLateCost(const int i) const {return jsp_latecost[i];}

  bool hasFloatCost() const {return jsp_floatcost != NULL; }
  double getJobFloatCost(const int i) const {return jsp_floatcost[i];}


  int getRankInJob(const int i, const int j=-1) const {return get_job_rank(i,(j==-1?getJob(i,0):j));}
  int getHeadInJob(const int i, const int j=-1) const {
    int rj = (j==-1?getJob(i,0):j);
    int rk = get_job_rank(i,rj);
    int head = 0;
    for(int k=0; k<rk; ++k) {
      head += getDuration(getJobTask(rj,k));
    }
    return head;
  }

  int getMakespanLowerBound();
  int getMakespanUpperBound(const int it=-1);

  int getEarlinessTardinessLowerBound(const int);
  int getEarlinessTardinessUpperBound(const int);

  void get_placement(const int i, const int j, std::vector<int>& intervals);
  void get_placement2(const int i, const int j);

  int nDisjuncts() const;
  int nPrecedences() const;

  // returns the real earliness/tardiness cost of finishing job i at time t
  double getRealCost(const int i, const int t) const;
  double getNormalizer() const;



  std::ostream& printStats(std::ostream& os);
};


double Instance::getRealCost(const int i, const int t) const {
  double cost = 0;
  if(hasEarlyCost() || hasLateCost()) {
    int dd = getJobDueDate(i);
    if(dd > t && hasLateCost()) {
      cost = (double)(dd - t) * getJobFloatCost(i);
    } else if(dd < t && hasEarlyCost()) {
      cost = (double)(t - dd) * getJobFloatCost(i);
    }
  }
  return cost;
}


Instance::Instance() {

  DBG("Build instance %s\n", params.data_file);

  dtp_nodes = 0;
  
  setup_time    = NULL;
  time_lag[0]   = NULL;
  time_lag[1]   = NULL;
  jsp_duedate   = NULL;
  jsp_latecost  = NULL;
  jsp_earlycost = NULL;
  jsp_floatcost = NULL;

  max_makespan  = INFTY;
}

Instance::~Instance() {
  int i, j;

  if(hasSetupTime()) {
    for(i=0; i<nMachines(); ++i) {
      for(j=0; j<nJobs(); ++j) {
	delete [] setup_time[i][j];
      }
      delete [] setup_time[i];
    }
    delete [] setup_time;
  }

  if(hasTimeLag()) {
    for(i=0; i<nJobs(); ++i) {
      delete [] time_lag[0][i];
      delete [] time_lag[1][i];
    }
    delete [] time_lag[0];
    delete [] time_lag[1];
  }

  if(hasJobDueDate()) {
    delete [] jsp_duedate;
  }

  if(hasLateCost()) {
    delete [] jsp_latecost;
  }

  if(hasEarlyCost()) {
    delete [] jsp_earlycost;
  }
}


int straight_compar(const void *x, const void *y) {
  int a = *(int*)x;
  int b = *(int*)y;
  return(a==b ? 0 : (a<b ? -1 : 1));
}

void Instance::get_placement(const int x, const int y, std::vector<int>& intervals) 
{  
  // we compute the allowed/forbidden start times of y relative to x
  // (like 'y in [x+30, x+70])

  // stores the start of each forbidden interval
  std::vector<int> forbidden_starts;
  // stores the end of each forbidden interval
  std::vector<int> forbidden_ends;

  // for each machine M the starting time of the task processed on M 
  // for job x
  int *timetable_x = new int[nMachines()];
  int *timetable_y = new int[nMachines()];

  // for each machine M the duration of the task processed on M
  // for job y
  int *duration_x = new int[nMachines()];
  int *duration_y = new int[nMachines()];


  int dur = 0;
  for(int k=0; k<nTasksInJob(x); ++k) {
    int t = getJobTask(x,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_x[m] = getDuration(t);
      timetable_x[m] = dur;
    }
    dur += getDuration(t);
  }
  dur = 0;
  for(int k=0; k<nTasksInJob(y); ++k) {
    int t = getJobTask(y,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_y[m] = getDuration(t);
      timetable_y[m] = dur;
    }
    dur += getDuration(t);
  }



  // compute all the forbidden intervals
  for(int k=0; k<nMachines(); ++k) {
    // forbidden interval for machine k
    // y cannot start in [st_x-duration_y+1, st_x+duration_x-1]
    // (since x is using machine k at this time)
    forbidden_starts.push_back(timetable_x[k]-timetable_y[k]-duration_y[k]+1);
    forbidden_ends.push_back(timetable_x[k]-timetable_y[k]+duration_x[k]);
  }


  // Now the cool part, we want to compute the 'real' intervals, since
  // some that we just computed can be merged if they overlap

  // we sort the intervals ends
  qsort(&(forbidden_starts[0]), forbidden_starts.size(), sizeof(int), straight_compar);
  qsort(&(forbidden_ends[0]), forbidden_ends.size(), sizeof(int), straight_compar);

  unsigned int i=0, j=0;
  int current = 0;    
  
  // now we can go forward from the earliest interval to latest
  while(i<forbidden_starts.size() && j<forbidden_ends.size()) {
    if(forbidden_starts[i]<forbidden_ends[j]) {
      ++i;

      // when we see the start of an interval, the current number
      // of machine that forbids this start time increases
      if(current++==0) {
	// if it goes from 0 to 1, it is the start of a forbidden interval
	intervals.push_back(forbidden_starts[i-1]-1);
      }
    } else if(forbidden_starts[i]>forbidden_ends[j]) {
      ++j;

      // when we see the end of an interval, the current number
      // of machine that forbids this start time decreases
      if(--current==0) {
	// if it goes from 1 to 0, it is the end of a forbidden interval
	intervals.push_back(forbidden_ends[j-1]);
      }
    } else {
      ++i;
      ++j;
    }
  }

  intervals.push_back(forbidden_ends.back());

//   int num__intervals = intervals.size/2+1;
//   int *max__intervals = new int[num__intervals];
//   int *min__intervals = new int[num__intervals];
//   min__intervals[0] = -NOVAL;
//   max__intervals[num__intervals-1] = NOVAL;
//   int k=0;
//   for(int i=0; i<num__intervals-1; ++i) {
//     max__intervals[i] = intervals[k++];
//     min__intervals[i+1] = intervals[k++];
//   }

}

void Instance::get_placement2(const int x, const int y) {
  // compute x's and y's lengths.
  int x_length=0, y_length=0;

  for(int k=0; k<nTasksInJob(x); ++k) {
    x_length += getDuration(getJobTask(x,k));
  }

  for(int k=0; k<nTasksInJob(y); ++k) {
    y_length += getDuration(getJobTask(y,k));
  }



  int *timetable_x_start = new int[nMachines()];
  int *timetable_x_end = new int[nMachines()];

  int *timetable_y_start = new int[nMachines()];
  int *timetable_y_end = new int[nMachines()];

  int *duration_x = new int[nMachines()];
  int *duration_y = new int[nMachines()];

  // the slack is the maximum right shift of y staying right of the given machine
  int *slack = new int[nMachines()];
  // the delay is the minimum right shift of y to get past the given machine
  int *delay = new int[nMachines()];

  //BitSet *timetable_y = new BitSet[nMachines()];

  
  int dur = 0;
  for(int k=0; k<nTasksInJob(x); ++k) {
    int t = getJobTask(x,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_x[m] = getDuration(t);
      timetable_x_start[m] = dur;
      timetable_x_end[m] = dur+getDuration(t);
    }
    dur += getDuration(t);
  }

  dur = -y_length;
  for(int k=0; k<nTasksInJob(y); ++k) {
    int t = getJobTask(y,k);
    for(int i=0; i<nMachines(t); ++i) {
      int m = getMachine(t,i);
      duration_y[m] = getDuration(t);
      timetable_y_start[m] = dur;
      timetable_y_end[m] = dur+getDuration(t);
    }
    dur += getDuration(t);
  }

  while(true) {
    int lb_shift=INFTY;
    int ub_shift=0;

    for(int i=0; i<nMachines(); ++i) {
      slack[i] = timetable_x_start[i]-timetable_y_end[i];
      delay[i] = timetable_x_end[i]-timetable_y_start[i];
      if(slack[i]<0) slack[i] = INFTY;

      if(slack[i]<lb_shift) lb_shift=slack[i];
      if(delay[i]>ub_shift) ub_shift=delay[i];
    }
  
    std::cout << "job" << x << ": " << x_length << " ";
    for(int k=0; k<nTasksInJob(x); ++k) {
      int t = getJobTask(x,k);
      for(int i=0; i<nMachines(t); ++i) {
	int m = getMachine(t,i);
	std::cout << "[" << std::setw(4) << timetable_x_start[m] << ":"
		  << std::setw(2) << getDuration(t) 
		  << " " << m
		  << "]";
      }
    }
    std::cout << std::endl;
    
    std::cout << "job" << y << ": " << y_length << " ";
    for(int k=0; k<nTasksInJob(y); ++k) {
      int t = getJobTask(y,k);
      for(int i=0; i<nMachines(t); ++i) {
	int m = getMachine(t,i);
	std::cout << "[" << std::setw(4) << timetable_y_start[m] << ":"
		  << std::setw(2) << getDuration(t) 
		  << " " << m
		  << "]";
      }
    }
    std::cout << std::endl;

    std::cout << "job" << y << ": " << y_length << " ";
    for(int k=0; k<nTasksInJob(y); ++k) {
      int t = getJobTask(y,k);
      for(int i=0; i<nMachines(t); ++i) {
	int m = getMachine(t,i);
	std::cout << "[" << std::setw(4) << slack[m] << "," << std::setw(4) << delay[m] << "]";
      }
    }
    std::cout << std::endl;

    //std::cout << min_forced_shift << " " << max_accepted_shift << " " << max_jump_shift << std::endl;
  

    
    // the shift should be, for each machine, either less than the slack, or more than the delay
    // therefore, we compute a forbidden interval, between the minimum slack, and the maximum delay
    if(lb_shift>=0) {
      std::cout << lb_shift << "]";
    }
    if(ub_shift>=0) {
      std::cout << "[" << ub_shift ;
    }

    std::cout << std::endl;

    
    for(int i=0; i<nMachines(); ++i) {
      timetable_y_start[i] += ub_shift;
      timetable_y_end[i] += ub_shift;
    }
  }

  exit(1);
}

int Instance::addTask(const int dur, const int job, const int machine) {
  int index = duration.size();
  if(job >= 0) addTaskToJob(index, job);
  if(machine >= 0) addTaskToMachine(index, machine);
  duration.push_back(dur);
  due_date.push_back(INFTY);
  release_date.push_back(0);
  
  return index;
}

void Instance::addTaskToJob(const unsigned int index, const unsigned int j) {
  if(tasks_in_job.size() <= j) tasks_in_job.resize(j+1);
  if(jobs_of_task.size() <= index) jobs_of_task.resize(index+1);
  if(task_rank_in_job.size() <= index) task_rank_in_job.resize(index+1);
  tasks_in_job[j].push_back(index);
  jobs_of_task[index].push_back(j);
  pair_ x(j, tasks_in_job[j].size()-1);
  task_rank_in_job[index].push_back(x);
}

void Instance::addTaskToMachine(const unsigned int index, const unsigned int j) {
  if(tasks_in_machine.size() <= j) tasks_in_machine.resize(j+1);
  if(machines_of_task.size() <= index) machines_of_task.resize(index+1);
  if(task_rank_in_machine.size() <= index) task_rank_in_machine.resize(index+1);
  tasks_in_machine[j].push_back(index);
  machines_of_task[index].push_back(j);
  pair_ x(j, tasks_in_machine[j].size()-1);
  task_rank_in_machine[index].push_back(x);
}

int Instance::getSetupTime(const int k, const int i, const int j) const {
  // get the rank of task i in machine k
  int ri = 0;
  for(unsigned int x=0; x<task_rank_in_machine[i].size(); ++x)
    if(task_rank_in_machine[i][x].element == k) {
      ri = task_rank_in_machine[i][x].rank;
      break;
    }

  int rj = 0;
  for(unsigned int x=0; x<task_rank_in_machine[j].size(); ++x)
    if(task_rank_in_machine[j][x].element == k) {
      rj = task_rank_in_machine[j][x].rank;
      break;
    }
  
  return setup_time[k][ri][rj];
}

std::ostream& Instance::print(std::ostream& os) {
  os << " c " << (nJobs()) << " jobs, " 
     << nMachines() << " machines ("
     << nTasks() << " tasks)" << std::endl;
  for(int i=0; i<nJobs(); ++i) {
    if(nTasksInJob(i) > 1) {
      os << " c ";
      for(int j=1; j<nTasksInJob(i); ++j)
	os << "  t" << tasks_in_job[i][j-1] << "+" << (duration[tasks_in_job[i][j-1]]) 
	   << " <= t" << tasks_in_job[i][j];
      if(hasJobDueDate()) {
	os << " dd=" << getJobDueDate(i);
      }
      os << std::endl;
    }
  }
  for(int i=0; i<nMachines(); ++i) {
    if(tasks_in_machine[i].size() > 0) {
      os << " c machine" << i << ": t" << tasks_in_machine[i][0];
      for(unsigned int j=1; j<tasks_in_machine[i].size(); ++j)
	os << ", t" << tasks_in_machine[i][j];
      os << std::endl;
    }
  }

//   for(int i=0; i<task_rank_in_machine.size(); ++i) {
//     std::cout << "task_" << i << " is"; 
//     for(int j=0; j<task_rank_in_machine[i].size(); ++j) {
//       pair_ p = task_rank_in_machine[i][j];
//       std::cout << " the " << p.rank << "/" << getMachine(i,j) << "th in machine_" << p.element;
//       if(task_rank_in_machine[i][j].element != machines_of_task[i][j]) {
// 	std::cout << "INCONSISTENCY" << std::endl;
// 	exit(1);
//       }
//     }
//     std::cout << std::endl;
//   }

  return os;
}

int Instance::nDisjuncts() const {
  int n_disjuncts = 0;
  for(int i=0; i<nMachines(); ++i) {
    n_disjuncts += (nTasksInMachine(i) * (nTasksInMachine(i)-1))/2;
  }
  return n_disjuncts;
}

int Instance::nPrecedences() const {
  int n_precedences = 0;
  for(int i=0; i<nJobs(); ++i) {
    n_precedences += (nTasksInJob(i)-1);
  }
  if(hasTimeLag()) {
    //n_precedences *= 2;
    for(int i=0; i<nJobs(); ++i) 
      for(int j=1; j<nTasksInJob(i); ++j) 
	if(getMaxLag(i,j-1) >= 0) 
	  ++n_precedences;
  }
  return n_precedences;
}

double Instance::getNormalizer() const {
  double normalizer=0.0, cost;
  int i, j, job_dur;
  for(i=0; i<nJobs(); ++i) {
    job_dur = 0;
    for(j=0; j<nTasksInJob(i); ++j) job_dur += getDuration(getJobTask(i,j));
    if(hasFloatCost())
      cost = getJobFloatCost(i);
    else if(hasEarlyCost() || hasLateCost()) {
      cost = 0.0;
      if(hasEarlyCost()) cost += getJobEarlyCost(i);
      if(hasLateCost()) cost += getJobLateCost(i);
    } else cost = 1.0;
    normalizer += ((double)job_dur * cost);
  }
  return normalizer;
}

std::ostream& Instance::printStats(std::ostream& os) {
  os << " c +===============[ instance ]================+" << std::endl
     << " d " << std::left << std::setw(25)  << " NUMTASKS "      << std::right << std::setw(19) << nTasks() << std::endl
     << " d " << std::left << std::setw(25)  << " NUMJOBS "       << std::right << std::setw(19) << nJobs() << std::endl
     << " d " << std::left << std::setw(25)  << " NUMMACHINES "   << std::right << std::setw(19) << nMachines() << std::endl
     << " d " << std::left << std::setw(25)  << " NUMDISJUNCTS "  << std::right << std::setw(19) << nDisjuncts() << std::endl
     << " d " << std::left << std::setw(25)  << " NUMPRECEDENCES "<< std::right << std::setw(19) << nPrecedences() << std::endl
//     << " d " << std::left << std::setw(25)  << "LBMAKESPAN "    << std::right << std::setw(19) << lb_C_max << std::endl
//     << " d " << std::left << std::setw(25)  << "UBMAKESPAN "    << std::right << std::setw(19) << ub_C_max << std::endl
     << " c +===============[ instance ]================+" << std::endl;
  return os;
}

int Instance::getMakespanLowerBound() {
  int mkp = 0, length;
  for(int i=0; i<nJobs(); ++i) {
    length = 0;
    for(int j=0; j<nTasksInJob(i); ++j)
      length += getDuration(getJobTask(i,j));
    if(mkp < length) mkp = length;
  }
  for(int i=0; i<nMachines(); ++i) {
    length = 0;
    for(int j=0; j<nTasksInMachine(i); ++j)
      length += getDuration(getMachineTask(i,j));
    if(mkp < length) mkp = length;
  }

  DBG("Get instance's makespan lb (%d)\n", mkp);

  return mkp;
}

int Instance::getMakespanUpperBound(const int iterations) {


  if(max_makespan < INFTY) return max_makespan;

  int best_makespan = INFTY;
  if(!hasTimeLag()) {
    int *current_job = new int[nJobs()];
    int *current_job_bound = new int[nJobs()];
    int *current_machine_bound = new int[nMachines()];
    int *current_machine_task = new int[nMachines()];
    int *ranks = new int[nTasks()];
    int *random_jobs = new int[nTasks()+1];
    int *m = new int[nMachines()];
    int k=0, i, j, t, n=0;

    std::fill(current_job, current_job+nJobs(), 0);
    while(k<nTasks()) {
      ranks[n++] = k;
      for(i=0; i<nJobs(); ++i)
	if(current_job[i]<nTasksInJob(i)) {
	  random_jobs[k++] = i;
	  ++current_job[i];
	  //std::cout << " " << i ;
	}
    }
    //std::cout << std::endl;
    ranks[n] = nTasks();

    int iter = iterations;
  
    while(iter--) {
      std::fill(current_job, current_job+nJobs(), 0);
      std::fill(current_job_bound, current_job_bound+nJobs(), 0);
      std::fill(current_machine_bound, current_machine_bound+nMachines(), 0);
      std::fill(current_machine_task, current_machine_task+nMachines(), -1);
  
      int makespan = 0;

      /*
      if(iter<iterations-1)
	for(i=0; i<n; ++i) {
	  for(t=ranks[i]; t<ranks[i+1]; ++t) {
	    j = randint(ranks[i+1]-t);
	    k = random_jobs[t];
	    random_jobs[t] = random_jobs[ranks[i]+j];
	    random_jobs[ranks[i]+j] = k;
	  }
	}
      */
      //      for(i=0; i<nTasks(); ++i) {
      //	j = randint(nTasks()-i);
      // 	k = random_jobs[i];
      // 	random_jobs[i] = random_jobs[i+j];
      // 	random_jobs[i+j] = k;
      //       }
    
      for(i=0; i<nTasks(); ++i) {

// 	for(int jj=0; jj<nJobs(); ++jj) {
// 	  std::cout << " " << std::setw(4) << current_job_bound[jj];
// 	}
// 	std::cout << std::endl;

// 	for(int jj=0; jj<nMachines(); ++jj) {
// 	  std::cout << " " << std::setw(4) << current_machine_bound[jj];
// 	}
// 	std::cout << std::endl;


	// pick the next job
	j = random_jobs[i];


	// find out which task is that
	t = getJobTask(j, current_job[j]);

	DBG("pick task t%d\n", t);

	// find out which machine is that
	for(k=0; k<nMachines(t); ++k) {
	  m[k] = getMachine(t,k);
	  DBG("  -> uses machine%d\n", m[k]);
	}
      
// 	std::cout << "pick task " << t << "(job=" << j 
// 		  << ", mach=" << m[0] << ")" << std::endl;


	// find the current timepoint for this job
	int j_mkp = current_job_bound[j];

	// find the current timepoint for this machine
	int m_mkp = current_machine_bound[m[0]];
	DBG("m%d = %d\n", m[0], current_machine_bound[m[0]]);
	for(k=1; k<nMachines(t); ++k)
	  if(m_mkp < current_machine_bound[m[k]]) {
	    m_mkp = current_machine_bound[m[k]];
	    DBG("m%d = %d\n", m[k], current_machine_bound[m[k]]);
	  }

	// get the start time for this task
	int timepoint = (j_mkp < m_mkp ? m_mkp : j_mkp);

	//	std::cout << "earliest start time = " << timepoint << std::endl;

	// check its release date
	if(getReleaseDate(t) > timepoint) timepoint = getReleaseDate(t);

	DBG("timepoint = %d\n", timepoint);

	// add setup time, if any
	if(hasSetupTime()) {
	  int setup = 0;
	  int setup_mk;
	  for(k=0; k<nMachines(t); ++k) {
	    if(current_machine_task[m[k]] >= 0) {
	      setup_mk = getSetupTime(m[k], current_machine_task[m[k]], t);
	      if(setup < setup_mk) setup = setup_mk;
	      DBG("setup = %d\n", setup_mk);
	    }
	  }
	  timepoint += setup;
	}

	// get the finish time for this task
	timepoint += getDuration(t);

	// update machin and job bounds
	for(k=0; k<nMachines(t); ++k) {
	  current_machine_bound[m[k]] = timepoint;
	  current_machine_task[m[k]] = t;
	}
	//current_machine_bound[m] = timepoint;
	current_job_bound[j] = timepoint;

	// get the final makespan
	if(makespan < timepoint) makespan = timepoint;

	++current_job[j];
      }
      if(best_makespan > makespan) best_makespan = makespan;

      //exit(1);
      //std::cout << "\t" << makespan << " " << best_makespan << std::endl;
    }

    delete [] current_job;
    delete [] current_job_bound;
    delete [] current_machine_bound;
    delete [] current_machine_task;
    delete [] ranks;
    delete [] random_jobs;
    delete [] m;

  } else {

    best_makespan = 0;
    for(int t=0; t<nTasks(); ++t) best_makespan += getDuration(t);

  }

  DBG("Get instance's makespan ub (%d)\n", best_makespan);

  return best_makespan;
}

int Instance::getEarlinessTardinessLowerBound(const int c_max) {
  return 0;
}
int Instance::getEarlinessTardinessUpperBound(const int c_max) {
  int i, ti, sum_ub = 0;

  if(hasEarlyCost()) 
    for(i=0; i<nJobs(); ++i) {
      ti = getLastTaskofJob(i);
      sum_ub += ((getJobDueDate(i) - (getReleaseDate(ti) + getDuration(ti)))*getJobEarlyCost(i));
    }
    
  if(hasLateCost()) 
    for(i=0; i<nJobs(); ++i) {
      ti = getLastTaskofJob(i);
      sum_ub += ((c_max - getJobDueDate(i))*getJobLateCost(i));
    }

  return sum_ub;
}


void Instance::osp_readData( const char* filename ) {

  DBG("Read (osp)%s\n", "");

  int opt, lb, i=0, j, k, nJobs, nMachines, dur, bufsize=1000;
   char *buf = new char[bufsize];
   std::ifstream infile( filename, std::ios_base::in );
	
   do {
     infile.getline( buf, bufsize, '\n' );
   } while( buf[0] == '#' );
	
   while( buf[i] != ' ' ) ++i;
   buf[i] = '\0';
   lb = atoi( buf );
	
   while( buf[i] == '\0' || buf[i] == ' ' || buf[i] == '*' ) ++i;
   j = i;
   while( buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\0' ) ++i;
   buf[i] = '\0';
   opt = atoi( &(buf[j]) );
	
   do {
     infile.get( buf[0] );
     if( buf[0] != '#' ) break;
     infile.getline( buf, bufsize, '\n' );
   } while( true );
   infile.unget();
	
   infile >> nJobs;
   infile >> nMachines;

   infile.getline( buf, bufsize, '\n' );
	
   do {
     infile.get( buf[0] );
     if( buf[0] != '#' ) break;
     infile.getline( buf, bufsize, '\n' );
   } while( true );
   infile.unget();
	
   k = 0;
   for(i=0; i<nJobs; ++i) {
     for(j=0; j<nMachines; ++j) {
       infile >> dur;
       addTask(dur, k, -1);

       addTaskToMachine(k, i);
       addTaskToMachine(k, nJobs+j);       

       ++k;
     }
   }

   delete [] buf;
}

void Instance::sds_readData( const char* filename ) {

  DBG("Read (sds)%s\n", "");

  int i, j, nJobs, nMachines, nFamilies, dur, mach;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> nMachines;
  infile >> nJobs;
  infile >> nFamilies;

  int **family_matrix = new int*[nFamilies+1];
  for(i=0; i<=nFamilies; ++i)
    family_matrix[i] = new int[nFamilies];
	
  setup_time = new int**[nMachines];
  for(i=0; i<nMachines; ++i) {
    setup_time[i] = new int*[nJobs];
    for(j=0; j<nJobs; ++j) {
      setup_time[i][j] = new int[nJobs];
    }
  }
	
  int **family = new int*[nJobs];
  for(i=0; i<nJobs; ++i) 
    family[i] = new int[nMachines];

  for(i=0; i<nJobs; ++i) {
		
    infile >> j;
    assert(j==nMachines);
		
    for(j=0; j<nMachines; ++j) {

      infile >> dur;
      infile >> mach;
      --mach;

      addTask(dur, i, mach);
      //addTaskToMachine(k++, mach);
			
      infile >> family[i][j];
      --family[i][j];
    }
  }
	
  for(i=0; i<=nFamilies; ++i)
    for(j=0; j<nFamilies; ++j)
      infile >> family_matrix[i][j];
	
  for(int k=0; k<nMachines; ++k) {
    for(i=0; i<nJobs; ++i) {
      for(j=0; j<nJobs; ++j) {
	setup_time[k][i][j] = family_matrix[1+family[i][k]][family[j][k]] ;
	//std::cout << " " << setup_time[k][i][j];
      }
      //std::cout << std::endl;
    }
    //std::cout << std::endl;
  }
  
  int k=0;
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      release_date[k++] = family_matrix[0][family[i][j]];
      //std::cout << release_date[k-1] << std::endl;
    }
  }	

}
void Instance::jtl_readData( const char* filename ) {

  DBG("Read (jtl)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines, opt;
  std::string tag;
  char c;
  std::ifstream infile( filename, std::ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {

      infile >> mach;
      infile >> dur;

      addTask(dur, i, mach);
    }
  }


  infile >> tag;
  infile >> opt;
  infile.ignore( 100, '\n' );
  
  infile.get(c);
  assert( c == 'T' );
  infile.get(c);
  assert( c == 'L' );
  infile.get(c);
  assert( c == '=' );

  
  time_lag[0] = new int*[nJobs];
  time_lag[1] = new int*[nJobs];

  for(i=0; i<nJobs; ++i) {
    time_lag[0][i] = new int[nMachines];
    time_lag[1][i] = new int[nMachines];
    for(j=0; j<nMachines; ++j) {
      infile >> time_lag[0][i][j];
      infile >> time_lag[1][i][j];
    }
  }

}
void Instance::now_readData( const char* filename ) {

  DBG("Read (now)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines;

  std::string tag;

  std::ifstream infile( filename, std::ios_base::in );

  infile >> nJobs;

  infile >> nMachines;

  for(i=0; i<nJobs; ++i) {

    for(j=0; j<nMachines; ++j) {

      infile >> mach;

      infile >> dur;

      addTask(dur, i, mach);

    }

  }

  time_lag[0] = new int*[nJobs];

  time_lag[1] = new int*[nJobs];

  for(i=0; i<nJobs; ++i) {

    time_lag[0][i] = new int[nMachines];

    time_lag[1][i] = new int[nMachines];

    for(j=0; j<nMachines; ++j) {

      time_lag[0][i][j] = 0;

      time_lag[1][i][j] = 0;

    }

  }

//   print(std::cout);
//   std::cout << std::endl;
 
//   Vector<int> intervals;

//   for(int i=0; i<nJobs; ++i) {
//     for(int j=i+1; j<nJobs; ++j) {
//       intervals.clear();
//       get_placement(i,j,intervals);
//       intervals.print(std::cout);
//     }
//     std::cout << std::endl;
//   }
//   //exit(1);
}
void Instance::jla_readData( const char* filename ) {

  DBG("Read (jla)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines;
  std::ifstream infile( filename, std::ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {

      infile >> mach;
      infile >> dur;

      addTask(dur, i, mach);
    }
  }

}
void Instance::tsp_readData( const char* filename ) {

  DBG("Read (tsp)%s\n", "");

}
void Instance::fsp_readData( const char* filename ) {

  DBG("Read (fsp)%s\n", "");

  int dur;
  std::vector<int> duration;
  double obj;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> obj;
  max_makespan = (int)obj;

  while( true ) {
    infile >> dur;
    if(!infile.good()) break;
    duration.push_back(dur);

    //std::cout << dur << std::endl;
  }

  int nJobs = duration.size()/2;

  for(int i=0; i<nJobs; ++i) {
    addTask(duration[i], i, 0);
    addTask(duration[i+nJobs], i, 1);
  }

  //this->print(std::cout);

  //exit(1);
}

void Instance::jsp_readData( const char* filename ) {

  DBG("Read (jsp)%s\n", "");

  int i, j, k, dur, mach;
  long int dump;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  int nJobs;
  int nMachines;
  infile >> nJobs;
  infile >> nMachines;
	
  infile >> dump;
  infile >> dump;
	
  infile >> dump;
  infile >> dump;
	
  infile >> tag;

  assert( tag == "Times" );

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      infile >> dur;
      addTask(dur, i, -1);
    }
  }
	
  infile >> tag;
  assert( tag == "Machines" );
  
  k = 0;
  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      infile >> mach;
      addTaskToMachine(k++, --mach);
    }
}

void Instance::jet_readData( const char* filename ) {

  DBG("Read (jet)%s\n", "");

  int i, j, dur, mach, nJobs, nMachines;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> nJobs;
  infile >> nMachines;
	
  jsp_duedate = new int[nJobs];
  jsp_earlycost = new int[nJobs];
  jsp_latecost = new int[nJobs];
	
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
			
      infile >> mach;
      infile >> dur;
      addTask(dur, i, mach);

    }
		
    infile >> jsp_duedate[i];
    infile >> jsp_earlycost[i];
    infile >> jsp_latecost[i];
  }

}

void Instance::dyn_readData( const char* filename, const int precision ) {

  DBG("Read (dyn)%s\n", "");

  int i, j, k, dur, mach, nJobs, nMachines;
  std::string tag;
  std::ifstream infile( filename, std::ios_base::in );
	
  infile >> nJobs;
  infile >> nMachines;
	
  jsp_duedate = new int[nJobs];
  jsp_earlycost = new int[nJobs];
  jsp_latecost = new int[nJobs];
  jsp_floatcost = new double[nJobs];

  for(i=0; i<nJobs; ++i) {
    k = release_date.size();

    for(j=0; j<nMachines; ++j) {
      
      infile >> mach;
      infile >> dur;
      
      if(mach != -1 && dur != -1) addTask(dur, i, mach);
    }		

    infile >> jsp_duedate[i];
    infile >> j;
    //infile >> jsp_latecost[i];
    infile >> jsp_floatcost[i];
    
    setReleaseDate(k, j);
    //jsp_earlycost[i] = jsp_latecost[i];
  }


  double min_cost = jsp_floatcost[0];
  double max_cost = jsp_floatcost[0];
  for(i=1; i<nJobs; ++i) {
    if(min_cost > jsp_floatcost[i]) min_cost = jsp_floatcost[i];
    if(min_cost < jsp_floatcost[i]) max_cost = jsp_floatcost[i];
  }
  for(i=0; i<nJobs; ++i) {
    int approx = (int)((jsp_floatcost[i]/max_cost)*precision);
    jsp_earlycost[i] = jsp_latecost[i] = approx;
    //std::cout << approx << std::endl;
  }

}

void Instance::dtp_readData( const char* filename ) {

  DBG("Read (dtp)%s\n", "");

  char comment = '#';
  std::ifstream infile( filename, std::ios_base::in );

  while(true) {
    comment = infile.get();
    if(comment == '#') {
      infile.ignore(1000, '\n');
    } else break;
  }
  infile.unget();
  
  std::string exp;
  int i, j, x, y, d;
  int step = 0;

  while(true) {
    infile >> exp;
    if(!infile.good())  break;

    if(step == 0) {
      std::vector< Term > clause;
      dtp_clauses.push_back(clause);
      dtp_weights.push_back(1); //randint(100)+1);
    }

    //

    if(step%2 == 0) {
      i = 0;
      while(exp[++i] != '-');
      x = atoi(exp.substr(1,i-1).c_str());
      j = i+2;
      while(exp[++j] != '<');
      y = atoi(exp.substr(i+2,j-i-2).c_str());
      d = atoi(exp.substr(j+2,10).c_str());

      Term t;
      t.X = x;
      t.Y = y;
      t.k = d;

      if(t.X > dtp_nodes) dtp_nodes = t.X;
      if(t.Y > dtp_nodes) dtp_nodes = t.Y;

      dtp_clauses.back().push_back(t);

      //std::cout << "t" << x << " - t" << y  << " <= " << d << std::endl << std::endl;
    }

    step = (step+1) % 4;
  }

    ++dtp_nodes;

  // for(int i=0; i<dtp_clauses.size(); ++i) {
  //   for(int j=0; j<dtp_clauses[i].size(); ++j) {
  //     std::cout << "t" << dtp_clauses[i][j].X << " - t" << dtp_clauses[i][j].Y << " <= " << dtp_clauses[i][j].k << " OR ";
  //   }
  //   std::cout << std::endl;
  // }


}

