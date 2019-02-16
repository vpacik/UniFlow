#include "Task.h"

Task::Task() :
  fNum{},
  fHist{}
{
  printf("Task::Task()");
}

Task::~Task()
{
  printf("Task::~Task()");
  if(fHist) delete fHist;
}
