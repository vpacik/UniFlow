#ifndef TASK_H
#define TASK_H

#include "TH1.h"

class Task
{
  public:
    Task();
    ~Task();
  private:
    Int_t fNum;
    TH1* fHist;

  ClassDef(Task,1);
};
#endif
