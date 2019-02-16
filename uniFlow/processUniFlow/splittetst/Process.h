#ifndef PROCESS_H
#define PROCESS_H

#include "TObject.h"
#include "Task.h"

class Process
{
  public:
    Process();
    ~Process();
  private:
    Double_t fVal;
    TObject fObject;
    Task fTask;

  ClassDef(Process,1);
};

#endif
