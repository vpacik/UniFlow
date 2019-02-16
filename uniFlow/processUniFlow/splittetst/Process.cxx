#include "Process.h"

ClassImp(Process);

Process::Process() :
  fVal(),
  fObject(),
  fTask()
{
  printf("Process::Process()");
}

Process::~Process()
{
  printf("Process::~Process()");
}
