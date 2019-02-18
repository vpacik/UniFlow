// -*- mode: C++ -*-
#ifndef CHILD1_H
#define CHILD1_H
#include "Base.h"

/** First child class */
class Child1 : public Base
{
public:
  /** Destructor */
  virtual ~Child1() {}
  /** Instantize */
  Child1() : Base() {}
  /** Function */
  virtual void func() {}
  // Dictonary and persistency
  ClassDef(Child1,1);
};
  
#endif 
//
// EOF
//
