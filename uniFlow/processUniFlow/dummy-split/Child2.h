// -*- mode: C++ -*-
#ifndef CHILD2_H
#define CHILD2_H
#include "Base.h"

/** First child class */
class Child2 : public Base
{
public:
  /** Destructor */
  virtual ~Child2() {}
  /** Instantize */
  Child2() : Base() {}
  /** Function */
  virtual void func();
  // Dictonary and persistency
  ClassDef(Child2,1);
};

#endif 
//
// EOF
//
