// -*- mode: C++ -*-
#ifndef GRANDCHILD2_H
#define GRANDCHILD2_H
#include "Child1.h"

/** First child class */
class GrandChild2 : public Child1
{
public:
  /** Destructor */
  virtual ~GrandChild2() {}
  /** Instantize */
  GrandChild2() : Child1() {}
  /** Function */
  virtual void func() { fBase /= 10; }
  // Dictonary and persistency
  ClassDef(GrandChild2,1);
};
  
#endif 
//
// EOF
//
