// -*- mode: C++ -*-
#ifndef GRANDCHILD1_H
#define GRANDCHILD1_H
#include "Child1.h"

/** First child class */
class GrandChild1 : public Child1
{
public:
  /** Destructor */
  virtual ~GrandChild1() {}
  /** Instantize */
  GrandChild1();
  /** Function */
  virtual void func();
  // Dictonary and persistency
  ClassDef(GrandChild1,1);
};

//
// Below is a hack for ROOT5
//
// We make sure that a dictionary (and compiled code) of Child1 is made.
// Child1 lives entirely in a header file, so we will not AcLIC it, but
// then we miss the binary load.  Hence we put the below statements in.
// Note, this means that GrandChild1.cxx _must_ be AcLIC'ed before
// GrandChild2.h
//
#ifdef __CINT__
#pragma link C++ class Child1+;
#endif

#endif 
//
// EOF
//
