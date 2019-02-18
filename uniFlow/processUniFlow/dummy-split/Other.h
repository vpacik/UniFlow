// -*- mode: C++ -*- 
#ifndef OTHER_H
#define OTHER_H
#include <TObject.h>

class Base;

/** Other class */
class Other : public TObject
{
public:
  Other();
  ~Other();
  void func(Base& b);
  // Dictonary and persistency
  ClassDef(Other,1);
};

#endif
