#include "Base.h"

Base::Base()
  : TObject(), fBase(42)
{}

int Base::value() const
{
  return fBase;
}

//
// EOF
//
