#include "Other.h"
#include "Base.h"
#include <iostream>

Other::Other()
  : TObject()
{}


Other::~Other()
{}

void
Other::func(Base& b)
{
  b.func();
  std::cout << "Value of object: " << b.value() << std::endl;
}

//
// EOF
//
