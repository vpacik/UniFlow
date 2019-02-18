// -*- mode: C++ -*- 
#ifndef BASE_H
#define BASE_H
#include <TObject.h>
/**
 * Base class 
 */
class Base : public TObject
{
public:
  /** Destructor */
  virtual ~Base() {}
  /** Function */
  virtual void func() = 0;
  /** Get value */
  virtual int value() const;
protected:
  /** Protected CTOR - cannot directly make objects */
  Base();
  int fBase;
  // Dictonary and persistency
  ClassDef(Base,1);
};

#endif
//
// EOF
//

