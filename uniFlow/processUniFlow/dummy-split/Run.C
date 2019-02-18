// _Only_ interpreted - assumes Load.C has been run _before_ for ROOT6

void DoIt(Other& o, Base& b, const char* which)
{
  std::cout << which << "\t" << std::flush;
  o.func(b);
}


void Run()
{

  gSystem->Load("VP.so");

  Other o;
  Child1 c1;
  Child2 c2;
  GrandChild1 g1;
  GrandChild2 g2;

  DoIt(o, c1, "Child1");
  DoIt(o, c2, "Child2");
  DoIt(o, g1, "GrandChild1");
  DoIt(o, g2, "GrandChild2");
}

//
// EOF
//
