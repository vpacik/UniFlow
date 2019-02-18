Bool_t
LoadOne(const char* clsname, const char* opt="", const char* suf="cxx")
{
  if (gROOT->GetClass(clsname) &&
      gROOT->GetClass(clsname)->HasDictionary())
    if (opt[0] == 'g' || opt[1] == 'g') 
      Info("LoadOne", "Class %s already known", clsname);

  int ret = gROOT->LoadMacro(Form("%s.%s+%s",clsname,suf,opt));
  if (ret != 0) 
    Warning("LoadOne", "Failed to load %s (%s, with +%s)",clsname,suf,opt);
  return ret == 0;
}
	
void
Load(const char* opt="")
{
  LoadOne("Base",        opt);
  LoadOne("Other",       opt);
  LoadOne("GrandChild1", opt);
  LoadOne("GrandChild2", opt, "h");
  LoadOne("Child2",      opt);
}

