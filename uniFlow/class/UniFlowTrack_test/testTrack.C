#include <vector>
#include <algorithm>
#include "AliAODTrack.h"

#include "UniFlowTrack.h"

bool sortPt(const UniFlowTrack& t1, const UniFlowTrack& t2) { return (t1.fPt < t2.fPt); }
bool sortEta(const UniFlowTrack& t1, const UniFlowTrack& t2) { return (t1.fEta < t2.fEta); }

void testTrack()
{

  AliAODTrack* aod = new AliAODTrack();
  aod->SetPt(1.2);

  // UniFlowTrack track = UniFlowTrack(aod);
  UniFlowTrack track = UniFlowTrack(aod);



  printf("%g\n", track.fPt);

  // std::vector<Double_t> vec = {1.2,3.4,5.,10.,0.2};
  // vec.push_back(3.4);
  //
  // printf("pre-sorting\n");
  // for(int i = 0; i < vec.size(); ++i) {
  //   printf("%g\n",vec[i]);
  // }
  //
  // std::sort(vec.begin(), vec.end());
  //
  // printf("post-sorting\n");
  // for(int i = 0; i < vec.size(); ++i) {
  //   printf("%g\n",vec[i]);
  // }

  // WORKS with Double_t (above)

  // std::vector<UniFlowTrack> vecT;
  // vecT.push_back(UniFlowTrack(0.2,1.0));
  // vecT.push_back(UniFlowTrack(10.2,2.3));
  // vecT.push_back(UniFlowTrack(10.,2.3));
  // vecT.push_back(UniFlowTrack(1.2,-1.));
  // vecT.push_back(UniFlowTrack(0.2,-1.1));
  //
  // printf("pre-sorting\n");
  // for(int i = 0; i < vecT.size(); ++i) {
  //   printf("%g | %g\n",vecT[i].fPt,vecT[i].fEta);
  // }
  //
  // // std::sort(vecT.begin(), vecT.end());
  // // std::sort(vecT.begin(), vecT.end(),sortPt);
  // std::sort(vecT.begin(), vecT.end(),sortEta);
  //
  // printf("post-sorting\n");
  // for(int i = 0; i < vecT.size(); ++i) {
  //   printf("%g | %g\n",vecT[i].fPt,vecT[i].fEta);
  // }


};
