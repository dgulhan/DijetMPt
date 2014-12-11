#ifndef getWeightedMean_h
#define getWeightedMean_h

#include <vector>

void getWeightedMean(std::vector<float>* inVal_p, std::vector<float>* inWeight_p, Int_t& outMult, Float_t& outMean, Float_t& outSig, Float_t& outWeight, Float_t& outWeight2)
{
  outMult = (Int_t)(inVal_p->size());
  outMean = 0;
  outSig = 0;
  outWeight = 0;
  outWeight2 = 0;

  if(outMult == 0) return;

  for(Int_t iter = 0; iter < outMult; iter++){
    outMean += inVal_p->at(iter)*inWeight_p->at(iter);
    outWeight += inWeight_p->at(iter);
    outWeight2 += inWeight_p->at(iter)*inWeight_p->at(iter);
  }

  outMean /= outWeight;

  for(Int_t iter = 0; iter < (Int_t)(inVal_p->size()); iter++){
    outSig += inWeight_p->at(iter)*(inVal_p->at(iter) - outMean)*(inVal_p->at(iter) - outMean);
  }

  outSig /= outWeight;

  return;
}

#endif
