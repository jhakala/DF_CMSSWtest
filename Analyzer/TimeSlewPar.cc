#include <iostream>
#include <cmath>
#include <climits>
#include "TimeSlewPar.h"
//#include "isFinite.h"


TimeSlewPar::TimeSlewPar() {
}

TimeSlewPar::~TimeSlewPar() { 
}

void TimeSlewPar::getParameters(const std::vector<double> & inputCharge, const std::vector<double> & inputPedestal, double & RatioTS54, double & TimeSlew, double & Pulse, TF1 *slewFit, PedestalSub pedSubFxn_) const
{

  std::vector<double> corrCharge = {};
  pedSubFxn_.Calculate(inputCharge, inputPedestal, corrCharge);
  //RatioTS54 = corrCharge[5] / corrCharge[4];
  RatioTS54 = (inputCharge[5]-inputPedestal[5]) / (inputCharge[4]-inputPedestal[4]);
  //Pulse = corrCharge[4];
  Pulse = (inputCharge[4]-inputPedestal[4]);
  TimeSlew = -slewFit->GetX(RatioTS54);


}
