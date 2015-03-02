#include "HcalTimeSlew.h"
#include <cmath>

static const double tzero[3]= {23.960177, 13.307784, 9.109694};
static const double slope[3] = {-3.178648,  -1.556668, -1.075824 };
static const double tmax[3] = {16.00, 10.00, 6.25 };

double HcalTimeSlew::delay(double fC, BiasSetting bias) {
  double rawDelay=tzero[bias]+slope[bias]*log(fC);
  return (rawDelay<0)?(0):((rawDelay>tmax[bias])?(tmax[bias]):(rawDelay));			   
}

double HcalTimeSlew::delay(double fC, ParaSource source, BiasSetting bias) {

  if (source==TestStand) {
    return HcalTimeSlew::delay(fC, bias);
  }
  else if (source==Data) {
    //from john 2/20 talk: indico.cern.ch/event/375365/contribution/9/material/slides/5.pdf
    return 13.98-3.20*log(fC+32)-2.82965+10;
  }
  else if (source==MC) {
    // from Xinmei
    return 10.491-2.25495*log(fC+7.95067);
  }
  
  std::cout << "rechit energy = " <<  std::endl;
  std::cout << "What are you doing for the time slew parameterization?" << std::endl;
  
  return 0;
  
}
