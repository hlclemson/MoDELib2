/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AnalyticalSolidSolutionWhiteNoise_cpp
#define model_AnalyticalSolidSolutionWhiteNoise_cpp

#include <AnalyticalSolidSolutionWhiteNoise.h>
//#include <complex.h>

namespace model
{

AnalyticalSolidSolutionWhiteNoise::AnalyticalSolidSolutionWhiteNoise(const std::string& type,
                                                                     const std::string& tag,
                                                                     const int& seed,
                                                                     const GridSizeType& gridSize,
                                                                     const GridSpacingType& gridSpacing,
                                                                     const Eigen::Matrix<double,2,2>& latticeBasis,
                                                                     //const double& a_in,
                                                                     const double& a_Cai_in,
                                                                     const double& dislocLength,
                                                                     const double& MSSS) :
  /* init */ GlidePlaneNoiseBase<2>(type,"AnalyticalSolidSolutionWhiteNoise"+tag,seed,gridSize,gridSpacing,latticeBasis)
  /* init */,a_cai(a_Cai_in)
  /* init */,uncorrKcoeff(MSSS*gridSpacing(0)/dislocLength)
{
}

std::array<AnalyticalSolidSolutionWhiteNoise::COMPLEX,2> AnalyticalSolidSolutionWhiteNoise::kCorrelations(const Eigen::Matrix<double,3,1>& kv,const Eigen::Matrix<int,3,1>& kvID) const
{
  // white noise implementation
  const double foldingFactor((kvID(2)==0 || kvID(2)==NZ/2)? 2.0 : 1.0);
  //std::array<AnalyticalSolidSolutionWhiteNoise::COMPLEX,2> temp { uncorrKcoeff*foldingFactor, uncorrKcoeff*foldingFactor };
  std::array<AnalyticalSolidSolutionWhiteNoise::COMPLEX,2> temp { uncorrKcoeff, uncorrKcoeff };
  if(a_cai>0.0)
  {
    const double wkc2(this->Wk_Cai_squared(kv(0),kv(1),kv(2), a_cai)); // using the square because this is before the square root
    temp[0]*=wkc2;
    temp[1]*=wkc2;
  }
  return temp; // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
}

}
#endif

