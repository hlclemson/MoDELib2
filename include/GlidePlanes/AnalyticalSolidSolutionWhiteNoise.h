/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AnalyticalSolidSolutionWhiteNoise_h
#define model_AnalyticalSolidSolutionWhiteNoise_h

#include <cmath>
#include <random>

#include <Eigen/Dense>

#include <PolycrystallineMaterialBase.h>
#include <NoiseTraits.h>
#include <GlidePlaneNoiseBase.h>

namespace model
{

struct AnalyticalSolidSolutionWhiteNoise: public GlidePlaneNoiseBase<2>
{
  typedef typename NoiseTraits<2>::REAL_SCALAR REAL_SCALAR;
  typedef typename NoiseTraits<2>::COMPLEX COMPLEX;
  typedef typename NoiseTraits<2>::GridSizeType GridSizeType;
  typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
  typedef typename NoiseTraits<2>::NoiseType NoiseType;
  typedef typename NoiseTraits<2>::NoiseContainerType NoiseContainerType;

  const REAL_SCALAR a_cai;
  const REAL_SCALAR uncorrKcoeff;

  AnalyticalSolidSolutionWhiteNoise(const std::string& tag, const int& seed,
                                    const GridSizeType& gridSize, const GridSpacingType& gridSpacing,
                                    const Eigen::Matrix<double,2,2>& latticeBasis,
                                    const double& a_Cai_in, const double& dislocLength_in, const double& MSSS,
                                  const double&  effsroAve_in);

  std::array<COMPLEX,2> kCorrelations(const Eigen::Matrix<double,3,1>& kv,const Eigen::Matrix<int,3,1>& kvID) const override;

};

}
#endif

