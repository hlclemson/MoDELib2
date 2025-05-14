/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MDStackingFaultNoise_h
#define model_MDStackingFaultNoise_h

#include <cmath>
#include <random>

#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>

#include <PolycrystallineMaterialBase.h>
#include <NoiseTraits.h>
#include <GlidePlaneNoiseBase.h>

namespace model
{
    struct MDStackingFaultNoise : public GlidePlaneNoiseBase<1>
    {
        typedef typename NoiseTraits<1>::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraits<1>::COMPLEX COMPLEX;
        typedef typename NoiseTraits<1>::GridSizeType GridSizeType;
        typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
        typedef typename NoiseTraits<1>::NoiseType NoiseType;
        typedef typename NoiseTraits<1>::NoiseContainerType NoiseContainerType;

        COMPLEX *Rk;
        const std::string correlationFile;

        MDStackingFaultNoise(const PolycrystallineMaterialBase& mat,
                             const std::string& tag,
                             const std::string& correlationFile_in,
                             const int& outputNoise,
                             const std::string& noiseFile,
                             const int& testNoiseSampling,
                             const int& seed,
                             const GridSizeType& gridSize,
                             const GridSpacingType& gridSpacing);

        std::array<COMPLEX,1> kCorrelations(const Eigen::Matrix<double, 3, 1> &kv, const Eigen::Matrix<int, 3, 1> &index) const override;
        Eigen::Matrix<double,2,2> invTransitionMatrix();
        Eigen::Matrix<double,2,2> nonOrthogonalBasisReader(const std::string& fileName_vtk) const;
        void StackingFaultCorrelationReader(const std::string &fileName_vtk, REAL_SCALAR *Rr, int NR);
        GridSizeType readVTKfileDimension(const char *fname);

        void sampleNoiseLocalInKspace(const PolycrystallineMaterialBase& mat, const int& localSeed, COMPLEX* localRk);
        void circularShift(REAL_SCALAR* rCorrelation);
        void sampleNoiseRepeatedly(const PolycrystallineMaterialBase& mat, const int& realizationNum);
        void Write_field_slice(const PolycrystallineMaterialBase& mat, const int& seed, const char *fname);
    };


}
#endif
