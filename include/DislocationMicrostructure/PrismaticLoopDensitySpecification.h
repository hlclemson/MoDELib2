/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PrismaticLoopDensitySpecification_H_
#define model_PrismaticLoopDensitySpecification_H_

#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{

struct PrismaticLoopDensitySpecification : public MicrostructureSpecificationBase
{
    double targetDensity;
    double radiusDistributionMean;
    double radiusDistributionStd;
    std::vector<int> allowedGrainIDs;
    std::vector<int> allowedSlipSystemIDs;

    PrismaticLoopDensitySpecification();
    PrismaticLoopDensitySpecification(const std::string& fileName);
};

}
#endif
