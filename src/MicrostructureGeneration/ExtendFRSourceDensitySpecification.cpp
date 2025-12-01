/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExtendFRSourceDensitySpecification_cpp_
#define model_ExtendFRSourceDensitySpecification_cpp_

#include <ExtendFRSourceDensitySpecification.h>

namespace model
{
    ExtendFRSourceDensitySpecification::ExtendFRSourceDensitySpecification():
    /* init */ MicrostructureSpecificationBase("ExtendFRSource","Density")
    /* init */,targetDensity(0.0)
    {
        
    }

    ExtendFRSourceDensitySpecification::ExtendFRSourceDensitySpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("ExtendFRSource","Density",fileName)
    /* init */,targetDensity(this->parser->readScalar<double>("targetDensity",true))
    {
        
    }
}
#endif
