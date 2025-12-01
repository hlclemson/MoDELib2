/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExtendFRSourceIndividualSpecification_cpp_
#define model_ExtendFRSourceIndividualSpecification_cpp_

#include <ExtendFRSourceIndividualSpecification.h>

namespace model
{
    ExtendFRSourceIndividualSpecification::ExtendFRSourceIndividualSpecification():
    /* init */ MicrostructureSpecificationBase("ExtendFRSource","Individual")
    {
        
    }

    ExtendFRSourceIndividualSpecification::ExtendFRSourceIndividualSpecification(const std::string& fileName):
    /* init */ MicrostructureSpecificationBase("ExtendFRSource","Individual",fileName)
    /* init */,slipSystemIDs(this->parser->readArray<int>("slipSystemIDs",true))
    /* init */,exitFaceIDs(this->parser->readArray<int>("exitFaceIDs",true))
    /* init */,ExtendFRSourceCenters(this->parser->readMatrix<double>("ExtendFRSourceCenters",slipSystemIDs.size(),3,true))
    /* init */,ExtendFRSourceHeights(this->parser->readArray<double>("ExtendFRSourceHeights",true))
    /* init */,nodesPerLine(this->parser->readArray<int>("nodesPerLine",true))
    /* init */,glideSteps(this->parser->readArray<double>("glideSteps",true))
    /* init */,ExtendFRSourceLengths(this->parser->readArray<double>("ExtendFRSourceLengths",true))

    {
        
    }
}
#endif
