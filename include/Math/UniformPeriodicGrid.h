/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UniformPeriodicGrid_h_
#define model_UniformPeriodicGrid_h_

#include <Eigen/Dense>
#include <CTM.h>


namespace model
{

template <int dim>
struct UniformPeriodicGrid
{

    typedef Eigen::Array<int   ,dim,1> ArrayDimI;
    typedef Eigen::Array<double,dim,1> ArrayDimD;

    // re-organized in the order of initialization
    const ArrayDimI gridSize_ss;
    const ArrayDimD gridSpacing_ss;
    const ArrayDimI gridSize_sf;
    const ArrayDimD gridSpacing_sf;

    // mutable to allow for dynamic selection 
    mutable ArrayDimI gridSize;
    mutable ArrayDimD gridSpacing;

    UniformPeriodicGrid(const ArrayDimI& gridSize_in_ss, const ArrayDimD& gridSpacing_in_ss,
                        const ArrayDimI& gridSize_in_sf, const ArrayDimD& gridSpacing_in_sf) :
        /* init */ gridSize_ss(gridSize_in_ss), gridSpacing_ss(gridSpacing_in_ss),
        /* init */ gridSize_sf(gridSize_in_sf), gridSpacing_sf(gridSpacing_in_sf),
        /* init */ gridSize(gridSize_in_ss), // default
        /* init */ gridSpacing(gridSpacing_in_ss) 
    {

    }

    // Dynamically choose grid size and spacing based on noise type
    void selectGridSizeAndSpace(const std::string& noiseType) const
    {
        if (noiseType == "stackingFault")
        {
            gridSize = gridSize_sf;
            gridSpacing = gridSpacing_sf;
        }
        else if (noiseType == "solidSolution")
        {
            gridSize = gridSize_ss;
            gridSpacing = gridSpacing_ss;
        }
        else
        {
            std::cout << "Error: noiseType not recognized" << std::endl;
        }
    }

    // UniformPeriodicGrid(const ArrayDimI& gridSize_in,const ArrayDimD& gridSpacing_in) :
    // /* init */ gridSize(gridSize_in)
    // /* init */,gridSpacing(gridSpacing_in)
    // {

    // }


    static std::array<ArrayDimI,CTM::pow(2,dim)> cornerIdx(const std::pair<ArrayDimI,ArrayDimI>& idx)
    {
        std::array<ArrayDimI,CTM::pow(2,dim)> temp;
        for(int p=0;p<std::pow(2,dim);++p)
        {
            for(int i=0;i<dim;++i)
            {
                const bool isSecond((p%CTM::pow(2,dim-i))/CTM::pow(2,dim-1-i));
                temp[p](i)= (isSecond? idx.second(i) : idx.first(i));
            }
        }
        return temp;
    }

    std::array<ArrayDimI,CTM::pow(2,dim)> posToCornerIdx(const ArrayDimD& localPos) const
    {
        return cornerIdx(posToIdx(localPos));
    }

    std::array<ArrayDimI,CTM::pow(2,dim)> posToPeriodicCornerIdx(const ArrayDimD& localPos) const
    {
        return cornerIdx(posToPeriodicIdx(localPos));
    }

    ArrayDimD idxToPos(const ArrayDimI& idx) const
    {
        return idx.template cast<double>()*gridSpacing;
    }

    std::array<double,CTM::pow(2,dim)> posToWeights(const ArrayDimD& localPos) const
    {
        // const auto cidx(posToCornerIdx(localPos));
        const auto cidx(cornerIdx(posToPeriodicIdx(localPos)));
        ArrayDimD gridTotalSize = gridSpacing * gridSize.template cast<double>();
        ArrayDimD PeriodicLocalPos = localPos - (localPos / gridTotalSize).floor() * gridTotalSize;
        const double vol(gridSpacing.prod());
        std::array<double,CTM::pow(2,dim)> temp;
        for(int p=0;p<CTM::pow(2,dim);++p)
        {
            temp[p]=(gridSpacing-(idxToPos(cidx[p])-PeriodicLocalPos).abs()).prod()/vol;
        }
        return temp;
    }

    std::pair<std::array<ArrayDimI,CTM::pow(2,dim)>,std::array<double,CTM::pow(2,dim)>> posToPeriodicCornerIdxAndWeights(const ArrayDimD& localPos) const
    {
        return std::make_pair(cornerIdx(posToPeriodicIdx(localPos)), posToWeights(localPos));
    }



    std::pair<ArrayDimI,ArrayDimI> posToIdx(const ArrayDimD& localPos) const
    {/*!\param[in] localPos the  position vector on the grid
          * \returns The grid indices (possibly outside the gridSize bounds) as the pair  (iLow,jLow),(iHIgh,jHigh)
          * ------------------
          * |   s2   |   s3     |
          * |         |(localPos)|
          * |____+_____ |
          * |        |          |
          * |        |          |
          * |   s0   |   s1     |
          * ------------------*
          * iLow,jLow
          *
          */
        const ArrayDimI highIndex((localPos/gridSpacing).ceil().template cast<int>());
        return std::make_pair(highIndex-1,highIndex);
    }

    std::pair<ArrayDimI,ArrayDimI> posToPeriodicIdx(const ArrayDimD& localPos) const
    {/*!\param[in] localPos the  position vector on the grid
          * \returns The grid index periodically wrapped within the gridSize bounds
          */
        return idxToPeriodicIdx(posToIdx(localPos));
    }

    ArrayDimI idxToPeriodicIdx(const ArrayDimI& gi) const
    {/*!\param[in] gi the  grid index possibly outside the gridSize bounds
          * \returns The grid index periodically wrapped within the gridSize bounds
          */
        const ArrayDimD gd(gi.template cast<double>());
        return gi-(gd/gridSize.template cast<double>()).floor().template cast<int>()*gridSize;
    }

    std::pair<ArrayDimI,ArrayDimI> idxToPeriodicIdx(const std::pair<ArrayDimI,ArrayDimI>& idx) const
    {/*!\param[in] gi the  grid index possibly outside the gridSize bounds
          * \returns The grid index periodically wrapped within the gridSize bounds
          */
        return std::make_pair(idxToPeriodicIdx(idx.first),idxToPeriodicIdx(idx.second));
    }



};
}


#endif
