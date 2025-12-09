/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneNoise_cpp
#define model_GlidePlaneNoise_cpp

#include <algorithm>
#include <filesystem>
#include <GlidePlaneNoise.h>

namespace model
{

GlidePlaneNoise::GlidePlaneNoise(const PolycrystallineMaterialBase& mat)
{
  const auto noiseFiles(TextFileParser(mat.materialFile).readStringVector("glidePlaneNoise"));
  for(const auto& pair : noiseFiles)
  {
    const auto noiseFileString(TextFileParser::removeSpaces(pair.first));
    if(noiseFileString!="none" && noiseFileString!="None")
    {
      const std::string noiseFileName(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+noiseFileString);
      TextFileParser parser(noiseFileName);
      const std::string type(TextFileParser::removeSpaces(parser.readString("type",true)));
      const std::string  tag(TextFileParser::removeSpaces(parser.readString("tag",true)));
      const int seed(parser.readScalar<int>("seed",true));
      const Eigen::Matrix<int,1,3> gridSize(parser.readMatrix<int,1,3>("gridSize",true));
      const Eigen::Matrix<double,1,3> gridSpacing(parser.readMatrix<double,1,3>("gridSpacing_SI",true)/mat.b_SI);

      if(type=="AnalyticalSolidSolutionWhiteNoise")
      {
        //const double a(parser.readScalar<double>("spreadLstress_SI",true)/mat.b_SI);      // spreading length for stresses [AA]
        const double a_Cai(parser.readScalar<double>("a_cai_SI",true)/mat.b_SI);
        const double dislocLength(parser.readScalar<double>("dislocation_length_SI",true)/mat.b_SI);
        const double MSSS(parser.readScalar<double>("MSSS_SI",true)/std::pow(mat.mu_SI,2));
        const auto success(solidSolutionNoise().emplace(tag,new AnalyticalSolidSolutionWhiteNoise(tag,seed,gridSize,gridSpacing,Eigen::Matrix<double,2,2>::Identity(),a_Cai,dislocLength,MSSS,0.0)));

        if(!success.second)
        {
          throw std::runtime_error("Could not insert noise "+tag);
        }

      }
      if(type=="AnalyticalSolidSolutionNoise")
      {
        //const bool isWhite(parser.readScalar<int>("white",true));
        const double a(parser.readScalar<double>("spreadLstress_SI",true)/mat.b_SI);      // spreading length for stresses [AA]
        const double a_Cai(parser.readScalar<double>("a_cai_SI",true)/mat.b_SI);
        const double MSSS(parser.readScalar<double>("MSSS_SI",true)/std::pow(mat.mu_SI,2));

        //const auto success(solidSolutionNoise().emplace(tag,new AnalyticalSolidSolutionNoise(tag,seed,gridSize,gridSpacing,Eigen::Matrix<double,2,2>::Identity(),isWhite,a,a_Cai,MSSS)));
        const auto success(solidSolutionNoise().emplace(tag,new AnalyticalSolidSolutionNoise(tag,seed,gridSize,gridSpacing,Eigen::Matrix<double,2,2>::Identity(),a,a_Cai,MSSS,0.0)));

        if(!success.second)
        {
          throw std::runtime_error("Could not insert noise "+tag);
        }
      }
      if(type=="MDSolidSolutionNoise")
      {
        const double a_Cai(parser.readScalar<double>("a_cai_SI",true)/mat.b_SI);
        const double MSSS_xz(parser.readScalar<double>("MSSS_xz_SI",true)); // Added
        const double MSSS_yz(parser.readScalar<double>("MSSS_yz_SI",true)); // Added

        // relative paths for correlation vtk files (easier to make everything self contained)
        const std::string correlationFile_xz(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+TextFileParser::removeSpaces(parser.readString("correlationFile_xz",true)));
        const std::string correlationFile_yz(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+TextFileParser::removeSpaces(parser.readString("correlationFile_yz",true)));

        const auto success(solidSolutionNoise().emplace(tag,new MDSolidSolutionNoise(mat,tag,correlationFile_xz,correlationFile_yz,seed,gridSize,gridSpacing,Eigen::Matrix<double,2,2>::Identity(),MSSS_xz,MSSS_yz,a_Cai,0.0)));
        if(!success.second)
        {
          throw std::runtime_error("Could not insert noise "+tag);
        }
      }
      if(type=="MDStackingFaultNoise")
      {
        // relative paths for correlation vtk files (easier to make everything self contained)
        const std::string correlationFile_stackingFault(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+TextFileParser::removeSpaces(parser.readString("correlationFile",true)));

        const auto success(stackingFaultNoise().emplace(tag,new MDStackingFaultNoise(mat,tag,correlationFile_stackingFault,seed,gridSize,gridSpacing,Eigen::Matrix<double,2,2>::Identity(),0.0)));

        if(!success.second)
        {
          throw std::runtime_error("Could not insert noise "+tag);
        }
      }
      if(type=="MDShortRangeOrderNoise")
      {
        // relative paths for correlation vtk files (easier to make everything self contained)
        const std::string correlationFile_ShortRangeOrder(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+TextFileParser::removeSpaces(parser.readString("correlationFile",true)));
        const double effsroAve(parser.readScalar<double>("effsroAve",true)/(mat.b_SI*mat.mu_SI));
        const auto success(shortRangeOrderNoise().emplace(tag, new MDShortRangeOrderNoise(mat,tag,correlationFile_ShortRangeOrder,seed,gridSize,gridSpacing,Eigen::Matrix<double,2,2>::Identity(),effsroAve)));

        if(!success.second)
        {
          throw std::runtime_error("Could not insert noise "+tag);
        }

      }
    }


  }
  for(auto& pair : solidSolutionNoise())
  {
    pair.second->computeRealNoise();
    //pair.second->computeRealNoiseStatistics(mat);
    pair.second->computeRealNoiseStatistics();
    pair.second->write_field_slice(mat);
  }

  for(auto& pair : stackingFaultNoise())
  {
    pair.second->computeRealNoise();
    //pair.second->computeRealNoiseStatistics(mat);
    pair.second->computeRealNoiseStatistics();
    pair.second->write_field_slice(mat);
  }

  for(auto& pair : shortRangeOrderNoise())
  {
    pair.second->computeRealNoise();
    //pair.second->computeRealNoiseStatistics(mat);
    pair.second->computeRealNoiseStatistics();
    pair.second->write_field_slice(mat);
  }

}

//std::tuple<double,double,double> GlidePlaneNoise::gridInterp(const Eigen::Matrix<double,2,1>& localPos) const
std::tuple<double,double,double,double,double> GlidePlaneNoise::gridInterp(const Eigen::Matrix<double,2,1>& localPos ) const
{   // Added by Hyunsoo (hyunsol@g.clemson.edu), Xin Liu (liuxin23@mails.tsinghua.edu.cn)
  double effsolNoiseXZ(0.0);
  double effsolNoiseYZ(0.0);
  // 1. Solid solution noise
  for(const auto& noise : solidSolutionNoise())
  {
    const auto idxAndWeights(noise.second->posToPeriodicCornerIdxAndWeights(localPos));
    for(size_t p=0;p<idxAndWeights.first.size();++p)
    {
      const int storageID(noise.second->storageIndex(idxAndWeights.first[p](0),idxAndWeights.first[p](1)));
      effsolNoiseXZ+=noise.second->operator[](storageID)(0)*idxAndWeights.second[p];
      effsolNoiseYZ+=noise.second->operator[](storageID)(1)*idxAndWeights.second[p];
    }
  }
  // 2. Stacking fault noise
  double effsfNoise(0.0);
  for(const auto& noise : stackingFaultNoise())
  {
    // rotate stacking fault noise field

    // transform the basis if it is non-orthogonal
    //const auto idxAndWeights(noise.second->posToPeriodicCornerIdxAndWeights(localPos));
    const auto idxAndWeights(noise.second->posToPeriodicCornerIdxAndWeights(noise.second->invTransposeLatticeBasis*localPos));

    for(size_t p=0;p<idxAndWeights.first.size();++p)
    {
      const int storageID(noise.second->storageIndex(idxAndWeights.first[p](0),idxAndWeights.first[p](1)));
      effsfNoise+=noise.second->operator[](storageID)*idxAndWeights.second[p];
    }
  }
   // 3. SRO noise
    double effsroNoise(0.0);
    double effsroAve(0.0);
    for(const auto& noise : shortRangeOrderNoise())
    {
      // rotate stacking fault noise field

      // transform the basis if it is non-orthogonal
      const auto idxAndWeights(noise.second->posToPeriodicCornerIdxAndWeights(noise.second->invTransposeLatticeBasis*localPos));
      effsroAve = noise.second->effsroAve;
      // std::cout<<"effsroAve="<<effsroAve<<std::endl;

      for(size_t p=0;p<idxAndWeights.first.size();++p)
      {
        const int storageID(noise.second->storageIndex(idxAndWeights.first[p](0),idxAndWeights.first[p](1)));
        effsroNoise+=noise.second->operator[](storageID)*idxAndWeights.second[p];
      }
  }
  return std::make_tuple(effsolNoiseXZ,effsolNoiseYZ,effsfNoise, effsroAve, effsroNoise);
}

std::tuple<double,double,double,double,double> GlidePlaneNoise::gridVal(const Eigen::Array<int,2,1>& idx) const
{   // Added by Hyunsoo (hyunsol@g.clemson.edu, Liuxin, liuxin23@mails.tsinghua.eud.cn)

  double effsolNoiseXZ(0.0);
  double effsolNoiseYZ(0.0);
  for(const auto& noise : solidSolutionNoise())
  {
    const Eigen::Array<int,2,1> pidx(noise.second->idxToPeriodicIdx(idx));
    const int storageID(noise.second->storageIndex(pidx(0),pidx(1)));
    effsolNoiseXZ+=noise.second->operator[](storageID)(0);
    effsolNoiseYZ+=noise.second->operator[](storageID)(1);
  }

  double effsfNoise(0.0);
  for(const auto& noise : stackingFaultNoise())
  {
    const Eigen::Array<int,2,1> pidx(noise.second->idxToPeriodicIdx(idx));
    const int storageID(noise.second->storageIndex(pidx(0),pidx(1)));
    effsfNoise+=noise.second->operator[](storageID);
  }

  double effsroNoise(0.0);
  double effsroAve(0.0);
  for(const auto& noise : shortRangeOrderNoise())
  {
    const double effsroAve = noise.second->effsroAve;

    const Eigen::Array<int,2,1> pidx(noise.second->idxToPeriodicIdx(idx));
    const int storageID(noise.second->storageIndex(pidx(0),pidx(1)));
    effsroNoise+=noise.second->operator[](storageID);
  }

  return std::make_tuple(effsolNoiseXZ,effsolNoiseYZ,effsfNoise,effsroAve, effsroNoise);

}

const typename GlidePlaneNoise::SolidSolutionNoiseContainer& GlidePlaneNoise::solidSolutionNoise() const
{
  return *this;
}

typename GlidePlaneNoise::SolidSolutionNoiseContainer& GlidePlaneNoise::solidSolutionNoise()
{
  return *this;
}

const typename GlidePlaneNoise::StackingFaultNoiseContainer& GlidePlaneNoise::stackingFaultNoise() const
{
  return *this;
}

typename GlidePlaneNoise::StackingFaultNoiseContainer& GlidePlaneNoise::stackingFaultNoise()
{
  return *this;
}

const typename GlidePlaneNoise::StackingFaultNoiseContainer& GlidePlaneNoise::shortRangeOrderNoise() const
{
  return *this;
}

typename GlidePlaneNoise::StackingFaultNoiseContainer& GlidePlaneNoise::shortRangeOrderNoise()
{
  return *this;
}

}
#endif

