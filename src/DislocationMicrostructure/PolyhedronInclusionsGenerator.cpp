/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionsGenerator_cpp_
#define model_InclusionsGenerator_cpp_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <PolyhedronInclusionsGenerator.h>
#include <PlaneLineIntersection.h>
#include <GmshReader.h>

namespace model
{

    bool PolyhedronInclusionsGenerator::generateSingle(MicrostructureGenerator& mg,const std::map<size_t,Eigen::Vector3d>& polyNodes,const std::map<size_t,std::vector<size_t>>& faceMap, const Eigen::Matrix<double,1,dim*dim>& eTrow, const double& vrc,const int&type)
    {
        Eigen::Matrix<double,dim,dim> eT(Eigen::Map<const Eigen::Matrix<double,dim,dim>>(eTrow.data(),dim,dim).transpose());
        mg.insertInclusion(polyNodes,faceMap, eT, vrc,type);
        return true;
    }

    PolyhedronInclusionsGenerator::PolyhedronInclusionsGenerator(const PolyhedronInclusionIndividualSpecification& spec,MicrostructureGenerator& mg) :
    /* init */ allowOverlap(false)
    {
        GmshReader mshReader(spec.mshFile);
        
//        std::cout<<spec.X0<<std::endl;
//        std::cout<<spec.F<<std::endl;
//        std::cout<<spec.eigenDistortion<<std::endl;
        
        Eigen::Vector3d nodeBary(Eigen::Vector3d::Zero());
        for(const auto& node : mshReader.nodes())
        {
            nodeBary+=node.second;
        }
        nodeBary/=mshReader.nodes().size();
        
        for(int j=0;j<spec.X0.rows();++j)
        {
            std::map<size_t,Eigen::Vector3d> scaledNodes;
            const Eigen::Vector3d x0(spec.X0.row(j));
            for(const auto& node : mshReader.nodes())
            {
                scaledNodes.emplace(node.first,spec.F*(node.second-nodeBary)+x0);
            }
            
            std::map<size_t,std::vector<size_t>> faces;
            size_t eleCounter(0);
            for(const auto& ele : mshReader.elements())
            {
                if(ele.second.type==2)
                {// 3-nodes triangle
                    faces.emplace(eleCounter,ele.second.nodeIDs);
                    eleCounter++;
                }
                
            }
            
            generateSingle(mg,scaledNodes,faces, spec.eigenDistortion, spec.velocityReductionFactor,spec.phaseID);
        }
    }

    //    void PeriodicDipoleGenerator::generate(MicrostructureGenerator& mg)
    //    {
    //        generateIndividual(mg);
    //    }


    //    void PolyhedronInclusionsGenerator::generateDensity(MicrostructureGenerator& mg)
    //    {
    ////        std::cout<<magentaBoldColor<<"Generating inclusions density"<<defaultColor<<std::endl;
    ////        const double targetInclusionDensity(this->parser.readScalar<double>("targetDensity",true));
    ////        if(targetInclusionDensity>0)
    ////        {
    //////            const double inclusionsDiameterLognormalDistribution_M(this->parser.readScalar<double>("diameterLognormalDistribution_M",true) );
    //////            const double inclusionsDiameterLognormalDistribution_S(this->parser.readScalar<double>("diameterLognormalDistribution_S",true));
    //////            const double inclusionsDiameterLognormalDistribution_A(this->parser.readScalar<double>("diameterLognormalDistribution_A",true));
    //////            const Eigen::Matrix<double,1,dim*dim> inclusionsTransformationEigenDistortion(this->parser.readMatrix<double>("transformationEigenDistortion",1,dim*dim,true));
    //////            const Eigen::Matrix<double,1,dim> patternVector(this->parser.readMatrix<double>("patternVector",1,dim,true)/mg.ddBase.poly.b_SI);
    //////            const double velocityReductionFactor(this->parser.readScalar<double>("velocityReductionFactor",true));
    //////            const int phaseIDs(this->parser.readScalar<int>("phaseID",true));
    //////
    //////
    //////            std::lognormal_distribution<double> distribution(log(inclusionsDiameterLognormalDistribution_M/inclusionsDiameterLognormalDistribution_A),inclusionsDiameterLognormalDistribution_S);
    //////
    //////            const double patternHeigth(patternVector.norm());
    //////            const bool applyPattern(patternHeigth>0.0);
    //////            const VectorDimD patternDir(applyPattern? (patternVector/patternHeigth).eval() : VectorDimD::Zero());
    //////
    //////            std::mt19937 generator;
    //////            double density(0.0);
    //////            while(density<targetInclusionDensity)
    //////            {
    //////
    //////                const double diameter = distribution(generator)*inclusionsDiameterLognormalDistribution_A/mg.ddBase.poly.b_SI;
    //////                const double radius(0.5*diameter);
    //////                std::pair<LatticeVector<dim>,int> pointPair=mg.ddBase.poly.randomLatticePointInMesh();
    //////                VectorDimD P=pointPair.first.cartesian();
    //////                const int& grainID(pointPair.second);
    //////
    //////                if(applyPattern)
    //////                {
    //////                    const VectorDimD globalVector(mg.ddBase.poly.grain(grainID).singleCrystal->C2G*patternVector.transpose());
    //////                    const VectorDimD globalDir(mg.ddBase.poly.grain(grainID).singleCrystal->C2G*patternDir);
    //////                    const long long pointHeigth=std::round(P.dot(globalDir)/patternHeigth);
    //////                    const VectorDimD O(pointHeigth*globalVector);
    //////                    P-=(P-O).dot(globalDir)*globalDir;
    //////                }
    //////                if(generateSingle(mg,P,radius,inclusionsTransformationEigenDistortion,velocityReductionFactor,phaseIDs))
    //////                {
    //////                    density+=1.0/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,3);
    //////                    std::cout<<"inclusion density="<<density<<std::endl;
    //////                }
    //////
    //////
    //////            }
    //////
    //////
    //////            //                            const double patternHeigth(currentPattern.norm());
    //////            //                            const bool applyPattern(patternHeigth>0.0);
    //////
    //////
    ////////            const double inclusionsDiameterLognormalDistribution_S(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_S",true) : std::vector<double>())
    ////////            const double inclusionsDiameterLognormalDistribution_A(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_A",true) : std::vector<double>())
    ////////            const double inclusionsTransformationStrains(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsTransformationStrains",targetInclusionDensities.size(),MicrostructureGenerator::dim*MicrostructureGenerator::dim,true) : Eigen::Matrix<double,Eigen::Dynamic,MicrostructureGenerator::dim*MicrostructureGenerator::dim>::Zero(1,MicrostructureGenerator::dim*MicrostructureGenerator::dim))
    ////////            const double patternVectors(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("patternVectors",targetInclusionDensities.size(),MicrostructureGenerator::dim,true) : Eigen::Matrix<double,Eigen::Dynamic,MicrostructureGenerator::dim>::Zero(1,MicrostructureGenerator::dim))
    ////
    ////        }
    //    }

    //    void PolyhedronInclusionsGenerator::generateIndividual(MicrostructureGenerator& mg)
    //    {
    //
    //        std::cout<<magentaBoldColor<<"Generating polyhedral individual inclusions"<<defaultColor<<std::endl;
    //
    //        const std::string mshFile=mg.traits().inputFilesFolder+"/"+this->parser.readString("mshFile",true);
    //        const Eigen::Matrix<double,3,Eigen::Dynamic> X0(this->parser.readMatrixCols<double>("x0",3,true).transpose());
    //        const Eigen::Matrix<double,3,3> A(this->parser.readMatrix<double>("A",3,3,true));
    //        const Eigen::Matrix<double,1,3*3> eT(this->parser.readMatrix<double>("inclusionsEigenDistortions",1,3*3,true));
    //        const double vrc(this->parser.readScalar<double>("inclusionVelocityReductionFactors",true));
    //        const int phaseIDs(this->parser.readScalar<int>("phaseIDs",true));
    //
    //        std::cout<<"mshFile="<<mshFile<<std::endl;
    //        GmshReader mshReader(mshFile);
    //
    //
    //        for(int j=0;j<X0.cols();++j)
    //        {
    //            std::map<size_t,Eigen::Vector3d> scaledNodes;
    //            const auto x0(X0.col(j));
    //            for(const auto& node : mshReader.nodes())
    //            {
    //                scaledNodes.emplace(node.first,A*(node.second-x0));
    //            }
    //
    //            std::map<size_t,std::vector<size_t>> faces;
    //            size_t eleCounter(0);
    //            for(const auto& ele : mshReader.elements())
    //            {
    //                if(ele.second.type==2)
    //                {// 3-nodes triangle
    //                    faces.emplace(eleCounter,ele.second.nodeIDs);
    //                    eleCounter++;
    //                }
    //
    //            }
    //
    //
    //            generateSingle(mg,scaledNodes,faces, eT, vrc,phaseIDs);
    //
    //        }
    //
    //    }



}
#endif
