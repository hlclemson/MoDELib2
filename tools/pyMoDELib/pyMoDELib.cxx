/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_pyMoDELib_cpp_
#define model_pyMoDELib_cpp_

#ifdef _MODEL_PYBIND11_
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#endif

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <set>
#include <memory>
#include <map>

#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <PolycrystallineMaterialBase.h>
#include <Polycrystal.h>
#include <MicrostructureBase.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <DefectiveCrystal.h>
#include <MicrostructureGenerator.h>
#include <MoDELib2Vtk.h>
#include <GlidePlaneNoiseBase.h>


using namespace model;
// https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
// https://pybind11.readthedocs.io/en/stable/advanced/classes.html

typedef Eigen::Matrix<double,3,1> VectorDim;
typedef Eigen::Matrix<double,3,3> MatrixDim;
typedef GlidePlane<3> GlidePlaneType;
typedef DislocationNetwork<3,0> DislocationNetworkType;
typedef typename TypeTraits<DislocationNetworkType>::LoopNodeType LoopNodeType;
typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
typedef typename TypeTraits<DislocationNetworkType>::NetworkNodeType NetworkNodeType;

#ifdef _MODEL_PYBIND11_

PYBIND11_MAKE_OPAQUE(std::map<typename LoopNodeType::KeyType,const std::weak_ptr<LoopNodeType>>);
PYBIND11_MAKE_OPAQUE(std::map<typename LoopType::KeyType,const std::weak_ptr<LoopType>>);
PYBIND11_MAKE_OPAQUE(std::vector<MeshedDislocationLoop>);
PYBIND11_MAKE_OPAQUE(GlidePlaneNoiseBase<1>); // opaque: pybind11 is not going to try to guess the data type for python
PYBIND11_MAKE_OPAQUE(GlidePlaneNoiseBase<2>);

PYBIND11_MODULE(pyMoDELib,m)
{
    namespace py=pybind11;

    // N=1 explicit
    py::class_<GlidePlaneNoiseBase<1>, std::shared_ptr<GlidePlaneNoiseBase<1>>>(m, "GlidePlaneNoiseBase1")
      .def(py::init<const std::string&, const int&,
                    const NoiseTraitsBase::GridSizeType&,
                    const NoiseTraitsBase::GridSpacingType&,
                    const Eigen::Matrix<double,2,2>&>())
      .def("averageNoiseCorrelation", &GlidePlaneNoiseBase<1>::averageNoiseCorrelation)
      .def("sampleAverageNoise", &GlidePlaneNoiseBase<1>::sampleAverageNoise)
    ;

    // Bind MDStackingFaultNoise
    py::class_<MDStackingFaultNoise, GlidePlaneNoiseBase<1>, std::shared_ptr<MDStackingFaultNoise>>(m, "MDStackingFaultNoise")
      .def(py::init<
          const model::PolycrystallineMaterialBase&,
          const std::string&,
          const std::string&,
          const int&,
          const model::NoiseTraitsBase::GridSizeType&,
          const model::NoiseTraitsBase::GridSpacingType&,
          const Eigen::Matrix<double, 2, 2>&
      >())
    ;

    // N=2 explicit
    py::class_<GlidePlaneNoiseBase<2>, std::shared_ptr<GlidePlaneNoiseBase<2>>>(m, "GlidePlaneNoiseBase2")
      .def(py::init<const std::string&, const int&,
                    const NoiseTraitsBase::GridSizeType&,
                    const NoiseTraitsBase::GridSpacingType&,
                    const Eigen::Matrix<double,2,2>&>())
      .def("averageNoiseCorrelation", &GlidePlaneNoiseBase<2>::averageNoiseCorrelation)
      .def("sampleAverageNoise", &GlidePlaneNoiseBase<2>::sampleAverageNoise)
    ;

    // Bind MDSolidSolutionNoise
    py::class_<MDSolidSolutionNoise, GlidePlaneNoiseBase<2>, std::shared_ptr<MDSolidSolutionNoise>>(m, "MDSolidSolutionNoise")
      .def(py::init<
          const model::PolycrystallineMaterialBase&,
          const std::string&,
          const std::string&,
          const std::string&,
          const int&,
          const model::NoiseTraitsBase::GridSizeType&,
          const model::NoiseTraitsBase::GridSpacingType&,
          const Eigen::Matrix<double, 2, 2>&,
          const double&
      >())
    ;

    // Bind AnalyticalSolidSolutionNoise
    py::class_<AnalyticalSolidSolutionNoise, GlidePlaneNoiseBase<2>,std::shared_ptr<AnalyticalSolidSolutionNoise>>(m, "AnalyticalSolidSolutionNoise")
      // Constructor
      .def(py::init<
          const std::string&,
          const int&,
          const model::NoiseTraitsBase::GridSizeType&,
          const model::NoiseTraitsBase::GridSpacingType&,
          const Eigen::Matrix<double, 2, 2>&,
          const double&,
          const double&,
          const double&
      >())
    ;

    py::class_<DDtraitsIO>(m,"DDtraitsIO")
        .def(py::init<const std::string&>())
        .def_readonly("simulationFolder", &DDtraitsIO::simulationFolder)
        .def_readonly("inputFilesFolder", &DDtraitsIO::inputFilesFolder)
        .def_readonly("evlFolder", &DDtraitsIO::evlFolder)
        .def_readonly("auxFolder", &DDtraitsIO::auxFolder)
        .def_readonly("fFolder", &DDtraitsIO::fFolder)
        .def_readonly("ddFile", &DDtraitsIO::ddFile)
        .def_readonly("fFile", &DDtraitsIO::fFile)
        .def_readonly("flabFile", &DDtraitsIO::flabFile)
        .def_readonly("polyFile", &DDtraitsIO::polyFile)
        .def_readonly("materialFile", &DDtraitsIO::materialFile)
        .def_readonly("microstructureFile", &DDtraitsIO::microstructureFile)
        .def_readonly("meshFile", &DDtraitsIO::meshFile)
    ;
    
    py::class_<DefectiveCrystalParameters>(m,"DefectiveCrystalParameters")
        .def(py::init<const std::string&>())
        .def_readwrite("runID", &DefectiveCrystalParameters::runID)
        .def_readonly("useFEM", &DefectiveCrystalParameters::useFEM)
        .def_readonly("traitsIO", &DefectiveCrystalParameters::traitsIO)
    ;
    
    py::class_<MeshRegionObserver<MeshRegion<3>>>(m,"MeshRegionObserver")
        .def(py::init<>())
    ;
    
    py::class_<SimplexReader<3>>(m,"SimplexReader")
        .def(py::init<>())
    ;
    
    py::class_<std::map<typename SimplexTraits<3,3>::SimplexIDType,const Simplex<3,3>>>(m,"SimplexIDMap")
        .def(py::init<>())
    ;
    
    py::class_<std::map<std::pair<size_t,size_t>,MeshRegionBoundary<3>>>(m,"MeshRegionBoundaryMap")
        .def(py::init<>())
    ;
    
    py::class_<SimplicialMesh<3>,
    /*      */ MeshRegionObserver<MeshRegion<3>>,
    /*      */ SimplexReader<3>,
    /*      */ std::map<typename SimplexTraits<3,3>::SimplexIDType,const Simplex<3,3>>,
    /*      */ std::map<std::pair<size_t,size_t>,MeshRegionBoundary<3>>>(m,"SimplicialMesh")
        .def(py::init<>())
//        .def(py::init<const std::string&,const Eigen::Matrix<double,3,3>&,const Eigen::Matrix<double,3,1>&,const std::set<int>&>())
        .def("xMin",static_cast<const Eigen::Matrix<double,3,1>& (SimplicialMesh<3>::*)() const>(&SimplicialMesh<3>::xMin))
        .def("xMax",static_cast<const Eigen::Matrix<double,3,1>& (SimplicialMesh<3>::*)() const>(&SimplicialMesh<3>::xMax))
        .def("volume",&SimplicialMesh<3>::volume)
        .def("xCenter",&SimplicialMesh<3>::xCenter)
    ;

    py::class_<PolycrystallineMaterialBase>(m,"PolycrystallineMaterialBase")
        .def(py::init<const std::string&,const double&>())
    ;

    py::class_<std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<3>>>>(m,"GrainBoundaryMap")
        .def(py::init<>())
    ;
    
    py::class_<Grain<3>,std::map<std::pair<size_t,size_t>,const std::shared_ptr<GrainBoundary<3>>>>(m,"Grain")
        .def(py::init<const MeshRegion<3>&,const PolycrystallineMaterialBase&,const std::string& >())
//            .def_readonly("grainID", &Grain<3>::grainID)
    ;
    
    py::class_<Polycrystal<3>,PolycrystallineMaterialBase>(m,"Polycrystal")
        .def(py::init<const std::string&,const SimplicialMesh<3>&>())
        .def("randomPoint", &Polycrystal<3>::randomPoint)
        .def_readonly("grains", &Polycrystal<3>::grains)
        .def("grain", &Polycrystal<3>::grain)
    ;
    
    py::class_<DislocationDynamicsBase<3>>(m,"DislocationDynamicsBase")
        .def(py::init<const std::string&>())
        .def_readonly("mesh", &DislocationDynamicsBase<3>::mesh)
        .def_readonly("poly", &DislocationDynamicsBase<3>::poly)
        .def_readonly("simulationParameters", &DislocationDynamicsBase<3>::simulationParameters)
    ;
    
    py::class_<MicrostructureBase<3>>(m,"MicrostructureBase")
        .def("displacement", static_cast<Eigen::Matrix<double,Eigen::Dynamic,3> (MicrostructureBase<3>::*)(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>) const>(&MicrostructureBase<3>::displacement))
        .def("stress", static_cast<std::vector<Eigen::Matrix<double,3,3>> (MicrostructureBase<3>::*)(Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>) const>(&MicrostructureBase<3>::stress))
    ;
    
    py::class_<MicrostructureContainer<3>,MicrostructureBase<3>>(m,"MicrostructureContainer",py::multiple_inheritance())
        .def(py::init<DislocationDynamicsBase<3>&>())
    ;

    //py::class_<GlidePlaneNoiseBase<1>>(m,"GlidePlaneNoiseBase")
    //    .def(py::init<const std::string&>())
    //;

    //py::class_<GlidePlaneNoiseBase<2>>(m,"GlidePlaneNoiseBase")
    //    .def(py::init<const std::string&>())
    //;

    py::class_<DD2OvitoVtk>(m,"DD2OvitoVtk")
        .def(py::init<const std::string&>())
    ;
    
    // LoopNode
    py::bind_map<std::map<typename LoopNodeType::KeyType,const std::weak_ptr<LoopNodeType>>>(m, "LoopNodeWeakPtrMap");

    py::class_<WeakPtrFactory<DislocationNetworkType,LoopNodeType>
    /*      */,std::map<typename LoopNodeType::KeyType,const std::weak_ptr<LoopNodeType>>
    /*      */>(m,"LoopNodeWeakPtrFactory")
        .def(py::init<>())
        .def("getRef",&WeakPtrFactory<DislocationNetworkType,LoopNodeType>::getRef,pybind11::return_value_policy::reference)
    ;
    
    py::class_<LoopNodeType
    /*      */ >(m,"LoopNode")
        .def(py::init<typename TypeTraits<DislocationNetworkType>::LoopNetworkType* const,
             const std::shared_ptr<LoopType>&,
             const std::shared_ptr<NetworkNodeType>&,
             const typename TypeTraits<DislocationNetworkType>::VectorDim&,
             const std::shared_ptr<PeriodicPlanePatch<3>>&,
             const std::pair<const std::shared_ptr<PeriodicPlaneEdge<3>>,const std::shared_ptr<PeriodicPlaneEdge<3>>>&>())
    ;
    
    //Loop
    py::bind_vector<std::vector<MeshedDislocationLoop>>(m, "MeshedDislocationLoopVector");

    py::bind_map<std::map<typename LoopType::KeyType,const std::weak_ptr<LoopType>>>(m, "LoopWeakPtrMap");

    py::class_<WeakPtrFactory<DislocationNetworkType,LoopType>
    /*      */,std::map<typename LoopType::KeyType,const std::weak_ptr<LoopType>>
    /*      */>(m,"LoopWeakPtrFactory")
        .def(py::init<>())
        .def("getRef",&WeakPtrFactory<DislocationNetworkType,LoopType>::getRef,pybind11::return_value_policy::reference)
    ;
    
    py::class_<LoopType
    /*      */ >(m,"Loop")
        .def(py::init<DislocationNetworkType* const,
             const VectorDim&,
             const std::shared_ptr<GlidePlaneType>&>())
        .def("solidAngle",&LoopType::solidAngle)
        .def("meshed",&LoopType::meshed)
    ;
    
    py::class_<MeshedDislocationLoop
    /*      */ >(m,"MeshedDislocationLoop")
//        .def(py::init<const VectorDim&,
//             const DislocationDynamicsBase<3>&,
//             const GlidePlane<3>&,const std::vector<Eigen::Matrix<double,3,1>>&,
//             const double&,const double&>())
        .def("solidAngle",&MeshedDislocationLoop::solidAngle)
        .def("plasticDisplacement",&MeshedDislocationLoop::plasticDisplacement)
        .def("plasticDisplacementKernel",&MeshedDislocationLoop::plasticDisplacementKernel)
        .def_readwrite("points", &MeshedDislocationLoop::points)
        .def_readwrite("triangles", &MeshedDislocationLoop::triangles)
    ;

    py::class_<LoopNetwork<DislocationNetworkType>
    /*      */,WeakPtrFactory<DislocationNetworkType,LoopType>
    /*      */,WeakPtrFactory<DislocationNetworkType,LoopNodeType>
    /*      */ >(m,"LoopNetwork",py::multiple_inheritance())
        .def(py::init<>())
        .def("loops", static_cast<const WeakPtrFactory<DislocationNetworkType,LoopType>& (LoopNetwork<DislocationNetworkType>::*)()const>(&LoopNetwork<DislocationNetworkType>::loops),pybind11::return_value_policy::reference)
        .def("loopNodes", static_cast<const WeakPtrFactory<DislocationNetworkType,LoopNodeType>& (LoopNetwork<DislocationNetworkType>::*)()const>(&LoopNetwork<DislocationNetworkType>::loopNodes),pybind11::return_value_policy::reference)
    ;
    
    py::class_<DislocationNetworkType
    /*      */,MicrostructureBase<3>
    /*      */,LoopNetwork<DislocationNetworkType>>(m,"DislocationNetwork")
        .def(py::init<MicrostructureContainer<3>&>())
    ;
    
    py::class_<DefectiveCrystal<3>
    /*      */,MicrostructureContainer<3>
    /*      */>(m,"DefectiveCrystal")
        .def(py::init<DislocationDynamicsBase<3>&>())
        .def("initializeConfiguration", static_cast<void (DefectiveCrystal<3>::*)(const DDconfigIO<3>&)>(&DefectiveCrystal<3>::initializeConfiguration))
        .def("dislocationNetwork", &DefectiveCrystal<3>::dislocationNetwork,pybind11::return_value_policy::reference)
        .def("runSingleStep",&DefectiveCrystal<3>::runSingleStep)
        .def("runSteps",&DefectiveCrystal<3>::runSteps)
    ;

    py::class_<DDconfigIO<3>
    /*      */>(m,"DDconfigIO")
        .def(py::init<const std::string&>())
    ;
    
    py::class_<MicrostructureGenerator
    /*      */>(m,"MicrostructureGenerator")
        .def(py::init<DislocationDynamicsBase<3>&>())
        .def_readonly("configIO", &MicrostructureGenerator::configIO)
        .def("addShearLoopDensity", &MicrostructureGenerator::addShearLoopDensity)
        .def("addShearLoopIndividual", &MicrostructureGenerator::addShearLoopIndividual)
        .def("addPeriodicDipoleDensity", &MicrostructureGenerator::addPeriodicDipoleDensity)
        .def("addPeriodicDipoleIndividual", &MicrostructureGenerator::addPeriodicDipoleIndividual)
        .def("addPrismaticLoopDensity", &MicrostructureGenerator::addPrismaticLoopDensity)
        .def("addPrismaticLoopIndividual", &MicrostructureGenerator::addPrismaticLoopIndividual)
        .def("addFrankLoopsDensity", &MicrostructureGenerator::addFrankLoopsDensity)
        .def("addFrankLoopsIndividual", &MicrostructureGenerator::addFrankLoopsIndividual)
        .def("addStackingFaultTetrahedraDensity", &MicrostructureGenerator::addStackingFaultTetrahedraDensity)
        .def("addStackingFaultTetrahedraIndividual", &MicrostructureGenerator::addStackingFaultTetrahedraIndividual)
        .def("addSphericalInclusionDensity", &MicrostructureGenerator::addSphericalInclusionDensity)
        .def("addSphericalInclusionIndividual", &MicrostructureGenerator::addSphericalInclusionIndividual)
        .def("addPolyhedronInclusionIndividual", &MicrostructureGenerator::addPolyhedronInclusionIndividual)
        .def("writeConfigFiles", &MicrostructureGenerator::writeConfigFiles)
    ;
        
    py::class_<ShearLoopDensitySpecification
    /*      */>(m,"ShearLoopDensitySpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("targetDensity", &ShearLoopDensitySpecification::targetDensity)
        .def_readwrite("numberOfSides", &ShearLoopDensitySpecification::numberOfSides)
        .def_readwrite("radiusDistributionMean", &ShearLoopDensitySpecification::radiusDistributionMean)
        .def_readwrite("radiusDistributionStd", &ShearLoopDensitySpecification::radiusDistributionStd)
        .def_readwrite("allowedGrainIDs", &ShearLoopDensitySpecification::allowedGrainIDs)
        .def_readwrite("allowedSlipSystemIDs", &ShearLoopDensitySpecification::allowedSlipSystemIDs)
    ;

    py::class_<ShearLoopIndividualSpecification
    /*      */>(m,"ShearLoopIndividualSpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("slipSystemIDs", &ShearLoopIndividualSpecification::slipSystemIDs)
        .def_readwrite("loopRadii", &ShearLoopIndividualSpecification::loopRadii)
        .def_property( "loopCenters",
                    [](const ShearLoopIndividualSpecification& self )
                    {// Getter
                        return self.loopCenters;
                    },
                    []( ShearLoopIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>& val )
                    {// Setter
                        self.loopCenters = val;
                    }
                )
        .def_readwrite("loopSides", &ShearLoopIndividualSpecification::loopSides)
    ;
    
    py::class_<PeriodicDipoleDensitySpecification
    /*      */>(m,"PeriodicDipoleDensitySpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("targetDensity", &PeriodicDipoleDensitySpecification::targetDensity)
    ;
    
    py::class_<PeriodicDipoleIndividualSpecification
    /*      */>(m,"PeriodicDipoleIndividualSpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("slipSystemIDs", &PeriodicDipoleIndividualSpecification::slipSystemIDs)
        .def_readwrite("exitFaceIDs", &PeriodicDipoleIndividualSpecification::exitFaceIDs)
        .def_property( "dipoleCenters",
                    [](const PeriodicDipoleIndividualSpecification& self )
                    {// Getter
                        return self.dipoleCenters;
                    },
                    []( PeriodicDipoleIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>& val )
                    {// Setter
                        self.dipoleCenters = val;
                    }
                )
        .def_readwrite("dipoleHeights", &PeriodicDipoleIndividualSpecification::dipoleHeights)
        .def_readwrite("nodesPerLine", &PeriodicDipoleIndividualSpecification::nodesPerLine)
        .def_readwrite("glideSteps", &PeriodicDipoleIndividualSpecification::glideSteps)
    ;
    
    py::class_<PrismaticLoopDensitySpecification
    /*      */>(m,"PrismaticLoopDensitySpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("targetDensity", &PrismaticLoopDensitySpecification::targetDensity)
        .def_readwrite("radiusDistributionMean", &PrismaticLoopDensitySpecification::radiusDistributionMean)
        .def_readwrite("radiusDistributionStd", &PrismaticLoopDensitySpecification::radiusDistributionStd)
        .def_readwrite("allowedGrainIDs", &PrismaticLoopDensitySpecification::allowedGrainIDs)
        .def_readwrite("allowedSlipSystemIDs", &PrismaticLoopDensitySpecification::allowedSlipSystemIDs)
    ;
    
    py::class_<PrismaticLoopIndividualSpecification
    /*      */>(m,"PrismaticLoopIndividualSpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("slipSystemIDs", &PrismaticLoopIndividualSpecification::slipSystemIDs)
        .def_readwrite("loopRadii", &PrismaticLoopIndividualSpecification::loopRadii)
        .def_property( "loopCenters",
                    [](const PrismaticLoopIndividualSpecification& self )
                    {// Getter
                        return self.loopCenters;
                    },
                    []( PrismaticLoopIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>& val )
                    {// Setter
                        self.loopCenters = val;
                    }
                )
        .def_readwrite("glideSteps", &PrismaticLoopIndividualSpecification::glideSteps)
    ;
    
    py::class_<FrankLoopsDensitySpecification
    /*      */>(m,"FrankLoopsDensitySpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("targetDensity", &FrankLoopsDensitySpecification::targetDensity)
        .def_readwrite("numberOfSides", &FrankLoopsDensitySpecification::numberOfSides)
        .def_readwrite("radiusDistributionMean", &FrankLoopsDensitySpecification::radiusDistributionMean)
        .def_readwrite("radiusDistributionStd", &FrankLoopsDensitySpecification::radiusDistributionStd)
        .def_readwrite("areVacancyLoops", &FrankLoopsDensitySpecification::areVacancyLoops)
    ;

    py::class_<FrankLoopsIndividualSpecification
    /*      */>(m,"FrankLoopsIndividualSpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("planeIDs", &FrankLoopsIndividualSpecification::planeIDs)
        .def_readwrite("loopRadii", &FrankLoopsIndividualSpecification::loopRadii)
        .def_property( "loopCenters",
                    [](const FrankLoopsIndividualSpecification& self )
                    {// Getter
                        return self.loopCenters;
                    },
                    []( FrankLoopsIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>& val )
                    {// Setter
                        self.loopCenters = val;
                    }
                )
        .def_readwrite("loopSides", &FrankLoopsIndividualSpecification::loopSides)
        .def_readwrite("isVacancyLoop", &FrankLoopsIndividualSpecification::isVacancyLoop)
    ;
    
    py::class_<StackingFaultTetrahedraDensitySpecification
    /*      */>(m,"StackingFaultTetrahedraDensitySpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("targetDensity", &StackingFaultTetrahedraDensitySpecification::targetDensity)
        .def_readwrite("sizeDistributionMean", &StackingFaultTetrahedraDensitySpecification::sizeDistributionMean)
        .def_readwrite("sizeDistributionStd", &StackingFaultTetrahedraDensitySpecification::sizeDistributionStd)
    ;

    py::class_<StackingFaultTetrahedraIndividualSpecification
    /*      */>(m,"StackingFaultTetrahedraIndividualSpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("planeIDs", &StackingFaultTetrahedraIndividualSpecification::planeIDs)
        .def_readwrite("areInverted", &StackingFaultTetrahedraIndividualSpecification::areInverted)
        .def_readwrite("sizes", &StackingFaultTetrahedraIndividualSpecification::sizes)
        .def_property( "basePoints",
                    [](const StackingFaultTetrahedraIndividualSpecification& self )
                    {// Getter
                        return self.basePoints;
                    },
                    []( StackingFaultTetrahedraIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>& val )
                    {// Setter
                        self.basePoints = val;
                    }
                )
    ;
    
    py::class_<SphericalInclusionDensitySpecification
    /*      */>(m,"SphericalInclusionDensitySpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("targetDensity", &SphericalInclusionDensitySpecification::targetDensity)
        .def_readwrite("diameterLognormalDistribution_M", &SphericalInclusionDensitySpecification::diameterLognormalDistribution_M)
        .def_readwrite("diameterLognormalDistribution_S", &SphericalInclusionDensitySpecification::diameterLognormalDistribution_S)
        .def_readwrite("diameterLognormalDistribution_A", &SphericalInclusionDensitySpecification::diameterLognormalDistribution_A)
        .def_property( "transformationEigenDistortion",
                    [](const SphericalInclusionDensitySpecification& self )
                    {// Getter
                        return self.transformationEigenDistortion;
                    },
                    []( SphericalInclusionDensitySpecification& self, const Eigen::Ref<const Eigen::Matrix<double,1,3*3>>& val )
                    {// Setter
                        self.transformationEigenDistortion = val;
                    }
                )
        .def_property( "patternVector_SI",
                    [](const SphericalInclusionDensitySpecification& self )
                    {// Getter
                        return self.patternVector_SI;
                    },
                    []( SphericalInclusionDensitySpecification& self, const Eigen::Ref<const Eigen::Matrix<double,1,3>>& val )
                    {// Setter
                        self.patternVector_SI = val;
                    }
                )
        .def_readwrite("allowOverlap", &SphericalInclusionDensitySpecification::allowOverlap)
        .def_readwrite("allowOutside", &SphericalInclusionDensitySpecification::allowOutside)
        .def_readwrite("velocityReductionFactor", &SphericalInclusionDensitySpecification::velocityReductionFactor)
        .def_readwrite("phaseID", &SphericalInclusionDensitySpecification::phaseID)
    ;
    
    py::class_<SphericalInclusionIndividualSpecification
    /*      */>(m,"SphericalInclusionIndividualSpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("radii_SI", &SphericalInclusionIndividualSpecification::radii_SI)
        .def_property( "centers",
                    [](const SphericalInclusionIndividualSpecification& self )
                    {// Getter
                        return self.centers;
                    },
                    []( SphericalInclusionIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>& val )
                    {// Setter
                        self.centers = val;
                    }
                )
        .def_property( "eigenDistortions",
                    [](const SphericalInclusionIndividualSpecification& self )
                    {// Getter
                        return self.eigenDistortions;
                    },
                    []( SphericalInclusionIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3*3>>& val )
                    {// Setter
                        self.eigenDistortions = val;
                    }
                )
        .def_readwrite("velocityReductionFactors", &SphericalInclusionIndividualSpecification::velocityReductionFactors)
        .def_readwrite("phaseIDs", &SphericalInclusionIndividualSpecification::phaseIDs)
        .def_readwrite("allowOverlap", &SphericalInclusionIndividualSpecification::allowOverlap)
        .def_readwrite("allowOutside", &SphericalInclusionIndividualSpecification::allowOutside)
    ;
    
    py::class_<PolyhedronInclusionIndividualSpecification
    /*      */>(m,"PolyhedronInclusionIndividualSpecification")
        .def(py::init<>())
        .def(py::init<const std::string&>())
        .def_readwrite("mshFile", &PolyhedronInclusionIndividualSpecification::mshFile)
        .def_property( "X0",
                    [](const PolyhedronInclusionIndividualSpecification& self )
                    {// Getter
                        return self.X0;
                    },
                    []( PolyhedronInclusionIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,Eigen::Dynamic,3>>& val )
                    {// Setter
                        self.X0 = val;
                    }
                )
        .def_property( "F",
                    [](const PolyhedronInclusionIndividualSpecification& self )
                    {// Getter
                        return self.F;
                    },
                    []( PolyhedronInclusionIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,3,3>>& val )
                    {// Setter
                        self.F = val;
                    }
                )
        .def_property( "eigenDistortion",
                    [](const PolyhedronInclusionIndividualSpecification& self )
                    {// Getter
                        return self.eigenDistortion;
                    },
                    []( PolyhedronInclusionIndividualSpecification& self, const Eigen::Ref<const Eigen::Matrix<double,1,3*3>>& val )
                    {// Setter
                        self.eigenDistortion = val;
                    }
                )
        .def_readwrite("velocityReductionFactor", &PolyhedronInclusionIndividualSpecification::velocityReductionFactor)
        .def_readwrite("phaseID", &PolyhedronInclusionIndividualSpecification::phaseID)
    ;
    
}
#endif


int main(int argc, char** argv)
{
    return 0;
}

#endif
