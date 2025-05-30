/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

//#ifndef model_DD2OvitoVtk_H_segmentMap
//#define model_DD2OvitoVtk_H_

#ifdef _MODEL_PYBIND11_
    #undef slots
    #include <pybind11/embed.h>
    #include <pybind11/eigen.h>
    #include <pybind11/stl.h>
    #define slots Q_SLOTS
#endif

#include <tuple>
#include <map>
#include <Eigen/Dense>

//#include <TextFileParser.h>
//#include <TerminalColors.h>
#include <DDconfigIO.h>
#include <PolycrystallineMaterialBase.h>
#include <Polycrystal.h>
#include <DDauxIO.h>
#include <filesystem>
#include <DDtraitsIO.h>

namespace model
{

class DD2OvitoVtk
{
  DDtraitsIO traitsIO;
  DDconfigIO<3> configIO;
  DDauxIO<3> auxIO;
  const std::set<int> periodicFaceIDs;

  // const std::string meshFilename;
  const SimplicialMesh<3> mesh;
  //        const double minSize;
  //        const double maxSize;
  Polycrystal<3> poly;

public:
  DD2OvitoVtk(const std::string &folderName) : 
    /* init*/ traitsIO(folderName)
    /* init*/,configIO(traitsIO.evlFolder)
    /* init*/,auxIO(traitsIO.auxFolder)
    /* init */,periodicFaceIDs(TextFileParser(traitsIO.polyFile).template readSet<int>("periodicFaceIDs", true))
    /* init */,mesh(traitsIO.meshFile, TextFileParser(traitsIO.polyFile).readMatrix<double>("F", 3, 3, true),TextFileParser(traitsIO.polyFile).readMatrix<double>("X0", 1, 3, true).transpose(), periodicFaceIDs)
    /* init*/,poly(traitsIO.polyFile, mesh)
  ///* init */ mat(TextFileParser("inputFiles/polycrystal.txt").readString("materialFile",true))
  ///* init*/ meshFilename(TextFileParser("./inputFiles/polycrystal.txt").readString("meshFile",true))
  ///* init*/,mesh(meshFilename,TextFileParser("./inputFiles/polycrystal.txt").readMatrix<double>("A",3,3,true),TextFileParser("./inputFiles/polycrystal.txt").readMatrix<double>("x0",1,3,true).transpose())
  ///* init*/,minSize(0.1*min(mesh.xMax(0)-mesh.xMin(0),min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
  ///* init*/,maxSize(max(mesh.xMax(0)-mesh.xMin(0),max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
  ///* init*/,poly("./inputFiles/polycrystal.txt",mesh)
  ///* init*/,glidePlaneFactory(poly)
  {
    meshBoundary2vtk();
    evl2vtkNodes(0);
    evl2vtkQuadraturePoints(traitsIO.evlFolder);
  }

  void meshBoundary2vtk()
  {

    size_t numFaces(0);
    size_t totFaceVerteces(0);

    std::set<const Simplex<3, 0> *> bndVertices;
    for (const auto &region : poly.mesh.regions())
    {
      for (const auto &face : region.second->faces())
      {
        for (const auto &vertex : face.second->convexHull())
        {
          bndVertices.insert(vertex);
          totFaceVerteces++;
        }
        totFaceVerteces += 2;
        numFaces++;
      }
    }

    std::ofstream ovitoFile(traitsIO.evlFolder + "/meshFaces.vtk");
    ovitoFile << "# vtk DataFile Version 3.0\n";
    ovitoFile << "# Mesh faces converted from MoDELib file\n";
    ovitoFile << "ASCII\n";
    ovitoFile << "DATASET UNSTRUCTURED_GRID\n";
    ovitoFile << "POINTS " + std::to_string(bndVertices.size()) + " double\n";

    for (const auto &v : bndVertices)
    {
      ovitoFile << std::setprecision(15) << std::scientific << v->P0.transpose() * poly.b_SI * 1.0e10 << "\n";
    }

    ovitoFile << "\nCELLS " + std::to_string(numFaces) + " " + std::to_string(totFaceVerteces) + "\n";
    for (const auto &region : poly.mesh.regions())
    {
      for (const auto &face : region.second->faces())
      {
        ovitoFile << face.second->convexHull().size() + 1 << " "; // not sure if this is the right number to write

        for (const auto &vertex : face.second->convexHull())
        {
          const auto vtxIter(bndVertices.find(vertex));
          const size_t vtxID(std::distance(bndVertices.begin(), vtxIter));
          ovitoFile << vtxID << " "; // not sure if this is the right number to write
        }
        const auto firstVertex(bndVertices.find(face.second->convexHull()[0]));
        const size_t firstVertexID(std::distance(bndVertices.begin(), firstVertex));
        ovitoFile << firstVertexID << "\n"; // not sure if this is the right number to write
      }
    }

    ovitoFile << "\nCELL_TYPES " + std::to_string(numFaces) + "\n";
    for (const auto &region : poly.mesh.regions())
    {
      for (const auto &face : region.second->faces())
      {
        ovitoFile << 7 << "\n"; // not sure if this is the right number to write
      }
    }
  }

  double latticeParameter() const
  {
    if (poly.crystalStructure == "BCC")
    {
      return 2.0 * poly.b_SI / sqrt(3.0);
    }
    else if (poly.crystalStructure == "FCC")
    {
      return 2.0 * poly.b_SI / sqrt(2.0);
    }
    //            else if(crystalStructure=="HCP")
    //            {
    //                return hcpLattice->reciprocalPlaneNormals(lat);
    //            }
    else
  {
      std::cout << "Unknown lattice parameter for " << poly.crystalStructure << "'. Exiting." << std::endl;
      exit(EXIT_FAILURE);
      return 0.0;
    }
  }

  void evl2vtkQuadraturePoints(const std::string &path)
  {

    const std::string subS1 = "evl_";
    const std::string subS2 = ".txt";
    const std::string subS3 = ".bin";

    for (const auto &p : std::filesystem::directory_iterator(path))
    {

      const size_t found1 = p.path().string().find(subS1);
      const size_t found2 = p.path().string().find(subS2);
      const size_t found3 = p.path().string().find(subS3);

      if (found1 != std::string::npos && (found2 != std::string::npos || found3 != std::string::npos))
      {

        const size_t foundExt(found3 != std::string::npos ? found3 : found2);

        const std::string stringID = p.path().string().substr(found1 + subS1.length(), foundExt - (found1 + subS1.length()));
        const size_t runID(std::stol(stringID));
        evl2vtkQuadraturePoints(runID);
      }
    }
  }

  void evl2vtkQuadraturePoints(const size_t &frameID)
  {
    if (configIO.isBinGood(frameID))
    {
      configIO.readBin(frameID);
      auxIO.readBin(frameID);
    }
    else
  {
      configIO.readTxt(frameID);
      auxIO.readTxt(frameID);
    }

    std::map<size_t, size_t> nodeIDmap; // key=nodeID, value=row index in file
    size_t k = 0;
    for (const auto &node : configIO.nodes())
    {
      if (node.meshLocation!=-2)
      {
        nodeIDmap.emplace(node.sID, k);
        k++;
      }

    }

    const auto segmentMap(configIO.segments());
    std::map<std::pair<size_t, size_t>, std::map<size_t, DislocationQuadraturePointIO<3>>> qPointsMap; // organized by segments
    std::map<std::pair<size_t, size_t>, DislocationSegmentIO<3>> qPointsBurgers;
    for (const auto &qp : auxIO.quadraturePoints())
    {
      const std::pair<size_t, size_t> outerKey(std::make_pair(qp.sourceID, qp.sinkID));
      const size_t innerKey(qp.qID);
      qPointsMap[outerKey].emplace(innerKey, qp);
      //                ovitoFile<<std::setprecision(15)<<std::scientific<<qp.r.transpose()*poly.b_SI*1.0e10<<"\n";
      // std::cout << "outerKey=" << outerKey.first << "," << outerKey.second << std::endl;
      // std::cout << "qp=" << qp.sourceID << "," << qp.sinkID << std::endl;
    }
    // for(const auto& outerPair : qPointsMap)
    //{
    //   qPointsBurgers[outerPair.first].b
    // }

    std::ofstream ovitoFile(traitsIO.evlFolder + "/quadrature_" + std::to_string(frameID) + ".vtk");
    ovitoFile << "# vtk DataFile Version 3.0\n";
    ovitoFile << "# Dislocation lines converted from MoDELib file\n";
    //            ovitoFile<<"# length unit = A\n";
    //            ovitoFile<<"# energy unit = eV\n";
    //            ovitoFile<<"# time unit = ps\n";
    //            ovitoFile<<"# stress unit = GPa";
    ovitoFile << "ASCII\n";
    ovitoFile << "DATASET UNSTRUCTURED_GRID\n";
    ovitoFile << "POINTS " + std::to_string(auxIO.quadraturePoints().size() + nodeIDmap.size()) + " double\n";
    //            for(const auto& qp : auxIO.quadraturePoints())
    //            {
    //                ovitoFile<<std::setprecision(15)<<std::scientific<<qp.r.transpose()*poly.b_SI*1.0e10<<"\n";
    //            }
    int nonGhostCount=0;
    for (const auto &node : configIO.nodes())
    { // write all nodes positions
      if (node.meshLocation!=-2)
      {
        ovitoFile << std::setprecision(15) << std::scientific << node.P.transpose() * poly.b_SI * 1.0e10 << "\n";
        nonGhostCount++;
      }
    }

    for (const auto &outerPair : qPointsMap)
    { // write all quadrature points positions
      for (const auto &innerPair : outerPair.second)
      {
        const auto &qp(innerPair.second);
        ovitoFile << std::setprecision(15) << std::scientific << qp.r.transpose() * poly.b_SI * 1.0e10 << "\n";
        //                    ovitoFile<<std::setprecision(15)<<std::scientific<<qp.r.transpose()<<"\n";
      }
    }
    // std::cout << "qPointsMap=" << qPointsMap.size() << std::endl;
    // std::cout << "qPointSize=" << auxIO.quadraturePoints().size() << std::endl;
    ovitoFile << "\nCELLS " + std::to_string(qPointsMap.size()) + " " + std::to_string(auxIO.quadraturePoints().size() + 3 * qPointsMap.size()) + "\n";
    size_t qpCounter(0);
    for (const auto &pair : qPointsMap)
    {

      const auto &sourceNodeID(pair.first.first);
      const auto &sinkNodeID(pair.first.second);

      const auto sourceNodeRow(nodeIDmap.at(sourceNodeID));
      const auto sinkNodeRow(nodeIDmap.at(sinkNodeID));

      ovitoFile << pair.second.size() + 2 << " " << sourceNodeRow << " ";
      //for (int k = 0; k < pair.second.size(); ++k)
      for (unsigned long k = 0; k < pair.second.size(); ++k)
      {
        ovitoFile << qpCounter + k + nodeIDmap.size() << " ";
      }
      ovitoFile << sinkNodeRow << "\n";
      qpCounter += pair.second.size();
    }

    ovitoFile << "\nCELL_TYPES " + std::to_string(qPointsMap.size()) + "\n";
    for (const auto &pair : qPointsMap)
    {
      ovitoFile << 4 << "\n"; // not sure if this is the right number to write
    }

    ovitoFile << "\nPOINT_DATA " + std::to_string(auxIO.quadraturePoints().size() + nodeIDmap.size()) + "\n";
    ovitoFile << "VECTORS PK_force double\n";
    //            for(const auto& qp : auxIO.quadraturePoints())
    //            {
    //                ovitoFile<<std::setprecision(15)<<std::scientific<<qp.pkForce.transpose()*poly.b_SI*poly.mu_SI*10.0/160.21766208<<"\n";
    //            }

    for (const auto &node : nodeIDmap)
    {
      ovitoFile << std::setprecision(15) << std::scientific << Eigen::Matrix<double, 1, 3>::Zero() << "\n";
    }

    for (const auto &outerPair : qPointsMap)
    {
      for (const auto &innerPair : outerPair.second)
      {
        const auto &qp(innerPair.second);
        ovitoFile << std::setprecision(15) << std::scientific << qp.pkForce.transpose() * poly.b_SI * poly.mu_SI * 10.0 / 160.21766208 << "\n";
        //                    ovitoFile<<std::setprecision(15)<<std::scientific<<qp.pkForce.transpose()<<"\n";
      }
    }

    ovitoFile << "\nVECTORS velocity double\n";
    //            for(const auto& qp : auxIO.quadraturePoints())
    //            {
    //                ovitoFile<<std::setprecision(15)<<std::scientific<<qp.glideVelocity.transpose()*poly.cs_SI*0.01<<"\n";
    //            }

    for (const auto &node : nodeIDmap)
    {
      //                ovitoFile<<std::setprecision(15)<<std::scientific<<node.V.transpose()*poly.cs_SI*0.01<<"\n";
      ovitoFile << std::setprecision(15) << std::scientific << Eigen::Matrix<double, 1, 3>::Zero() << "\n";
    }

    for (const auto &outerPair : qPointsMap)
    {
      for (const auto &innerPair : outerPair.second)
      {
        const auto &qp(innerPair.second);
        ovitoFile << std::setprecision(15) << std::scientific << qp.velocity.transpose() * poly.cs_SI * 0.01 << "\n";
      }
    }

    // do I need cell types to load the force?
    //ovitoFile << "\nCELL_TYPES " + std::to_string(qPointsMap.size()) + "\n";
    //for (const auto &pair : qPointsMap)
    //{
    //    ovitoFile << 4 << "\n"; // not sure if this is the right number to write
    //}
    //// Add stacking fault force in vector format
    //ovitoFile << "\nPOINT_DATA " + std::to_string(auxIO.quadraturePoints().size() + nodeIDmap.size()) + "\n";
    ovitoFile << "VECTORS SF_force double\n";
    for (const auto &node : nodeIDmap)
    {
      ovitoFile << std::setprecision(15) << std::scientific << Eigen::Matrix<double, 1, 3>::Zero() << "\n";
    }

    for (const auto &outerPair : qPointsMap)
    {
      for (const auto &innerPair : outerPair.second)
      {
        const auto &qp(innerPair.second);
        ovitoFile << std::setprecision(15) << std::scientific << qp.stackingFaultForce.transpose() * poly.b_SI * poly.mu_SI << "\n";
      }
    }


    ovitoFile << "\nTENSORS stress double\n";
    //            for(const auto& qp : auxIO.quadraturePoints())
    //            {
    //                ovitoFile<<std::setprecision(15)<<std::scientific<<qp.stress*poly.mu_SI*1e-9/160.21766208<<"\n";
    //            }

    for (const auto &node : nodeIDmap)
    {
      ovitoFile << std::setprecision(15) << Eigen::Matrix<double, 3, 3>::Zero() << "\n\n";
    }

    for (const auto &outerPair : qPointsMap)
    {
      for (const auto &innerPair : outerPair.second)
      {
        const auto &qp(innerPair.second);
        ovitoFile << std::setprecision(15) << std::scientific << qp.stress * poly.mu_SI * 1e-9 << "\n\n";
        // ovitoFile<<std::setprecision(15)<<std::scientific<<qp.stress*1e-9<<"\n\n";
      }
    }

    ovitoFile << "\nCELL_DATA " + std::to_string(qPointsMap.size()) + "\n";
    ovitoFile << "SCALARS dislocation_index int\n";
    ovitoFile << "LOOKUP_TABLE default\n";
    for (size_t k = 0; k < qPointsMap.size(); ++k)
    {
      ovitoFile << k << "\n"; // not sure if this is the right number to write
    }

    ovitoFile << "\nVECTORS burgers_vector_global double\n";

    for (const auto &outerPair : qPointsMap)
    {
      // std::cout<<"BVType:"<<std::__cmp_cat::type(segmentMap.at(outerPair.first).b.transpose())<<std::endl;
      ovitoFile << segmentMap.at(outerPair.first).b.transpose() * poly.b_SI * 1.0e10 << "\n";
    }
    std::cout<<"Poly.b.si="<<poly.b_SI<<std::endl;
    ovitoFile << "\nVECTORS burgers_vector_local double\n";
    const auto &grain(poly.grains.begin()->second);
    const auto G2C(poly.grains.begin()->second.singleCrystal->C2G.inverse());
    std::cout<< "G2C= " <<G2C<<std::endl;
    // const auto G2C(poly.grain.singleCrystal->C2G.inverse());

    for (const auto &outerPair : qPointsMap)
    {
      ovitoFile << (G2C * segmentMap.at(outerPair.first).b).transpose() * poly.b_SI / latticeParameter() << "\n";
    }
    // for(const auto& loop : configIO.loops())
    // {
    //   ovitoFile<<(G2C*loop.B).transpose()*poly.b_SI/latticeParameter()<<"\n"; // not sure if this is the right number to write
    // }

    // ovitoFile<<"\nVECTORS burgers_vector_world double\n";
    // for(const auto& loop : configIO.loops())
    // {
    //   ovitoFile<<loop.B.transpose()*poly.b_SI*1.0e10<<"\n"; // not sure if this is the right number to write
    // }
  }

  void evl2vtkNodes(const size_t &frameID)
  {
    if (configIO.isBinGood(frameID))
    {
      configIO.readBin(frameID);
    }
    else
  {
      configIO.readTxt(frameID);
    }

    //DDconfigIO<3> configIONoGhost;

    // for (const auto &node : configIO.nodes())
    // {
    //     if (node.meshLocation != -2)
    //     {
    //         configIONoGhost.nodes().push_back(node);
    //     }
    // }
    // for (const auto &loopLink : configIO.loopLinks())
    // {   
    //     configIONoGhost.loopLinks().push_back(loopLink);
    // }
    // for (const auto &loop : configIO.loops())
    // {
    //     configIONoGhost.loops().push_back(loop);
    // }
    // for (const auto &loopNode : configIO.loopNodes())
    // {
    //     configIONoGhost.loopNodes().push_back(loopNode);
    //}


    //            for(const auto& gp : quadraturePoints())
    //            {evl2vtkNodes.se
    //
    //                ovitoFile<<std::setprecision(15)<<std::scientific<<gp.r.*poly.b_SI*1.0e10<<"\n";
    //
    //            }
    int nonGhostCount=0;
    for (const auto &node : configIO.nodes())
    {
      if (node.meshLocation!=-2)
      {
        nonGhostCount++;
      }

    }
    std::ofstream ovitoFile(traitsIO.evlFolder + "/config_" + std::to_string(frameID) + ".vtk");
    ovitoFile << "# vtk DataFile Version 3.0\n";
    ovitoFile << "# Dislocation lines converted from MoDELib file\n";
    ovitoFile << "ASCII\n";
    ovitoFile << "DATASET UNSTRUCTURED_GRID\n";
    ovitoFile << "POINTS " + std::to_string(nonGhostCount) + " double\n";
    for (const auto &node : configIO.nodes())
    {
      if (node.meshLocation!=-2)
      {
        ovitoFile << std::setprecision(15) << std::scientific << node.P.transpose() * poly.b_SI * 1.0e10 << "\n";
      }

    }

    std::map<size_t, std::map<size_t, size_t>> linkMap; // loopID,sourceID,sinkID
    for (const auto &link : configIO.loopLinks())
    {
      linkMap[link.loopID].emplace(link.sourceID, link.sinkID);
    }

    ovitoFile << "\nCELLS " + std::to_string(linkMap.size()) + " " + std::to_string(nonGhostCount + 2 * configIO.loops().size()) + "\n";
    for (const auto &pair : linkMap)
    {
      ovitoFile << pair.second.size() + 1 << " ";

      auto edgeIter(pair.second.begin());
      size_t currentSource(edgeIter->first);
      size_t currentSink(edgeIter->second);
      for (unsigned long k = 0; k < pair.second.size(); ++k)
      {
        ovitoFile << currentSource << " ";
        edgeIter = pair.second.find(currentSink);
        currentSource = edgeIter->first;
        currentSink = edgeIter->second;
      }
      ovitoFile << currentSource << "\n";
    }

    ovitoFile << "\nCELL_TYPES " + std::to_string(configIO.loops().size()) + "\n";
    for (const auto &loop : configIO.loops())
    {
      ovitoFile << configIO.loops().size() << "\n"; // not sure if this is the right number to write
    }

    ovitoFile << "\nCELL_DATA " + std::to_string(configIO.loops().size()) + "\n";
    ovitoFile << "SCALARS dislocation_index int\n";
    ovitoFile << "LOOKUP_TABLE default\n";
    for (size_t k = 0; k < configIO.loops().size(); ++k)
    {
      ovitoFile << k << "\n"; // not sure if this is the right number to write
    }

    ovitoFile << "\nVECTORS burgers_vector_local double\n";
    const auto &grain(poly.grains.begin()->second);
    const auto G2C(poly.grains.begin()->second.singleCrystal->C2G.inverse());
    for (const auto &loop : configIO.loops())
    {
      ovitoFile << (G2C * loop.B).transpose() * poly.b_SI / latticeParameter() << "\n"; // not sure if this is the right number to write
      //std::cout<<"loop.b.transpose: "<<loop.B.transpose()<<std::endl;
    }
    //std::cout<<"C2G: "<<poly.grains.begin()->second.singleCrystal->C2G<< " \nG2C: "<<G2C<<std::endl;

    ovitoFile << "\nVECTORS burgers_vector_world double\n";
    for (const auto &loop : configIO.loops())
    {
      ovitoFile << loop.B.transpose() * poly.b_SI * 1.0e10 << "\n"; // not sure if this is the right number to write
    }
  }
};

}
