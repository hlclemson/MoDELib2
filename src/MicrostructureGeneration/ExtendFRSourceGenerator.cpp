/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExtendFRSourceGenerator_cpp_
#define model_ExtendFRSourceGenerator_cpp_


#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>
#include <random>
#include <iomanip>

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
#include <Line.h>
#include <MicrostructureGenerator.h>
#include <ExtendFRSourceGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

    ExtendFRSourceGenerator::ExtendFRSourceGenerator(const ExtendFRSourceIndividualSpecification& spec,MicrostructureGenerator& mg)
    {
        std::cout<<magentaBoldColor<<"Generating individual ExtendFRSources"<<defaultColor<<std::endl;
        if(spec.slipSystemIDs.size())
        {
            if(spec.slipSystemIDs.size()!=spec.exitFaceIDs.size())
            {
                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO exitFaceIDs.size()="+std::to_string(spec.exitFaceIDs.size()));
            }
            if(int(spec.slipSystemIDs.size())!=spec.ExtendFRSourceCenters.rows())
            {
                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO ExtendFRSourceCenters.rows()="+std::to_string(spec.ExtendFRSourceCenters.rows()));
            }
            if(spec.slipSystemIDs.size()!=spec.ExtendFRSourceHeights.size())
            {
                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO ExtendFRSourceHeights.size()="+std::to_string(spec.ExtendFRSourceHeights.size()));
            }
            if(spec.slipSystemIDs.size()!=spec.nodesPerLine.size())
            {
                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO nodesPerLine.size()="+std::to_string(spec.nodesPerLine.size()));
            }
            if(spec.slipSystemIDs.size()!=spec.glideSteps.size())
            {
                throw std::runtime_error("slipSystemIDs.size()="+std::to_string(spec.slipSystemIDs.size())+" NOT EQUAL TO ExtendFRSourceGlideSteps.size()="+std::to_string(spec.glideSteps.size()));
            }
            
            for(size_t k=0;k<spec.slipSystemIDs.size();++k)
            {
                generateSingle(mg,spec.slipSystemIDs[k],spec.ExtendFRSourceCenters.row(k),spec.exitFaceIDs[k],spec.ExtendFRSourceHeights[k],spec.nodesPerLine[k],spec.glideSteps[k],spec.ExtendFRSourceLengths[k]); // changed
            }
        }
        
    }

    ExtendFRSourceGenerator::ExtendFRSourceGenerator(const ExtendFRSourceDensitySpecification& spec,MicrostructureGenerator& mg)
    {
        std::cout<<magentaBoldColor<<"Generating ExtendFRSource density"<<defaultColor<<std::endl;
    //    const double targetDensity(this->parser.readScalar<double>("targetDensity",true));
        if(spec.targetDensity>0.0)
        {
            std::mt19937 generator;
            double density=0.0;
            while(density<spec.targetDensity)
            {
                const std::pair<LatticeVector<3>, int> rp(mg.ddBase.poly.randomLatticePointInMesh());
                const LatticeVector<3> L0=rp.first;
                const size_t grainID=rp.second;
                std::uniform_int_distribution<> ssDist(0,mg.ddBase.poly.grain(grainID)->slipSystems().size()-1);
                const int rSS(ssDist(generator)); // a random SlipSystem
                std::uniform_int_distribution<> fDist(0,mg.ddBase.poly.grain(grainID)->region.faces().size()-1);
                const int rF(fDist(generator)); // a random face
                auto faceIter(mg.ddBase.poly.grain(grainID)->region.faces().begin());
                std::advance(faceIter,rF);

                try
                {
                    generateSingle(mg,rSS,L0.cartesian(),faceIter->first,100,0,20.0,40.0);
                    density+=2.0*faceIter->second->periodicFacePair.first.norm()/mg.ddBase.mesh.volume()/std::pow(mg.ddBase.poly.b_SI,2);
                    std::cout<<"Extend FR Source density="<<density<<std::endl;
                }
                catch(const std::exception& e)
                {
                    
                }
            }
        }
    }

    bool ExtendFRSourceGenerator::generateSingle(MicrostructureGenerator& mg,const int& rSS,const VectorDimD& ExtendFRSourcePoint,const int& exitFaceID,const int& ExtendFRSourceHeight,const int& ExtendFRSourceNodes, double glideStep, double ExtendFRSourceLength)
    {
        
        
        if(rSS>=0)
        {
            std::pair<bool,const Simplex<3,3>*> found(mg.ddBase.mesh.search(ExtendFRSourcePoint));
            if(!found.first)
            {
                std::cout<<"Point "<<ExtendFRSourcePoint.transpose()<<" is outside mesh. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            const int grainID(found.second->region->regionID);
            assert(mg.ddBase.poly.grains.size()==1 && "Extended FRSource(periodic) only supported for single crystals");
            const auto& grain(mg.ddBase.poly.grain(grainID));

            if(rSS<int(grain->slipSystems().size()))
            {
                const auto periodicFaceIter(grain->region.faces().find(exitFaceID));
                if(periodicFaceIter!=grain->region.faces().end())
                {
                    const auto periodicFaceA(periodicFaceIter->second);
                    const auto periodicFaceB(periodicFaceA->periodicFacePair.second);

                    if(periodicFaceB!=nullptr)
                    {

                        const auto& faceAshift(periodicFaceA->periodicFacePair.first);
                        const auto faceAlatticeShift(grain->latticeVector(faceAshift));
                        
                        const auto& slipSystem(*grain->slipSystems().at(rSS));
                        
                        if(slipSystem.n.dot(faceAlatticeShift)==0)
                        {
                            
                            const VectorDimD planePoint(ExtendFRSourcePoint-0.5*ExtendFRSourceHeight*slipSystem.n.interplaneVector());
                            
                            // Construct lower glide plane
                            const long int planeIndex(slipSystem.n.closestPlaneIndexOfPoint(planePoint));
                            GlidePlaneKey<3> glidePlaneKey(planeIndex, slipSystem.n);
                            std::shared_ptr<PeriodicGlidePlane<3>> glidePlane(mg.ddBase.periodicGlidePlaneFactory.get(glidePlaneKey));

                             const VectorDimD P0(grain->snapToLattice(planePoint).cartesian()); // WARNING: this may shift the point compared to the input.
                            
                            // Construct Prismatic Plane
                            const auto prismaticPlaneNormal(grain->reciprocalLatticeDirection(glidePlane->referencePlane->unitNormal.cross(faceAshift)));
                            const long int prismaticPlaneIndex(prismaticPlaneNormal.closestPlaneIndexOfPoint(P0));
                            GlidePlaneKey<3> prismaticPlaneKey(prismaticPlaneIndex, prismaticPlaneNormal);
                            std::shared_ptr<PeriodicGlidePlane<3>> prismaticGlidePlane(mg.ddBase.periodicGlidePlaneFactory.get(prismaticPlaneKey));

                            PlanePlaneIntersection<3> ppi(*glidePlane->referencePlane,*prismaticGlidePlane->referencePlane);
                            Line<3> ppLine(ppi.P,ppi.d);
                            //const VectorDimD P0(ppi.P+(planePoint-ppi.P).dot(ppi.d)*ppi.d);
                            // const VectorDimD P0(ppLine.snap(planePoint));

                            
                            PlaneLineIntersection<3> pliA(periodicFaceA->center(),periodicFaceA->outNormal(),P0,faceAshift);
                            PlaneLineIntersection<3> pliB(periodicFaceB->center(),periodicFaceB->outNormal(),P0,faceAshift);
                            
                            if(pliA.type==PlaneLineIntersection<3>::INCIDENT && pliB.type==PlaneLineIntersection<3>::INCIDENT)
                            {
                                const VectorDimD AB(pliB.P-pliA.P);
                                const VectorDimD AB_dir = AB/AB.norm();

                                if((AB-faceAshift).norm()<FLT_EPSILON)
                                {

                                    GlidePlaneKey<3> parallelGlidePlaneKey(planeIndex+ExtendFRSourceHeight, slipSystem.n);
                                    std::shared_ptr<PeriodicGlidePlane<3>> parallelglidePlane(mg.ddBase.periodicGlidePlaneFactory.get(parallelGlidePlaneKey));


                                    if(parallelglidePlane && prismaticGlidePlane)
                                    {
                                        
                                        std::map<VectorDimD,size_t,CompareVectorsByComponent<double,3,float>> uniqueNetworkNodeMap; // networkNodePosition->networkNodeID
                                        // The prismatic loop
                                        const int nShift(ExtendFRSourceLength / 2.0);
                                        // const int nShift(2);
                                        const VectorDimD lShift(-nShift*AB_dir);
                                        const VectorDimD rShift(nShift*AB_dir);
                                        const VectorDimD startPoint(0.5*(pliB.P+pliA.P));
                                        std::vector<VectorDimD> prismaticNodePos;
                                        // left and right node                                       
                                        const VectorDimD leftStartPoint(startPoint+lShift);
                                        const VectorDimD rightStartPoint(startPoint+rShift);
                                        prismaticNodePos.push_back(leftStartPoint); //HERE ADD N*AB
                                        prismaticNodePos.push_back(rightStartPoint);
                                        prismaticNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(rightStartPoint));
                                        prismaticNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(leftStartPoint));
                                        
                                        const VectorDimD lineP0 = 0.5*(startPoint+lShift+parallelglidePlane->referencePlane->snapToPlane(startPoint+lShift));
                                        const VectorDimD lineP1 = 0.5*(startPoint+rShift+parallelglidePlane->referencePlane->snapToPlane(startPoint+rShift));
                                        
                                        mg.insertJunctionLoop(prismaticNodePos,prismaticGlidePlane,
                                                           slipSystem.s.cartesian(),prismaticGlidePlane->referencePlane->unitNormal,
                                                           P0,grainID,DislocationLoopIO<3>::SESSILELOOP);
                                        
                                        if(std::fabs(glideStep)>FLT_EPSILON)
                                        {
                                            const int halfNodeNumber(ExtendFRSourceNodes/2);
                                            if(halfNodeNumber < 0)
                                            {
                                                throw std::runtime_error("ExtendFRSourceNodes="+std::to_string(halfNodeNumber)+"is smaller than 0");
                                            }
                                            // const double halfDipoleLength(AB.norm()*0.5);
                                            // const VectorDimD dipoleDir(AB/AB.norm());
                                            // const VectorDimD shiftNodeLength(halfDipoleLength/(1.0*(halfNodeNumber+1))*dipoleDir);
                                            
                                            // First glide loop
                                            //                                    const double glideStep=50.0;
                                            std::vector<VectorDimD> firstNodePos;
                                            firstNodePos.push_back(leftStartPoint);
                                            firstNodePos.push_back(rightStartPoint);
                                            if((rShift-lShift).cross(glideStep*prismaticGlidePlane->referencePlane->unitNormal).dot(glidePlane->referencePlane->unitNormal)>0.0)
                                            {// loop is right-handed to glidePlane->referencePlane->unitNormal
                                                glideStep*=-1;
                                            }

                                            const VectorDimD middlePoint(startPoint+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                            const VectorDimD leftdirVec(middlePoint-leftStartPoint);
                                            const VectorDimD leftdir = leftdirVec/leftdirVec.norm();

                                            const VectorDimD rightdirVec(rightStartPoint-middlePoint);
                                            const VectorDimD rightdir = rightdirVec/rightdirVec.norm();

                                            const double halfExtendedFRSLength = sqrt((ExtendFRSourceLength * 0.5) * (ExtendFRSourceLength * 0.5) + glideStep * glideStep);
                                            const VectorDimD shiftNodeLength_l(halfExtendedFRSLength/(1.0*(halfNodeNumber+1))*leftdir);
                                            const VectorDimD shiftNodeLength_r(halfExtendedFRSLength/(1.0*(halfNodeNumber+1))*rightdir);

                                            firstNodePos.push_back(startPoint+glideStep*prismaticGlidePlane->referencePlane->unitNormal);

                                            for(int k=1; k<=halfNodeNumber; k++)
                                            { // Add nodes on rigth arm
                                                firstNodePos.push_back(startPoint+rShift-k*shiftNodeLength_r);
                                            }
                                            for(int k=halfNodeNumber+1; k>0; k--)
                                            { // Add nodes on left arm
                                                // firstNodePos.push_back(startPoint+lShift+k*shiftNodeLength+glideStep*prismaticGlidePlane->referencePlane->unitNormal);
                                                firstNodePos.push_back(startPoint+lShift+k*shiftNodeLength_l);

                                            }
                                            
                                            mg.insertJunctionLoop(firstNodePos,glidePlane,
                                                                  -slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
                                                                  P0,grainID,DislocationLoopIO<3>::GLISSILELOOP);
                                            
                                            
                                            FiniteLineSegment<3> mirrowLine(lineP0,lineP1);
                                            
                                            // Second glide loop
                                            std::vector<VectorDimD> secondNodePos;
                                            for(const auto& pos : firstNodePos)
                                            {
                                                VectorDimD c=mirrowLine.snapToInfiniteLine(pos);
                                                secondNodePos.push_back(2.0*c-pos);
                                                //                                        secondNodePos.push_back(parallelglidePlane->referencePlane->snapToPlane(pos));
                                            }
                                            
                                            
                                            mg.insertJunctionLoop(secondNodePos,parallelglidePlane,
                                                                  slipSystem.s.cartesian(),parallelglidePlane->referencePlane->unitNormal,
                                                                  parallelglidePlane->referencePlane->snapToPlane(P0),grainID,DislocationLoopIO<3>::GLISSILELOOP);
                                            return true;
                                        }
                                        
                                    }
                                    else
                                    {
                                        std::cout<<"Cannot create glide planes"<<std::endl;
                                    }
                                }
                                else
                                {
                                    throw std::runtime_error("periodic line intersection with faces is not the face shift vector");
                                }
                            }
                            else
                            {
                                throw std::runtime_error("periodic line direction does not form an incident intersecitn with periodic faces");
                            }
                        }
                        else
                        {
                            throw std::runtime_error("planeNormal of slipSystem "+std::to_string(rSS)+" is not othogonal to facelatticeShift "+std::to_string(exitFaceID));
                        }
                    }
                    else
                    {
                        throw std::runtime_error("Mesh face "+std::to_string(exitFaceID)+" is not a periodic face.");
                    }
                }
                else
                {
                    throw std::runtime_error("Mesh face "+std::to_string(exitFaceID)+" not found.");
                }
            }
            else
            {
                throw std::runtime_error("slipSystem "+std::to_string(rSS)+" not found, skipping.");
            }
        }
        else
        {
            std::cout<<"Skipping slip system "<<rSS<<std::endl;
        }
        return false;
    }

}
#endif
