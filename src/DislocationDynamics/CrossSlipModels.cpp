/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CrossSlipModels_cpp
#define model_CrossSlipModels_cpp

#include <DislocationNetwork.h>
#include <CrossSlipModels.h>

namespace model
{

    template <typename DislocationNetworkType>
    BaseCrossSlipModel<DislocationNetworkType>::BaseCrossSlipModel(const DDtraitsIO& traitsIO) :
    /* init */ crossSlipDeg(TextFileParser(traitsIO.ddFile).readScalar<double>("crossSlipAngle_deg",true))
    /* init */,sinCrossSlip(std::sin(crossSlipDeg*std::numbers::pi/180.0))
    {
        
        if(crossSlipDeg<0.0 || crossSlipDeg>90.0)
        {
            throw std::runtime_error("crossSlipDeg must be 0.0<= crossSlipDeg <= 90.0");
        }
        
    }

    template <typename DislocationNetworkType>
    std::pair<bool,int> BaseCrossSlipModel<DislocationNetworkType>::isBaseCrossSlipLink(const NetworkLinkType& link) const
    {
        if(   !link.isBoundarySegment()
           && !link.isGrainBoundarySegment()
           && link.loopLinks().size()==1
           && link.unitDirection().cross(link.burgers().normalized()).norm()<=sinCrossSlip
           && link.chordLength()>2.0*link.network().networkRemesher.Lmin
           && link.grains().size()==1
           )
        {
            if((*link.loopLinks().begin())->loop->slipSystem())
            {
                const auto& grain(**link.grains().begin());
                return std::make_pair(true,grain.grainID);
            }
        }
        return std::make_pair(false,-1);
    }

    template <typename DislocationNetworkType>
    AthermalCrossSlipModel<DislocationNetworkType>::AthermalCrossSlipModel(const DDtraitsIO& traitsIO) :
    /* init */ BaseCrossSlipModel<DislocationNetworkType>(traitsIO)
    {
        std::cout<<greenBoldColor<<"Creating AthermalCrossSlipModel"<<defaultColor<<std::endl;
    }

    template <typename DislocationNetworkType>
    std::pair<bool,std::pair<int,int>> AthermalCrossSlipModel<DislocationNetworkType>::isCrossSlipLink(const NetworkLinkType& link)
    {
        
        const auto isBaseCSLink(this->isBaseCrossSlipLink(link));

        if(isBaseCSLink.first)
        {
            // Find slip system on which screw segment has highest glide force
            const auto loop((*link.loopLinks().begin())->loop);
            const auto slipSystem(loop->slipSystem());
            if(slipSystem)
            {
                const auto& grain(**link.grains().begin());
                const auto pki=link.pkIntegral();
                const double pkGlide=std::fabs((pki-pki.dot(slipSystem->unitNormal)*slipSystem->unitNormal).norm());

                std::map<double,int> ssMap;
                for(size_t s=0;s<grain.singleCrystal->slipSystems().size();++s)
                {
                    const auto& crossSlipSystem(grain.singleCrystal->slipSystems()[s]);
                    if((loop->flow()-crossSlipSystem->s).squaredNorm()==0)
                    {
                        const double pkCross=std::fabs((pki-pki.dot(crossSlipSystem->unitNormal)*crossSlipSystem->unitNormal).norm());
                        ssMap.emplace(pkCross,s);
                    }
                }
                
                if(ssMap.size())
                {
                    const auto& crossSlipSystem(grain.singleCrystal->slipSystems()[ssMap.rbegin()->second]);
                    if(slipSystem->unitNormal.cross(crossSlipSystem->unitNormal).squaredNorm()>FLT_EPSILON
                       && ssMap.rbegin()->first>1.1*pkGlide)
                    {// highest glide PK-force is not on current slip system
                        return std::make_pair(true,std::make_pair(grain.grainID,ssMap.rbegin()->second));
                    }
                }
            }
        }
        return std::make_pair(false,std::make_pair(-1,-1));
    }

    template <typename DislocationNetworkType>
    void AthermalCrossSlipModel<DislocationNetworkType>::addToCrossSlip(const NetworkLinkType& link,
                                                                        CrossSlipContainerType& crossSlipDeq)
    {
        const auto isCSLink(isCrossSlipLink(link));
        if(isCSLink.first)
        {
            crossSlipDeq.emplace_back(link.source,
                                      link.sink,
                                      isCSLink.second.first,
                                      isCSLink.second.second);
        }
    }

//// Added for FCC
//template <typename DislocationNetworkType>
//EscaigCrossSlipModelFCC<DislocationNetworkType>::EscaigCrossSlipModelFCC(const PolycrystallineMaterialBase& material,const DDtraitsIO& traitsIO) :
///* init */ BaseCrossSlipModel<DislocationNetworkType>(traitsIO)
///* init */,gen(rd())
///* init */,dist(0.0,1.0)
///* init */,kB_DD(TextFileParser(material.materialFile).readScalar<double>("kB_DD",true) ) // /material.b_SI
///* init */,lengthC(TextFileParser(material.materialFile).readScalar<double>("lengthC",true)) // *material.eV2J/(material.mu_SI*std::pow(material.b_SI,3)) )
///* init */,activeVol(TextFileParser(material.materialFile).readScalar<double>("activeVol",true)) // Units is b /std::pow(material.b_SI,3) )
///*init*/,Tcs(TextFileParser(material.materialFile).readScalar<double>("Tcs",true)) // Units in burgers and K
///*init*/,omegaA(TextFileParser(material.materialFile).readScalar<double>("omegaA",true)) // Units in burgers and K
///*init*/,Ea(TextFileParser(material.materialFile).readScalar<double>("Ea_SI",true)*material.eV2J/(material.mu_SI*std::pow(material.b_SI,3))) // input units is modelib2
///*init*/,Tau_3(TextFileParser(material.materialFile).readScalar<double>("Tau_3",true)) // input units is modelib2
//{
//    std::cout<<greenBoldColor<<"Creating EscaigCrossSlipModelFCC"<<defaultColor<<std::endl;
//}
//
//// Added for FCC
//template <typename DislocationNetworkType>
//double EscaigCrossSlipModelFCC<DislocationNetworkType>::crossSlipRateKernel(const int& k,const NetworkLinkType& link,const SlipSystem& css,
//                                                                            const double& w,
//                                                                            const double& L,
//                                                                            const double& dE,
//                                                                            const double& Veff
//                                                                           ) const
//{ /*!@param[in] k the current quadrature point
//   *
//   */
//    const double eg(link.slipSystem()->unitNormal.transpose()*link.quadraturePoint(k).stress*link.slipSystem()->unitNormal.cross(link.slipSystem()->s.cartesian().normalized()));
//    const double ss(css.unitNormal.transpose()*link.quadraturePoint(k).stress*css.s.cartesian().normalized());
//    const double ec(css.unitNormal.transpose()*link.quadraturePoint(k).stress*css.unitNormal.cross(css.s.cartesian().normalized()));
//    // const double dG(dE-Veff*(-eg + ec)); // deltaSigma
//    const double dG(dE-std::fabs(ss)*Veff+eg*Veff-ec*Veff);
//
//    return w*std::exp(-dG/link.network().poly.kB/link.network().poly.T)*link.quadraturePoint(k).j/L;
//}
//
//// Added for FCC
//template <typename DislocationNetworkType>
//std::pair<bool,std::pair<int,int>> EscaigCrossSlipModelFCC<DislocationNetworkType>::isCrossSlipLink(const NetworkLinkType& link)
//{
//
//    const auto isBaseCSLink(this->isBaseCrossSlipLink(link));
//    if(isBaseCSLink.first)
//    {
//        const auto& grain(**link.grains().begin());
//        std::map<double,int> ssMap;
//        for(size_t s=0;s<grain.singleCrystal->slipSystems().size();++s)
//        {
//            const auto& slipSystem(grain.singleCrystal->slipSystems()[s]);
//            if(((*link.loopLinks().begin())->loop->slipSystem()->s-slipSystem->s).squaredNorm()==0
//               && slipSystem->unitNormal.cross(link.glidePlaneNormal()).squaredNorm()>FLT_EPSILON) // not current slip system
//            {
//                if(link.slipSystem()->mobility && slipSystem->mobility)
//                {
//                    
//                    double csRate(0.0);
//                    NetworkLinkType::QuadratureDynamicType::integrate(link.quadraturePoints().size(),this,csRate,&EscaigCrossSlipModelFCC<DislocationNetworkType>::crossSlipRateKernel,link,*slipSystem,
//                                                                          omegaA,lengthC,Ea,activeVol);
//                       
//                    const double csProb(1.0-std::exp(-csRate*link.network().simulationParameters.dt));
//                    ssMap.emplace(csProb,s);
//                }
//            }
//        }
//        
//        if(ssMap.size())
//        {
//            const double r(dist(gen));
//            if(ssMap.rbegin()->first > r)
//            {
//                return std::make_pair(true,std::make_pair(grain.grainID,ssMap.rbegin()->second));
//            }
//        }
//    }
//    return std::make_pair(false,std::make_pair(-1,-1));
//}
//
//
//template <typename DislocationNetworkType>
//void EscaigCrossSlipModelFCC<DislocationNetworkType>::addToCrossSlip(const NetworkLinkType& link,
//                                                                     CrossSlipContainerType& crossSlipDeq)
//{
//    const auto isCSLink(isCrossSlipLink(link));
//    if(isCSLink.first)
//    {
//        crossSlipDeq.emplace_back(link.source,
//                                  link.sink,
//                                  isCSLink.second.first,
//                                  isCSLink.second.second);
//    }
//}

    template <typename DislocationNetworkType>
    EscaigCrossSlipModelHEX<DislocationNetworkType>::EscaigCrossSlipModelHEX(const PolycrystallineMaterialBase& material,const DDtraitsIO& traitsIO) :
    /* init */ BaseCrossSlipModel<DislocationNetworkType>(traitsIO)
    /* init */,gen(rd())
    /* init */,dist(0.0,1.0)
    /* init */,enableBasalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enableBasalSlipSystems",false))
    /* init */,enablePrismaticSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePrismaticSlipSystems",false))
    /* init */,enablePyramidalSlipSystems(TextFileParser(material.materialFile).readScalar<int>("enablePyramidalSlipSystems",false))
    /* init */,csw(TextFileParser(material.materialFile).readScalar<double>("csw_SI",true)*material.b_SI/material.cs_SI)
    /* init */,csLb(enableBasalSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csLb_SI",true)/material.b_SI : 0.0)
    /* init */,csLp(enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csLp_SI",true)/material.b_SI : 0.0)
    /* init */,csb2pE(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pE_eV",true)*material.eV2J/(material.mu_SI*std::pow(material.b_SI,3)) : 0.0)
    /* init */,csb2pVeg(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pVeg_SI",true)/std::pow(material.b_SI,3) : 0.0)
    /* init */,csb2pVsc(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pVsc_SI",true)/std::pow(material.b_SI,3) : 0.0)
    /* init */,csb2pVec(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csb2pVec_SI",true)/std::pow(material.b_SI,3) : 0.0)
    /* init */,csp2bE(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bE_eV",true)*material.eV2J/(material.mu_SI*std::pow(material.b_SI,3)) : 0.0)
    /* init */,csp2bVeg(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bVeg_SI",true)/std::pow(material.b_SI,3) : 0.0)
    /* init */,csp2bVsc(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bVsc_SI",true)/std::pow(material.b_SI,3) : 0.0)
    /* init */,csp2bVec(enableBasalSlipSystems && enablePrismaticSlipSystems? TextFileParser(material.materialFile).readScalar<double>("csp2bVec_SI",true)/std::pow(material.b_SI,3) : 0.0)
    {
        std::cout<<greenBoldColor<<"Creating EscaigCrossSlipModelHEX"<<defaultColor<<std::endl;
    }

    template <typename DislocationNetworkType>
    double EscaigCrossSlipModelHEX<DislocationNetworkType>::crossSlipRateKernel(const int& k,const NetworkLinkType& link,const SlipSystem& css,
                                                                                const double& w,
                                                                                const double& L,
                                                                                const double& dE,
                                                                                const double& Veg,
                                                                                const double& Vsc,
                                                                                const double& Vec) const
    { /*!@param[in] k the current quadrature point
       *
       */
        const double eg(link.slipSystem()->unitNormal.transpose()*link.quadraturePoint(k).stress*link.slipSystem()->unitNormal.cross(link.slipSystem()->s.cartesian().normalized()));
        const double ss(css.unitNormal.transpose()*link.quadraturePoint(k).stress*css.s.cartesian().normalized());
        const double ec(css.unitNormal.transpose()*link.quadraturePoint(k).stress*css.unitNormal.cross(css.s.cartesian().normalized()));
        const double dG(dE-std::fabs(ss)*Vsc+eg*Veg-ec*Vec);
        return w*std::exp(-dG/link.network().ddBase.poly.kB/link.network().ddBase.poly.T)*link.quadraturePoint(k).j/L;
    }

    template <typename DislocationNetworkType>
    std::pair<bool,std::pair<int,int>> EscaigCrossSlipModelHEX<DislocationNetworkType>::isCrossSlipLink(const NetworkLinkType& link)
    {
        const auto isBaseCSLink(this->isBaseCrossSlipLink(link));
        if(isBaseCSLink.first)
        {
            const auto& grain(**link.grains().begin());
            std::map<double,int> ssMap;
            for(size_t s=0;s<grain.singleCrystal->slipSystems().size();++s)
            {
                const auto& slipSystem(grain.singleCrystal->slipSystems()[s]);
                if(((*link.loopLinks().begin())->loop->slipSystem()->s-slipSystem->s).squaredNorm()==0
                   && slipSystem->unitNormal.cross(link.glidePlaneNormal()).squaredNorm()>FLT_EPSILON) // not current slip system
                {
                    
                    if(link.slipSystem()->mobility && slipSystem->mobility)
                    {
                        
                        double csRate(0.0);
                        if(   link.slipSystem()->mobility->name.find("basal") != std::string::npos
                           && slipSystem->mobility->name.find("prismatic") != std::string::npos)
                        {// basal to prismatic crossSlip
                            NetworkLinkType::QuadratureDynamicType::integrate(link.quadraturePoints().size(),this,csRate,&EscaigCrossSlipModelHEX<DislocationNetworkType>::crossSlipRateKernel,link,*slipSystem,
                                                                              csw,csLb,csb2pE,csb2pVeg,csb2pVsc,csb2pVec);
                        }
                        else if(   link.slipSystem()->mobility->name.find("prismatic") != std::string::npos
                                && slipSystem->mobility->name.find("basal") != std::string::npos)
                        {// prismatic to basal  crossSlip
                            NetworkLinkType::QuadratureDynamicType::integrate(link.quadraturePoints().size(),this,csRate,&EscaigCrossSlipModelHEX<DislocationNetworkType>::crossSlipRateKernel,link,*slipSystem,
                                                                              csw,csLp,csp2bE,csp2bVeg,csp2bVsc,csp2bVec);
                        }
                        const double csProb(1.0-std::exp(-csRate*link.network().ddBase.simulationParameters.dt));
                        ssMap.emplace(csProb,s);
                    }
                }
            }
            
            if(ssMap.size())
            {
                const double r(dist(gen));
                if(ssMap.rbegin()->first > r)
                {
                    return std::make_pair(true,std::make_pair(grain.grainID,ssMap.rbegin()->second));
                }
            }
        }
        return std::make_pair(false,std::make_pair(-1,-1));
    }

    template <typename DislocationNetworkType>
    void EscaigCrossSlipModelHEX<DislocationNetworkType>::addToCrossSlip(const NetworkLinkType& link,
                                                                         CrossSlipContainerType& crossSlipDeq)
    {
        
        const auto isCSLink(isCrossSlipLink(link));
        if(isCSLink.first)
        {
            crossSlipDeq.emplace_back(link.source,
                                      link.sink,
                                      isCSLink.second.first,
                                      isCSLink.second.second);
        }
    }

    template struct AthermalCrossSlipModel<DislocationNetwork<3,0>>;
    template struct EscaigCrossSlipModelHEX<DislocationNetwork<3,0>>;

}
#endif
