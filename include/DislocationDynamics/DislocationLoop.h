/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_h_
#define model_DislocationLoop_h_


#include <DislocationDynamicsModule.h>
#include <MeshedDislocationLoop.h>


#ifndef NDEBUG
#define VerboseDislocationLoop(N,x) if(verboseDislocationLoop>=N){std::cout<<x;}
#else
#define VerboseDislocationLoop(N,x)
#endif

namespace model
{
    template <int _dim, short unsigned int corder>
    class DislocationLoop : public Loop<DislocationLoop<_dim,corder>>
    {

    public:
        
        typedef TypeTraits<DislocationLoop<_dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::MatrixDim MatrixDim;
        typedef typename TraitsType::GrainType GrainType;
        typedef typename TraitsType::GlidePlaneType GlidePlaneType;
        typedef typename TraitsType::PeriodicGlidePlaneType PeriodicGlidePlaneType;
        typedef typename TraitsType::ReciprocalLatticeDirectionType ReciprocalLatticeDirectionType;
        static int verboseDislocationLoop;

        const std::shared_ptr<GlidePlaneType> glidePlane;
        const std::shared_ptr<PeriodicGlidePlaneType> periodicGlidePlane;
        DislocationLoopPatches<_dim> _patches;
        const GrainType& grain;
        const int loopType;

    private:
        
        VectorDim nA;
        VectorDim nAR;
        double _slippedArea;
        double _slippedAreaRate;
        VectorDim _rightHandedUnitNormal;
        VectorDim _rightHandedUnitNormal_old;
        ReciprocalLatticeDirectionType _rightHandedNormal;
        std::shared_ptr<SlipSystem> _slipSystem;

    public:
        
        DislocationLoop(LoopNetworkType* const,
                        const VectorDim&,
                        const std::shared_ptr<GlidePlaneType>& glidePlane_in);
        ~DislocationLoop();
        std::shared_ptr<LoopType> clone() const;
        const double& slippedArea() const;
        const double& slippedAreaRate() const;
        const VectorDim& rightHandedUnitNormal() const;
        const VectorDim& rightHandedUnitNormalOld() const;
        const ReciprocalLatticeDirectionType& rightHandedNormal() const;
        std::tuple<double,double,double,double> loopLength() const;
        std::shared_ptr<SlipSystem> searchSlipSystem() const;
        void updateSlipSystem();
        void updateGeometry();
        void updateRates();
        MatrixDim averagePlasticDistortion() const;
        MatrixDim averagePlasticDistortionRate() const;
        VectorDim burgers() const;
        const std::shared_ptr<SlipSystem>&  slipSystem() const;
        bool isVirtualBoundaryLoop() const;
        double solidAngle(const VectorDim& x) const;
        void crossSlipBranches(std::deque<std::pair<std::deque<std::shared_ptr<LoopNodeType>>,int>>& csNodes) const;
        const DislocationLoopPatches<_dim>& patches() const;
        std::vector<MeshedDislocationLoop> meshed(const double& meshSize,const double& localMeshSize) const;
        static void initFromFile(const std::string&);
        static double planarSolidAngle(const VectorDim& x,const VectorDim& planePoint,const VectorDim& rhN,const std::vector<std::pair<VectorDim,VectorDim>>& polygonSegments);
        template <typename T> static int sgn(const T& val);
   };
    
}
#endif
