/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNode_h_
#define model_DislocationNode_h_

//#include <PlanarDislocationNode.h>

#include <DislocationDynamicsModule.h>
#include <PlaneLineIntersection.h>
#include <DislocationNetwork.h>
#ifndef NDEBUG
#define VerboseDislocationNode(N,x) if(this->network().verboseDislocationNode>=N){std::cout<<x;}
#else
#define VerboseDislocationNode(N,x)
#endif

namespace model
{
    template <int dim, short unsigned int corder>
    class DislocationNode : public NetworkNode<DislocationNode<dim,corder>>
    /*                   */,public SplineNode<DislocationNode<dim,corder>,dim,corder,Hermite>
    {
        
        public:
        typedef TypeTraits<DislocationLoop<dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef SplineNode<DislocationNode<dim,corder>,dim,corder,Hermite> SplineNodeType;
//        typedef ConfinedDislocationObject<dim> ConfinedDislocationObjectType;
        typedef typename TypeTraits<NetworkNodeType>::MeshLocation MeshLocation;
        typedef std::vector<VectorDim> VectorOfNormalsType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef std::set<const GlidePlaneType *> GlidePlaneContainerType;
        typedef PlanarMeshFace<dim> PlanarMeshFaceType;
        typedef std::set<const PlanarMeshFaceType *> PlanarMeshFaceContainerType;
        typedef Grain<dim> GrainType;
        typedef std::set<const GrainType*> GrainContainerType;
        typedef Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClimbVelocityScalarType;

        const Simplex<dim,dim>* p_Simplex;
        ClimbVelocityScalarType climbVelocityScalar;
        VectorDim velocity;
        VectorDim vOld;         //! The previous velocity vector of *this PlanarDislocationNode
        double velocityReductionCoeff;

        bool setNodalVelocityToZeroInX;
        bool setNodalVelocityToZeroInY;
        bool setNodalVelocityToZeroInZ;

//        std::shared_ptr<NetworkNodeType> virtualNode;
//        NetworkNodeType* const masterNode;
        static int totalCappedNodes; //gives the fraction of the capped velocities
        
        DislocationNode(LoopNetworkType* const,const VectorDim&,const VectorDim&,const ClimbVelocityScalarType&,const double&);
        DislocationNode(LoopNetworkType* const,const VectorDim&);
        ~DislocationNode();
        std::shared_ptr<DislocationNode> clone() const;
        const Simplex<dim,dim>* get_includingSimplex(const VectorDim&,const Simplex<dim,dim>* const) const;
        const Simplex<dim,dim>* includingSimplex() const;
        bool isMovableTo(const VectorDim&) const;
        const double& velocityReduction() const;
//        bool isVirtualBoundaryNode() const __attribute__ ((deprecated));
        bool isBoundaryNode() const;
        bool isGrainBoundaryNode() const;
        void updateGeometry();
        bool set_P(const VectorDim&);
        bool trySet_P(const VectorDim&);
        const VectorDim& get_V() const;
        std::set<LoopType*> sessileLoops() const;
        MeshLocation meshLocation() const;
        void set_V(const VectorDim& vNew,const bool& isClimbingStep);
        void projectVelocity(const bool& isClimbingStep);
        VectorDim invariantDirectionOfMotion() const;
//        const std::shared_ptr<NetworkNodeType>& virtualBoundaryNode() const;
        GlidePlaneContainerType glidePlanes() const;
        PlanarMeshFaceContainerType meshFaces() const;
        bool isOnExternalBoundary() const;
        bool isOnInternalBoundary() const;
        bool isOnBoundary() const;
        VectorDim bndNormal() const;
        VectorDim snapToGlidePlanesinPeriodic(const VectorDim& P) const;
        GrainContainerType grains() const;
        VectorDim climbDirection() const;
    };
    
}
#endif
