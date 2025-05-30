/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Loop_H_
#define model_Loop_H_

#include <iostream>
#include <list>
#include <deque>
#include <set>
#include <memory>
#include <iterator>
#include <algorithm>
#include <StaticID.h>
#include <CRTP.h>
#include <NetworkBase.h>

#ifndef NDEBUG
#define VerboseLoop(N,x) if(verboseLevel>=N){std::cout<<cyanColor<<x<<defaultColor;}
#else
#define VerboseLoop(N,x)
#endif

namespace model
{
    
    
    template<typename Derived>
    class Loop : public CRTP<Derived>
    /*        */,public StaticID<Derived>
    /*        */,public NetworkBase<Derived,size_t>
    /*        */,private std::set<typename TypeTraits<Derived>::LoopLinkType*>
    {
        
    public:

        
        typedef typename TypeTraits<Derived>::LoopNetworkType LoopNetworkType;
        typedef NetworkBase<Derived,size_t> NetworkBaseType;
        typedef typename TypeTraits<Derived>::LoopNodeType LoopNodeType;
        typedef typename TypeTraits<Derived>::LoopLinkType LoopLinkType;
        typedef std::set<LoopLinkType*> LoopLinkContainerType;
        typedef std::deque<const LoopLinkType*> LoopLinkSequenceType;
        typedef std::deque<std::pair<std::shared_ptr<LoopNodeType>,std::shared_ptr<LoopNodeType>>> LoopNodeSequenceType;
        typedef typename TypeTraits<LoopLinkType>::FlowType FlowType;
        
    private:
        
        FlowType _flow;

    public:
        
        static int verboseLevel;


        
        Loop(LoopNetworkType* const loopNetwork_in,
             const FlowType& f) :
        /* init */ NetworkBaseType(loopNetwork_in,&loopNetwork_in->loops(),this->sID)
        /* init */,_flow(f)
        {
            VerboseLoop(1,"Creating Loop "<<tag()<<std::endl);
        }
        
        ~Loop()
        {
            VerboseLoop(1,"Destroying Loop "<<tag()<<std::endl);
        }
        
        void update()
        {
            
        }
        
        const FlowType& flow() const
        {
            return _flow;
        }
        
        LoopLinkContainerType& loopLinks()
        {
            return *this;
        }
        
        const LoopLinkContainerType& loopLinks() const
        {
            return *this;
        }
        
        void addLoopLink(LoopLinkType* const pL)
        {
            VerboseLoop(2,"Loop "<<tag()<<" addLoopLink "<<pL->tag()<<std::endl);
            const bool success=loopLinks().insert(pL).second;
            if(!success)
            {
                std::cout<<"DislocationLoop "<<this->sID<<" cannot add LoopLink "<<pL->tag()<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        void removeLoopLink(LoopLinkType* const pL)
        {
            VerboseLoop(2,"Loop "<<tag()<<" removeLoopLink "<<pL->tag()<<std::endl);
            const size_t erased=loopLinks().erase(pL);
            if(erased!=1)
            {
                throw std::runtime_error("Cannot remove link "+pL->tag()+" from loop "+tag());
            }
        }
        
        LoopLinkSequenceType linkSequence() const
        {
            //RECODE THIS USING prev/next
            //typename LoopLinkContainerType::const_iterator iter;
            VerboseLoop(3,"Loop "<<tag()<<" loopLinks().size()= "<<loopLinks().size()<<std::endl);
            LoopLinkSequenceType temp;
            if(loopLinks().size())
            {
                //Start at the link starting at the minimum Node
                std::map<size_t,const LoopLinkType*> orderedNodeIDandLinks;
                for (const auto& link : loopLinks())
                {
                    orderedNodeIDandLinks.emplace(link->source->sID,link);
                }
                const LoopLinkType* pL=orderedNodeIDandLinks.begin()->second;
                for(size_t k=0;k<loopLinks().size();++k)
                {
                    if(pL)
                    {
                        if(pL->next)
                        {
                            temp.push_back(pL);
                            pL=pL->next;
                        }
                    }
                }
            }
            VerboseLoop(3,"Loop "<<tag()<<" linkSequence.size()= "<<temp.size()<<std::endl);
            return temp;
        }
        
        LoopNodeSequenceType nodeSequence() const
        {
            LoopNodeSequenceType temp;
            for(const auto& link : linkSequence())
            {
                temp.emplace_back(link->source,link->sink);
            }
            return temp;
        }
        
        std::deque<size_t> nodeIDSequence() const
        {
            std::deque<size_t> temp;
            for(const auto& link : linkSequence())
            {
                temp.emplace_back(link->source->sID);
            }
            return temp;
        }
        
        std::pair<bool,LoopLinkType*> linkStartingAt(const std::shared_ptr<LoopNodeType>& Ni) const
        {
            
            std::pair<bool,LoopLinkType*> temp=std::make_pair(false,static_cast<LoopLinkType*>(nullptr));
            for(const auto& link : loopLinks())
            {
                if(link->source==Ni)
                {
                    temp=std::make_pair(true,link);
                    break;
                }
            }
            return temp;
        }

        bool isLoop() const
        {
            const LoopLinkSequenceType linkSeq(linkSequence());
            bool temp(linkSeq.size()>=3);
            for(typename LoopLinkSequenceType::const_iterator iter=linkSeq.begin();iter!=linkSeq.end();++iter)
            {
                auto next=std::next(iter,1);
                if(next==linkSeq.end())
                {
                    next=linkSeq.begin();
                }
                
                temp=(temp && ((*iter)->sink->sID==(*next)->source->sID));
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }

        bool isIsolated() const
        {
            for(const auto& link : loopLinks())
            {
                if(link->networkLink())
                {
                    if(link->networkLink()->loopLinks().size()!=1)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        
        std::string tag() const
        {/*!\returns the string "i->j" where i is source()->sID and j=sink()->sID
          */
            return std::to_string(this->sID);
        }
        
        void printLoop() const
        {
            std::cout<<"Loop "<<this->sID<<std::endl;
            const LoopLinkSequenceType linkSeq(linkSequence());
            for(const auto& link : linkSeq)
            {
                std::cout<<"    "<<link->tag()
                <<" (prev "<<link->prev->tag()<<")"
                <<" (next "<<link->next->tag()<<")"<<std::endl;
            }
        }
        
    };
    
    template<typename Derived>
    int Loop<Derived>::verboseLevel=0;

}
#endif
