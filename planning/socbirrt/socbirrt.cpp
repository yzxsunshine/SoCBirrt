/*
   This code is very early version, by directly and ugly modify the source code from CoMPs
*/
/* Copyright (c) 2010 Carnegie Mellon University and Intel Corporation
   Author: Dmitry Berenson <dberenso@cs.cmu.edu>

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
     * Neither the name of Intel Corporation nor Carnegie Mellon University,
       nor the names of their contributors, may be used to endorse or
       promote products derived from this software without specific prior
       written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL INTEL CORPORATION OR CARNEGIE MELLON
   UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
   OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/** \file cbirrt.cpp
    \brief Implements the cbirrt planner class.
 */
#include "stdafx.h"

const dReal STEP_LENGTH        = 0.05; /// step length in radians
const int   MAX_NUM_ATTACHEDIKSOLS = 30; /// if attached IK solver gives too many solutions for a single query, sample this many from the set

const int   UNDEFINED          = -1;
const dReal INF                = 1000000.0;
const dReal TINY               = 0.001; // used when checking equality of floats/doubles

//the following are defaults, they can be ignored depending on external input
const int   NUM_OPT_ITERATIONS = 300;  // optimization iterations
const double MAX_TIME           = 25.0; // in seconds
const double MAX_FIRSTIK_TIME   = 5.0; // in seconds

//#define TRACK_COLLISIONS //comment this in if you want to create a file with collision information
//#define RECORD_TIMES //comment this in if you want to record how much time the components of the algorithm take

void  PrintMatrix(dReal* pMatrix, int numrows, int numcols, const char * statement)
{
    stringstream s;
    s.setf(ios::fixed,ios::floatfield);
    s.precision(5);
    s << "\n"; 
    if(statement != NULL)
        s << statement <<"\n";
    for(int i = 0; i < numrows; i++)
    {
        for(int j = 0; j < numcols; j++)
        {
            s << pMatrix[i*numcols + j] << " \t";
            if(fabs(pMatrix[i*numcols + j]) > 10000)
                RAVELOG_INFO("WARNING: %d %d is huge!!!\n",i,j);
        }
        s << endl;
    }
    RAVELOG_INFO(s.str().c_str());

}


std::vector<RrtNode>* SoCBirrtPlanner::treenodes;

SoCBirrtPlanner::~SoCBirrtPlanner()
{
  
}

//clean up and return
OpenRAVE::PlannerStatus SoCBirrtPlanner::CleanUpReturn(bool retval)
{
    if(bdelete_distmetric)
        delete pdistmetric;

    delete _pForwardTree->_pMakeNext;
    delete _pBackwardTree->_pMakeNext;

    _pForwardTree->DeleteNodes();
    _pBackwardTree->DeleteNodes();

#ifdef TRACK_COLLISIONS
     colfile.close();
#endif
    
    if(retval)
        return PS_HasSolution;
    else
        return PS_Failed;

}


bool SoCBirrtPlanner::InitPlan(RobotBasePtr  pbase, PlannerParametersConstPtr pparams)
{
#ifdef TRACK_COLLISIONS
        colfile.open("cols.txt",ios::app);
#endif

    RAVELOG_DEBUG("Initializing Planner\n");
    if( pparams.get() != NULL )
    {
        _parameters.reset(new SoCBirrtParameters());
        _parameters->copy(pparams);
        //_parameters = pparams;
    }
    else
    {
        RAVELOG_DEBUG("CBirrt::InitPlan - Error: No parameters passed in to initialization");
        _outputstream << "Error: No parameters passed in to initialization\n";
        return false;
    }
    
    _pRobot = pbase;
    for(int i = 0; i < _parameters->vTSRChains.size(); i++)
        _parameters->vTSRChains[i].Initialize(GetEnv());

    RAVELOG_DEBUG("RrtPlanner: Using Default Distance Function\n");
    bdelete_distmetric = true;
    pdistmetric = new SoCBirrtDistanceMetric(this);
    

    if( _pRobot.get() == NULL || pdistmetric == NULL ) {
        _pRobot.reset();
        return false;
    }

    max_planning_time = MAX_TIME;
    max_firstik_time = MAX_FIRSTIK_TIME;
    if(_parameters->timelimit > 0)
    {
        max_planning_time = _parameters->timelimit;
        if(max_planning_time < max_firstik_time)
            max_firstik_time = max_planning_time;
    }
   
    
    _numdofs = _pRobot->GetActiveDOF();
    // initialize the jointResolutions
    _pRobot->GetActiveDOFResolutions(_jointResolution);
    //get the limits
    _pRobot->GetActiveDOFLimits(_lowerLimit,_upperLimit);
    
    
    _viscircular.resize(GetNumDOF(),0);
    for(size_t i = 0; i < _pRobot->GetActiveDOFIndices().size(); ++i) {
       int dof = _pRobot->GetActiveDOFIndices().at(i);
       KinBody::JointPtr pjoint = _pRobot->GetJointFromDOFIndex(dof);
       _viscircular[i] = pjoint->IsCircular(dof-pjoint->GetDOFIndex());
    }

    //need to know mapping between real joint indices and active indices in case there is an attached IK solver
    vJointIndToActiveInd.resize(_pRobot->GetDOF(),-1);
    for(int i = 0; i < _pRobot->GetActiveDOF(); i++)
    {
        //check for affine dofs
        if(_pRobot->GetActiveDOFIndices()[i] >= 0)
            vJointIndToActiveInd[_pRobot->GetActiveDOFIndices()[i]] = i;
    }

    int numTSRMimicDOF = 0;
    int old_numTSRMimicDOF = 0;
    int numTSRTotalDOF = 0;
    int old_numTSRTotalDOF = 0;

    vTSRChainIndexToConfigurationOffset.resize(_parameters->vTSRChains.size(),0);
    vTSRChainIndexToGuessOffset.resize(_parameters->vTSRChains.size(),0);

    vvManipToConstraintTSRChainind.resize(_pRobot->GetManipulators().size());
    vvManipToGoalSamplingTSRChainind.resize(_pRobot->GetManipulators().size());
    vvManipToStartSamplingTSRChainind.resize(_pRobot->GetManipulators().size());

    vMimicBodies.clear();

    bool bDoingConstraint = false;
    bool bDoingSampling = false;
    RobotBasePtr probottemp;
    printf("Chain Started\n");
    for(int i = 0; i < _parameters->vTSRChains.size();i++)
    {
        if(!_parameters->vTSRChains[i].RobotizeTSRChain(GetEnv(),probottemp))
        {
            RAVELOG_FATAL("Unable to robotize TSR for TSR Chain %d, cannot initialize planner.\n",i);
            _outputstream << "Unable to robotize TSR for TSR Chain " << i << ", cannot initialize planner\n";
            return false;
        }

        if(probottemp.get() != NULL)
            sIgnoredBodies.push_back(probottemp);//ignore TSR robot for collision checks

        if(_parameters->vTSRChains[i].IsForConstraint())
        {
            if(_parameters->vTSRChains[i].GetMimicBody() != NULL)
                vMimicBodies.push_back(_parameters->vTSRChains[i].GetMimicBody());

            vvManipToConstraintTSRChainind[_parameters->vTSRChains[i].GetManipInd()].push_back(i);
            old_numTSRMimicDOF = numTSRMimicDOF;
            numTSRMimicDOF += _parameters->vTSRChains[i].GetNumMimicDOF();

            if(numTSRMimicDOF > old_numTSRMimicDOF)
            {
                _numdofs = _pRobot->GetActiveDOF()+numTSRMimicDOF;
                _lowerLimit.resize(GetNumDOF(),0);
                _upperLimit.resize(GetNumDOF(),0);

                vTSRChainIndexToConfigurationOffset[i] =_pRobot->GetActiveDOF()+old_numTSRMimicDOF;

                _parameters->vTSRChains[i].GetChainJointLimits(&_lowerLimit[_pRobot->GetActiveDOF()+old_numTSRMimicDOF],&_upperLimit[_pRobot->GetActiveDOF()+old_numTSRMimicDOF]);
            }
            else
            {
                //this is a TSR chain with 0 dofs, i.e. a point tsr, so set this to point to nowhere
                vTSRChainIndexToConfigurationOffset[i] = 0;
            }


            old_numTSRTotalDOF = numTSRTotalDOF;
            numTSRTotalDOF += _parameters->vTSRChains[i].GetNumDOF();
    
            if(numTSRTotalDOF > old_numTSRTotalDOF)
            {
                vTSRChainIndexToGuessOffset[i] = old_numTSRTotalDOF;
            }
            else
            {
                //this is a TSR chain with 0 dofs, i.e. a point tsr, so set this to point to nowhere
                vTSRChainIndexToGuessOffset[i] = 0;
            }
            bDoingConstraint = true;
        }

        if(_parameters->vTSRChains[i].IsForGoalSampling())
        {
            vvManipToGoalSamplingTSRChainind[_parameters->vTSRChains[i].GetManipInd()].push_back(i);
            bDoingSampling = true;
        }

        if(_parameters->vTSRChains[i].IsForStartSampling())
        {
            vvManipToStartSamplingTSRChainind[_parameters->vTSRChains[i].GetManipInd()].push_back(i);
            bDoingSampling = true;
        }


    }
    _jointResolution.resize(GetNumDOF(),1);
    vTSRChainValues_temp.resize(numTSRTotalDOF);

    RAVELOG_INFO("numdofs: %d\n",_numdofs);
    printf("numdofs: %d\n",_numdofs);
    // invert for speed
    _jointResolutionInv.resize(_jointResolution.size());
    for (int i = 0; i < GetNumDOF(); i++) {
        _jointResolutionInv[i] = (dReal)1.0/_jointResolution[i];
        RAVELOG_DEBUG("   joint Resolution %d : %f   joint Resolution inv: %f\n", i, _jointResolution[i], _jointResolutionInv[i]);
    }

    P_SAMPLE_IK = _parameters->Psample;
    if(P_SAMPLE_IK == 0 && bDoingSampling)
    {
        RAVELOG_INFO("ERROR: there are chains defined for sampling but psample is 0.\n");
        _outputstream << "ERROR: there are chains defined for sampling but psample is 0\n";
        return false;
    }
    else if(P_SAMPLE_IK > 0 && !bDoingSampling)
    {
        RAVELOG_INFO("ERROR: psample is set to >0 but there are no chains defined for sampling.\n");
        _outputstream << "ERROR: psample is set to >0 but there are no chains defined for sampling\n";
        return false;
    }

    if(_parameters->bgrabbed)
    {
        std::vector<KinBodyPtr> vbodies;
        _pRobot->GetGrabbed(vbodies);
        pheldobject = vbodies[0];
    }

    checkConfig.resize(GetNumDOF());

    _randomConfig.resize(GetNumDOF());
    _validRange.resize(GetNumDOF());
    _jointIncrement.resize(GetNumDOF());

  
    for (int i = 0; i < GetNumDOF(); i++) {
        _validRange[i] = _upperLimit[i] - _lowerLimit[i];
        assert(_validRange[i] > 0);
    }

    _pForwardTree->_pMakeNext = new MakeNext(_pForwardTree->GetFromGoal(),GetNumDOF(),_pRobot,this);
    _pBackwardTree->_pMakeNext = new MakeNext(_pBackwardTree->GetFromGoal(),GetNumDOF(),_pRobot,this);

    if(treenodes != NULL && treenodes->size() > 0) {
    	_pForwardTree->AddNodes(treenodes);
    }
    printf("Initial nodes: t-%d; f-%d; b-%d\n", treenodes->size(), _pForwardTree->GetSize(), _pBackwardTree->GetSize());
    RAVELOG_INFO("grabbed: %d\n",_parameters->bgrabbed);

    bSmoothPath = _parameters->bsmoothpath;


    _pIkSolver = RaveCreateIkSolver(GetEnv(),"GeneralIK");
    _pIkSolver->Init(_pRobot->GetActiveManipulator());

    RAVELOG_INFO("Psample: %f\n",_parameters->Psample);


    _pInitConfig.resize(GetNumDOF());
    std::vector<dReal> tempjointvals;
    std::vector<dReal> tempconfig;
    bool bProjected;
    int start_index = 0;
    int num_starts = 0;
    int projectednode_id =0;
    if(_parameters->vinitialconfig.size() != 0)
    {
        while(1)
        {
            for(int i = 0 ; i < GetNumDOF()- numTSRMimicDOF; i++)
            {
                if(start_index < (int)_parameters->vinitialconfig.size())
                {
                  _pInitConfig[i] = _parameters->vinitialconfig[start_index];
                }
                else
                {
                    RAVELOG_INFO("CBirrtPlanner::InitPlan - Error: Starts are improperly specified.\n");
                    _outputstream << "Error: Starts are improperly specified\n";
                    return false;
                }

                start_index++;
            }



            stringstream s;
            tempjointvals = _pInitConfig;
            ClearTSRChainValues();
            tempconfig = _pInitConfig;
            bool bProjected = false;
            //WARNING: DISTANCE CHECKING FOR STARTCONFIG PROJECTION IS CURRENTLY DISABLED
            //This is because there is always some numerical error so can't meet lower-dimensional constraints exactly
            if(!_pForwardTree->_pMakeNext->ConstrainedNewConfig(_pInitConfig,tempjointvals,vTSRChainValues_temp,false))
            {
                s << "Start " << num_starts << ": ";
                RAVELOG_INFO("CBirrtPlanner::InitPlan - Error: Start configuration %d does not meet constraints:\n",num_starts);
                for(int j = 0 ; j < GetNumDOF(); j++)
                    s<< _pInitConfig[j] << " ";
                s << endl;
                RAVELOG_INFO(s.str().c_str());
                _outputstream << "Error: Start configuration "<< num_starts << " does not meet constraints\n";

                SetDOF(_pInitConfig);
                if(!_pForwardTree->_pMakeNext->CheckSupport())
                {
                    RAVELOG_INFO("CBirrtPlanner::InitPlan - Error: Start configuration %d is not in balance.\n",num_starts);
                    _outputstream << "Error: Start configuration "<< num_starts << " is not in balance\n";
                }

                return false;
            }

            if(_CheckCollision(_pInitConfig))
            {
                s << "Start " << num_starts << ": ";
                RAVELOG_INFO("CBirrtPlanner::InitPlan - Error: Start configuration in collision:\n");
                for(int j = 0 ; j < GetNumDOF(); j++)
                    s<< _pInitConfig[j] << " ";
                s << endl;
                RAVELOG_INFO(s.str().c_str());
                _outputstream << "Error: Start configuration in collision\n";
                return false;
            }


            // set up the start states
            _pActiveNode = new RrtNode(GetNumDOF(),false,num_starts+projectednode_id);
            _pActiveNode->SetConfig(&tempconfig[0]);
            _pActiveNode->SetTSRChainValues(vTSRChainValues_temp);
            _pForwardTree->AddNode(*_pActiveNode);
	    delete _pActiveNode;
	    RAVELOG_INFO("Checking for start projection...\n");
		for(int i = 0; i < GetNumDOF() - numTSRMimicDOF; i++)
		{
		   if(_pInitConfig[i] != tempconfig[i])
		   {
                RAVELOG_INFO("CBirrtPlanner::InitPlan - WARNING: Start configuration %d was projected to meet constraint:\n",num_starts);
                        for(int j = 0 ; j < GetNumDOF(); j++)
                            s<< _pInitConfig[j] << " ";
                        s << endl;
                        RAVELOG_INFO(s.str().c_str());
                        _outputstream << "WARNING: Start configuration "<< num_starts << " was projected to meet constraint:\n"<< s.str();
                bProjected = true;
                projectednode_id++;
			    _pActiveNode = new RrtNode(GetNumDOF(),false,num_starts+projectednode_id);
			    _pActiveNode->SetConfig(&_pInitConfig[0]);
			    _pActiveNode->SetTSRChainValues(vTSRChainValues_temp);
			    _pActiveNode->SetParent(num_starts+projectednode_id-1);
			    printf("[socbirrt.cpp-SoCBirrtPlanner::InitPlan] Add node to parent's child list\n");
			    _pForwardTree->GetNode(num_starts+projectednode_id-1)->AddChild(num_starts+projectednode_id);
			    printf("[socbirrt.cpp-SoCBirrtPlanner::InitPlan] Add current node to tree list\n");
			    _pForwardTree->AddNode(*_pActiveNode);
			    delete _pActiveNode;

			break;
		   }
		}



	    num_starts++;


            if(start_index == (int)_parameters->vinitialconfig.size())
                break;

        }
    }

    RAVELOG_INFO("Start Node(s) Created\n");

    // reset all RRT parameters
    _pConnectNode = NULL;

    if(_parameters->vgoalconfig.size() != 0)
    {
        _pGoalConfig.resize(GetNumDOF());
        RAVELOG_DEBUG("Num Active DOFs: %d\n",GetNumDOF());
        //read in all goals
        int goal_index = 0;
        int num_goals = 0;
        while(1)
        {
            for(int i = 0 ; i < GetNumDOF(); i++)
            {
                if(goal_index < (int)_parameters->vgoalconfig.size())
                    _pGoalConfig[i] = _parameters->vgoalconfig[goal_index];
                else
                {
                    RAVELOG_INFO("CBirrtPlanner::InitPlan - Error: goals are improperly specified.\n");
                    _outputstream << "Error: goals are improperly specified\n";
                    return false;
                }
                goal_index++;
            }


            stringstream s;


            tempjointvals = _pGoalConfig;
            ClearTSRChainValues();
            //don't care about distance check for goal config
            if(!_pBackwardTree->_pMakeNext->ConstrainedNewConfig(_pGoalConfig,tempjointvals,vTSRChainValues_temp,false))
            {
                s << "Goal " << num_goals << ": ";
                RAVELOG_INFO("CBirrtPlanner::InitPlan - Error: Goal configuration does not meet constraints:\n");
                for(int j = 0 ; j < GetNumDOF(); j++)
                    s<< _pGoalConfig[j] << " ";
                s << endl;
                RAVELOG_INFO(s.str().c_str());
                _outputstream << "Error: Goal configuration does not meet constraints:\n" << s.str();

                SetDOF(_pGoalConfig);
                if(!_pForwardTree->_pMakeNext->CheckSupport())
                {
                    RAVELOG_INFO("CBirrtPlanner::InitPlan - Error: Goal configuration %d is not in balance.\n",num_goals);
                    _outputstream << "Error: Goal configuration "<< num_goals << " is not in balance\n";
                }

                return false;
            }

            if(_CheckCollision(_pGoalConfig))
            {
                s << "Skipping Goal " << num_goals << ": ";
                for(int j = 0 ; j < GetNumDOF(); j++)
                    s<< _pGoalConfig[j] << " ";
                s << endl;
                RAVELOG_INFO(s.str().c_str());
                RAVELOG_INFO(" because it is in collision:\n");
                _outputstream << s.str() << " because it is in collision:\n";
                if(goal_index == (int)_parameters->vgoalconfig.size())
                  break;
                else
                  continue;
                //return false;
            }

            // set up the goal states
            _pActiveNode = new RrtNode(GetNumDOF(),true,num_goals);
            _pActiveNode->SetConfig(&_pGoalConfig[0]);
            _pActiveNode->SetTSRChainValues(vTSRChainValues_temp);
            num_goals++;

            _pBackwardTree->AddNode(*_pActiveNode);

            delete _pActiveNode;

            if(goal_index == (int)_parameters->vgoalconfig.size())
                break;

        }

        if (num_goals !=0)
          RAVELOG_DEBUG("Goal State Node(s) Created\n");
        else
        {
          RAVELOG_DEBUG("Could not initialize Goal Configuration\n");
          _outputstream << "Could not initialize Goal Configuration\n";
          return false;
        }
    }

    RAVELOG_DEBUG("Initializaing Start State\n");

    //don't touch ikguess if one given, if not use start configuration, if no start use random (new random each time)
    if(_parameters->vikguess.size() == 0 && _pForwardTree->GetSize() > 0)
    {
        _parameters->vikguess = *_pForwardTree->GetNode(0)->GetDataVector();
        RAVELOG_DEBUG("Using first start configuration as IK guess.\n");
    }


    bInit = true;
    RAVELOG_DEBUG("RrtPlanner::InitPlan - RRT Planner Initialized\n");
    return true;
}


OpenRAVE::PlannerStatus SoCBirrtPlanner::PlanPath(TrajectoryBasePtr ptraj)
{
    if(!bInit)
    {
        RAVELOG_INFO("SoCBirrtPlanner::PlanPath - Error, planner not initialized\n");
        _outputstream << "Error, planner not initialized\n";
        return CleanUpReturn(false);
    }
    RAVELOG_DEBUG("Starting PlanPath\n");


    bool bContinuePlan = true;

    NodeTree* TreeA = _pForwardTree;
    NodeTree* TreeB = _pBackwardTree;
    NodeTree* ptreetemp = NULL;
    int iclosest = -1;
    int iConnectedA = 1;
    int iConnectedB = -1;
    std::vector<RrtNode> NewNodes;
    NewNodes.reserve(512);

    nntime = 0;
    projectiontime = 0;
    collisionchecktime = 0;
    nncalls = 0;
    collisioncheckcalls = 0;
    projectioncalls = 0;
    printf("[socbirrt.cpp-PlanPath-276] Start Path Planning\n");
    double starttime = timeGetThreadTime();
    int numitr = 0;

    // the main planning loop
    while(bContinuePlan)
    {
        //stop planning if you see PS_PlanFailed
        if(*(_parameters->pplannerstate) == PS_PlanFailed)
        {
            RAVELOG_INFO("Planning terminated from external function.\n");
            _outputstream << "Planning terminated from external function\n";
            return CleanUpReturn(false);
        }

        if(_pForwardTree->GetSize() == 0)
        {
            RAVELOG_DEBUG("Forward Tree is empty\n");
            printf("[socbirrt.cpp-PlanPath-294] Forward Tree is empty\n");
            while(1)
            {
                if(_parameters->vikguess.size() == 0)
                {
                    _PickRandomConfig();
                    //printf("[socbirrt.cpp-PlanPath-300] _randomConfig\n");
					//for (int v=0; v<_randomConfig.size(); v++) {
					//	printf(" %f", _randomConfig[v]);
					//}
					//printf("\n");
                    if(_pForwardTree->_pMakeNext->AddRootConfiguration(_pForwardTree,_randomConfig) )
                        break;
                }
                else
                {
                	//printf("[socbirrt.cpp-PlanPath-305 vikguess\n");
                	//for (int v=0; v<_parameters->vikguess.size(); v++) {
                	//	printf(" %f", _parameters->vikguess[v]);
                	//}
                	//printf("\n");
                	if(_pForwardTree->_pMakeNext->AddRootConfiguration(_pForwardTree,_parameters->vikguess) )
                        break;
                }
                if(timeGetThreadTime() - starttime > max_firstik_time)
                {
                    RAVELOG_FATAL("Unable to find a start IK solution in %f seconds, planner failed\n",(dReal)max_firstik_time);
                    _outputstream << "Unable to find a start IK solution in " << (dReal)max_firstik_time << " seconds, planner failed\n";
                    return CleanUpReturn(false);
                }
            }
            RAVELOG_INFO("Got first start ik solution!\n");
            printf("[socbirrt.cpp-PlanPath-276] Got first start ik sollution\n");
        }


        if(_pBackwardTree->GetSize() == 0)
        {
            RAVELOG_DEBUG("Backward Tree is empty\n");
            printf("[socbirrt.cpp-PlanPath-333] Backward Tree is empty\n");
            while(1)
            {
            	//printf("[socbirrt.cpp-PlanPath-633] AddRootConfiguration\n");
                if(_pBackwardTree->_pMakeNext->AddRootConfiguration(_pBackwardTree,_parameters->vikguess) )
                    break;
                if(timeGetThreadTime() - starttime > max_firstik_time)
                {
                    RAVELOG_FATAL("Unable to find a goal IK solution in %f seconds, planner failed\n",(dReal)max_firstik_time);
                    _outputstream << "Unable to find a start IK solution in " << (dReal)max_firstik_time << " seconds, planner failed\n";
                    return CleanUpReturn(false);
                }
            }
            RAVELOG_INFO("Got first goal ik solution!\n");
            printf("[socbirrt.cpp-PlanPath-644] Got first goal ik solution\n");
        }

        if(RANDOM_FLOAT() < P_SAMPLE_IK)
        {
            RAVELOG_DEBUG("Sampling IK...\n");
            if(_parameters->bsamplingstart)
                _pForwardTree->_pMakeNext->AddRootConfiguration(_pForwardTree,_parameters->vikguess);
            if(_parameters->bsamplinggoal)
                _pBackwardTree->_pMakeNext->AddRootConfiguration(_pBackwardTree,_parameters->vikguess);
            
            continue;
        }

        //RAVELOG_INFO("TreeA: %d    TreeB: %d\n",TreeB->GetSize(),  TreeB->GetSize());

        RAVELOG_DEBUG("Sampling Random Config\n");
        _PickRandomConfig();
        //printf("[socbirrt.cpp-PlanPath-364] _randomConfig\n");
		//for (int v=0; v<_randomConfig.size(); v++) {
		//	printf(" %f", _randomConfig[v]);
		//}
		//printf("\n");
        //Tree A
        iclosest = _FindClosest(TreeA);
        RAVELOG_DEBUG("Extending TreeA toward random config\n");

        TreeA->_pMakeNext->MakeNodes(TreeA->GetSize(),TreeA->GetNode(iclosest),&_randomConfig,&NewNodes,this,false);
        
        if(NewNodes.size() > 0)
        {
            iConnectedA = NewNodes.back().GetID();
            RAVELOG_DEBUG("TreeA MakeNodes Has Generated: \n");
            for(unsigned int i = 0; i < NewNodes.size(); i++)
            {            
                NewNodes[i].Print();
            }
            TreeA->AddNodes(&NewNodes);
            
            //set up randomConfig for treeB
            _randomConfig = *NewNodes.back().GetDataVector();

        }
        else
        {
            RAVELOG_DEBUG("TreeA MakeNodes Has Generated No New Nodes\n");
            RAVELOG_DEBUG("Using Node %d as target\n",iclosest);
            
            _randomConfig = *TreeA->GetNode(iclosest)->GetDataVector();
            iConnectedA = iclosest;

        }
        //Tree B
        iclosest = _FindClosest(TreeB);
        RAVELOG_DEBUG("Extending TreeB toward random config\n");
        if(TreeB->_pMakeNext->MakeNodes(TreeB->GetSize(),TreeB->GetNode(iclosest),&_randomConfig,&NewNodes,this,false))
        {
            RAVELOG_DEBUG("TreeB MakeNodes Has Generated: \n");
            for(unsigned int i = 0; i < NewNodes.size(); i++)
            {
                NewNodes[i].Print();
            }
            if( NewNodes.size() > 0 )
                iConnectedB = NewNodes.back().GetID();
            else
                iConnectedB = -1;
        }
        else
        {
            RAVELOG_DEBUG("TreeB Has Not Reached Target\n");
            iConnectedB = -1;
        }



        TreeB->AddNodes(&NewNodes);
        
        if(iConnectedB != -1)
            bContinuePlan = false;
        else
        {
            ptreetemp = TreeA;        
            TreeA = TreeB;
            TreeB = ptreetemp;
        }
        numitr++;
        
        

        if(timeGetThreadTime()-starttime > max_planning_time)
        {
#ifdef RECORD_TIMES
            RAVELOG_FATAL("Component       \t# calls\t total time (s)\t avg. time (s)\n");
            RAVELOG_FATAL("-------------------------------------------------------------\n");
            RAVELOG_FATAL("Projection      \t %d \t %f \t %f\n", projectioncalls, projectiontime, projectiontime/projectioncalls);
            RAVELOG_FATAL("CollisionCheck  \t %d \t %f \t %f\n", collisioncheckcalls, collisionchecktime, collisionchecktime/collisioncheckcalls);
            RAVELOG_FATAL("NearestNeighbor \t %d \t %f \t %f\n", nncalls, nntime, nntime/nncalls);
#endif
            RAVELOG_FATAL("Planner Failed: %fs time limit reached\n",max_planning_time);
            _outputstream << "Planner Failed: " << max_planning_time << "s time limit reached\n";
            return CleanUpReturn(false);
        }

    }

    //figure out which tree is which
    if(TreeB->GetFromGoal())
    {
        _pForwardTree->_iConnected = iConnectedA;
        _pBackwardTree->_iConnected = iConnectedB;
    }
    else
    {
        _pForwardTree->_iConnected = iConnectedB;
        _pBackwardTree->_iConnected = iConnectedA;

    }

#ifdef RECORD_TIMES
    RAVELOG_FATAL("Component       \t# calls\t total time (s)\t avg. time (s)\n");
    RAVELOG_FATAL("-------------------------------------------------------------\n");
    RAVELOG_FATAL("Projection      \t %d \t %f \t %f\n", projectioncalls, projectiontime, projectiontime/projectioncalls);
    RAVELOG_FATAL("CollisionCheck  \t %d \t %f \t %f\n", collisioncheckcalls, collisionchecktime, collisionchecktime/collisioncheckcalls);
    RAVELOG_FATAL("NearestNeighbor \t %d \t %f \t %f\n", nncalls, nntime, nntime/nncalls);

    //this is a more compact way to dispay the information if you're doing lots of trials
    //RAVELOG_FATAL("%f \t %f \t %f \t %f \t %f \t %f\n",projectiontime,collisionchecktime,nntime,projectiontime/projectioncalls,collisionchecktime/collisioncheckcalls,nntime/nncalls);
#endif
    RAVELOG_FATAL("Planning time: %fs\n",timeGetThreadTime()-starttime);
    
    //construct optimized trajectory
    _JoinPathHalves();
    //if(_parameters->bmergetree) {
    	MergeTree();
    //}
    bool bTerminated;
    _OptimizePath(bTerminated, starttime);
    if(bTerminated)
    {
        _outputstream << "Smoothing terminated from external function\n";
        return CleanUpReturn(false);
    }

    _CreateTraj(ptraj);

    return CleanUpReturn(true);
}


int SoCBirrtPlanner::MakeNext::SampleTSRChainIndex(const std::vector<int>& TSRChaininds)
{
    if(TSRChaininds.size() == 0)
        return -1;

    if(TSRChaininds.size() == 1)
        return TSRChaininds[0];


    //sampling weighted by sumbounds
    //calculate the sum of the sumbounds

    dReal sumsumbounds = 0;
    for(int i = 0; i < TSRChaininds.size(); i++)
        sumsumbounds += _planner->_parameters->vTSRChains[TSRChaininds[i]].GetSumOfBounds();

    
    //calculate the cumsum
    std::vector<dReal> _cumsum(TSRChaininds.size()+1);
    _cumsum[0] = 0;


    for(int i = 0; i < TSRChaininds.size(); i++)
        _cumsum[i+1] = _cumsum[i] + _planner->_parameters->vTSRChains[TSRChaininds[i]].GetSumOfBounds()/sumsumbounds;

    
    //see which bin the random number falls into
    dReal r = RANDOM_FLOAT();
    //RAVELOG_INFO("r: %f\n",r);
    for(int i = 0; i < _cumsum.size()-1; i++)
    {
        //RAVELOG_INFO("Check bin: %f %f\n",_scorecumsum[i],_scorecumsum[i+1]);
        if( r >= _cumsum[i] && r <= _cumsum[i+1])
            return TSRChaininds[i];
    }
    return -1;

}

/// create a random permutation of the numbers 0 through N
void rand_perm(int N, std::vector<int>& p)
{
    p.resize(N,0);
    for (int i = 0; i < N; ++i) {
      int j = rand() % (i + 1);
      p[i] = p[j];
      p[j] = i;
    }
}

bool SoCBirrtPlanner::MakeNext::AddRootConfiguration(SoCBirrtPlanner::NodeTree* ptree, std::vector<dReal>& guess)
{
  //RAVELOG_INFO("starting addroot\n");

    /// _planner->_parameters->bikfastsinglesolution controls whether to use the attached IK solver for a single IK solution with collision checking(1) or many solutions with no collision checking(0), if 0 these solutions will be col-checked later

    std::vector<Transform> vTtarg(_planner->_pRobot->GetManipulators().size());
    std::vector<dReal> q_near;//not used, just need to fill an argument

    _planner->_randomConfig = guess;
    std::vector<dReal> q_temp;
    int activeind;

    std::vector< std::vector<dReal> > qGuesses;

    //pick the TSRChaininds corresponding to this tree
    std::vector<std::vector<int > >* pvvManipToSamplingTSRChainind;

    if(ptree->GetFromGoal())
    {
        pvvManipToSamplingTSRChainind = &_planner->vvManipToGoalSamplingTSRChainind;
    }
    else
    {
        pvvManipToSamplingTSRChainind = &_planner->vvManipToStartSamplingTSRChainind;
    }

    // generate sample and do ik
    for(int i = 0; i < pvvManipToSamplingTSRChainind->size(); i++)
    {
        if(pvvManipToSamplingTSRChainind->at(i).size() == 0)
            continue;

        int sampleind = SampleTSRChainIndex(pvvManipToSamplingTSRChainind->at(i));
        //RAVELOG_INFO("tsrchain ind: %d\n",sampleind);
        vTtarg[i] = _planner->_parameters->vTSRChains[sampleind].GenerateSample();

        //if this manipulator has an attached IK solver, use it for initial guess(es)
        if(!!(_probot->GetManipulators()[i]->GetIkSolver()))
        {
            _planner->SetDOF(_planner->_randomConfig);

            std::vector< std::vector<dReal> > qSolutions;
            std::vector< dReal > qSolution;

            if(_planner->_parameters->bikfastsinglesolution && _probot->GetManipulators()[i]->FindIKSolution(_planner->_parameters->Tattachedik_0[i]*vTtarg[i], qSolution, true))
            {
                RAVELOG_INFO("Attached Solver Found a Solution!\n");
                qGuesses.resize(1);
                qGuesses[0].resize(_planner->GetNumDOF());
                for(int k = 0; k < _probot->GetManipulators()[i]->GetArmIndices().size(); k++)
                {
                    activeind = _planner->vJointIndToActiveInd[_probot->GetManipulators()[i]->GetArmIndices()[k]];
                    if(activeind == -1)
                    {
                        RAVELOG_INFO("ERROR: An index used by the attached IK solver of manipulator %d is not active, terminating AddRootConfiguration!\n",i);
                        return false;
                    }
                    qGuesses[0][activeind] = qSolution[k];
                }
            }
            else if(!_planner->_parameters->bikfastsinglesolution && _probot->GetManipulators()[i]->FindIKSolutions(_planner->_parameters->Tattachedik_0[i]*vTtarg[i], qSolutions, false))
            {
                RAVELOG_INFO("Attached Solver Found %d Solutions!\n",qSolutions.size());
                //sometimes we get too many solutions, so sample some number from this set
                if(qSolutions.size() > MAX_NUM_ATTACHEDIKSOLS)
                {
                    RAVELOG_INFO("Too many IK solutions, picking %d.\n",MAX_NUM_ATTACHEDIKSOLS);
                    std::vector< std::vector<dReal> > qSolutions2(MAX_NUM_ATTACHEDIKSOLS);
                    std::vector<int> perm;
                    rand_perm(qSolutions.size(), perm);
                    //take first MAX_NUM_ATTACHEDIKSOLS from random permutation
                    for(int j = 0; j < qSolutions2.size(); j++)
                        qSolutions2[j] = qSolutions[perm[j]];
                    qSolutions = qSolutions2;
                }

                //TODO: Currently, IK solutions for multiple manipulators are paired up randomly, might be better to collision check each one by itself then use all col-free combinations
//                 if(qGuesses.size() != 0 && qGuesses.size() < qSolutions.size())
//                 {
//                     RAVELOG_INFO("One manipulator has more IK solutions than another, duplicated some.\n");
//                     std::vector< std::vector<dReal> > qSolutions2(MAX_NUM_ATTACHEDIKSOLS);
//                     std::vector<int> perm;
//                     rand_perm(qSolutions.size(), perm);
//                 }
                qGuesses.resize(qSolutions.size());
                for(int j = 0; j < qGuesses.size(); j++)
                {
                    qGuesses[j].resize(_planner->GetNumDOF());
                    for(int k = 0; k < _probot->GetManipulators()[i]->GetArmIndices().size(); k++)
                    {
                        activeind = _planner->vJointIndToActiveInd[_probot->GetManipulators()[i]->GetArmIndices()[k]];
                        if(activeind == -1)
                        {
                            RAVELOG_INFO("ERROR: An index used by the attached IK solver of manipulator %d is not active, terminating AddRootConfiguration!\n",i);
                            return false;
                        }
                        qGuesses[j][activeind] = qSolutions[j][k];
                    }
                }
            }
            else
            {
              //RAVELOG_INFO("Attached IK solver could not find a solution\n");
                return false;
            }
        }
    }

    //if no guesses, put in the default guess
    if(qGuesses.size() == 0)
        qGuesses.push_back(_planner->_randomConfig);

    bool bGotSolution = false;
    for(int i = 0; i < qGuesses.size(); i++)
    {
        _planner->ClearTSRChainValues();
        // if a guess is valid, this will confirm it and return quickly, if not, it will try to use generalik to get an answer
        if(ConstrainedNewConfig(qGuesses[i], q_near, _planner->vTSRChainValues_temp,false, &vTtarg, !ptree->GetFromGoal()) && !(_planner->_CheckCollision(qGuesses[i]) ))
        {
            
            RrtNode * pikNode = new RrtNode(_planner->GetNumDOF(),ptree->GetFromGoal(),ptree->GetSize());
            pikNode->SetConfig(&qGuesses[i][0]);
            pikNode->SetTSRChainValues(_planner->vTSRChainValues_temp); 

            ptree->AddNode(*pikNode); //this will automatically increment the size
            pikNode->Print();
            delete pikNode;
            //return true;
            bGotSolution = true;
        }
    }
    if(bGotSolution)
        RAVELOG_INFO("Constrained IK solution(s) found!\n");
    else
        RAVELOG_INFO("No Constrained IK solution found\n");

    return bGotSolution;
}

bool SoCBirrtPlanner::_CheckCollision(std::vector<dReal>& pConfig)
{
    collisioncheckcalls++;
#ifdef RECORD_TIMES
    double starttime = timeGetThreadTime();
#endif

    SetDOF(pConfig);

    CollisionReportPtr preport(new CollisionReport());

    
#ifdef TRACK_COLLISIONS
    GetEnv()->SetCollisionOptions(CO_Contacts);
#endif


    if(GetEnv()->CheckCollision(KinBodyConstPtr(_pRobot), sIgnoredBodies,sIgnoredLinks, preport) ) {
        RAVELOG_DEBUG("SoCBirrt collision: %s:%s with %s:%s\n", preport->plink1->GetParent()->GetName().c_str(), preport->plink1->GetName().c_str(), preport->plink2->GetParent()->GetName().c_str(), preport->plink2->GetName().c_str());

#ifdef TRACK_COLLISIONS
        Transform linktm = preport->plink2->GetTransform();

        //RAVELOG_INFO("num contacts: %d\n",preport->contacts.size());
        
        //RAVELOG_INFO("%f %f %f\n",preport->contacts[0].pos.x,preport->contacts[0].pos.y,preport->contacts[0].pos.z);
        for(int i = 0; i < preport->contacts.size(); i++)
            colfile << preport->contacts[i].pos.x << " " << preport->contacts[i].pos.y << " " << preport->contacts[i].pos.z << endl;

        //colfile << linktm.trans.x << " " << linktm.trans.y << " " << linktm.trans.z << endl;
        
#endif

#ifdef RECORD_TIMES
        collisionchecktime += timeGetThreadTime() - starttime;
#endif
        return true;
    }

    //MIMIC BODY COLLISION CHECKING IS CURRENTLY DISABLED (to allow door opening, etc.)
    for(int i = 0; i < vMimicBodies.size(); i++)
    {
        GetEnv()->CheckCollision(KinBodyConstPtr(vMimicBodies[i]), sIgnoredBodies, sIgnoredLinks, preport);
        RAVELOG_DEBUG("Mimic body is in collision\n");
        //RAVELOG_DEBUG("SoCBirrt collision: %S:%S with %S:%S\n", preport->plink1->GetParent()->GetName(), preport->plink1->GetName(), preport->plink2->GetParent()->GetName(), preport->plink2->GetName());
    }

    //GetEnv()->CheckSelfCollision(KinBodyConstPtr(_pRobot), preport) // doesn't work for grabbed bodies
    if(_pRobot->CheckSelfCollision(preport))
    {
        RAVELOG_DEBUG("SoCBirrt self-collision: %s:%s with %s:%s\n", preport->plink1->GetParent()->GetName().c_str(), preport->plink1->GetName().c_str(), preport->plink2->GetParent()->GetName().c_str(), preport->plink2->GetName().c_str());

#ifdef RECORD_TIMES
        collisionchecktime += timeGetThreadTime() - starttime;
#endif
        return true;
    }

#ifdef RECORD_TIMES
    collisionchecktime += timeGetThreadTime() - starttime;
#endif
    return false;
}


bool SoCBirrtPlanner::_CheckCollision(std::vector<dReal>& pQ0, std::vector<dReal>&  pQ1, IntervalType interval)
{
  // set the bounds based on the interval type
  int start;
  bool bCheckEnd;
  switch (interval) {
  case OPEN:
    start = 1;  bCheckEnd = false;
    break;
  case OPEN_START:
    start = 1;  bCheckEnd = true;
    break;
  case OPEN_END:
    start = 0;  bCheckEnd = false;
    break;
  case CLOSED:
    start = 0;  bCheckEnd = true;
    break;
  default:
    cerr << "RrtPlanner: ERROR - unknown interval type." << endl;
  }

  // first make sure the end is free
  if (bCheckEnd)
    if (_CheckCollision(pQ1))
      return true;

  // compute  the discretization
  int i, numSteps = 1;
  for (i = 0; i < GetNumDOF(); i++) {
    int steps = (int)(fabs(pQ1[i] - pQ0[i]) * _jointResolutionInv[i]);
    if (steps > numSteps)
      numSteps = steps;
  }
  //cerr << "CheckCollision: number of steps: " << numSteps << endl;

  // compute joint increments
  for (i = 0; i < GetNumDOF(); i++)
    _jointIncrement[i] = (pQ1[i] - pQ0[i])/((dReal)numSteps);

  // check for collision along the straight-line path
  // NOTE: this does not check the end config, and may or may
  // not check the start based on the value of 'start'
  for (int f = start; f < numSteps; f++) {
    for (i = 0; i < GetNumDOF(); i++)
      checkConfig[i] = pQ0[i] + (_jointIncrement[i] * f);
    
    if (_CheckCollision(checkConfig))
      return true;
  }

  return false;
}


void SoCBirrtPlanner::_PickRandomConfig()
{

    for (int i = 0; i < GetNumDOF(); i++)
    {
        if( _viscircular[i]) {
            _randomConfig[i] = -PI + 2*PI*RANDOM_FLOAT();
        }
        else{
            _randomConfig[i] = _lowerLimit[i] + RANDOM_FLOAT(_validRange[i]);

        }
    }
}


int SoCBirrtPlanner::_FindClosest(NodeTree* pNodeTree)
{
  nncalls++;

#ifdef RECORD_TIMES
  double starttime = timeGetThreadTime();
#endif

  int iBest = UNDEFINED;
  dReal minDist = INF;
  dReal dist;

  for (int i = 0; i < pNodeTree->GetSize(); i++) 
  {
      dist = pdistmetric->Eval(&_randomConfig[0],pNodeTree->GetNode(i)->GetData());
      //RAVELOG_DEBUG("dist: %f\n",dist);
      if (dist < minDist) {
        minDist = dist;
        iBest = i;
      }
  }
  
  assert(iBest != UNDEFINED);

#ifdef RECORD_TIMES
  nntime += timeGetThreadTime() - starttime;
#endif

  return iBest;
}

bool SoCBirrtPlanner::SetVecPathFromTraj(TrajectoryBasePtr ptraj_in)
{
    //TODO: Change this to use configurationspecs
    std::vector<int> activedofinds = _pRobot->GetActiveDOFIndices();

    vecpath.clear();
    if(GetNumDOF() != activedofinds.size())
    {
        RAVELOG_FATAL("SoCBirrt::SetVecPathFromTraj ERROR: GetNumDOF (%d) != activedofinds.size() (%d).\n",GetNumDOF(),activedofinds.size());
        return false;
    }

    //vTSRChainValues_temp.resize(numTSRTotalDOF);
    _randomConfig.resize(GetNumDOF());
    std::vector<dReal> vtemp;
    for(int i =0; i < ptraj_in->GetNumWaypoints(); i++)
    {
        for(int j = 0; j < activedofinds.size(); j++)
        {
            //_randomConfig[j] =  ptraj_in->GetPoints()[i].q[activedofinds[j]];
            ptraj_in->GetWaypoint(i,vtemp);
            _randomConfig[j] = vtemp[activedofinds[j]];
        }

        _pActiveNode = new RrtNode(GetNumDOF(),false,i);
        _pActiveNode->SetConfig(&_randomConfig[0]);
        _pActiveNode->SetTSRChainValues(vTSRChainValues_temp);
        _pActiveNode->SetParent(i-1);
        vecpath.push_back(*_pActiveNode);
        delete _pActiveNode;
        RAVELOG_DEBUG("%d\n",i);

    }

    return true;
}


bool SoCBirrtPlanner::_OptimizePath(bool &bTerminated, double starttime)
{

    //ofstream outfile("pathlength.txt",ios::out);
    bTerminated = false;
    bool bResult = false;
    int numSuccessful = 0;
    int origLength = (int)vecpath.size();
    
    
    std::vector<RrtNode> NewNodes;
    std::vector<RrtNode> temp;
    temp.reserve(origLength);
    int iStart=-1;
    int iGoal=-1;
    int j;
    int iLimit;

    if(_parameters->smoothingitrs  == -1)
        iLimit = NUM_OPT_ITERATIONS;
    else
        iLimit = _parameters->smoothingitrs;

    dReal origpathlength;
    dReal straightpathlength;
    RAVELOG_INFO("Smoothing...\n");
    for (int i = 0; i < iLimit  ; i++)
    {
        if(timeGetThreadTime()-starttime > max_planning_time)
        {
            // This doesn't count as terminated because we just ran out of time, path is still valid
            RAVELOG_INFO("Smoothing terminated because %dms time limit reached\n",max_planning_time);
            return bResult;
        }

        //terminate smoothing if PS_PlanFailed was set
        if(*(_parameters->pplannerstate) == PS_PlanFailed)
        {
            RAVELOG_INFO("Smoothing terminated from external function.\n");
            bTerminated = true;
            return bResult;
        }
        //outfile << GetPathLength(vecpath,0,vecpath.size()) << endl;

        temp.clear();
        if((int)vecpath.size() <= 2)
            break;
        iGoal = 1+RANDOM_INT((int)vecpath.size()-1);
        iStart = RANDOM_INT(iGoal);
        //RAVELOG_INFO("Start: %d   Goal: %d  Size: %d TreeSize: %d\n",iStart,iGoal,vecpath.size(),_pForwardTree->GetSize());
        
        origpathlength = GetPathLength(vecpath,iStart,iGoal+1);
        straightpathlength = pdistmetric->Eval(vecpath[iStart].GetData(), vecpath[iGoal].GetData());
        //check if we should try to smooth this segment
        RAVELOG_DEBUG("straight: %f current: %f\n",straightpathlength,origpathlength);
        if( 1.1*straightpathlength >  origpathlength)
        {
            RAVELOG_DEBUG("saved a smoothing iteration\n");
            continue;
        }


        if (_pForwardTree->_pMakeNext->MakeNodes(0,&vecpath[iStart],vecpath[iGoal].GetDataVector(),&NewNodes,this,false)) 
        {
            if(pdistmetric->Eval(vecpath[iStart].GetData(), NewNodes[0].GetData()) + GetPathLength(NewNodes,0,NewNodes.size()) < origpathlength)
            {
                for(j = 0; j <= iStart; j++)
                    temp.push_back(vecpath[j]);
                for(j = 0; j < (int)NewNodes.size(); j++)
                    temp.push_back(NewNodes[j]);
                for(j = iGoal+1; j < (int)vecpath.size(); j++)
                    temp.push_back(vecpath[j]);

                vecpath = temp;

                numSuccessful++;
                bResult = true;
            }
        }
    }
    // statistics
    RAVELOG_DEBUG("Number of Path Points:  original=%d  final=%d\n", origLength, vecpath.size());
    RAVELOG_INFO("Done\n");
    //outfile.close();
    return bResult;
}


dReal SoCBirrtPlanner::GetPathLength(std::vector<RrtNode>& path, int indexstart, int indexend)
{
    dReal length = 0;
    for(int i = indexstart; i < indexend-1; i++)
        length += pdistmetric->Eval(path[i].GetData(), path[i+1].GetData());

    return length;
}


bool SoCBirrtPlanner::_CreateTraj(TrajectoryBasePtr ptraj)
{

    ConfigurationSpecification planningspec =  _pRobot->GetActiveConfigurationSpecification();


    RobotBasePtr tempbody;
    for(int i = 0; i < _parameters->vTSRChains.size(); i++)
    {
        tempbody = _parameters->vTSRChains[i].GetMimicBody();
        if(tempbody == NULL)
            continue;

        tempbody->SetActiveDOFs(_parameters->vTSRChains[i].GetMimicDOFInds());
        planningspec = planningspec + tempbody->GetActiveConfigurationSpecification();
    }


    ptraj->Init(planningspec);


    vector<dReal> pvectemp;
    vector<dReal> vtraj_data;

    for (unsigned int f = 0; f < vecpath.size(); f++) {
        pvectemp = *vecpath[f].GetDataVector();
        for (int i = 0; i < _pRobot->GetActiveDOF(); i++)
        {
            vtraj_data.push_back(pvectemp[i]);
        }

        for(int i = 0; i < _parameters->vTSRChains.size(); i++)
        {
            tempbody = _parameters->vTSRChains[i].GetMimicBody();
            if(tempbody == NULL)
                continue;
            //RAVELOG_INFO("Processing TSR Chain with Mimic Body %d\n", i);

            for (int k = 0; k < _parameters->vTSRChains[i].GetNumMimicDOF(); k++)
                vtraj_data.push_back(*(vecpath[f].GetData()+vTSRChainIndexToConfigurationOffset[i] + k));
        }

        
    }
    ptraj->Insert(ptraj->GetNumWaypoints(),vtraj_data,planningspec);
    _parameters->_sPostProcessingPlanner = "lineartrajectoryretimer";
    _parameters->_configurationspecification = planningspec;
    _parameters->_vConfigLowerLimit = _lowerLimit;
    _parameters->_vConfigUpperLimit = _upperLimit;

    std::vector<dReal> velocitylimits(_parameters->_vConfigLowerLimit.size(),1);
    _parameters->_vConfigVelocityLimit = velocitylimits;
    _parameters->_vConfigAccelerationLimit = velocitylimits;
    _parameters->_vConfigResolution = velocitylimits;

    std::vector<dReal> initconfig;
    ptraj->GetWaypoints(0,0,initconfig);
    _parameters->vinitialconfig = initconfig;

    _ProcessPostPlanners(_pRobot,ptraj);


    if(ptraj->GetNumWaypoints() != vecpath.size())
    {
        RAVELOG_INFO("Error, trajectory timer changed the number of points: %d before vs. %d after\n",vecpath.size(),ptraj->GetNumWaypoints());
        return false;
    }

/*
    //OLD way to generate traj for mimic bodies
    //if using chains and they have mimic bodies, make a trajectory for the mimic bodies and print it out
    //calc timing for the robot's traj to ensure syncing
    std::vector<dReal> mimicbodyvals;
    RobotBasePtr tempbody;
    for(int i = 0; i < _parameters->vTSRChains.size(); i++)
    {
        tempbody = _parameters->vTSRChains[i].GetMimicBody();
        if(tempbody == NULL)
            continue;

        vtraj_data.clear();

        TrajectoryBasePtr pbodytraj = RaveCreateTrajectory(GetEnv(),"");//tempbody->GetDOF());
        pbodytraj->Init(tempbody->GetActiveConfigurationSpecification());


        for (unsigned int j = 0; j < ptraj->GetNumWaypoints(); j++)
        {
            _parameters->vTSRChains[i].MimicValuesToFullMimicBodyValues(vecpath[j].GetData()+vTSRChainIndexToConfigurationOffset[i], mimicbodyvals);
            for (int k = 0; k < tempbody->GetDOF(); k++)
                    vtraj_data.push_back(mimicbodyvals[k]);
                    //p.q[k] = mimicbodyvals[k];

            //p.time = pfulltraj->GetPoints()[j].time;
            //p.trans = tempbody->GetTransform();

            //pbodytraj->AddPoint(p);
    
        }

        pbodytraj->Insert(pbodytraj->GetNumWaypoints(),vtraj_data);

        OpenRAVE::planningutils::RetimeActiveDOFTrajectory(pbodytraj,tempbody,false,1,"LinearTrajectoryRetimer");
        // save the constrained trajectory
        char trajfilename[128];
        sprintf(trajfilename,"%s_traj.txt",tempbody->GetName().c_str());
        ofstream outfile(trajfilename,ios::out);
        //pbodytraj->Write(outfile, Trajectory::TO_IncludeTimestamps|Trajectory::TO_IncludeBaseTransformation);
        pbodytraj->serialize(outfile);
        chmod(trajfilename, S_IRWXG | S_IRWXO | S_IRWXU); //chmod 777
        outfile.close();
    }
    */

    return true;
}


bool SoCBirrtPlanner::_JoinPathHalves()
{
    vector<RrtNode> vectemp;
    //get forward path
    int index = -1;

    
    _pActiveNode = _pForwardTree->GetNode(_pForwardTree->_iConnected);
    while(true)
    {
        //RAVELOG_DEBUG("Node Added\n");
        //_pActiveNode->Print();
        vectemp.push_back(*_pActiveNode);
        index = _pActiveNode->GetParent();
        if(index == -1)
            break;
        _pActiveNode = _pForwardTree->GetNode(index);
    }
    


    for(int i = (int)vectemp.size()-1; i >= 0; i--)
        vecpath.push_back(vectemp[i]);

    
    
    if(_pBackwardTree->GetNode(_pBackwardTree->_iConnected)->GetParent() == -1)
        _pActiveNode = _pBackwardTree->GetNode(_pBackwardTree->_iConnected);
    else
    {
        _pActiveNode = _pBackwardTree->GetNode(_pBackwardTree->GetNode(_pBackwardTree->_iConnected)->GetParent());
        
        bool bduplicate = true;
        while(true)
        {
            //RAVELOG_DEBUG("Node Added\n");
            //_pActiveNode->Print();
            vecpath.push_back(*_pActiveNode);
            index = _pActiveNode->GetParent();
            if(index == -1)
                break;
            _pActiveNode = _pBackwardTree->GetNode(index);
        }
    }

    RAVELOG_DEBUG("Joined Halves\n");
    return true;
}

dReal SoCBirrtPlanner::SoCBirrtDistanceMetric::Eval(const void* c0, const void* c1)
{
    return sqrt(OpenRAVE::mathextra::lengthsqr((dReal*)c0, (dReal*)c1,_pplanner->GetNumDOF()));
}


SoCBirrtPlanner::MakeNext::MakeNext(bool bFromGoal, int numdof, RobotBasePtr  robot, SoCBirrtPlanner * planner)
{
    _bFromGoal = bFromGoal;
    _probot = robot; 
    _probot->GetActiveDOFLimits(_lowerLimit,_upperLimit);
    _planner = planner;    

    diff.resize(numdof);
    dist.resize(numdof);
    newConfig.resize(numdof);
    oldConfig.resize(numdof);
    oldoldConfig.resize(numdof);
    ancientConfig.resize(numdof);
    lastcollchecked.resize(numdof);
}

bool SoCBirrtPlanner::MakeNext::MakeNodes(int TreeSize, RrtNode* StartNode, std::vector<dReal>* targetConfig,std::vector<RrtNode>* voutput,SoCBirrtPlanner* planner,bool bAllowGoPast)
{
   
    assert(StartNode != NULL);
    RAVELOG_DEBUG("In MakeNodes\n");
    bExtendConnect = true;

    startConfig = StartNode->GetData();

    numDOF = (unsigned int)targetConfig->size(); 
    
    success = false;
    normDiff = 0;
    configDiff = 0;
    normDist =0;
    targetReached = false;

    unsigned int i;

    voutput->clear();
    
    // calculate diffence vector
    for (i = 0; i < numDOF; i++)
    {
        lastcollchecked[i] = oldConfig[i] = startConfig[i];
        diff[i] = targetConfig->at(i) - oldConfig[i];
        normDiff += diff[i]*diff[i];
    }
    
    normDiff = sqrt(normDiff);
    //printf("Diff: \n");
    for (i = 0; i < numDOF; i++)
    {
        diff[i] = diff[i]/normDiff;
        //printf("%f ", diff[i]);
    }
    //printf("\n");

    _planner->vTSRChainValues_temp = StartNode->GetTSRChainValues();

    dReal oldnormDist;
    dReal ancientnorm;
    dReal colnorm;
    int numsteps = 0;
    normDist = INF;

    // when bExtendConnect is false, this loop executes only once
    do
    {
        if(numsteps > 2)
        {
            ancientnorm = 0;
            for (i = 0; i < numDOF; i++)
            {        
                ancientnorm += (ancientConfig[i]- oldConfig[i])*(ancientConfig[i] - oldConfig[i]);
            }
            //PrintMatrix(voutput->at(voutput->size() - 2).GetData(),numDOF,1," oldold cfg: ");
            //PrintMatrix(oldConfig,numDOF,1," cur cfg: ");
            ancientnorm = sqrt(ancientnorm);
            //RAVELOG_INFO("ancientnorm: %f\n",ancientnorm);
        }

        oldnormDist = normDist;
        normDist =0;

        //printf("Dist: \n");
        for (i = 0; i < numDOF; i++)
        {
            dist[i] = targetConfig->at(i) - oldConfig[i]; 
            //printf("%f ", dist[i]);
            normDist += dist[i]*dist[i]; 
        }
        //printf("\n");
        normDist = sqrt(normDist);

        //get new config
        if(normDist < STEP_LENGTH)
        {
            targetReached = true;
            for (i = 0; i < numDOF; i++)
                newConfig[i] = targetConfig->at(i);
        }
        else if(!bAllowGoPast && normDist > oldnormDist)
        {   
            RAVELOG_DEBUG("Extension is going past target\n");
            break;
        }
        else if(normDist > oldnormDist - TINY && normDist < oldnormDist + TINY)
        {   
            RAVELOG_DEBUG("Extension has stalled\n");
            break;
        }
        //this is to prevent bouncing between two configurations
        else if(numsteps > 2 && ancientnorm < TINY)
        {
            RAVELOG_DEBUG("Extension is bouncing between two configurations\n");
            break;

        }
        else
            for (i = 0; i < numDOF; i++)
                newConfig[i] = oldConfig[i] + (STEP_LENGTH * dist[i]/normDist);

        // code for handling constraints
        std::vector<dReal> q_near(numDOF);
        std::vector<dReal> q_s(numDOF);
        for(i = 0; i < numDOF; i++)
        {
            q_s[i] = newConfig[i];
            q_near[i] = oldConfig[i];
        }
           
        // constrain the new configuration
        if(ConstrainedNewConfig(q_s,q_near,_planner->vTSRChainValues_temp))
        {
            for(i = 0; i < numDOF; i++)
                newConfig[i] = q_s[i];
        }
        else
        {
            return false;
        }

        // keep track of how much you've moved since last collision check
        colnorm = 0;
        for(i = 0; i < numDOF; i++)
            colnorm += (lastcollchecked[i]-newConfig[i])*(lastcollchecked[i]-newConfig[i]);
        colnorm = sqrt(colnorm);

        //RAVELOG_INFO("Colnorm: %f\n", colnorm);

        bool bForwardBack = StartNode->GetFromGoal();


        if(colnorm >= STEP_LENGTH - TINY || targetReached)
        {
            if(planner->_CheckCollision(lastcollchecked,newConfig,OPEN_START))
            {
                RAVELOG_DEBUG("Collision!");
                break;
            }
            for (i = 0; i < numDOF; i++)
                lastcollchecked[i] = newConfig[i];
        }

        if(!bExtendConnect)
        {
            normDist =0;

            for (i = 0; i < numDOF; i++)
            {
                dist[i] = targetConfig->at(i) - newConfig[i]; 
                normDist += dist[i]*dist[i]; 
            }
            RAVELOG_DEBUG("Dist: %f\n",normDist);
            normDist = sqrt(normDist);

            if(normDist < STEP_LENGTH)
            {
                targetReached = true;
                for (i = 0; i < numDOF; i++)
                    newConfig[i] = targetConfig->at(i);
            }
        
        }

        configDiff = 0.0f;
        for (i = 0; i < numDOF; i++)
        {
            configDiff += fabs(oldConfig[i]-newConfig[i]);
            ancientConfig[i] = oldoldConfig[i];
            oldoldConfig[i] = oldConfig[i];
            oldConfig[i] = newConfig[i];
            
        }


        //if there is a significant configuration difference
        if(colnorm >= STEP_LENGTH - TINY || targetReached)
        {

            //create a node with this configuration
            RrtNode* pnewnode = new RrtNode(numDOF,StartNode->GetFromGoal(),TreeSize++);
            if(voutput->size() == 0) {
                pnewnode->SetParent(StartNode->GetID());
                //printf("[socbirrt.cpp-SoCBirrtPlanner::MakeNext::MakeNodes] Add node to parent (id %d) 's child list\n", TreeSize-1);
                StartNode->AddChild(TreeSize-1);
            }
            else {
                pnewnode->SetParent(voutput->back().GetID());
                //printf("[socbirrt.cpp-SoCBirrtPlanner::MakeNext::MakeNodes] Add node to parent (id %d) 's child list\n", TreeSize-1);
                voutput->back().AddChild(TreeSize-1);
            }
            pnewnode->SetConfig(&newConfig[0]);
            pnewnode->SetTSRChainValues(_planner->vTSRChainValues_temp); 
            //printf("[socbirrt.cpp-SoCBirrtPlanner::MakeNext::MakeNodes] push current node to tree list list\n");
            voutput->push_back(*pnewnode);

            delete pnewnode;
        }
        else
            RAVELOG_DEBUG("Configdiff is small\n");

        if(targetReached)
        {
            success = true;
            break;
        }

        numsteps++;
        //RAVELOG_DEBUG("numsteps: %d\n",numsteps);

    }while(bExtendConnect);  

    return success;
}


bool SoCBirrtPlanner::MakeNext::ConstrainedNewConfig(std::vector<dReal>& q_s, std::vector<dReal>& q_near, std::vector<dReal>& vTSRChainValues, bool bCheckDistance, std::vector<Transform>* pvTtarg, bool bStartvsGoalSampling)
{
    _planner->projectioncalls++;
    std::vector<dReal> q0 = q_s;

    //PrintMatrix(&q_s[0], 1, q_s.size(), "q_sin");
    
    _planner->SetDOF(q_s);
    std::vector<dReal> ikparams;  
    dReal totaltempdist = 0;

    Transform Toffset;
    KinBody::LinkPtr  plink;
    
    //make a reference just for convenience
    std::vector<std::vector<int> >& pvvManipToConstraintChain = _planner->vvManipToConstraintTSRChainind;
  
    //pick the samplingTSRinds corresponding to starts or goals
    std::vector<std::vector<int> >* pvvManipToSamplingChain;

    if(bStartvsGoalSampling)
        pvvManipToSamplingChain = &_planner->vvManipToStartSamplingTSRChainind;
    else
        pvvManipToSamplingChain = &_planner->vvManipToGoalSamplingTSRChainind;

    int nummanipsused = 0;
    for(int i = 0; i < pvvManipToConstraintChain.size(); i++)
    {
        if(pvvManipToConstraintChain[i].size() != 0 || (pvvManipToSamplingChain->at(i).size() != 0 && pvTtarg != NULL))
            nummanipsused++;
    }

    ikparams.push_back(nummanipsused);

    Transform EEtm;
    int numtargets = 0;
    int minind;
    dReal mindist;
    dReal tempdist;
    Transform minTtarg;
    Transform tempTtarg;
    TaskSpaceRegionChain* minpTSRChain;
    TaskSpaceRegionChain* temppTSRChain;
    std::vector<dReal> tempTSRvals;
    std::vector<dReal> minTSRvals;
    for(int i = 0; i < pvvManipToConstraintChain.size(); i++)
    {
        //if we have to deal with any constraints for this manipulator
        if(pvvManipToConstraintChain[i].size() != 0)
        {
            if(pvTtarg == NULL)
                EEtm = _probot->GetManipulators()[i]->GetEndEffectorTransform();
            else
                EEtm = pvTtarg->at(i);

            mindist = INF;
            //PrintMatrix(&q_s[0], 1, q_s.size(), "q_s before");
            //PrintMatrix(&vTSRChainValues[0], 1, vTSRChainValues.size(), "vTSRChainValues_temp before");
            for(int j = 0; j < pvvManipToConstraintChain[i].size(); j++)
            {
                temppTSRChain = &(_planner->_parameters->vTSRChains[pvvManipToConstraintChain[i][j]]);
                tempTSRvals.resize(temppTSRChain->GetNumDOF());
                for(int k =0; k < temppTSRChain->GetNumDOF(); k++)
                    tempTSRvals[k] = vTSRChainValues[_planner->vTSRChainIndexToGuessOffset[pvvManipToConstraintChain[i][j]] + k];

#ifdef RECORD_TIMES
                double starttime = timeGetThreadTime();
#endif

                tempdist = temppTSRChain->GetClosestTransform(EEtm,&tempTSRvals[0],tempTtarg);

#ifdef RECORD_TIMES
                _planner->projectiontime += timeGetThreadTime() - starttime;
#endif
                //RAVELOG_INFO("tempdist: %f\n",tempdist);
                if(tempdist >= mindist)
                    continue;
                
                minind = j;
                mindist = tempdist;
                minTtarg = tempTtarg;
                minTSRvals = tempTSRvals;
                minpTSRChain = temppTSRChain;
            }
            //RAVELOG_INFO("manip %d minind %d\n",i,minind);

            if(pvTtarg != NULL && mindist > TINY)
            {
              //RAVELOG_INFO("ERROR: No constraint chain was able to achieve the sampled target for manipulator %d (dist: %f), are the constraint TSR chains defined to be a superset of the sampling TSR chains?\n",i,mindist);
                return false;
            }

            //RAVELOG_INFO("guess offset: %d\n",_planner->vTSRChainIndexToGuessOffset[pvvManipToConstraintChain[i][minind]] );
            for(int k =0; k < minpTSRChain->GetNumDOF(); k++)
                vTSRChainValues[_planner->vTSRChainIndexToGuessOffset[pvvManipToConstraintChain[i][minind]] + k] = minTSRvals[k];
           
            //PrintMatrix(&vTSRChainValues[0], 1, vTSRChainValues.size(), "vTSRChainValues_temp after");
            //RAVELOG_INFO("extracting to: %d\n",_planner->vTSRChainIndexToConfigurationOffset[pvvManipToConstraintChain[i][minind] ]);
            minpTSRChain->ExtractMimicDOFValues(&minTSRvals[0], &q_s[_planner->vTSRChainIndexToConfigurationOffset[pvvManipToConstraintChain[i][minind] ] ]);
            _planner->SetDOF(q_s); //mimic body dofs will be set by SetDOF
        }
        else if(pvvManipToConstraintChain[i].size() == 0 && pvvManipToSamplingChain->at(i).size() != 0 && pvTtarg != NULL) //if only sampling for this manipulator
        {
            minTtarg = pvTtarg->at(i);
        }
        else //otherwise we don't care about this manipulator at all
        {
            continue;
        }
        //PrintMatrix(&q_s[0], 1, q_s.size(), "q_s after");
        ikparams.push_back(i);
        ikparams.push_back(minTtarg.rot.x); ikparams.push_back(minTtarg.rot.y); ikparams.push_back(minTtarg.rot.z); ikparams.push_back(minTtarg.rot.w);
        ikparams.push_back(minTtarg.trans.x); ikparams.push_back(minTtarg.trans.y); ikparams.push_back(minTtarg.trans.z);
        numtargets++;
    }


    ikparams.push_back(0); //don't do any balancing inside the IK solver
    ikparams.push_back(0); //select the mode
    ikparams.push_back(0); //do rotation
    

    boost::shared_ptr<std::vector<dReal> > pq_s(new std::vector<dReal>(q_s) );
    //PrintMatrix(&ikparams[0], 1, ikparams.size(), "ikparams");
    //PrintMatrix(&q0[0], 1, q0.size(), "q0");
#ifdef RECORD_TIMES
    double starttime = timeGetThreadTime();
#endif

    if(numtargets == 0 || _planner->_pIkSolver->Solve(IkParameterization(), q0, ikparams, false, pq_s) )
    {   
#ifdef RECORD_TIMES
        _planner->projectiontime += timeGetThreadTime() - starttime;
#endif
        q_s = *pq_s.get();
        //RAVELOG_INFO("robot ik found\n");
        if(bCheckDistance)
        {
            
            totaltempdist = 0;
            for (int i = 0; i < q_s.size(); i++)
            {
                tempdist = q_near[i] - q_s[i]; 
                //RAVELOG_INFO("tempdist: %d %f\n",i,tempdist);
                totaltempdist += tempdist*tempdist; 
            }

            if(sqrt(totaltempdist) > STEP_LENGTH*2)
            {   
                //PrintMatrix(&q_near[0], 1, q_near.size(), "q_near");        
                //PrintMatrix(&q_s[0], 1, q_s.size(), "q_s");
                RAVELOG_DEBUG("distance: %f\n",sqrt(totaltempdist));            
                return false;
            }
            //RAVELOG_INFO("distance: %f\n",sqrt(totaltempdist));            
        }
        _planner->SetDOF(q_s);

        return CheckSupport();
    }
    else
    {
#ifdef RECORD_TIMES
        _planner->projectiontime += timeGetThreadTime() - starttime;
#endif
        q_s = *pq_s.get();
        //RAVELOG_INFO("Projection failed\n");
        return false;
    }
    return false;
}


bool SoCBirrtPlanner::MakeNext::CheckSupport(bool bDraw)
{
    bDraw = false;
    if(_planner->_parameters->vsupportpolyx.size() == 0)
        return true;

    Vector center;  
    dReal fTotalMass = 0;
    std::vector<KinBody::LinkPtr >::const_iterator itlink;
    std::vector<KinBody::LinkPtr>vlinks = _probot->GetLinks();

    if(_planner->_parameters->bgrabbed)
    {
        for(int i = 0; i < _planner->pheldobject->GetLinks().size(); i++)
            vlinks.push_back(_planner->pheldobject->GetLinks()[i]);
    }

    FORIT(itlink, vlinks) {
        //RAVELOG_INFO("comoffset: %f %f %f\n", (*itlink)->GetCOMOffset().x,(*itlink)->GetCOMOffset().y,(*itlink)->GetCOMOffset().z);
        //GetEnv()->plot3((*itlink)->GetTransform() * (*itlink)->GetCOMOffset(), 1, 0, 10, Vector(0,0,1) );
        center += ((*itlink)->GetTransform() * (*itlink)->GetCOMOffset() * (*itlink)->GetMass());
        fTotalMass += (*itlink)->GetMass();
    }

    if( fTotalMass > 0 )
        center /= fTotalMass;
    //RAVELOG_INFO("\nmass: %f\ncog: %f %f %f\n",fTotalMass,center.x,center.y,center.z);


    dReal testx = center.x;
    dReal testy = center.y;
    dReal * vertx = &_planner->_parameters->vsupportpolyx[0];
    dReal * verty = &_planner->_parameters->vsupportpolyy[0];
    int nvert = _planner->_parameters->vsupportpolyx.size();

    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
    }
   
    center.z = 0;
    if(c)    
    {
        if(bDraw)
        {
            _planner->figure = GetEnv()->plot3(&(DoubleVectorToFloatVector(center)[0]), 1, 0, 10, Vector(0,1,0) );
            //_planner->figure = GetEnv()->plot3(&vpoints[0].x,vpoints.size(),sizeof(vpoints[0]),16.0,&vcolors[0],0);
            graphptrs.push_back(_planner->figure);
        }

        return true;
    }
    else
    {
        if(bDraw)
        {
            _planner->figure = GetEnv()->plot3(&(DoubleVectorToFloatVector(center)[0]), 1, 0, 10, Vector(1,0,0) );
            //_planner->figure = GetEnv()->plot3(&vpoints[0].x,vpoints.size(),sizeof(vpoints[0]),16.0,&vcolors[0],0);
            graphptrs.push_back(_planner->figure);
        }

        RAVELOG_DEBUG("Robot unbalanced\n");
        return false;
    }
}

bool SoCBirrtPlanner::SetForwardTree(NodeTree* tree) {
	if (tree != NULL && _pForwardTree != tree) {
		if (_pForwardTree != NULL) {
			delete _pForwardTree;
		}
		_pForwardTree = tree;
	}
	return true;
}
bool SoCBirrtPlanner::MergeTree(void) {
	printf("[socbirrt.cpp-1944] Merge Tree.\n");
	if(_pForwardTree->_iConnected < 0 || _pBackwardTree->_iConnected < 0)
		return false;
	RrtNode* _pForwardTreeNode = _pForwardTree->GetNode(_pForwardTree->_iConnected);
	int parentId = -1;
	/*//printf("Reverse forward nodes: ");
	// reverse the path of forward tree
	if(_pForwardTree->GetNode(_pForwardTree->_iConnected)->GetParent() == -1) {
		_pActiveNode = _pForwardTree->GetNode(_pForwardTree->_iConnected);
		// no need reverse as the connect point (next init) is the same with current init.
	}
	else
	{
		int curID = _pForwardTreeNode->GetID();
		parentId = _pForwardTreeNode->GetParent();
		_pForwardTreeNode->DeleteChild(parentId);
		_pForwardTreeNode->SetParent(-1);
		//printf("[%d - %d] ", curID, parentId);
		bool bduplicate = true;
		while(true)
		{
			_pActiveNode = _pForwardTree->GetNode(parentId);

			parentId = _pActiveNode->GetParent();
			_pActiveNode->AddChild(parentId);
			_pActiveNode->DeleteChild(curID);
			_pActiveNode->SetParent(curID);
			curID = _pActiveNode->GetID();
			//printf("[%d - %d] ", curID, parentId);
			if(parentId == -1 || curID == -1)
				break;
		}
	}
	//printf("\n");*/
	printf("Reverse backward nodes: ");
	// find root first and then start from the root, push them into the new tree.
	// add backward tree to path
	/*RrtNode root(_pForwardTreeNode->GetDataVector()->size(), false, 0);
	root.SetConfig(*(_pForwardTreeNode->GetDataVector()));
	root.SetTSRChainValues(_pForwardTreeNode->GetTSRChainValues());
	root.SetParent(-1);

	treenodes->push_back(root);*/
	treenodes->clear();
	treenodes->reserve(_pBackwardTree->GetSize());
	parentId = 0;	// current configuration
	int rootId = 0;
	_pActiveNode = _pBackwardTree->GetNode(_pBackwardTree->_iConnected);
	int duplicateID = _pActiveNode->GetID();
	bool bduplicate = true;
	// find root
	while(true)
	{
		rootId = _pActiveNode->GetID();
		int index = _pActiveNode->GetParent();
		if(index == -1)
			break;
		_pActiveNode = _pBackwardTree->GetNode(index);
	}

	_pActiveNode = _pBackwardTree->GetNode(rootId);
	AddNodeToTreeNodeList(_pBackwardTree, _pActiveNode, -1, -1);	// add child nodes
				/*RrtNode pNode(_pActiveNode->GetDataVector()->size(), false, _pForwardTree->GetSize());
				pNode.SetConfig(*(_pActiveNode->GetDataVector()));
				pNode.SetTSRChainValues(_pActiveNode->GetTSRChainValues());
				pNode.SetParent(parentId);
				(*treenodes)[parentId].AddChild(_pForwardTree->GetSize());
				//printf("%d - ", parentId);
				parentId = treenodes->size();
				treenodes->push_back(pNode);
				AddNodeToTreeNodeList(_pBackwardTree, _pActiveNode, duplicateID, parentId);	// add child nodes*/

	printf("\n");
	//printf("Start Copy");
	//(*treenodes) = _pForwardTree->GetNodes();
	//for (int i=0; i<_pForwardTree->GetNodes().size(); i++) {
	//	treenodes->push_back(_pForwardTree->GetNodes()[i]);
	//}
	//std::copy(_pForwardTree->GetNodes().begin(), _pForwardTree->GetNodes().end(), treenodes->begin());
	//printf("Tree size %d \n", treenodes->size());
	return true;
}
// recursively add back tree nodes to node list
bool SoCBirrtPlanner::AddNodeToTreeNodeList(NodeTree* tree, RrtNode* node, int duplicateID, int parentID) {
	//printf("%d - ", parentID);
	RrtNode pNode(node->GetDataVector()->size(), false, treenodes->size());
	pNode.SetConfig(*(node->GetDataVector()));
	std::vector<dReal> tsrChainValues = node->GetTSRChainValues();
	pNode.SetTSRChainValues(tsrChainValues);
	pNode.SetParent(parentID);
	parentID = treenodes->size();
	//_pForwardTree->AddNode(pNode);
	treenodes->push_back(pNode);
	std::vector<int> cids = node->GetChildren();	// children ids
	for (int i=0; i<cids.size(); i++) {	// add child trees to new tree
		if(cids[i] >= 0 && tree->GetNode(cids[i])->GetID() != duplicateID) {
			AddNodeToTreeNodeList(tree, tree->GetNode(cids[i]), -1, parentID);
		}
	}
	return true;
}
