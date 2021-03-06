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
/** \file cbirrtproblem.cpp
    \brief Implements the cbirrt problem class.
 */
#include "stdafx.h"
#include <queue>
//#define SAVE_DATA


std::vector<GraphHandlePtr> graphptrs;



SoCBirrtProblem::SoCBirrtProblem(EnvironmentBasePtr penv) : ProblemInstance(penv)
{
    __description = ":Interface Author: Dmitry Berenson\nInterface to SoCBirrt planner that parses input and passes it to the planner. Can also call the GeneralIK plugin. \n\n`C++ Documentation <http://automation.berkeley.edu/~berenson/docs/SoCBirrt/index.html>`_";
    RegisterCommand("GrabBody",boost::bind(&SoCBirrtProblem::GrabBody,this,_1,_2),
                    "Robot calls ::Grab on a body with its current manipulator");

    RegisterCommand("RunSoCBirrt",boost::bind(&SoCBirrtProblem::RunSoCBirrt,this,_1,_2),
                    "Run the SoCBirrt Planner");

    RegisterCommand("DoGeneralIK",boost::bind(&SoCBirrtProblem::DoGeneralIK,this,_1,_2),
                    "Calls the General IK solver");

    RegisterCommand("CheckSupport",boost::bind(&SoCBirrtProblem::CheckSupport,this,_1,_2),
                    "Checks whether the cg of the robot is over the given support polygon");

    RegisterCommand("CheckSelfCollision",boost::bind(&SoCBirrtProblem::CheckSelfCollision,this,_1,_2),
                    "Checks whether the robot is in self-collision");

    RegisterCommand("GetJointAxis",boost::bind(&SoCBirrtProblem::GetJointAxis,this,_1,_2),
                    "Returns the joint axis of a given kinbody");

    RegisterCommand("GetJointTransform",boost::bind(&SoCBirrtProblem::GetJointTransform,this,_1,_2),
                    "Returns the joint transform of a given kinbody");

    RegisterCommand("GetCamView",boost::bind(&SoCBirrtProblem::GetCamView,this,_1,_2),
                    "Get the currnet camera view");

    RegisterCommand("SetCamView",boost::bind(&SoCBirrtProblem::SetCamView,this,_1,_2),
                    "Set the currnet camera view");

    RegisterCommand("Traj",boost::bind(&SoCBirrtProblem::Traj,this,_1,_2),
                    "Execute a trajectory from a file on the local filesystem");

    RegisterCommand("GetPlannerState",boost::bind(&SoCBirrtProblem::GetPlannerState,this,_1,_2),
                    "Returns the current state of the planner");

    RegisterCommand("StopPlanner",boost::bind(&SoCBirrtProblem::StopPlanner,this,_1,_2),
                    "Stops the planner if it is running");

    RegisterCommand("ClearDrawn",boost::bind(&SoCBirrtProblem::ClearDrawn,this,_1,_2),
                    "Clears objects drawn by SoCBirrt planner and problem");


    _reusePlanner = false;
    _plannerState = PS_Idle;
    _plannerThread.reset();

}

void SoCBirrtProblem::Destroy()
{

}

SoCBirrtProblem::~SoCBirrtProblem()
{
    Destroy();
}

int SoCBirrtProblem::main(const std::string& cmd)
{
   
    RAVELOG_DEBUG("env: %s\n", cmd.c_str());

    const char* delim = " \r\n\t";
    string mycmd = cmd;
    char* p = strtok(&mycmd[0], delim);
    if( p != NULL )
        _strRobotName = p;

    std::vector<RobotBasePtr> robots;
    GetEnv()->GetRobots(robots);
    SetActiveRobots(robots);

     _pIkSolver = RaveCreateIkSolver(GetEnv(),"GeneralIK");
    _pIkSolver->Init(robot->GetActiveManipulator());
    RAVELOG_INFO("IKsolver initialized\n");
    
    return 0;
}

void SoCBirrtProblem::SetActiveRobots(const std::vector<RobotBasePtr >& robots)
{


    if( robots.size() == 0 ) {
        RAVELOG_WARNA("No robots to plan for\n");
        return;
    }

    vector<RobotBasePtr >::const_iterator itrobot;
    FORIT(itrobot, robots) {
        if( strcmp((*itrobot)->GetName().c_str(), _strRobotName.c_str() ) == 0  ) {
            robot = *itrobot;
            break;
        }
    }

    if( robot == NULL ) {
        RAVELOG_ERRORA("Failed to find %S\n", _strRobotName.c_str());
        return;
    }
}


bool SoCBirrtProblem::SendCommand(std::ostream& sout, std::istream& sinput)
{
    ProblemInstance::SendCommand(sout,sinput);
    return true;
}

bool SoCBirrtProblem::ClearDrawn(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());
    graphptrs.clear();
    return true;
}


int SoCBirrtProblem::CheckSelfCollision(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    if(robot->CheckSelfCollision())
        sout << "1";
    else
        sout << "0";
    
    return true;
}

bool SoCBirrtProblem::DoGeneralIK(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    RAVELOG_DEBUG("Starting General IK...\n");

    TransformMatrix tmtarg;
    TransformMatrix tmrobot;
    Transform Ttarg;

    Vector cogtarg;
    bool bExecute = false;
    bool bBalance = false;
    bool bGetTime = false;
    bool bNoRot = false;
    bool bReturnClosest = false;
    std::vector<string> supportlinks;
    std::vector<dReal> ikparams;
    std::vector<dReal> supportpolyx;
    std::vector<dReal> supportpolyy;
    dReal temp;
    string cmd;
    Vector polyscale(1.0,1.0,1.0);
    Vector polytrans(0,0,0);
    while(!sinput.eof()) {
        sinput >> cmd;
        if( !sinput )
            break;

       if( stricmp(cmd.c_str(), "nummanips") == 0 ) {
            sinput >> temp;
            ikparams.push_back(temp);
       }
       else if( stricmp(cmd.c_str(), "exec") == 0 ) {
            bExecute = true;

       }
       else if( stricmp(cmd.c_str(), "gettime") == 0 ) {
            bGetTime = true;
       }
       else if( stricmp(cmd.c_str(), "norot") == 0 ) {
            bNoRot = true;
       }
       else if( stricmp(cmd.c_str(), "returnclosest") == 0 ) {
            bReturnClosest = true;
       }
       else if( stricmp(cmd.c_str(), "robottm") == 0 ) {
            sinput >> tmrobot.m[0];
            sinput >> tmrobot.m[4];
            sinput >> tmrobot.m[8];
            sinput >> tmrobot.m[1];
            sinput >> tmrobot.m[5];
            sinput >> tmrobot.m[9];
            sinput >> tmrobot.m[2];
            sinput >> tmrobot.m[6];
            sinput >> tmrobot.m[10];
            sinput >> tmrobot.trans.x;
            sinput >> tmrobot.trans.y;
            sinput >> tmrobot.trans.z;
            robot->SetTransform(tmrobot);
       }
       else if( stricmp(cmd.c_str(), "maniptm") == 0 ) {
            sinput >> temp;
            ikparams.push_back(temp);

            sinput >> tmtarg.m[0];
            sinput >> tmtarg.m[4];
            sinput >> tmtarg.m[8];
            sinput >> tmtarg.m[1];
            sinput >> tmtarg.m[5];
            sinput >> tmtarg.m[9];
            sinput >> tmtarg.m[2];
            sinput >> tmtarg.m[6];
            sinput >> tmtarg.m[10];
            sinput >> tmtarg.trans.x;
            sinput >> tmtarg.trans.y;
            sinput >> tmtarg.trans.z;
            Ttarg = Transform(tmtarg);

            ikparams.push_back(Ttarg.rot.x); ikparams.push_back(Ttarg.rot.y); ikparams.push_back(Ttarg.rot.z); ikparams.push_back(Ttarg.rot.w);
            ikparams.push_back(Ttarg.trans.x); ikparams.push_back(Ttarg.trans.y); ikparams.push_back(Ttarg.trans.z);
        }
        else if( stricmp(cmd.c_str(), "movecog") == 0 ) {
            sinput >> cogtarg.x;
            sinput >> cogtarg.y;
            sinput >> cogtarg.z;
            bBalance = true;

        }
        else if(stricmp(cmd.c_str(), "supportlinks") == 0 ){
            int numsupportlinks;
            string tempstring;
            sinput >> numsupportlinks;
            for(int i =0; i < numsupportlinks; i++)
            {
                sinput >> tempstring;
                supportlinks.push_back(tempstring);
            }
        }
        else if(stricmp(cmd.c_str(), "polyscale") == 0 ){
            sinput >> polyscale.x;
            sinput >> polyscale.y;
            sinput >> polyscale.z;
        }
        else if(stricmp(cmd.c_str(), "polytrans") == 0 ){
            sinput >> polytrans.x;
            sinput >> polytrans.y;
            sinput >> polytrans.z;
        }
        else break;


        if( !sinput ) {
            RAVELOG_DEBUG("failed\n");
            return false;
        }
    }
    //don't do cog movement
    if(bBalance == false)
    {
        ikparams.push_back(0);

    }
    else
    {
        ikparams.push_back(1);
        ikparams.push_back(cogtarg.x);
        ikparams.push_back(cogtarg.y);
        ikparams.push_back(cogtarg.z);
            

        if(supportlinks.size() == 0)
        {
            RAVELOG_INFO("ERROR: Must specify support links to do balancing\n");
            return false;
        }
            

        GetSupportPolygon(supportlinks,supportpolyx,supportpolyy,polyscale,polytrans);

        ikparams.push_back(supportpolyx.size());

        for(int i = 0 ; i < supportpolyx.size(); i ++)
            ikparams.push_back(supportpolyx[i]);

        for(int i = 0 ; i < supportpolyy.size(); i ++)
            ikparams.push_back(supportpolyy[i]);

    }

    std::vector<dReal> q0(robot->GetActiveDOF());
    robot->GetActiveDOFValues(q0);
    
    std::vector<dReal> qResult(robot->GetActiveDOF());

    ikparams.push_back(0);//select the mode
    if(bNoRot)
        ikparams.push_back(1);//don't do rotation
    else
        ikparams.push_back(0);//do rotation

    unsigned long starttime = timeGetTime();


    boost::shared_ptr<std::vector<dReal> > pqResult(new std::vector<dReal> );
    //PrintMatrix(&ikparams[0], 1, ikparams.size(), "ikparams");

    if(_pIkSolver->Solve(IkParameterization(), q0, ikparams, false, pqResult) )
    {

        qResult = *pqResult.get();
        int timetaken = timeGetTime() - starttime;
        for(int i = 0; i < qResult.size(); i++)
        {
            sout << qResult[i] << " ";

        }
        if(bGetTime)
            sout << timetaken << " ";

        RAVELOG_INFO("Solution Found! (%d ms)\n",timetaken);
        if(bExecute)
            robot->SetActiveDOFValues(qResult);
        return true;
    }
    else
    {
        int timetaken = timeGetTime() - starttime;
        if(bReturnClosest)
        {
            qResult = *pqResult.get();
            for(int i = 0; i < qResult.size(); i++)
            {
                sout << qResult[i] << " ";
            }
        }

        if(bGetTime)
            sout << timetaken << " ";
        RAVELOG_INFO("No IK Solution Found (%d ms)\n",timetaken);
        return true;
    }
}


bool SoCBirrtProblem::CheckSupport(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    RAVELOG_DEBUG("Starting Check Support...\n");


    std::vector<dReal> polyx;
    std::vector<dReal> polyy;
    int nvert = 0;
    int numsupportlinks = 0;
    string cmd;
    std::vector<string> supportlinks;
    KinBodyPtr pheld;
    string tempstring;
    dReal temp;
    bool bdraw = false;
    Vector polyscale(1.0,1.0,1.0);
    Vector polytrans(0,0,0);
    while(!sinput.eof()) {
        sinput >> cmd;
        if( !sinput )
            break;

        if(stricmp(cmd.c_str(), "supportlinks") == 0 ){
            sinput >> numsupportlinks;
            for(int i =0; i < numsupportlinks; i++)
            {
                sinput >> tempstring;
                supportlinks.push_back(tempstring);
            }

        }
        else if( stricmp(cmd.c_str(), "heldobject") == 0 ) {
            sinput >> tempstring;
            pheld = GetEnv()->GetKinBody(tempstring.c_str());
            if(pheld == NULL)
            {
                RAVELOG_INFO("Error: could not find the specified held object\n");
                return false;
            }
        }
        else if(stricmp(cmd.c_str(), "polyscale") == 0 ){
            sinput >> polyscale.x;
            sinput >> polyscale.y;
            sinput >> polyscale.z;
        }
        else if(stricmp(cmd.c_str(), "polytrans") == 0 ){
            sinput >> polytrans.x;
            sinput >> polytrans.y;
            sinput >> polytrans.z;
        }
        else if(stricmp(cmd.c_str(), "draw") == 0 ){
            bdraw = true;
        }
        else break;


        if( !sinput ) {
            RAVELOG_DEBUG("failed\n");
            return false;
        }
    }
    
    GetSupportPolygon(supportlinks,polyx,polyy,polyscale,polytrans);
    int numPointsOut = polyx.size();

    Vector center;  
    dReal fTotalMass = 0;
    std::vector<KinBody::LinkPtr>::const_iterator itlink;

    std::vector<KinBody::LinkPtr> vlinks = robot->GetLinks();
    if(pheld != NULL)
    {
        for(int i = 0; i < pheld->GetLinks().size(); i++)
            vlinks.push_back(pheld->GetLinks()[i]);
    }

    FORIT(itlink, vlinks) {
        //RAVELOG_INFO("Name %s comoffset: %f %f %f\n", (*itlink)->GetName().c_str(), (*itlink)->GetCOMOffset().x,(*itlink)->GetCOMOffset().y,(*itlink)->GetCOMOffset().z);
        if(bdraw)
        {
            GetEnv()->plot3(&(DoubleVectorToFloatVector((*itlink)->GetTransform() * (*itlink)->GetCOMOffset())[0]), 1, 0, (*itlink)->GetMass()/50, Vector(0,0,1),1 );
        }
        center += ((*itlink)->GetTransform() * (*itlink)->GetCOMOffset() * (*itlink)->GetMass());
        fTotalMass += (*itlink)->GetMass();
    }

    if( fTotalMass > 0 )
        center /= fTotalMass;
    RAVELOG_INFO("\nmass: %f\ncog: %f %f %f\n",fTotalMass,center.x,center.y,center.z);


    dReal testx = center.x;
    dReal testy = center.y;
    dReal * vertx = &polyx[0];
    dReal * verty = &polyy[0];
    nvert = numPointsOut;

    int i, j, c = 0;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
    }
   
    Vector cgline[2];
    cgline[0] = center;
    cgline[1] = center;
    cgline[1].z = 0;


    std::vector<RaveVector<float> > fcgline;
    for(int i =0; i < 2; i++)
        fcgline.push_back(RaveVector<float>(cgline[i].x,cgline[i].y,cgline[i].z));

    if(c)    
    {
        //GetEnv()->plot3(center, 1, 0, 0.03, Vector(0,1,0),1 );
        if(bdraw)
        {
            GetEnv()->drawlinestrip(&(fcgline[0].x),2,sizeof(RaveVector<float>(0, 1, 0, 0)),5, RaveVector<float>(0, 1, 0, 0));
        }
        RAVELOG_INFO("Supported\n");
        sout << "1";
    }
    else
    {
        if(bdraw)
        {
            GetEnv()->drawlinestrip(&(fcgline[0].x),2,sizeof(RaveVector<float>(0, 1, 0, 0)),5, RaveVector<float>(1, 0, 0, 0));
        }
        //GetEnv()->plot3(center, 1, 0, 0.03, Vector(1,0,0),1 );
        RAVELOG_INFO("Not Supported\n");
        sout << "0";
    }

    return true;


}


void SoCBirrtProblem::PlannerWorker(PlannerBasePtr _pTCplanner, TrajectoryBasePtr ptraj, string filename)
{
    RAVELOG_INFO("Creating lock in new thread...\n");
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());
    RAVELOG_INFO("Lock created in new thread.\n");
    _plannerState = PS_Planning;
    OpenRAVE::PlannerStatus bSuccess = _pTCplanner->PlanPath(ptraj);
    //NOTE: _plannerState is a different class than OpenRAVE::PlannerStatus (includes a "Planning" state)
    if(bSuccess == PS_HasSolution)
    {
        WriteTraj(ptraj, filename);
        _plannerState = PS_PlanSucceeded;
    }
    else
    {
        _plannerState = PS_PlanFailed;
    }

    /* Reset joint limits */
    for (int j=0; j<_limadj_joints.size(); j++) 
    { 
        KinBody::JointPtr joint;
        vector<dReal> j_lower;
        vector<dReal> j_upper;
        joint = _limadj_joints[j];
        j_lower = _limadj_lowers[j];
        j_upper = _limadj_uppers[j];
        joint->SetLimits(j_lower, j_upper); 
    }

    RAVELOG_INFO("Planner worker thread terminating.\n");
}

void SoCBirrtProblem::WriteTraj(TrajectoryBasePtr ptraj, string filename)
{
    //TrajectoryBasePtr pfulltraj = RaveCreateTrajectory(GetEnv(),robot->GetDOF());
    //robot->GetFullTrajectoryFromActive(pfulltraj, ptraj);

    // save the constrained trajectory

    //pfulltraj->CalcTrajTiming(robot, pfulltraj->GetInterpMethod(), true, false);
    //OpenRAVE::planningutils::RetimeActiveDOFTrajectory(ptraj, robot,false,1,"LinearTrajectoryRetimer");

    //print out groups
//    int numgroups = ptraj->GetConfigurationSpecification()._vgroups.size();
//    RAVELOG_INFO("num groups in traj: %d\n",numgroups);
//    for(int i = 0; i < numgroups; i++)
//    {
//        RAVELOG_INFO("name %d: %s\n",i,ptraj->GetConfigurationSpecification()._vgroups[i].name.c_str());
//    }



    //this is the equivalent of GetFullTrajectoryFromActive
    ConfigurationSpecification activespec = ptraj->GetConfigurationSpecification();
    std::set<KinBodyPtr> sbodies;
    FOREACH(itgroup, activespec._vgroups) {
      stringstream ss(itgroup->name);
      string type, bodyname;
      ss >> type;
      ss >> bodyname;
      //RAVELOG_INFO("Bodyname: %s\n",bodyname.c_str());
      if(GetEnv()->GetKinBody(bodyname))
      {
        sbodies.insert(GetEnv()->GetKinBody(bodyname));
      }
    }

    ConfigurationSpecification spec;
    set <KinBodyPtr>::iterator si;
    for (si=sbodies.begin(); si!=sbodies.end(); si++)
    {
        std::vector<int> indices((*si)->GetDOF());
        for(int i = 0; i < (*si)->GetDOF(); ++i) {
            indices[i] = i;
        }
        ConfigurationSpecification tempspec = (*si)->GetConfigurationSpecificationIndices(indices, "linear");
        tempspec.AddDerivativeGroups(1,true); // add velocity + timestamp
        spec += tempspec;
    }

    OpenRAVE::planningutils::ConvertTrajectorySpecification(ptraj, spec);


    ofstream outfile(filename.c_str(),ios::out);
    outfile.precision(16); 
    //pfulltraj->Write(outfile, Trajectory::TO_IncludeTimestamps|Trajectory::TO_IncludeBaseTransformation);
    ptraj->serialize(outfile);
    outfile.close();
    chmod(filename.c_str(), S_IRWXG | S_IRWXO | S_IRWXU); //chmod 777
}


int SoCBirrtProblem::GetPlannerState(ostream& sout, istream& sinput)
{
    PlannerState _currentState = _plannerState; //this should be atomic, so don't need a mutex

    switch (_currentState) {
    case PS_Idle:
        sout << "idle";
        break;
    case PS_Planning:
        sout << "planning";
        break;
    case PS_PlanSucceeded:
        sout << "plansucceeded";
        break;
    case PS_PlanFailed:
        sout << "planfailed";
        break;
    default:
        sout << "ERROR: UNKNOWN PLANNER STATE!";
        return 0;
    }

    return 1;
}

int SoCBirrtProblem::StopPlanner(ostream& sout, istream& sinput)
{

    PlannerState _currentState = _plannerState; //this should be atomic, so don't need a mutex

    switch (_currentState) {
    case PS_Idle:
        sout << "ERROR: PLANNER IS NOT CURRENTLY PLANNING!";
        return 0;
    case PS_Planning:
        //terminate the planner
        _plannerState = PS_PlanFailed; //this should be atomic, so don't need a mutex
        sout << "planner stopped";
        break;
    case PS_PlanSucceeded:
        sout << "planner already finished";
        break;
    case PS_PlanFailed:
        sout << "planner already finished";
        break;
    default:
        sout << "ERROR: UNKNOWN PLANNER STATE!";
        return 0;
    }

    return 1;
}

int SoCBirrtProblem::RunSoCBirrt(ostream& sout, istream& sinput)
{
    //lock mutex

    boost::shared_ptr<SoCBirrtParameters> params;
    params.reset(new SoCBirrtParameters());
    std::vector<dReal> goals;
    std::vector<dReal> starts;

    bool bPlanInNewThread = false;
    bool bSmoothTrajOnly = false;
    int nvert = 0;
    int numsupportlinks = 0;
    std::vector<string> supportlinks;
    KinBodyPtr pheld;
    string tempstring;
    Vector polyscale(1.0,1.0,1.0);
    Vector polytrans(0,0,0);
    string cmd;
    bool bAllowLimAdj = false;
    KinBody::JointPtr joint;
    vector<dReal> j_lower;
    vector<dReal> j_upper;
    
    double startTime = 0.0;
    double endTime = 0.0;
    double lastTime = 0.0;
    double sensorDeltaTime = 0.030; // 30 milisecond per frame
    double velMag = 0.001;	// magnitude of velocity
    double goalthreshold = 0.15;	// palm size, safe area for the goal

    string filename = "cmovetraj.txt";
    string smoothtrajfilename;
    params->Tattachedik_0.resize(robot->GetManipulators().size());

    string goalObjectName;
    int isMergeTree = 0;
    _planning_theta.resize(4, 0);
    _planning_theta[0] = 2;
    _execution_theta.resize(4, 0);
    _execution_theta[0] = 0.11;
    int isScreenShot = 0;
    dReal thickness = 0.05;	// size of box for sampling
    while(!sinput.eof()) {
        sinput >> cmd;
        if( !sinput )
            break;

        if( stricmp(cmd.c_str(), "jointgoals") == 0 ) {
            //note that this appends to goals, does not overwrite them
            int temp;
            sinput >> temp;           
            int oldsize = goals.size();
            goals.resize(oldsize+temp);
            for(size_t i = oldsize; i < oldsize+temp; i++)
            {            
                sinput >> goals[i];
            }   
            
        }
        else if( stricmp(cmd.c_str(), "jointstarts") == 0 ) {
            //note that this appends to starts, does not overwrite them
            int temp;
            sinput >> temp;
            int oldsize = starts.size();
            starts.resize(oldsize+temp);
            for(size_t i = oldsize; i < oldsize+temp; i++)
            {            
                sinput >> starts[i];
            }   
            
        }
        else if( stricmp(cmd.c_str(), "ikguess") == 0 ) {
            int temp;
            sinput >> temp;           
            params->vikguess.resize(temp);
            for(size_t i = 0; i < temp; i++)
            {            
                sinput >> params->vikguess[i];
            }   
            
        }
        else if( stricmp(cmd.c_str(), "Tattachedik_0") == 0) {
            int manipind;
            sinput >> manipind;
            TransformMatrix tmattachedik_0;
            sinput >> tmattachedik_0.m[0];
            sinput >> tmattachedik_0.m[4];
            sinput >> tmattachedik_0.m[8];
            sinput >> tmattachedik_0.m[1];
            sinput >> tmattachedik_0.m[5];
            sinput >> tmattachedik_0.m[9];
            sinput >> tmattachedik_0.m[2];
            sinput >> tmattachedik_0.m[6];
            sinput >> tmattachedik_0.m[10];
            sinput >> tmattachedik_0.trans.x;
            sinput >> tmattachedik_0.trans.y;
            sinput >> tmattachedik_0.trans.z;
            params->Tattachedik_0[manipind] = Transform(tmattachedik_0);
        }
        else if( stricmp(cmd.c_str(), "smoothingitrs") == 0) {
            sinput >> params->smoothingitrs;
        }
        else if( stricmp(cmd.c_str(), "filename") == 0) {
            sinput >> filename;
        }
        else if( stricmp(cmd.c_str(), "timelimit") == 0) {
            sinput >> params->timelimit;
        }
        else if( stricmp(cmd.c_str(), "TSRChain") == 0) {
            TaskSpaceRegionChain temp;
            temp.deserialize_from_matlab(robot,GetEnv(),sinput);
            params->vTSRChains.push_back(temp);
            if(temp.IsForStartSampling())
            {
                RAVELOG_INFO("Doing Start sampling\n");
                params->bsamplingstart = true;
            }
            if(temp.IsForGoalSampling())
            {
                RAVELOG_INFO("Doing Goal sampling\n");
                params->bsamplinggoal = true;
            }
        }
        else if(stricmp(cmd.c_str(), "supportlinks") == 0 ){
            sinput >> numsupportlinks;
            for(int i =0; i < numsupportlinks; i++)
            {
                sinput >> tempstring;
                supportlinks.push_back(tempstring);
            }
        }
        else if(stricmp(cmd.c_str(), "polyscale") == 0 ){
            sinput >> polyscale.x;
            sinput >> polyscale.y;
            sinput >> polyscale.z;
        }
        else if(stricmp(cmd.c_str(), "polytrans") == 0 ){
            sinput >> polytrans.x;
            sinput >> polytrans.y;
            sinput >> polytrans.z;
        }
        else if( stricmp(cmd.c_str(), "heldobject") == 0 ) {
            sinput >> tempstring;
            pheld = GetEnv()->GetKinBody(tempstring.c_str());
            if(pheld == NULL)
            {
                RAVELOG_INFO("Error: could not find the specified held object\n");
                sout << 0;
                return -1;
            }
            robot->Grab(pheld);
            robot->CheckSelfCollision();
            params->bgrabbed = true;
        }
        else if( stricmp(cmd.c_str(), "psample") == 0)
        {
            sinput >> params->Psample;

        }
        else if( stricmp(cmd.c_str(), "bikfastsinglesolution") == 0)
        {
               sinput >> params->bikfastsinglesolution;
        }
        else if( stricmp(cmd.c_str(), "allowlimadj") == 0 )
        { 
               sinput >> bAllowLimAdj; 
        }
        else if( stricmp(cmd.c_str(), "sensorfps") == 0 )
		{
			   sinput >> sensorDeltaTime;
		}
        else if( stricmp(cmd.c_str(), "goalvelocitymagnitude") == 0 )
		{
			   sinput >> velMag;
		}
        else if( stricmp(cmd.c_str(), "goalobject") == 0 )
		{
			   sinput >> goalObjectName;
		}
        else if( stricmp(cmd.c_str(), "planning_theta") == 0 )
		{
        	int num = 4;
			   sinput >> num;
			   cout << "planning params: " << num;
			   _planning_theta.resize(num);
			   for (int i=0; i<num; i++) {
				   sinput >> _planning_theta[i];
				   cout << " " << _planning_theta[i];
			   }
			   cout << std::endl;
		}
        else if( stricmp(cmd.c_str(), "execution_theta") == 0 )
		{
			int num = 4;
			   sinput >> num;
			   cout << "execution params: " << num;
			   _execution_theta.resize(num);
			   for (int i=0; i<num; i++) {
				   sinput >> _execution_theta[i];
				   cout << " " << _execution_theta[i];
			   }
			   cout << std::endl;
		}
        else if( stricmp(cmd.c_str(), "goalthreshold") == 0 )
		{
			   sinput >> goalthreshold;
		}
        else if( stricmp(cmd.c_str(), "mergetree") == 0 )
		{
        	printf("[*********]merge tree started\n");
			   sinput >> isMergeTree;
			   if(isMergeTree > 0) {
				   params->bmergetree = true;
				   printf("[*********]merge tree is true\n");
			   }
			   else
				   params->bmergetree = false;
		}
        else if( stricmp(cmd.c_str(), "screenshot") == 0 ){
        	sinput >> isScreenShot;
        }
        else if( stricmp(cmd.c_str(), "thickness") == 0 ){
			sinput >> thickness;
		}
        else break;
        if( !sinput ) {
            RAVELOG_DEBUG("failed\n");
            sout << 0;
            return -1;
        }
    }

    if(supportlinks.size() != 0)
        GetSupportPolygon(supportlinks,params->vsupportpolyx,params->vsupportpolyy,polyscale,polytrans);


    if(starts.size() == 0 && !params->bsamplingstart)
    {
        RAVELOG_INFO("Setting default init config to current config\n");
        //default case: when not sampling starts and no starts specified, use current config as start
        params->vinitialconfig.resize(robot->GetActiveDOF());
        robot->GetActiveDOFValues(params->vinitialconfig);  
    }
    else
    {
        //add any starts if they were specified
        for(int i = 0; i < starts.size(); i++)
            params->vinitialconfig.push_back(starts[i]);
    }


    //add any goals if they were specified
    for(int i = 0; i < goals.size(); i++)
        params->vgoalconfig.push_back(goals[i]);


    //boost::shared_ptr<SoCBirrtPlanner> _pTCplanner;

        
    OpenRAVE::PlannerStatus bSuccess = PS_Failed;
    u32 startplan = timeGetTime();
    RAVELOG_DEBUG("starting planning\n");
    
    //for error reporting
    stringstream outputstream;
    stringstream command;
    command << "GetOutputMessage";

    /* Send over pointer to our planner state */
    params->pplannerstate = &_plannerState;


    /* Fix joint limits if current values are outside! */ 
    _limadj_joints.clear();
    _limadj_lowers.clear(); 
    _limadj_uppers.clear();
    if (bAllowLimAdj)
    {
        vector<dReal> start;
        vector<dReal> lower;
        vector<dReal> upper;
        vector<int> inds;
        vector<dReal> j_start;
        
        robot->GetActiveDOFValues(start);
        robot->GetActiveDOFLimits(lower, upper);
        inds = robot->GetActiveDOFIndices();
        
        /* Save all current joint values */
        for (int i=0; i<start.size(); i++) if (start[i] < lower[i] || start[i] > upper[i])
        {
            RAVELOG_INFO("Temporarily adjusting joint limit [%d].\n", i);
            joint = robot->GetJointFromDOFIndex(inds[i]);
            joint->GetLimits(j_lower, j_upper);
            _limadj_joints.push_back(joint);
            _limadj_lowers.push_back(j_lower);
            _limadj_uppers.push_back(j_upper);
        }
        for (int i=0; i<start.size(); i++) if (start[i] < lower[i] || start[i] > upper[i])
        {
            joint = robot->GetJointFromDOFIndex(inds[i]);
            joint->GetValues(j_start);
            joint->GetLimits(j_lower, j_upper);
            for (int j=0; j<j_start.size(); j++)
            { 
                if (j_start[j] < j_lower[j]) j_lower[j] = j_start[j];
                if (j_start[j] > j_upper[j]) j_upper[j] = j_start[j];
            }
            joint->SetLimits(j_lower, j_upper);
        }
    }

    _plannerState = PS_Idle;
    startTime = timeGetThreadTime();
    endTime = startTime;
    lastTime = endTime;


    double sigma = 0.0;  //
    SensorConfiguration sc(GetEnv(), sensorDeltaTime, velMag, goalObjectName);
    std::vector<dReal> curConfig, goalConfig;
    robot->GetActiveDOFValues(curConfig);

    std::vector<dReal> dofVals;
    std::vector<TaskSpaceRegionChain> goalChains = params->vTSRChains;
    int curGoalChainId = 0;
    RaveVector<dReal> goal = goalChains[0].TSRChain[0].T0_w.trans;
    sc.SetInitGoal(goal);
    SoCBirrtPlanner::treenodes = new std::vector<RrtNode>();
    OpenRAVE::KinBodyPtr ball;
    OpenRAVE::KinBodyPtr boxes;
    std::vector<OBB> obbs;
    bool reachGoal = false;
    InitializeRotationMatrixes();
#ifdef SAVE_DATA
    std::ofstream fs;
    fs.open("data.txt", std::ofstream::out | std::ofstream::app);
    if(!fs) {	// file not exist, create one
    	fs.open("data.txt", std::ofstream::out);
    }
#endif
    int iterCounter = 0;
    while (endTime - startTime < params->timelimit && !reachGoal) {
    	printf("[socbirrtproblem.cpp] %d th iteration started\n", iterCounter);
    	// Generate new goal
    	RaveVector<dReal> newGoal;
    	dReal transDist;
    	std::vector<int> topNNodeIds;
		std::vector<dReal> topNNodeDists;
#ifdef SAVE_DATA
    	dReal planningStartTime, planningEndTime, executeStartTime, executeEndTime;
    	planningStartTime = timeGetThreadTime();
#endif
    	dReal sigma;
    	{
    	    EnvironmentMutex::scoped_lock envlock(GetEnv()->GetMutex());
			// Estimate planning & execution time
    	    if(ball) {
    	    	GetEnv()->Remove(ball);
    	    	ball->Destroy();
    	    }
    	    if(boxes) {
				GetEnv()->Remove(boxes);
				boxes->Destroy();
    	    }
    	    newGoal = sc.ReadSensorData(endTime);
    	    TaskSpaceRegionChain tsrc;
			tsrc.SetParameters(false, true, false);
			tsrc.SetMimicBodyName("NULL");
			TaskSpaceRegion tsr;
			tsr.T0_w.trans = newGoal;
			tsr.T0_w.rot.Set4(0, 0, 0, 1);
			tsr.Bw[3][0] = goalChains[0].TSRChain[0].Bw[3][0];//-PI;
			tsr.Bw[3][1] = goalChains[0].TSRChain[0].Bw[3][1];//PI;
			tsr.Bw[4][0] = goalChains[0].TSRChain[0].Bw[4][0];//-PI;
			tsr.Bw[4][1] = goalChains[0].TSRChain[0].Bw[4][1];//PI;
			tsr.Bw[5][0] = goalChains[0].TSRChain[0].Bw[5][0];//-PI;
			tsr.Bw[5][1] = goalChains[0].TSRChain[0].Bw[5][1];//PI;

			robot->GetActiveDOFValues(curConfig);

			double estimated_time = 0.0;
			// Calculate sigma
			// first get transform of end effector
	    	Transform T0_s = robot->GetActiveManipulator()->GetEndEffectorTransform();
			std::vector<dReal> dx;
			dReal distMag = tsr.DistanceToTSR(T0_s, dx);
			transDist = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
			topNNodeIds.resize(5, 0);
			topNNodeDists.resize(5, -1.0);
			GenerateNodesWorkSpaceData(newGoal, goalthreshold, topNNodeIds, topNNodeDists);	// generate work space value for tree nodes
			// find nearest points
#ifdef SAVE_DATA
			fs << transDist << "\t" << SoCBirrtPlanner::treenodes->size() << "\t";
			for (int ndi=0; ndi<topNNodeDists.size(); ndi++)
				fs << topNNodeDists[ndi] << "\t";
#endif
			dReal planningTime = EstimatePlanningTime(transDist, topNNodeDists[0], curConfig);
			dReal executionTime = EstimateExecutionTime(transDist, topNNodeDists[0], curConfig);
			sigma = (planningTime + executionTime) / sensorDeltaTime * sc.GetVelocityMagitude();
			printf("[socbirrtproblem.cpp-987] sigma: %f, dist:%f\n", sigma, transDist);

			if(sigma < goalthreshold) {
				tsr.Tw_e.trans.Set3(0, 0, 0.05);

				tsr.Bw[0][0] = -sigma;
				tsr.Bw[0][1] = sigma;
				tsr.Bw[1][0] = -sigma;
				tsr.Bw[1][1] = sigma;
				tsr.Bw[2][0] = -sigma;
				tsr.Bw[2][1] = sigma;
				reachGoal = true;
#ifdef SAVE_DATA
				fs << 4 / 3 * PI * sigma * sigma * sigma << "\t";	// volume
#endif
			}
			else {
				RaveVector<dReal> dir, xaxis, axis;
				dir.Set3(dx[0]/transDist, dx[1]/transDist, dx[2]/transDist);
				xaxis.Set3(1, 0, 0);
				// Set Bound, may change ball to a surface
				tsr.Tw_e.trans.Set3(dir.x * (sigma + thickness*0.25), dir.y * (sigma + thickness), dir.z * (sigma + thickness) + 0.05);
				dReal theta = acos(dir.dot3(xaxis));
				axis = dir.Cross(xaxis);
				axis = axis.normalize3() * sin(theta/2);


				tsr.Tw_e.rot.Set4(axis.x, axis.y, axis.z, cos(theta/2));
				tsr.Bw[0][0] = -thickness;
				tsr.Bw[0][1] = thickness;
				tsr.Bw[1][0] = -0.5 * sigma;
				tsr.Bw[1][1] = 0.5 * sigma;
				tsr.Bw[2][0] = -0.5 * sigma;
				tsr.Bw[2][1] = 0.5 * sigma;
			}
			tsr.relativelinkname = "NULL";
			tsr.relativebodyname = "NULL";
			tsr.manipind = 0;
			tsr.Initialize(GetEnv());


			tsrc.AddTSR(tsr);
			tsrc.Initialize(GetEnv());

			std::vector<TaskSpaceRegionChain> tsrChains;
			tsrChains.push_back(tsrc);
			if(sigma >= goalthreshold) {
				for(int ti=0; ti<4; ti++) {
					TaskSpaceRegionChain tmptsrc;
					tmptsrc.SetParameters(false, true, false);
					tmptsrc.SetMimicBodyName("NULL");
					TaskSpaceRegion tmptsr;
					tmptsr.T0_w.trans = newGoal;
					tmptsr.T0_w.rot.Set4(0, 0, 0, 1);
					tmptsr.Bw[3][0] = goalChains[0].TSRChain[0].Bw[3][0];//-PI;
					tmptsr.Bw[3][1] = goalChains[0].TSRChain[0].Bw[3][1];//PI;
					tmptsr.Bw[4][0] = goalChains[0].TSRChain[0].Bw[4][0];//-PI;
					tmptsr.Bw[4][1] = goalChains[0].TSRChain[0].Bw[4][1];//PI;
					tmptsr.Bw[5][0] = goalChains[0].TSRChain[0].Bw[5][0];//-PI;
					tmptsr.Bw[5][1] = goalChains[0].TSRChain[0].Bw[5][1];//PI;
					RaveVector<dReal> dir, xaxis, axis;
					dir.Set3(dx[0]/transDist, dx[1]/transDist, dx[2]/transDist);
					dir = MultiplyMatrix(ti, dir);
					xaxis.Set3(1, 0, 0);
					tmptsr.Tw_e.trans.Set3(dir.x * (sigma + thickness*0.25), dir.y * (sigma + thickness), dir.z * (sigma + thickness) + 0.05);
					dReal theta = acos(dir.dot3(xaxis));
					axis = dir.Cross(xaxis);
					axis = axis.normalize3() * sin(theta/2);


					tmptsr.Tw_e.rot.Set4(axis.x, axis.y, axis.z, cos(theta/2));
					tmptsr.Bw[0][0] = -thickness;
					tmptsr.Bw[0][1] = thickness;
					tmptsr.Bw[1][0] = -0.5 * sigma;
					tmptsr.Bw[1][1] = 0.5 * sigma;
					tmptsr.Bw[2][0] = -0.5 * sigma;
					tmptsr.Bw[2][1] = 0.5 * sigma;
					tmptsr.relativelinkname = "NULL";
					tmptsr.relativebodyname = "NULL";
					tmptsr.manipind = 0;
					tmptsr.Initialize(GetEnv());
					tmptsrc.AddTSR(tmptsr);
					tsrChains.push_back(tmptsrc);
				}
#ifdef SAVE_DATA
				fs << 5 * sigma * sigma * thickness * 2 << "\t";
#endif
			}


			//tsrChains.push_back(goalChains[curGoalChainId]);
			//sc.SetGoalTranform(goalChains[curGoalChainId].TSRChain[0].T0_w);
			curGoalChainId++;

			bool doSampling = false;
			//printf("[socbirrtproblem.cpp-RunSoCBiRRT-1029] Start Config");
			//for (int p = 0; p<curConfig.size(); p++) {
			//	printf(" %f", curConfig[p]);
			//}
			//printf("\n");
			params->vTSRChains = tsrChains;
			params->vinitialconfig = curConfig;
			params->vgoalconfig = goalConfig;
			params->vikguess.clear();
#ifdef SAVE_DATA
			params->timelimit = 10;
#else
			params->timelimit = planningTime;
#endif
			printf("[socbirrtproblem.cpp-RunSoCBiRRT-1034] Planning Started\n");
			PlannerBasePtr _pTCplanner;
			_pTCplanner = RaveCreatePlanner(GetEnv(),"SoCBirrt");
			if( _pTCplanner == NULL ) {
				RAVELOG_INFO("failed to create planner\n");
				sout << 0;
				return -1;
			}

			if( !_pTCplanner->InitPlan(robot, params) ) {
				RAVELOG_INFO("InitPlan failed\n");
				_pTCplanner->SendCommand(outputstream,command);
				sout << 0 << " InitPlan failed\n" << outputstream.str();
				reachGoal = false;
				continue;
			}

			TrajectoryBasePtr ptraj = RaveCreateTrajectory(GetEnv(),"");
			ptraj->Init(robot->GetActiveConfigurationSpecification());

			bSuccess = _pTCplanner->PlanPath(ptraj);
			printf("[socbirrtproblem.cpp-RunSoCBiRRT-1046] Planning Finished\n");

			if(bSuccess == PS_HasSolution)
				_plannerState = PS_PlanSucceeded;
			else {
				_plannerState = PS_PlanFailed;
				reachGoal = false;
			}
			_plannerState = PS_Idle;

			/* Reset joint limits */
			printf("[socbirrtproblem.cpp-RunSoCBiRRT-1101] Reset joint limits\n");
			for (int j=0; j<_limadj_joints.size(); j++)
			{
				joint = _limadj_joints[j];
				j_lower = _limadj_lowers[j];
				j_upper = _limadj_uppers[j];
				joint->SetLimits(j_lower, j_upper);
			}

			if(bSuccess != PS_HasSolution)
			{
				//_pTCplanner->SendCommand(outputstream,command);
				printf("[socbirrtproblem.cpp-RunSoCBiRRT-1120] No Solution\n");
				//sout << 0 << " " << outputstream.str();
			}
			else {
				printf("[socbirrtproblem.cpp-RunSoCBiRRT-1120] Set Path Start\n");
				robot->GetController()->SetPath(ptraj);
				printf("[socbirrtproblem.cpp-RunSoCBiRRT-1120] Set Path End\n");
			}
/*
			obbs.clear();
			for(int ti=0; ti<tsrc.TSRChain.size(); ti++) {
				OBB obb;
				obb.pos = tsrc.TSRChain[ti].Tw_e.trans;
				obb.extents.Set3(tsrc.TSRChain[ti].Bw[0][0], tsrc.TSRChain[ti].Bw[1][0], tsrc.TSRChain[ti].Bw[2][0]);
				geometry::RaveTransformMatrix<dReal> mat = geometry::matrixFromQuat(tsrc.TSRChain[ti].Tw_e.rot);
				obb.right.Set3(mat.m[0], mat.m[1], mat.m[2]);
				obb.dir.Set3(mat.m[4], mat.m[5], mat.m[6]);
				obb.up.Set3(mat.m[8], mat.m[9], mat.m[10]);
				obbs.push_back(obb);
			}
*/
    	}

#ifdef SAVE_DATA
    	planningEndTime = timeGetThreadTime();
    	executeStartTime = planningEndTime;
#endif

    	ball = RaveCreateKinBody(GetEnv(),"");
		std::vector<Vector> spheres;
		spheres.push_back(Vector(newGoal.x, newGoal.y, newGoal.z, sigma));
		ball->InitFromSpheres(spheres, true);
		ball->SetName("ball");

		//KinBody::LinkPtr link = ball->GetLinks()[0];
		//link->GetGeometry(0)->SetDiffuseColor(Vector(0.0,1.0,0.0));
		GetEnv()->Add(ball);
/*		boxes = RaveCreateKinBody(GetEnv(),"");
		boxes->InitFromBoxes(obbs, true);
		boxes->SetName("boxes");
		GetEnv()->Add(boxes);
*/
    	while(bSuccess == PS_HasSolution && !robot->GetController()->IsDone()) {
    		 boost::this_thread::sleep(boost::posix_time::milliseconds(1));
    	}
#ifdef SAVE_DATA
    	executeEndTime = timeGetThreadTime();
    	fs << planningEndTime - planningStartTime << "\t" << executeEndTime - executeStartTime << std::endl;
#endif
    	endTime = timeGetThreadTime();
    	double deltTime = endTime - lastTime;
    	lastTime = endTime;

    	printf("[socbirrtproblem.cpp] %d th iteration end\n", iterCounter);
    	iterCounter++;
    	if(isScreenShot)
    		getchar();
    	//endTime = timeGetThreadTime();
    	//lastTime = endTime;
    }
    delete SoCBirrtPlanner::treenodes;
#ifdef SAVE_DATA
    fs.clear();
    fs.close();
#endif

    //WriteTraj(ptraj, filename);
    //_pTCplanner->SendCommand(outputstream,command);
    //sout << 1 << " " << outputstream.str();
    return 1;
}



string getfilename_withseparator(istream& sinput, char separator)
{
    string filename;
    if( !getline(sinput, filename, separator) ) {
        // just input directly
        RAVELOG_ERRORA("filename not terminated with ';'\n");
        sinput >> filename;
    }

    // trim leading spaces
    size_t startpos = filename.find_first_not_of(" \t");
    size_t endpos = filename.find_last_not_of(" \t");

    // if all spaces or empty return an empty string  
    if(( string::npos == startpos ) || ( string::npos == endpos))
        return "";

    filename = filename.substr( startpos, endpos-startpos+1 );
    return filename;
}

bool SoCBirrtProblem::GetJointAxis(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    RAVELOG_DEBUG("Starting GetJointAxis...\n");

    KinBodyPtr pobject;
    int jointind = -1;

    dReal temp;
    string cmd, tempstring;

    while(!sinput.eof()) {
        sinput >> cmd;
        if( !sinput )
            break;

        if( stricmp(cmd.c_str(), "name") == 0 ) {
            sinput >> tempstring;
            pobject = GetEnv()->GetKinBody(tempstring.c_str());
            if(pobject == NULL)
            {
                RAVELOG_INFO("Error: could not find the specified object to get axis\n");
                return false;
            }
        }
        else if( stricmp(cmd.c_str(), "jointind") == 0 ) {
            sinput >> jointind;
        }
        else break;

        if( !sinput ) {
            RAVELOG_DEBUG("failed\n");
            return false;
        }
    }


    std::vector<KinBody::JointPtr> _vecjoints = pobject->GetJoints();

    if(jointind >= _vecjoints.size() || jointind < 0)
    {   
        RAVELOG_INFO("Joint index out of bounds or not specified\n");
        return false;
    }

    Vector axis = _vecjoints[jointind]->GetAxis(0);

    sout << axis.x << " " << axis.y << " " << axis.z;

    return true;
}

bool SoCBirrtProblem::GetJointTransform(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    RAVELOG_DEBUG("Starting GetJointTransform...\n");

    KinBodyPtr pobject;
    int jointind = -1;

    dReal temp;
    string cmd, tempstring;

    while(!sinput.eof()) {
        sinput >> cmd;
        if( !sinput )
            break;

        if( stricmp(cmd.c_str(), "name") == 0 ) {
            sinput >> tempstring;
            pobject = GetEnv()->GetKinBody(tempstring.c_str());
            if(pobject == NULL)
            {
                RAVELOG_INFO("Error: could not find the specified object to get transform\n");
                return false;
            }
        }
        else if( stricmp(cmd.c_str(), "jointind") == 0 ) {
            sinput >> jointind;
        }
        else break;

        if( !sinput ) {
            RAVELOG_DEBUG("failed\n");
            return false;
        }
    }


    std::vector<KinBody::JointPtr> _vecjoints = pobject->GetJoints();

    if(jointind >= _vecjoints.size() || jointind < 0)
    {   
        RAVELOG_INFO("Joint index out of bounds or not specified\n");
        return false;
    }

    TransformMatrix jointtm = _vecjoints[jointind]->GetFirstAttached()->GetTransform();
    jointtm.trans = _vecjoints[jointind]->GetAnchor();
    sout << jointtm.m[0] << " ";
    sout << jointtm.m[4] << " ";
    sout << jointtm.m[8] << " ";
    sout << jointtm.m[1] << " ";
    sout << jointtm.m[5] << " ";
    sout << jointtm.m[9] << " ";
    sout << jointtm.m[2] << " ";
    sout << jointtm.m[6] << " ";
    sout << jointtm.m[10] << " ";
    sout << jointtm.trans.x << " ";
    sout << jointtm.trans.y << " ";
    sout << jointtm.trans.z << " ";

    return true;
}


bool SoCBirrtProblem::GrabBody(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    RAVELOG_DEBUG("Starting GrabBody...\n");

    KinBodyPtr ptarget;

    string cmd;
    while(!sinput.eof()) {
        sinput >> cmd;
        if( !sinput )
            break;
        
        if( stricmp(cmd.c_str(), "name") == 0 ) {
            string name;
            sinput >> name;
            ptarget = GetEnv()->GetKinBody(name.c_str());
        }
        else break;

        if( !sinput ) {
            RAVELOG_DEBUG("failed\n");
            return false;
        }
    }

    if(ptarget == NULL) {
        RAVELOG_INFO("ERROR Manipulation::GrabBody - Invalid body name.\n");
        return false;
    }
    
    robot->Grab(ptarget);
    return true;
}

int SoCBirrtProblem::convexHull2D(coordT* pointsIn, int numPointsIn, coordT** pointsOut, int* numPointsOut) {
 
    char flags[250];
    int exitcode;
    facetT *facet, *newFacet;
    int curlong, totlong;
    vertexT *vertex, *vertexA, *vertexB;
    int j;


    sprintf (flags, "qhull QJ Pp s Tc ");
    //FILE* junk = fopen("qhullout.txt","w");

    exitcode= qh_new_qhull (2, numPointsIn, pointsIn, false,
                                      flags, NULL, stderr);
    //fclose(junk);
    *numPointsOut = qh num_vertices;
    *pointsOut = (coordT *)malloc(sizeof(coordT)*(*numPointsOut)*2);

    FORALLfacets {
      facet->seen = 0;
    }

    FORALLvertices {
      vertex->seen = 0;
    }

    facet=qh facet_list;
    j=0;

    while(1) {
    if (facet==NULL) {
        // empty hull
        break;
    }
            vertexA = (vertexT*)facet->vertices->e[0].p;
            vertexB = (vertexT*)facet->vertices->e[1].p;
            if (vertexA->seen==0) {
                    vertexA->seen = 1;
                    (*pointsOut)[j++] = vertexA->point[0];
                    (*pointsOut)[j++] = vertexA->point[1];
            }
            if (vertexB->seen==0) {
                    vertexB->seen = 1;
                    (*pointsOut)[j++] = vertexB->point[0];
                    (*pointsOut)[j++] = vertexB->point[1];
            }


            //qh_printfacet(stderr, facet);
            facet->seen = 1;
            newFacet = (facetT*)facet->neighbors->e[0].p;
            if (newFacet->seen==1) newFacet = (facetT*)facet->neighbors->e[1].p;
            if (newFacet->seen==1) { break; }
            facet = newFacet;
    }

    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);

    return exitcode;
}


void SoCBirrtProblem::compute_convex_hull(void)
{  

    int dim =2;  	              /* dimension of points */
    int numpoints = 4;            /* number of points */
    coordT points[8] = {0, 0,  1, 0,  1, 1,  0, 1};           /* array of coordinates for each point */
    boolT ismalloc = false;           /* True if qhull should free points in qh_freeqhull() or reallocation */
    char flags[]= "qhull Tv"; /* option flags for qhull, see qh_opt.htm */
    FILE *outfile= stdout;    /* output from qh_produce_output()
                                 use NULL to skip qh_produce_output() */
    FILE *errfile= stderr;    /* error messages from qhull code */
    int exitcode;             /* 0 if no error from qhull */
    facetT *facet;	          /* set by FORALLfacets */
    int curlong, totlong;	  /* memory remaining after qh_memfreeshort */

    /* initialize dim, numpoints, points[], ismalloc here */
    exitcode= qh_new_qhull (dim, numpoints, points, ismalloc,
                                                    flags, outfile, errfile);
    if (!exitcode) { /* if no error */
            /* 'qh facet_list' contains the convex hull */
            FORALLfacets {
                    /* ... your code ... */
            }
    }
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
    if (curlong || totlong)
            fprintf (errfile, "qhull internal warning (main): did not free %d bytes of long memory (%d pieces)\n",
                         totlong, curlong);

    RAVELOG_INFO("qhull done\n");
}


void SoCBirrtProblem::GetSupportPolygon(std::vector<string>& supportlinks, std::vector<dReal>& polyx, std::vector<dReal>& polyy, Vector polyscale, Vector polytrans)
{
    int numsupportlinks = supportlinks.size();

    std::vector<boost::shared_ptr<KinBody::Link::Geometry> >_listGeomProperties;
    std::vector<Vector> points;
    AABB bounds;
    //get points on trimeshes of support links
    vector<KinBody::LinkPtr> vlinks = robot->GetLinks();
    for(int i = 0 ; i < numsupportlinks; i++)
    {   
        for(int j =0; j < vlinks.size(); j++)
        {
            if(strcmp(supportlinks[i].c_str(), vlinks[j]->GetName().c_str()) == 0 )
            {
                RAVELOG_DEBUG("Found match!\n");
                _listGeomProperties = vlinks[j]->GetGeometries();

                //compute AABBs for the link at identity
                for(int k = 0; k < _listGeomProperties.size(); k++)
                {
                //if( _listGeomProperties.size() == 1){
                    Transform _t = _listGeomProperties[k]->GetTransform().inverse();
                    //bounds = _listGeomProperties[k]->ComputeAABB(_t);
                    bounds = _listGeomProperties[k]->ComputeAABB(_listGeomProperties[k]->GetTransform());
                    Transform offset = vlinks[j]->GetTransform()*_t;
                    points.push_back(offset*Vector(bounds.pos.x + bounds.extents.x,bounds.pos.y + bounds.extents.y,bounds.pos.z - bounds.extents.z));
                    points.push_back(offset*Vector(bounds.pos.x - bounds.extents.x,bounds.pos.y + bounds.extents.y,bounds.pos.z - bounds.extents.z));
                    points.push_back(offset*Vector(bounds.pos.x - bounds.extents.x,bounds.pos.y - bounds.extents.y,bounds.pos.z - bounds.extents.z));
                    points.push_back(offset*Vector(bounds.pos.x + bounds.extents.x,bounds.pos.y - bounds.extents.y,bounds.pos.z - bounds.extents.z));
                    
                    points.push_back(offset*Vector(bounds.pos.x + bounds.extents.x,bounds.pos.y + bounds.extents.y,bounds.pos.z + bounds.extents.z));
                    points.push_back(offset*Vector(bounds.pos.x - bounds.extents.x,bounds.pos.y + bounds.extents.y,bounds.pos.z + bounds.extents.z));
                    points.push_back(offset*Vector(bounds.pos.x - bounds.extents.x,bounds.pos.y - bounds.extents.y,bounds.pos.z + bounds.extents.z));
                    points.push_back(offset*Vector(bounds.pos.x + bounds.extents.x,bounds.pos.y - bounds.extents.y,bounds.pos.z + bounds.extents.z));
                }
                break;
            }


        }
    }
    RAVELOG_INFO("Num support points in to qhull: %d\n",points.size());
    std::vector<coordT> pointsin(points.size()*2);
    std::vector<RaveVector<float> > plotvecs(points.size());
    for(int i = 0; i < points.size();i++)
    {
        pointsin[i*2 + 0] = points[i].x;
        pointsin[i*2 + 1] = points[i].y;
        plotvecs[i] = RaveVector<float>(points[i].x,points[i].y,points[i].z);
    }
    //GraphHandlePtr graphptr1 = GetEnv()->plot3(&plotvecs[0].x,plotvecs.size(),sizeof(plotvecs[0]),5, RaveVector<float>(1,0, 0, 1));
    //graphptrs.push_back(graphptr1);

    coordT* pointsOut = NULL;

    int numPointsOut = 0;

    convexHull2D(&pointsin[0], points.size(), &pointsOut, &numPointsOut);
    RAVELOG_INFO("Num support points out of qhull:: %d\n",numPointsOut);
    
    std::vector<RaveVector<float> > tempvecs(numPointsOut +1);
    polyx.resize(numPointsOut);
    polyy.resize(numPointsOut);

    Point2D polygon[numPointsOut];
    dReal centerx = 0;
    dReal centery = 0;
    for(int i =0; i < numPointsOut; i++)
    {
        polyx[i] = pointsOut[(i)*2 + 0];       
        polyy[i] = pointsOut[(i)*2 + 1];
        polygon[i].x = polyx[i];
        polygon[i].y = polyy[i];
        //centerx += polyx[i];
        //centery += polyy[i];
    }

    size_t vertexCount = sizeof(polygon) / sizeof(polygon[0]);
    Point2D centroid = compute2DPolygonCentroid(polygon, vertexCount);
    //std::cout << "Centroid is (" << centroid.x << ", " << centroid.y << ")\n";


    centerx = centroid.x;//centerx/numPointsOut;
    centery = centroid.y;//centery/numPointsOut;
    //RAVELOG_INFO("center %f %f\n",centerx, centery);

    for(int i =0; i < numPointsOut; i++)
    {
        polyx[i] = polyscale.x*(polyx[i] - centerx) + centerx + polytrans.x;       
        polyy[i] = polyscale.y*(polyy[i] - centery) + centery + polytrans.y;
        tempvecs[i] = RaveVector<float>(polyx[i],polyy[i],0);
    }

    

    //close the polygon
    tempvecs[tempvecs.size()-1] = RaveVector<float>(polyx[0],polyy[0],0);
    GraphHandlePtr graphptr = GetEnv()->drawlinestrip(&tempvecs[0].x,tempvecs.size(),sizeof(tempvecs[0]),5, RaveVector<float>(0, 1, 1, 1));
    graphptrs.push_back(graphptr);


    free(pointsOut);
}


bool SoCBirrtProblem::Traj(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    string filename; sinput >> filename;
    if( !sinput )
        return false;

    TrajectoryBasePtr ptraj = RaveCreateTrajectory(GetEnv(), "");
    //ptraj->Init(robot->GetConfigurationSpecification());

    ifstream infile(filename.c_str(),ios::in);
    infile.precision(16); 
    ptraj->deserialize(infile);

    RAVELOG_VERBOSE("executing traj\n");
    robot->GetController()->SetPath(ptraj);
    sout << "1";
    return true;

}




SoCBirrtProblem::Point2D SoCBirrtProblem::compute2DPolygonCentroid(const SoCBirrtProblem::Point2D* vertices, int vertexCount)
{
    Point2D centroid = {0, 0};
    double signedArea = 0.0;
    double x0 = 0.0; // Current vertex X
    double y0 = 0.0; // Current vertex Y
    double x1 = 0.0; // Next vertex X
    double y1 = 0.0; // Next vertex Y
    double a = 0.0;  // Partial signed area

    // For all vertices except last
    int i=0;
    for (i=0; i<vertexCount-1; ++i)
    {
        x0 = vertices[i].x;
        y0 = vertices[i].y;
        x1 = vertices[i+1].x;
        y1 = vertices[i+1].y;
        a = x0*y1 - x1*y0;
        signedArea += a;
        centroid.x += (x0 + x1)*a;
        centroid.y += (y0 + y1)*a;
    }

    // Do last vertex
    x0 = vertices[i].x;
    y0 = vertices[i].y;
    x1 = vertices[0].x;
    y1 = vertices[0].y;
    a = x0*y1 - x1*y0;
    signedArea += a;
    centroid.x += (x0 + x1)*a;
    centroid.y += (y0 + y1)*a;

    signedArea *= 0.5;
    centroid.x /= (6.0*signedArea);
    centroid.y /= (6.0*signedArea);

    return centroid;
}



void SoCBirrtProblem::GetDistanceFromLineSegment(dReal cx, dReal cy, dReal ax, dReal ay,  dReal bx, dReal by, dReal& distanceSegment, dReal& xout, dReal& yout)
{

    dReal distanceLine;

    dReal r_numerator = (cx-ax)*(bx-ax) + (cy-ay)*(by-ay);
    dReal r_denomenator = (bx-ax)*(bx-ax) + (by-ay)*(by-ay);
    dReal r = r_numerator / r_denomenator;

    dReal px = ax + r*(bx-ax);
    dReal py = ay + r*(by-ay);
     
    dReal s =  ((ay-cy)*(bx-ax)-(ax-cx)*(by-ay) ) / r_denomenator;

	distanceLine = fabs(s)*sqrt(r_denomenator);

    dReal xx = px;
    dReal yy = py;

    if ( (r >= 0) && (r <= 1) )
    {
            distanceSegment = distanceLine;
    }
    else
    {

            dReal dist1 = (cx-ax)*(cx-ax) + (cy-ay)*(cy-ay);
            dReal dist2 = (cx-bx)*(cx-bx) + (cy-by)*(cy-by);
            if (dist1 < dist2)
            {
                    xx = ax;
                    yy = ay;
                    distanceSegment = sqrt(dist1);
            }
            else
            {
                    xx = bx;
                    yy = by;
                    distanceSegment = sqrt(dist2);
            }


    }
    xout = xx;
    yout = yy;
}

bool SoCBirrtProblem::GetCamView(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    Transform camtm = GetEnv()->GetViewer()->GetCameraTransform();
    sout << camtm.rot.x << " " << camtm.rot.y << " " << camtm.rot.z << " " << camtm.rot.w << " " << camtm.trans.x << " " << camtm.trans.y << " " << camtm.trans.z;
    return true;
}

bool SoCBirrtProblem::SetCamView(ostream& sout, istream& sinput)
{
    //lock mutex
    EnvironmentMutex::scoped_lock lockenv(GetEnv()->GetMutex());

    Transform camtm;
    sinput >> camtm.rot.x;
    sinput >> camtm.rot.y;
    sinput >> camtm.rot.z;
    sinput >> camtm.rot.w;
    sinput >> camtm.trans.x;
    sinput >> camtm.trans.y;
    sinput >> camtm.trans.z;
    
    GetEnv()->GetViewer()->SetCamera(camtm);
    return true;
}

dReal SoCBirrtProblem::EstimatePlanningTime(dReal WSDist, dReal closestDist, std::vector<dReal>& q) {
	dReal planningTime = _planning_theta[0];
	planningTime += _planning_theta[1] * WSDist;
	planningTime += _planning_theta[2] * SoCBirrtPlanner::treenodes->size();
	planningTime += _planning_theta[3] * closestDist;

	return planningTime;
}

dReal SoCBirrtProblem::EstimateExecutionTime(dReal WSDist, dReal closestDist, std::vector<dReal>& q) {
	dReal planningTime = _execution_theta[0];
	planningTime += _execution_theta[1] * WSDist;
	planningTime += _execution_theta[2] * SoCBirrtPlanner::treenodes->size();
	planningTime += _execution_theta[3] * closestDist;
	return planningTime;
}

struct NodeDist {
	int nodeID;
	dReal dist;
};

struct comparator {
	bool operator() (NodeDist& n1, NodeDist& n2) {
		if(n1.dist < n2.dist)
			return true;
		return false;
	}
};


dReal EuclideanDistance(const RaveVector<dReal>& a, const RaveVector<dReal>& b) {
	return sqrt((a-b).dot(a-b));	// for translation, as w is always 1, then, it will not affact the distance
}

void SoCBirrtProblem::GenerateNodesWorkSpaceData(const RaveVector<dReal>& goal, dReal palmSize, std::vector<int>& topNNodeIds, std::vector<dReal>& topNNodeDists) {
	std::vector<dReal> curConfig;
	int topN = topNNodeIds.size();
	robot->GetActiveDOFValues(curConfig);
	std::vector<RrtNode>::iterator iter = SoCBirrtPlanner::treenodes->begin();
	std::priority_queue<NodeDist, std::vector<NodeDist>, comparator> maxHeap;
	for (int i=0; iter != SoCBirrtPlanner::treenodes->end(); iter++, i++) {
		robot->SetActiveDOFValues((*iter->GetDataVector()));
		iter->wsTransform = robot->GetActiveManipulator()->GetEndEffectorTransform();
		NodeDist nd;
		nd.nodeID = i;
		nd.dist = EuclideanDistance(iter->wsTransform.trans, goal);
		if(maxHeap.size() < topN)
			maxHeap.push(nd);
		else if(maxHeap.top().dist > nd.dist) {
			maxHeap.push(nd);
			maxHeap.pop();	// remove current top and reorganize the heap
		}
	}
	if(maxHeap.size() == topN) {
		for(int i=0; i<topN; i++) {
			topNNodeIds[topN-i-1] = maxHeap.top().nodeID;
			topNNodeDists[topN-i-1] = maxHeap.top().dist;
			maxHeap.pop();
		}
	}
	robot->SetActiveDOFValues(curConfig);	// reset to init state
}

void SoCBirrtProblem::InitializeRotationMatrixes() {
	//rotMats.clear();
	rotMats[0][0][0] = cos(PI/3);	rotMats[0][0][1] = 0; rotMats[0][0][2] = -sin(PI/3);
	rotMats[0][1][0] = 0;			rotMats[0][1][1] = 1; rotMats[0][1][2] = 0;
	rotMats[0][2][0] = sin(PI/3);	rotMats[0][2][1] = 0; rotMats[0][2][2] = cos(PI/3);

	rotMats[1][0][0] = cos(-PI/3);	rotMats[1][0][1] = 0; rotMats[1][0][2] = -sin(-PI/3);
	rotMats[1][1][0] = 0;			rotMats[1][1][1] = 1; rotMats[1][1][2] = 0;
	rotMats[1][2][0] = sin(-PI/3);	rotMats[1][2][1] = 0; rotMats[1][2][2] = cos(-PI/3);

	rotMats[2][0][0] = cos(PI/3);	rotMats[2][0][1] = sin(PI/3); rotMats[2][0][2] = 0;
	rotMats[2][1][0] = -sin(PI/3);	rotMats[2][1][1] = cos(PI/3); rotMats[2][1][2] = 0;
	rotMats[2][2][0] = 0;			rotMats[2][2][1] = 0; rotMats[2][2][2] = 1;

	rotMats[3][0][0] = cos(-PI/3);	rotMats[3][0][1] = sin(-PI/3); rotMats[3][0][2] = 0;
	rotMats[3][1][0] = -sin(-PI/3);	rotMats[3][1][1] = cos(-PI/3); rotMats[3][1][2] = 0;
	rotMats[3][2][0] = 0;			rotMats[3][2][1] = 0; rotMats[3][2][2] = 1;
}

RaveVector<dReal> SoCBirrtProblem::MultiplyMatrix(int matID, RaveVector<dReal> originVec) {
	RaveVector<dReal> vec;
	vec.x = rotMats[matID][0][0]*originVec.x + rotMats[matID][0][1]*originVec.y + rotMats[matID][0][2]*originVec.z;
	vec.y = rotMats[matID][1][0]*originVec.x + rotMats[matID][1][1]*originVec.y + rotMats[matID][1][2]*originVec.z;
	vec.z = rotMats[matID][2][0]*originVec.x + rotMats[matID][2][1]*originVec.y + rotMats[matID][2][2]*originVec.z;
	return vec;
}
