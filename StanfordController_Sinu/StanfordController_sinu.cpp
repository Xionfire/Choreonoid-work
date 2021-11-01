/**
   Simple controller for the Stanford robot
   @author Rafael Cisneros
 */

#include <cnoid/SimpleController>
#include <cnoid/SimulatorItem>
#include <cnoid/Item>
#include <vector>
#include <iostream>
#include <fstream>


using namespace cnoid;

const double pgain[] = {
  200000, 200000, 200000, 500, 500, 500,
  5000, 5000 };

const double dgain[] = {
  500, 5000, 5000, 5, 5, 5,
  150, 150 };
 
int i=1;


class StanfordController_sinu : public SimpleController
{
  Body* ioBody;
  double dt;
  std::vector<double> qref;
  std::vector<double> qold;
  std::vector<double> sin_param;
  int id;
  double Am,F;
  double qi,q1,q2,q3,qd;
  Link* C_joint;

  //Matrix for compute Qd
  Eigen::MatrixXf M(8,8);
  Eigen::MatrixXf minv(8,8);
  Eigen::VectorXf Q(8);
  Eigen::VectorXf a(8);

  //time var
  double t0,t1,t2,t3,T,currentime;


  std::vector<double> Qdout,Tout;

public:
	
   void writetxt(std::string fname,std::vector<double> varpush){
   	std::ofstream sin("/home/xionfire/Bureau/Japan_Project/Controller/StanfordController_Sinu/"+ fname, std::ios::app);  
   	if(sin.is_open()){
        	for (int i=0;i<varpush.size();i++){
        		sin << std::to_string(varpush[i]) << std::endl;
        	}	
        }
    sin.close();
   }	
   

		
  virtual bool initialize(SimpleControllerIO* io) override
  {
  
    ioBody = io->body();
    t0 = io->currentTime();
    dt = io->timeStep();   
    
    //reading of our sinusoidal parameter 
    std::ifstream file ("/home/xionfire/Bureau/Japan_Project/Controller/StanfordController_Sinu/sin_param.txt");
    if (file.is_open()) {
    	std::string line;
    	while (std::getline(file, line)) {
    	    if((i%2)==0) {
            	sin_param.push_back(std::stod(line.c_str()));
            }    
            i=i+1;  
    	}
    	file.close();
    } 
    //-----------------------------------------------------------------------------------------------------------
    //Compute var 
    id = sin_param[0];
    Am = sin_param[1];
    F = sin_param[2];
    T=1/F;
    t1=(T/3);
    t2=(2*T/3);
    t3=T;
    C_joint=ioBody->joint(id);
    qi=C_joint->q();
    q1=(Am/2)+qi;
    q2=qi-(Am/2);
    q3=qi;
    std::cout << "qi = " <<qi << std::endl;
    std::cout << "q1 = " <<q1 << std::endl;
    std::cout << "q2 = " <<q2 << std::endl;
    
    //Compute a
    M<<1,t0,pow(t0,2),pow(t0,3),pow(t0,4),pow(t0,5),pow(t0,6),pow(t0,7),
       0,1,2*t0,3*pow(t0,2),4*pow(t0,3),5*pow(t0,4),6*pow(t0,5),7*pow(t0,6),
       0,0,2,6*t0,12*pow(t0,2),20*pow(t0,3),30*pow(t0,4),42*pow(t0,5),
       1,t1,pow(t1,2),pow(t1,3),pow(t1,4),pow(t1,5),pow(t1,6),pow(t1,7),
       1,t2,pow(t2,2),pow(t2,3),pow(t2,4),pow(t2,5),pow(t2,6),pow(t2,7),
       1,t3,pow(t3,2),pow(t3,3),pow(t3,4),pow(t3,5),pow(t3,6),pow(t3,7),
       0,1,2*t3,3*pow(t3,2),4*pow(t3,3),5*pow(t3,4),6*pow(t3,5),7*pow(t3,6),
       0,0,2,6*t3,12*pow(t3,2),20*pow(t3,3),30*pow(t3,4),42*pow(t3,5);
    Q<<qi,
       0,
       0,
       q1,
       q2,
       q3,
       0,
       0;
     
    minv = M.completeOrthogonalDecomposition().pseudoInverse();
    a=minv*Q;
    
    //Get the ref pos 
    for (int i = 0; i < ioBody->numJoints(); ++i) {
      Link* joint = ioBody->joint(i);
      joint->setActuationMode(Link::JOINT_TORQUE);
      io->enableIO(joint);
      qref.push_back(joint->q());
    }
    qold=qref;
    return true; 
  }
  
  
  
  
  
  virtual bool control() override
  { 
    //Compute qd
    currentime = ioBody->currentTime();
    qd=a(1,0)+a(1,0)*currentime+a(2,0)*pow(currentime,2)+a(3,0)*pow(currentime,3)+a(4,0)*pow(currentime,4)+a(5,0)*pow(currentime,5)+a(6,0)*pow(currentime,6)+a(7,0)*pow(currentime,7);
    qref[id]=qd;
    Qdout.push_back(qd);
    Tout.push_back(currentime);
   
    
    for (int i = 0; i < ioBody->numJoints(); ++i) {
      Link* joint = ioBody->joint(i);
      double q = joint->q();
      double dq = (q - qold[i]) / dt;
      double u = (qref[i] - q) * pgain[i] + (0.0 - dq) * dgain[i];
      qold[i] = q;
      joint->u() = u;
    }
    return true;
  }
  
  
  virtual void stop (){
   //push the time and the result of qd on a .txt for print our sinusoidal angle trajectory    
    writetxt("sinus.txt",Qdout);
    writetxt("time.txt",Tout);
   // -------------------------------------------------------------------------------------------------------------------------
  }
  
  
};

CNOID_IMPLEMENT_SIMPLE_CONTROLLER_FACTORY(StanfordController_sinu)
