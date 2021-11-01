/**
   Simple controller for the Stanford robot
   @author Rafael Cisneros
 */

#include <cnoid/SimpleController>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace cnoid;

const double pgain[] = {
  200000, 200000, 200000, 500, 500, 500,
  5000, 5000 };

const double dgain[] = {
  500, 5000, 5000, 5, 5, 5,
  150, 150 };

double d[] ={0.2,0,0.089};
int i=0;
double tf=5.0;



class StanfordController_3D : public SimpleController
{
  Body* ioBody;
  double dt,ct;
  std::vector<double> qref,qold,qac;
  Eigen::Vector3d z0,z1,z2,z3,z4,z5;
  Eigen::Matrix<double,1,3>o6,o3;
  Eigen::Matrix<double, 6, 1>j1,j2,j3,j4,j5,j6;
  Eigen::Matrix<double, 3, 1>Xi,Xref;
  Eigen::Transform<double,3,Eigen::Affine> A1,A2,A3,A4,A5,A6;
  Eigen::Matrix3d K;

public:


  void HomogeneousM(std::vector<double> qac)
  {
    A1.linear()<<cos(qac[0]),0,-sin(qac[0]),
    		 sin(qac[0]),0,cos(qac[0]),
    		 0,-1,0; 
    		 
    A2.linear()<<cos(qac[1]),0,sin(qac[1]),
    		 sin(qac[1]),0,-cos(qac[1]),
    		 0,1,0;
    A2.translation()<<0,
                      0,
                      d[0];
                      	 
    A3.linear()<<1,0,0,
                 0,1,0,
                 0,0,1;  
    A3.translation()<<0,
                      0,
                      d[1];
                      
    A4.linear()<<cos(qac[3]),0,-sin(qac[3]),
    		 sin(qac[3]),0,cos(qac[3]),
    		 0,-1,0; 
    		 
    A5.linear()<<cos(qac[4]),0,sin(qac[4]),
    		 sin(qac[4]),0,-cos(qac[4]),
    		 0,-1,0;
    A6.linear()<<cos(qac[5]),-sin(qac[5]),0,
    		 sin(qac[5]),cos(qac[5]),0,
    		 0,0,1;
    A6.translation()<<0,
                      0,
                      d[2];
  } 	
  Eigen::Matrix<double, 6, 6> Jacobien(std::vector<double> qac){
	o3 <<cos(qac[0])*sin(qac[1])*d[1]-sin(qac[0])*d[0],
	     sin(qac[0])*sin(qac[1])*d[1]+cos(qac[0])*d[0],
	     cos(qac[1])*d[1];
  	o6 <<cos(qac[0])*sin(qac[1])*d[1]-sin(qac[0])*d[0]+d[2]*(cos(qac[1])*cos(qac[0])*cos(qac[3])*sin(qac[4])+cos(qac[4])*cos(qac[0])*sin(qac[1])-sin(qac[0])*sin(qac[3]*sin(qac[4]))),
  	     sin(qac[0])*sin(qac[1])*d[1]-cos(qac[0])*d[0]+d[2]*(cos(qac[0])*sin(qac[3]*sin(qac[4])+cos(qac[1])*sin(qac[0])*cos(qac[3])*sin(qac[4])+cos(qac[4])*sin(qac[0])*sin(qac[1]))),
  	     cos(qac[1])*d[1]+d[2]*(cos(qac[1])*cos(qac[4]-cos(qac[3])*sin(qac[1])*sin(qac[4])));
  	z0<<0,
  	    0,
  	    1;
  	z1<<-sin(qac[0]),
  	    cos(qac[0]),
  	    0;
  	z2<<cos(qac[0])*sin(qac[1]),
  	    sin(qac[0])*sin(qac[1]),
  	    cos(qac[1]);
  	z3=z2;
  	z4<<-cos(qac[0])*cos(qac[1])*sin(qac[3])-sin(qac[0])*cos(qac[3]),   
  	    -sin(qac[0])*cos(qac[1])*sin(qac[3])+cos(qac[0])*cos(qac[3]),
  	    sin(qac[1])*sin(qac[3]);
  	    
  	z5<<cos(qac[0])*cos(qac[1])*cos(qac[3])*sin(qac[4])-sin(qac[0])*sin(qac[4])*sin(qac[3])+sin(qac[1])*cos(qac[0])*cos(qac[4]),
  	    sin(qac[0])*cos(qac[1])*cos(qac[3])*sin(qac[4])+sin(qac[3])*cos(qac[0])*sin(qac[4])+sin(qac[0])*cos(qac[4])*sin(qac[1]),
  	    -sin(qac[1])*cos(qac[3])*sin(qac[4])+cos(qac[1])*cos(qac[4]);    
  	    
  	j1<<z0.cross(o6),z0;
  	j2<<z1.cross(o6),z1;
  	j3<<z2,0,0,0;
  	j4<<z3.cross(o6-o3),z3;
  	j5<<z4.cross(o6-o3),z4;
  	j6<<z5.cross(o6-o3),z5;
  	Eigen::Matrix<double, 6, 6> J;
  	J<<j1,j2,j3,j4,j5,j6;
  	
  	
  return J;
  }
  
   Eigen::Matrix<double, 3, 1> Trajectory( Eigen::Matrix<double, 3, 1>Xref,Eigen::Matrix<double, 3, 1>Xi,double tf,double ct){
   	double r=10*pow(ct/tf,3)-15*pow(ct/tf,4)+6*pow(ct/tf,5);
   	double r_point=30.0*(pow(ct,2)/pow(tf,3)) - 60.0*(pow(ct,3)/pow(tf,4)) + 30.0*(pow(ct,4)/pow(tf,5));
  	Eigen::Matrix<double, 3, 1>D=Xref-Xi;
   	Eigen::Matrix<double, 3, 1> Xd=Xi+r*D;
   	Eigen::Matrix<double, 3, 1> X_pointd=r_point*D;
   return Xd,X_pointd;
   }


  virtual bool initialize(SimpleControllerIO* io) override
  {
    K << 10,0,0,
         0,10,0,
         0,0,10;
    ioBody = io->body();
    dt = io->timeStep();
    for (int i = 0; i < ioBody->numJoints(); ++i) {
      Link* joint = ioBody->joint(i);
      joint->setActuationMode(Link::JOINT_TORQUE);
      io->enableIO(joint);
      qold.push_back(joint->q());
    }
    qref=qold;

    
    // Get Xref from the file
    std::ifstream file ("/home/xionfire/Bureau/Japan_Project/Controller/StanfordController_3D/pos.txt");
    if (file.is_open()) {
    	std::string line;
    	while (std::getline(file, line)) {
    	    if(i!=0) {
            	Xref(i-1)=std::stod(line.c_str()); 
            }  
            i=i+1;;  
    	}
    	file.close();
    }
    //----------------------------------------------------
    for (int i = 0; i < ioBody->numJoints(); ++i) {
     Link* joint = ioBody->joint(i);
     if (i==2){
      	qac.push_back(joint->q());
      	d[1]=joint->q();
      	}
     else{qac.push_back((180*joint->q())/M_PI);}
    } 
    //compute the initial position
    HomogeneousM(qac);
    Eigen::Transform<double,3,Eigen::Affine> T=A1*A2*A3*A4*A5*A6;
    Xi=T.translation();	
    return true; 
  }






  virtual bool control() override
  {
  
   ct = ioBody->currentTime();
   for (int i = 0; i < ioBody->numJoints(); ++i) {
      Link* joint = ioBody->joint(i);
      if (i==2){
      		qac.push_back(joint->q());
      		d[1]=joint->q();
      		}
      else{qac.push_back((180*joint->q())/M_PI);}
    } 
    
    
    HomogeneousM(qac);
    Eigen::Matrix<double, 6, 6> J=Jacobien(qac);
    Eigen::MatrixXd Jinv=J.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::Transform<double,3,Eigen::Affine> T=A1*A2*A3*A4*A5*A6;
    Eigen::Matrix<double, 3, 1>X=T.translation();	
    Eigen::Matrix<double, 3, 1> Xd,X_pointd=Trajectory(Xref,Xi,tf,ct);
    Eigen::Matrix<double, 3, 1> Ep=K*(Xd-X)+X_pointd;
    Eigen::Matrix<double, 3, 6> Jinv_resize=Jinv.block(0, 0, 3, 6);
    Eigen::Matrix<double, 6, 1> Qpoint=Jinv_resize.transpose()*Ep;
    
    for (int i = 0; i < ioBody->numJoints()-2; ++i) {
     	 double qd = Qpoint[i]*dt+qold[i];
         qref[i]=qd;
      }
      
      
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
    HomogeneousM(qac);
    Eigen::Transform<double,3,Eigen::Affine> T=A1*A2*A3*A4*A5*A6;
    Eigen::Matrix<double, 3, 1>X=T.translation();
    std::cout<< X << std::endl;
 }
 
};

CNOID_IMPLEMENT_SIMPLE_CONTROLLER_FACTORY(StanfordController_3D)
