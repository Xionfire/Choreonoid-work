/**
   Simple controller for the Stanford robot
   @author Rafael Cisneros
 */

#include <cnoid/SimpleController>
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
  

class StanfordController : public SimpleController
{
  Body* ioBody;
  double dt;
  std::vector<double> qref;
  std::vector<double> qold;

public:
  	

  virtual bool initialize(SimpleControllerIO* io) override
  {
    ioBody = io->body();
    dt = io->timeStep();
   //Get Qref from a file
    std::ifstream file ("/home/xionfire/Bureau/Japan_Project/Controller/StanfordController/pos.txt");
    if (file.is_open()) {
    	std::string line;
    	while (std::getline(file, line)) {
            qref.push_back(std::stod(line.c_str()));       
    	}
    	file.close();
    }


    for (int i = 0; i < ioBody->numJoints(); ++i) {
      Link* joint = ioBody->joint(i);
      joint->setActuationMode(Link::JOINT_TORQUE);
      if (i==2){joint->setActuationMode(Link::JOINT_ANGLE);
      }
      std::cout<<joint->actuationMode()<< std::endl;
      io->enableIO(joint);
      qold.push_back(joint->q());
    }
        
   /*other Qref tested
    qref={-0.00969635,0.921505,0.295909,0.00646624,0.0193816,0.0452513,-7.34932e-05,7.34932e-05}
          1.69445,-1.37405,0.294984,0.00646679,0.0192387,0.0452518,1.6052e-05,-1.6052e-05*/


    return true;
  }

  virtual bool control() override
  {
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
};

CNOID_IMPLEMENT_SIMPLE_CONTROLLER_FACTORY(StanfordController)
