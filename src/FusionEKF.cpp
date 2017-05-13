#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.01, 0,
        0, 0.01;
 //0225
  //measurement covariance matrix - radar
  R_radar_ << 0.009, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.009;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

	  float rho;
	  float phi;
	  float rhodot;
	  float px;
	  float py;
	  float vx;
	  float vy;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    //state covariance matrix P
   ekf_.P_ = MatrixXd(4, 4);
   ekf_.P_ <<  1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 1000, 0,
			   0, 0, 0, 1000;

   //the initial transition matrix F_
   ekf_.F_ = MatrixXd(4, 4);

   // inital Q matrix gets overwritten below
   ekf_.Q_ = MatrixXd(4, 4);


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        rho = measurement_pack.raw_measurements_(0);
        phi = measurement_pack.raw_measurements_(1);
        rhodot = measurement_pack.raw_measurements_(2);

        px = rho * cos(phi);
        py = rho * sin(phi);
        vx = rhodot * cos(phi);
        vy = rhodot * sin(phi);

        ekf_.x_ << px,py,vx,vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
     ekf_.x_ << measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1),0.0,0.0;

    }
    // time starts from zero seconds ant not from Jan 1,1970
    previous_timestamp_ = measurement_pack.timestamp_;


    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt =  (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    float dt_4_4 = dt_4 / 4;
    float dt_3_2 = dt_3 / 2;
    float noise_ax = 4;
    float noise_ay = 4;


    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ <<  1, 0, dt, 0,
                0, 1, 0, dt,
                0, 0, 1,  0,
                0, 0, 0,  1;

    ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
  	         0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
  	         dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
   	         0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;

    ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	    // Radar updates
	    Tools tools;
		VectorXd z = VectorXd(3);
		z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],
				measurement_pack.raw_measurements_[2];
	    Hj_ = tools.CalculateJacobian(ekf_.x_);
	    ekf_.H_ = Hj_;
	    ekf_.R_ = R_radar_;
	    // if we have an initialised Jacobian update
		if (!Hj_.isZero(0)){
		    ekf_.UpdateEKF(z);
		}
	  } else {
	    // Laser updates
	    VectorXd zL(2);
	    zL << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1];
	    ekf_.H_ = H_laser_;
	    ekf_.R_ = R_laser_;
	    ekf_.Update(zL);

	  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
