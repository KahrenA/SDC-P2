#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//************************************************************************
// UKF
//
// Initializes Unscented Kalman filter
//************************************************************************
UKF::UKF() 
{

	// Init 
	is_initialized_ = false; 

	// if this is false, laser measurements will be ignored,except during init
	use_laser_ = false;     // true;

	// if this is false, radar measurements will be ignored,except during init
  	use_radar_ = true;

  	// initial state vector
  	x_ = VectorXd(5);

  	// initial covariance matrix
  	P_ = MatrixXd(5, 5);

	P_ <<  	1.0, 0.0, 0.0, 0.0, 0.0,
		  	0.0, 1.0, 0.0, 0.0, 0.0, 
			0.0, 0.0, 1.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 0.0, 1.0;


  	// Process noise standard deviation longitudinal acceleration in m/s^2
  	std_a_ = 0.3; 		// 30;

  	// Process noise standard deviation yaw acceleration in rad/s^2
  	std_yawdd_ = 0.3;   // 30;

  	// Laser measurement noise standard deviation position1 in m
  	std_laspx_ = 0.15;

  	// Laser measurement noise standard deviation position2 in m
  	std_laspy_ = 0.15;

  	// Radar measurement noise standard deviation radius in m
  	std_radr_ = 0.3;

  	// Radar measurement noise standard deviation angle in rad
  	std_radphi_ = 0.03;

  	// Radar measurement noise standard deviation radius change in m/s
  	std_radrd_ = 0.3;

  	// TODO:
  	// Complete the initialization. See ukf.h for other member properties.
	n_x_ = 5;	
	n_aug_ = 7; 
	n_z_= 3;
 	lambda_ = 3 - n_aug_;

	x_aug = VectorXd(n_aug_);

  	//create augmented covariance matrix
	Q = MatrixXd(2, 2);
	Q << (std_a_ * std_a_), 	0, 
			0, 				(std_yawdd_ * std_yawdd_);

	P_aug = MatrixXd(n_aug_, n_aug_);

	
	// -------Sigma points matrix------------
	Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
	Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

	// Measurement 
	z = VectorXd(n_z_);
	Zsig = MatrixXd(n_z_, 2*n_aug_+1);

//	H_laser_ = MatrixXd(2,4);
//	H_laser_ << 1, 0, 0, 0,
//	            0, 1, 0, 0;

	// -------------set weights---------------
	weights_ = VectorXd(2*n_aug_+1);
  	double weight_0 = lambda_/(lambda_+n_aug_);	
  	weights_(0) = weight_0;

  	for (int i=1; i < 2*n_aug_+1; i++) 
	{  
    	double weight = 0.5/(n_aug_+ lambda_);
    	weights_(i) = weight;
  	}
	//cout << "created weights_(i) = " << weights_ << "\n\n";

  	// Hint: one or more values initialized above might be wildly off...
	
}
//************************************************************************
// ~UKF
//
//************************************************************************
UKF::~UKF() {}



//************************************************************************
// ProcessMeasurement
//
// @param {MeasurementPackage} meas_package The latest measurement data of
// either radar or laser.
//************************************************************************
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  //  TODO:
  // Complete this function! Make sure you switch between lidar and radar
  // measurements.
 
	//-------------------
	//  Initialization
	//-------------------
  	if (!is_initialized_) 
  	{
  
		// first measurement
	    cout << "UKF: Init" << endl;
	    x_ << 1, 1, 1, 1, 1;

	    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) 
	    {
			//	Convert radar from polar to cartesian coordinates and initialize state
			// px = py = phi; vx = vy = 0 per forum notes 
			x_ << meas_package.raw_measurements_[1], 
						meas_package.raw_measurements_[1], 0, 0, 0;  

			cout << "x_ raw measurements = " << x_ << "\n\n";	
	
			time_us_ = meas_package.timestamp_;
 		
    	}
    	else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
 		{
 			// Initialize state.
 			// set the state with the initial location and zero velocity
			x_ << meas_package.raw_measurements_[0], 
						meas_package.raw_measurements_[1], 0, 0, 0;

			cout << "Read laser data\n"; 

			time_us_ = meas_package.timestamp_;
 		}

		cout << "read x_ successfully \n\n";

		//=============================
		GenerateAugmentedSigmaPoints();


		// done initializing, no need to predict or update
    	is_initialized_ = true;
		cout << "Init Done \n\n";
    	return;
  	}

	

  	//==============================
  	// Prediction
  	//==============================
	//cout << "EKF Predict \n\n";
	
  	//**TODO:
  	// Time is measured in seconds.
  	//
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	
	time_us_ = meas_package.timestamp_;
		
	
	//---------- Predict ---------------
 	Prediction(dt);
	
	
  	//=================
  	// Update
  	//=================

  	// TODO:
  	// Use the sensor type to perform the update step.
  	// Update the state and covariance matrices.
  	

  	//-------------------------------------------------------------
  	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) 
	{
		//------- Radar updates ---------
		UpdateRadar(meas_package);

		// Generate Sigma Points based on last measured x & P
		GenerateAugmentedSigmaPoints();

  	} 
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true)
	{
		// --------Laser updates----------
 		UpdateLidar(meas_package);
	
		// Generate Sigma Points based on last measured x & P
		GenerateAugmentedSigmaPoints();
  	
	}
	else
	{
		cout << "Either Use flags are false or data corruption!!! \n\n";
	}
	
	return;
}

//************************************************************************
// Prediction
// 
//Predicts sigma points, the state, and the state covariance matrix.
// @param {double} delta_t the change in time (in seconds) between the last
//measurement and this one.
//************************************************************************
void UKF::Prediction(double delta_t) 
{

//	cout << "In Prediction \n\n";
//  cout << "delta_t = " << delta_t << "\n\n";

  	// TODO: Complete this function! 
	// Estimate the object's location. Modify the state vector x_
  	// Predict sigma points, the state, and the state covariance matrix.

	//----------------------
	// predict sigma points
	//----------------------
	Xsig_pred_.fill(0.0);
  	for (int i = 0; i< 2*n_aug_+1; i++)
	{
		//extract values for better readability
		double p_x      = Xsig_aug(0,i);
		double p_y      = Xsig_aug(1,i);
		double v        = Xsig_aug(2,i);
		double yaw      = Xsig_aug(3,i);
		double yawd     = Xsig_aug(4,i);
		double nu_a     = Xsig_aug(5,i);
		double nu_yawdd = Xsig_aug(6,i);

		//cout << "i = " << i << "\n";
		//cout << "p_x = " << p_x << "\t" << "p_y = " << p_y << "\t" << "v = " << v << "\t" << "yaw = " << yaw << "\t" << "yawd = " << yawd << 
		//															"\t" << "nu_a = " << nu_a << "\t" << "nu_yawdd = " << nu_yawdd << "\n";

		//predicted state values
		double px_p, py_p;

    	//avoid division by zero
    	if (fabs(yawd) > 0.001)
		{
        	px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        	py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd*delta_t) );
    	}
    	else 
		{
	        px_p = p_x + v*delta_t*cos(yaw);
        	py_p = p_y + v*delta_t*sin(yaw);
    	}
		//cout << "px_p = " << px_p << "\t" << "py_p = " << py_p << "\n";

    	double v_p = v;
    	double yaw_p = yaw + yawd*delta_t;
    	double yawd_p = yawd;

    	//add noise
    	px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    	py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    	v_p = v_p + nu_a * delta_t;
		//cout << "Noise added: px_p = " << px_p << "\t" << "py_p = " << py_p << "\t" << "v_p = " << v_p << "\n";

    	yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    	yawd_p = yawd_p + nu_yawdd * delta_t;
		//cout << "yaw_p = " << yaw_p << "\t" << "yawd_p = " << yawd_p << "\n\n";

    	//write predicted sigma point into right column
    	Xsig_pred_(0,i) = px_p;
    	Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
    	Xsig_pred_(4,i) = yawd_p;

	} // for
	//cout << "Done with Predicting Sigma Points\n";


	//-------------------------------
	//predicted state mean
  	//------------------------------
  	x_.fill(0.0);

	// iterate over sigma points
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{  
	    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  	}
	
	cout << "Predicted x_ = \n" << x_ << "\n\n";

	//---------------------------------
  	//predicted state covariance matrix
	//---------------------------------
  	// P_.fill(0.0);

	// iterate over sigma points
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{  

    	// state difference
    	VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
		//angle normalization
    	while (x_diff(3) > M_PI) 
				x_diff(3) -= 2.*M_PI;

    	while (x_diff(3) < -M_PI) 
				x_diff(3) += 2.*M_PI;

    	P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

 	} // for

//	cout << "Predicted P_ = " << "\n" << P_ << "\n\n";
  
}

//*****************************************************************************
// UpdateLidar
//
// Updates the state and the state covariance matrix using a laser measurement.
// @param {MeasurementPackage} meas_package
//*****************************************************************************
void UKF::UpdateLidar(MeasurementPackage meas_package) 
{
	// TODO:
  	// Complete this function! Use lidar data to update the belief about the object's
  	// position. Modify the state vector, x_, and covariance, P_.


}

//*****************************************************************************
// UpdateRadar
//
// Updates the state and the state covariance matrix using a radar measurement.
// @param {MeasurementPackage} meas_package
//*****************************************************************************
void UKF::UpdateRadar(MeasurementPackage meas_package) 
{
	//cout << "In UpdateRadar\n\n";

  // TODO:
  // Complete this function! Use radar data to update the belief about the object's
  // position. Modify the state vector, x_, and covariance, P_.
 
	//transform sigma points into measurement space
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{  
    	// extract values for better readibility
		double p_x = Xsig_pred_(0,i);
		double p_y = Xsig_pred_(1,i);
		double v   = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		if(fabs(Zsig(0,i)) < 0.001)
			Zsig(0,i) = 0.001;	

		Zsig(1,i) = atan2(p_y,p_x);                                 //phi
		if(fabs(Zsig(1,i)) < 0.001)
			Zsig(1,i) = 0.001;	

		Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
		if(fabs(Zsig(2,i)) < 0.001)
			Zsig(2,i) = 0.001;	
	}

	// mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	z_pred.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++) 
	{
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}
	cout << "UpdateRadar: z_pred = \n" << z_pred << "\n\n";
	
	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{  
    	VectorXd z_diff = Zsig.col(i) - z_pred;		//residual

    	//angle normalization
    	while (z_diff(1) > M_PI) 
			z_diff(1) -= 2.*M_PI;
    	
		while (z_diff(1) < -M_PI) 
			z_diff(1) += 2.*M_PI;

    	S = S + weights_(i) * z_diff * z_diff.transpose();
  	}

  	//add measurement noise covariance matrix
  	MatrixXd R = MatrixXd(n_z_, n_z_);
  	R <<    std_radr_ * std_radr_, 	0, 							0,
    	    0, 						std_radphi_*std_radphi_, 	0,
    	    0, 						0,							std_radrd_*std_radrd_;

	S = S + R;
	//cout << "updateRadar: S = \n" << S << "\n\n";


	//-----------------------
	//create matrix for cross correlation Tc
  	MatrixXd Tc = MatrixXd(n_x_, n_z_);

	//calculate cross correlation matrix
  	Tc.fill(0.0);
  	for (int i = 0; i < 2 * n_aug_ + 1; i++) 
	{ 
	    //residual
	    VectorXd z_diff = Zsig.col(i) - z_pred;

	    //angle normalization
	    while (z_diff(1) > M_PI) 
			z_diff(1) -= 2.*M_PI;
	    while (z_diff(1) < -M_PI) 
			z_diff(1) += 2.*M_PI;

    	// state difference
	    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
		//angle normalization
    	while (x_diff(3) > M_PI) 
			x_diff(3) -= 2.*M_PI;

	    while (x_diff(3) < -M_PI) 
			x_diff(3) += 2.*M_PI;

	    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

  	//Kalman gain K;
  	MatrixXd K = Tc * S.inverse();

  	//residual
	z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
	cout << "UpdateRadar: measured z = \n" << z << "\n\n";
  	
	VectorXd z_diff = z - z_pred;

  	//angle normalization
  	while (z_diff(1) > M_PI) 
			z_diff(1) -= 2.*M_PI;
  	while (z_diff(1) < -M_PI) 
		z_diff(1) += 2.*M_PI;

  	//update state mean and covariance matrix
  	x_ = x_ + K * z_diff;
	cout << "UpdateRadar: x_ = \n" << x_ << "\n\n";

  	P_ = P_ - K * S * K.transpose();
//	cout << "UpdateRadar: P_ = \n" << P_ << "\n\n"; 


 
  	// You'll also need to calculate the radar NIS.
	// E = NIS_radar_ = (zk+1 - zk+1|k)T * Sinv k+1|k * ( zk+1 - zk+1|k)
	NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z-z_pred);
	cout << "NIS_radar_ = " << NIS_radar_ << "\n\n";
	
	//if( NIS_radar_ > 7.815 ) ... adjust something 
 
}


//**************************************************************
//GenerateAugmentedSigmaPoints
// 
// Based on last P_ and measured x_, generated sigma points
//**************************************************************
void UKF::GenerateAugmentedSigmaPoints()
{

	cout << "In GenerateAugmentedSigmaPoints\n\n";

	// GENERATE AUGMENTED SIGMA POINTS 
	//-----------------------------------------------------------
	//create augmented mean state
	//cout << "n_x_ = " << n_x_ << "\n\n";

	for (int i=0; i < n_x_; i++)
	{
//		cout << "x_(i) = " << x_(i) << "\n";
	    x_aug(i) = x_(i);
	}

	x_aug(5) = x_aug(6) = 0;
//	std::cout << "\nx_aug = \n" << x_aug << "\n\n";

	//----------------------------------------------------------
	//create square root matrix
	for (int i = 0; i < n_x_; i++)
	{
    	for (int j = 0; j < n_x_; j++)
   		{
        	P_aug(i,j) = P_(i,j);
    	}
	}	
	for (int i = 0; i < 2; i++)
	{
    	for (int j = 0; j < 2; j++)
    	{
		    P_aug(n_x_+i, n_x_+j) = Q(i, j); 
    	}
	}
//	std::cout << "P_aug = " << "\n" << P_aug << "\n\n";

	MatrixXd A = P_aug.llt().matrixL();
//	std::cout << "A = " << A << "\n\n";

	// create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for( int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}

//	cout << "Xsig_aug = " << "\n" << Xsig_aug << "\n\n";

}
