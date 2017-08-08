#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

//******************************************************************
// 
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
 	 
	//cout << "In Tools:CalculateRMSE\n\n";

	VectorXd rmse(4);
	rmse << 0,0,0,0;
	
    // check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0)
	{
		cout << "ground_truth vector size is 0 or not equal! \n";
		return rmse;
	}

	//cout << "Tools CalculateRMSE : estimations.size() = " << estimations.size() << "\n\n";
	
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i)
	{
        VectorXd diff = estimations[i] - ground_truth[i];
		//cout << "diff1 = " << diff << "\n\n\n";

		// coefficient-wise multiplication
		diff = diff.array() * diff.array();
		//cout << "diff2 = " << diff << "\n\n\n";

		rmse += diff;
	}

	//calculate the mean
	rmse = rmse / estimations.size();
	
	//calculate the squared root
	rmse = rmse.array().sqrt();
	
	cout << "RMSE = " << rmse(0) << "\t" << rmse(1) << "\t" << rmse(2) << "\t" << rmse(3) << "\n\n";
	cout << "----------------------------------------------\n\n";
	
	//return the result
	return rmse;


}
