/*
 * SensorConfiguration.h
 *
 *  Created on: Apr 7, 2015
 *      Author: jason
 */

#ifndef SENSORCONFIGURATION_H_
#define SENSORCONFIGURATION_H_

#include "stdafx.h"

class SensorConfiguration {
public:
	SensorConfiguration(EnvironmentBasePtr penv, double deltaTime, double velMag, string goalName);
	virtual ~SensorConfiguration() {};

private:
	double m_prevTime;
	double m_deltaTime;
	double m_velMag;	// velocity magnitude
	EnvironmentBasePtr m_env;
	RaveVector<dReal> m_goal;
	KinBodyPtr m_goalObject;
public:
	void SetInitGoal(RaveVector<dReal> goal) {
		m_goal = goal;
	}

	dReal GetVelocityMagitude() {
		return m_velMag;
	}

	RaveVector<dReal> ReadSensorData(double curTime);

	void SetGoalTranform(Transform Tg);
};

#endif /* SENSORCONFIGURATION_H_ */
