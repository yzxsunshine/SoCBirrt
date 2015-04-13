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
	SensorConfiguration(EnvironmentBasePtr penv, double deltaTime, double velMag);
	virtual ~SensorConfiguration() {};

private:
	double m_prevTime;
	double m_deltaTime;
	double m_velMag;	// velocity magnitude
	EnvironmentBasePtr m_env;
	RaveVector<dReal> m_goal;
public:
	void SetInitGoal(RaveVector<dReal> goal) {
		m_goal = goal;
	}

	dReal GetVelocityMagitude() {
		return m_velMag;
	}

	RaveVector<dReal> ReadSensorData(double curTime);
};

#endif /* SENSORCONFIGURATION_H_ */
