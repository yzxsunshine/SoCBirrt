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
	SensorConfiguration(EnvironmentBasePtr penv, double deltaTime) : m_prevTime(-1.0), m_deltaTime(deltaTime), m_env(penv) {

	}
	virtual ~SensorConfiguration();

private:
	double m_prevTime;
	double m_deltaTime;
	double m_velMag;	// velocity magnitude
	EnvironmentBasePtr m_env;
	RaveVector<dReal> m_goal;
public:
	void SetInitTSR(RaveVector<dReal> goal) {
		m_goal = goal;
	}
	RaveVector<dReal> ReadSensorData(double curTime) {
#ifdef REAL_SENSOR

#else
		if (m_prevTime > 0) {
			// first time to generate goal
			double passTime = curTime - m_prevTime;
			int moveCount = passTime / m_deltaTime;
			for (int i=0; i<moveCount; i++) {
				RaveVector<dReal> vel;
				vel.Set3(RANDOM_FLOAT(), RANDOM_FLOAT(), RANDOM_FLOAT());
				vel.normalize3();
				vel *= m_velMag;
				m_goal += vel;
			}
		}
		m_prevTime = curTime;
		return m_goal;
	}
#endif
};

#endif /* SENSORCONFIGURATION_H_ */
