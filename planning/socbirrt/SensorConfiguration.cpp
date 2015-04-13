/*
 * SensorConfiguration.cpp
 *
 *  Created on: Apr 13, 2015
 *      Author: jason
 */
#include "SensorConfiguration.h"

SensorConfiguration::SensorConfiguration(EnvironmentBasePtr penv, double deltaTime, double velMag) : m_prevTime(-1.0), m_deltaTime(deltaTime), m_velMag(velMag), m_env(penv) {

}

RaveVector<dReal> SensorConfiguration::ReadSensorData(double curTime) {
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
#endif
}




