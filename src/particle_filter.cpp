/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#include <map>

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

num_particles = 50;

default_random_engine gen;
normal_distribution<double> dist_x(x, std[0]);
normal_distribution<double> dist_y(y, std[1]);
normal_distribution<double> dist_theta(theta, std[2]);

for (int ii=0; ii<num_particles; ++ii){
  
  weights.push_back(1.0);
  
  Particle myparticle;
  myparticle.id = ii;
  myparticle.x = dist_x(gen);
  myparticle.y = dist_y(gen);
  myparticle.theta = dist_theta(gen);
  myparticle.weight = 1.0;
  particles.push_back(myparticle);
  
}

is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

default_random_engine gen;
normal_distribution<double> dist_x(0.0, std_pos[0]);
normal_distribution<double> dist_y(0.0, std_pos[1]);
normal_distribution<double> dist_theta(0.0, std_pos[2]);

/*
double x, y, theta;

if (fabs(yaw_rate) > 1.0E-6) {
  double vr = velocity / yaw_rate;
  double angle = yaw_rate * delta_t;
  for (std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
    x = it->x;
    y = it->y;
    theta = it->theta;
    it->x = x + vr*(sin(theta + angle) - sin(theta)) + dist_x(gen);
    it->y = y + vr*(cos(theta) - cos(theta + angle)) + dist_y(gen);
    theta = theta + angle + dist_theta(gen);
    if (theta >= 0) {
      theta = theta - 2*M_PI*floor(theta/(2*M_PI));
    }else{
      theta = theta + 2*M_PI*floor(fabs(theta)/(2*M_PI));
    }
    it->theta = theta;
  }
}else{
  for (std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
    x = it->x;
    y = it->y;
    theta = it->theta;
    it->x = x + velocity * cos(theta) + dist_x(gen);
    it->y = y + velocity * sin(theta) + dist_y(gen);
    theta = theta + dist_theta(gen);
    if (theta >= 0) {
      theta = theta - 2*M_PI*floor(theta/(2*M_PI));
    }else{
      theta = theta + 2*M_PI*floor(fabs(theta)/(2*M_PI));
    }
    it->theta = theta;
  }
}

*/


for (int ii=0; ii<num_particles; ++ii) {
  if (yaw_rate != 0){
    particles[ii].x = particles[ii].x + velocity / yaw_rate*(sin(particles[ii].theta + yaw_rate * delta_t) - sin(particles[ii].theta)) + dist_x(gen);
    particles[ii].y = particles[ii].y + velocity / yaw_rate*(cos(particles[ii].theta) - cos(particles[ii].theta + yaw_rate * delta_t)) + dist_y(gen);
    particles[ii].theta = particles[ii].theta + yaw_rate * delta_t + dist_theta(gen);
  } else {
    particles[ii].x = particles[ii].x + velocity * cos(particles[ii].theta) + dist_x(gen);
    particles[ii].y = particles[ii].y + velocity * sin(particles[ii].theta) + dist_y(gen);
    particles[ii].theta = particles[ii].theta + dist_theta(gen);
  } 
}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> landmarks, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  
  double distance_map_min = 3.42; //half min distance between any pair of landmarks in given map. used to reducee search time. it can be computed when a new map is loaded.

  int count_outofrange = 0;
  for (int ii=0; ii<observations.size(); ++ii) {

    double distance_min = std::numeric_limits<unsigned int>::max();
    double obs_x = observations[ii].x;
    double obs_y = observations[ii].y;
    int landmark_id;
    
    for (int jj=0; jj<landmarks.size(); ++jj) {
      double landmark_x, landmark_y;
      landmark_x = landmarks[jj].x;
      landmark_y = landmarks[jj].y;
      double distance = dist(obs_x, obs_y, landmark_x, landmark_y);
      if (distance < distance_min) {
          distance_min = distance;
          landmark_id = landmarks[jj].id; 
      }
    } //for jj

/*
    //look for nearest landmark to each observation. 
    int jj=0;
    bool found = false;
    while (!found && jj<landmarks.size()){
      double landmark_x, landmark_y;
      landmark_x = landmarks[jj].x;
      landmark_y = landmarks[jj].y;
      double distance = dist(obs_x, obs_y, landmark_x, landmark_y);
      if (distance < distance_min) {
        distance_min = distance;
        landmark_id = landmarks[jj].id; 
      }
      if (distance_min < distance_map_min) {
        found = true;
      }
      jj++;
    }//while

*/
    observations[ii].id = landmark_id;
    
  }//for ii

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

std::vector<int> associations;
std::vector<double> sense_x;
std::vector<double> sense_y;

long double factor = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);

//init landmarks to be used by dataAssociation()
std::map<int, Map::single_landmark_s> landmarks_dictionary;
std::vector<LandmarkObs> landmarks_vector;
for (int jj=0; jj<map_landmarks.landmark_list.size(); ++jj) {
  LandmarkObs landmark;
  landmark.id = map_landmarks.landmark_list[jj].id_i;
  landmark.x = map_landmarks.landmark_list[jj].x_f;
  landmark.y = map_landmarks.landmark_list[jj].y_f;
  landmarks_vector.push_back(landmark);
  landmarks_dictionary[landmark.id] = map_landmarks.landmark_list[jj];
}

weights.clear();
for (int ip=0; ip<num_particles; ++ip) {

  associations.clear();
  sense_x.clear();
  sense_y.clear();

  //transform coordinates
  std::vector<LandmarkObs> observations_trans;
  double x, y, theta;
  x = particles[ip].x;
  y = particles[ip].y;
  theta = particles[ip].theta;
  for (int ii=0; ii<observations.size(); ++ii) {
    LandmarkObs obs_trans, obs; 
    obs = observations[ii];
    obs_trans.x = x + obs.x*cos(theta) - obs.y*sin(theta);
    obs_trans.y = y + obs.x*sin(theta) + obs.y*cos(theta);
    observations_trans.push_back(obs_trans);
  }//ii

  //look for the nearest landmarks to a setof observations. 
  dataAssociation(landmarks_vector, observations_trans);

  for (int ii=0; ii<observations_trans.size(); ++ii) {
    double obs_x = observations_trans[ii].x;
    double obs_y = observations_trans[ii].y;
    double landmark_x = landmarks_dictionary[observations_trans[ii].id].x_f;
    double landmark_y = landmarks_dictionary[observations_trans[ii].id].y_f;
    double distance = dist(obs_x, obs_y, landmark_x, landmark_y);
    if (distance < sensor_range) {
      associations.push_back(observations_trans[ii].id);
      sense_x.push_back(observations_trans[ii].x);
      sense_y.push_back(observations_trans[ii].y);
    }
  }//ii

  SetAssociations(particles[ip], associations, sense_x, sense_y);
  
  //update weigth using multivariate gaussian probability distribution
  
  long double new_weight = 1.0;
  for (int ii=0; ii<associations.size(); ++ii) {
    int landmark_id = associations[ii];
    double obs_x = sense_x[ii];
    double obs_y = sense_y[ii];
    double landmark_x = landmarks_dictionary[landmark_id].x_f;
    double landmark_y = landmarks_dictionary[landmark_id].y_f;
    double dx = pow((obs_x-landmark_x)/std_landmark[0],2);
    double dy = pow((obs_y-landmark_y)/std_landmark[1],2);
    double multiplier = factor*exp(-0.5*(dx+dy));
    if (multiplier > 0) {
      new_weight *= multiplier;
    }
  }//ii
  
  particles[ip].weight = new_weight;
  weights.push_back(new_weight);


}//for ip


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

default_random_engine gen;
std::discrete_distribution<> distribution(weights.begin(), weights.end());
std::vector<Particle> particles_new;

for (int ii=0; ii<num_particles; ++ii) {
  particles_new.push_back(particles[distribution(gen)]);
}

particles = particles_new;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

