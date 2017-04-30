/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <iostream>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 10;
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2] ; // Standard deviations for x, y, and theta
	std::default_random_engine gen;

	// This line creates a normal (Gaussian) distribution for x
	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

	for(unsigned i = 0 ; i < num_particles;i++){
		Particle particle = Particle();
		particle.id = i ;
		particle.theta = dist_theta(gen);
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
        particle.weight = 1.0 ;
		particles.push_back(particle);
        weights.push_back(particle.weight);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2] ; // Standard deviations for x, y, and theta
  std::default_random_engine gen;

  if(fabs(yaw_rate)<0.001){
      for(unsigned i = 0 ; i < num_particles;i++){

        Particle& particle = particles.at(i);
        // This line creates a normal (Gaussian) distribution for x

        particle.x = particle.x + delta_t*velocity*cos(particle.theta);
        particle.y = particle.y + delta_t*velocity*sin(particle.theta);
        particle.theta = particle.theta;

        std::normal_distribution<double> dist_x(particle.x, std_x);
        std::normal_distribution<double> dist_y(particle.y, std_y);
        std::normal_distribution<double> dist_theta(particle.theta, std_theta);

        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);


      }
    }else{
      for(unsigned i = 0 ; i < num_particles;i++){
        Particle& particle = particles.at(i);
        // This line creates a normal (Gaussian) distribution for x

        particle.x = particle.x + (velocity/yaw_rate)*(sin(particle.theta + yaw_rate*delta_t)-sin(particle.theta));
        particle.y = particle.y + (velocity/yaw_rate)*(cos(particle.theta)-cos(particle.theta+yaw_rate*delta_t));
        particle.theta = particle.theta + delta_t*yaw_rate;

        std::normal_distribution<double> dist_x(particle.x, std_x);
        std::normal_distribution<double> dist_y(particle.y, std_y);
        std::normal_distribution<double> dist_theta(particle.theta, std_theta);

        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
      }
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  // Iterate through observations
  for(LandmarkObs& currentObservation:observations){

    double min_distance = INFINITY;
    int min_id = 0 ; 
    for(LandmarkObs& currentPred:predicted){
      double distance = dist(currentPred.x, currentPred.y, currentObservation.x, currentObservation.y);
      if(distance < min_distance){
        min_distance = distance;
        min_id = currentPred.id;
      }
      currentObservation.id=min_id;
    }
  }

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
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

  for (unsigned i = 0; i < particles.size(); i++)
  {
    Particle& current_particle = particles.at(i);

    /**
     * Translated observations
     */
    std::vector<LandmarkObs> translated_observations;
    /**
     * Translate vehicle coordinates into map coordinates
     */
    for (LandmarkObs& cur_obs : observations)
    {
      LandmarkObs landmark_t;
      landmark_t.id = cur_obs.id;
      landmark_t.x = cur_obs.x * cos(current_particle.theta) - cur_obs.y * sin(current_particle.theta) + current_particle.x;
      landmark_t.y = cur_obs.x * sin(current_particle.theta) + cur_obs.y * cos(current_particle.theta) + current_particle.y;
      translated_observations.push_back(landmark_t);
    }

    std::vector<LandmarkObs> valid_landmarks;

    for (Map::single_landmark_s landmark : map_landmarks.landmark_list) {
      /**
       * If landmark is within sensor range of the particle,
       * keep it else discard it from the list of landmarks for the given particle
       */
      if (dist(landmark.x_f, landmark.y_f, current_particle.x, current_particle.y) <= sensor_range) {
        LandmarkObs lm  = LandmarkObs();
        lm.id = landmark.id_i;
        lm.x = landmark.x_f;
        lm.y = landmark.y_f;
        valid_landmarks.push_back(lm);
      }
    }

    /**
     * Associate each observation with the corresponding landmark
     */
    dataAssociation(valid_landmarks, translated_observations);


    /**
     * Update weights
     */
    double new_weight = 1;
    for (LandmarkObs& current_observation : translated_observations)
    {
      Map::single_landmark_s cur_pred = map_landmarks.landmark_list[current_observation.id-1];

      double dx = current_observation.x - cur_pred.x_f;
      double dy = current_observation.y - cur_pred.y_f;

      double updated_weight = 1 / (M_PI * 2 * std_landmark[0] * std_landmark[1]) *
          std::exp(-1 * (((dx*dx) / (std_landmark[0]*std_landmark[0])) + ((dy*dy) / (std_landmark[1]*std_landmark[1]))));

      new_weight *= updated_weight;
    }

    current_particle.weight = new_weight;
    weights.at(i) = new_weight;
  }


  /**
   * Normalize weights
   */
  double weight_sum = 0 ;

  for(double weight:weights){
    weight_sum+=weight;
  }

  for(double& weight:weights){
    weight = weight/weight_sum;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::default_random_engine gen;
  std::discrete_distribution<int> dist {weights.begin(), weights.end()};

  std::vector<Particle> updated_particles;

  for (unsigned i = 0; i < num_particles; i++) {
    int current_particle_index = dist(gen);
    Particle cur_particle = particles[current_particle_index];
    updated_particles.push_back(cur_particle);
  }
  particles = updated_particles;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
