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
#include <unordered_map>

#include "particle_filter.h"

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1000;	

	normal_distribution<double> dist_x(0, std[0]);
	normal_distribution<double> dist_y(0, std[1]);
	normal_distribution<double> dist_theta(0, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		particle.id = i;
		particle.x = x + dist_x(gen);
		particle.y = y + dist_y(gen);
		particle.theta = theta + dist_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	for (int i = 0; i < particles.size(); ++i) {
		if (fabs(yaw_rate) < 0.001) {
			particles[i].x += velocity * cos(particles[i].theta) * delta_t + dist_x(gen);
			particles[i].y += velocity * sin(particles[i].theta) * delta_t + dist_y(gen);
			particles[i].theta += dist_theta(gen);
		}
		else {
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(gen);
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) + dist_y(gen);
			particles[i].theta += yaw_rate * delta_t + dist_theta(gen);	
		}
	}
}	

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// for each obsereved landmarks
	for (int i = 0; i < observations.size(); ++i) {
		// observed position
		double x = observations[i].x;
		double y = observations[i].y;
		
		// minimum distance and ID to be set
		double minDist = numeric_limits<double>::max();
		int assignID = -1;

		// get minimum distance and associated ID
		for (auto landmark : predicted) {
			double curDist = dist(x, y, landmark.x, landmark.y);
			if (curDist < minDist) {
				minDist = curDist;
				assignID = landmark.id; 
			}
		}

		// assign id to the observation
		observations[i].id = assignID;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	unordered_map<int, pair<float, float>> hash_landmarks;
	for (auto map_landmark : map_landmarks.landmark_list) { 
		hash_landmarks[map_landmark.id_i] = pair<float, float>(map_landmark.x_f, map_landmark.y_f);
	}

	// for each particle
	for (int i = 0; i < particles.size(); ++i) {
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		// vector to collect all the landmarks within sensor range
		vector<LandmarkObs> predicted;

		// collect landmarkers within particle's sensor range
		for (auto map_landmark : map_landmarks.landmark_list) {
			double curDist = dist(x, y, map_landmark.x_f, map_landmark.y_f);
			if (curDist <= sensor_range) {
				predicted.push_back(LandmarkObs{map_landmark.id_i, map_landmark.x_f, map_landmark.y_f});
			}
		}

		// transform observation from current particle's frame to world frame
		vector<LandmarkObs> transObs;
		for (auto observation : observations) {
			double wordlX = cos(theta) * observation.x - sin(theta) * observation.y + x;
			double wordlY = sin(theta) * observation.x + cos(theta) * observation.y + y;
			transObs.push_back(LandmarkObs{observation.id, wordlX, wordlY});
		}

		// data association for current observation
		dataAssociation(predicted, transObs);

		// update particle's wheight
		double weight = 1.0;
		for (auto transOb : transObs) {
			// define inputs
			double sig_x= std_landmark[0];
			double sig_y= std_landmark[1];
			double xObs = transOb.x; 
			double yObs = transOb.y;
			double mu_x = hash_landmarks[transOb.id].first;
			double mu_y = hash_landmarks[transOb.id].second;

			// calculate normalization term
			double obsW = (1/(2 * M_PI * sig_x * sig_y)) * exp(- 0.5 * (pow(xObs - mu_x, 2) / pow(sig_x, 2) + pow(yObs - mu_y, 2) / pow(sig_y, 2)));
			weight *= obsW;
		}
		particles[i].weight = weight;
	}
}

void ParticleFilter::resample() {
	vector<Particle> newParticles;

	// get max weight
	double maxW = 0.0;
	for (auto particle : particles) {
		maxW = max(maxW, particle.weight);
	}

	// starting index distribution
	uniform_int_distribution<int> dist_index(0, num_particles - 1);
	int index = dist_index(gen);
	// beta distribution
	uniform_real_distribution<double> dist_beta(0, 2 * maxW);
	double beta = 0.0;
	for (int i = 0; i < num_particles; ++i) {
		beta += dist_beta(gen);
		while (beta > particles[index].weight) {
			beta -= particles[index].weight;
			index = (index + 1) % num_particles;
		}
		newParticles.push_back(particles[index]);
	}

	particles = newParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
