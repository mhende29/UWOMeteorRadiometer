from __future__ import print_function, division, absolute_import

import argparse
import os
import sys
import math
import csv
import numpy as np
import scipy.linalg
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def sphere2cartesian(theta, phi):
	""" Convert the polar angle and azimuthal angle to cartesian coordinates. Assumes the radius r to be unity. 
	Arguments:
	    theta: [float] polar angle (90 - elevation)
	    phi: [float] azimuthal angle
	Return:
	    (x, y, z): [tuple of floats] Cartesian coordinates
	"""
	theta = np.deg2rad(theta)
	phi = np.deg2rad(phi)

	x = np.cos(theta)*np.cos(phi)
	y = np.cos(theta)*np.sin(phi)
	z = np.sin(theta)

	return x, y, z

def cartesian2sphere(x, y, z):
	""" Convert the cartesian coordinates to spherical coordinates. Assumes the radius r to be unity and therefore doesn't return it. 
	Arguments:
	    x: [float] distance along the x axis
	    y: [float] distance along the y axis
	    z: [float] distance along the z axis 
	Return:
	    (theta, phi): [tuple of floats] Spherical coordinates
	"""
	r = np.sqrt(x**2+y**2+z**2)
	theta = np.arccos(z/r)
	phi = np.arctan2(y,x)

	return theta, phi

def getAngle(p1,p2,p0=None):
	""" Computes the angle between 2 points in space with respect to a point. 
	Arguments:
	    p1: [1D ndarray] Coordinates of the first point
	    p2: [1D ndarray] Coordinates of the second point
	    p0: [1D ndarray] Coordinates of the origin point, zeros by default
	Return:
	    angle: [float] Angle in radians between the points
	"""	
	if(p0 is None):
		return np.arccos(np.dot(p1,p2)/(np.linalg.norm(p1)*np.linalg.norm(p2)))
	else:
		return np.arccos(np.dot(p1-p0,p2-p0)/(np.linalg.norm(p1-p0)*np.linalg.norm(p2-p0)))


def greatCircle(t, theta0, phi0):
	""" 
	Calculates the point on a great circle defined my theta0 and phi0 in Cartesian coordinates. 

	Sources:
		- http://demonstrations.wolfram.com/ParametricEquationOfACircleIn3D/
	Arguments:
		t: [float or 1D ndarray] phase angle of the point in the great circle
		theta0: [float] inclination of the great circle
		phi0: [float] nodal angle of the great circle
	Return:
		[tuple or 2D ndarray] a tuple of (X, Y, Z) coordinates in 3D space (becomes a 2D ndarray if the input
			parameter t is also a ndarray)
	"""
	# Calculate individual cartesian components of the great circle points
	x = -np.cos(t)*np.sin(phi0) + np.sin(t)*np.cos(theta0)*np.cos(phi0)
	y =  np.cos(t)*np.cos(phi0) + np.sin(t)*np.cos(theta0)*np.sin(phi0)
	z =  np.sin(t)*np.sin(theta0)

	return x, y, z

def fitGreatCircle(x, y, z):
	""" Fits a great circle to points in 3D space. 
	Arguments:
	"""

	# Add (0, 0, 0) to the data, as the great circle should go through the origin
	x = np.append(x, 0)
	y = np.append(y, 0)
	z = np.append(z, 0)

	# Fit a linear plane through the data points
	A = np.c_[x, y, np.ones(x.shape[0])]
	C,_,_,_ = scipy.linalg.lstsq(A, z)

	# Calculate the great circle parameters
	z2 = C[0]**2 + C[1]**2

	theta0 = np.arcsin(z2/np.sqrt(z2 + z2**2))
	phi0 = np.arctan2(C[1], C[0])

	return C, theta0, phi0

def getTraj2sphere(filename):
	traj_data = []
	with open(filename,"r") as file:
		data = file.readlines()
		for row in data:
			if("#" not in row or "unix" in row):
				traj_data.append(row)
	time_start_index = traj_data[0].find(": ") + len(": ")
	time_string = traj_data[0][time_start_index:].replace("\n","")
	time = float(time_string)
	times = []
	polar_angle = []
	elevation = []
	azimuth_north = []
	for row in traj_data:
		if ("#" not in row):
			time_start_index = 13
			elevation_index = 54
			azimuth_north_index = 63
			
			offset_len = 8
			elevation_len = 8
			azimuth_len = 8

			times.append(time + float(row[time_start_index - offset_len:time_start_index]))
			polar_angle.append(float(row[elevation_index - elevation_len:elevation_index]))
			elevation.append(90 - float(row[elevation_index - elevation_len:elevation_index]))
			azimuth = 90 - float(row[azimuth_north_index - azimuth_len:azimuth_north_index])
			# Radiometer azimuth is 90 ahead 90-azimuth
			if(azimuth>360):
				azimuth-=360
			elif(azimuth<0):
				azimuth+=360
			azimuth_north.append(azimuth)

	return np.array(times), np.array(polar_angle), np.array(elevation), np.array(azimuth_north)

def sensitivityCorrection():
	pass

if __name__ == "__main__":

	sps = 2100

	path = "/home/michael/Radiometer Trajectories/ev_20180815_073510A_02A.txt"
	time, given_polar_angle, given_elevation, given_azimuth = getTraj2sphere(path)
	x,y,z = sphere2cartesian(given_polar_angle, given_azimuth)
	C, theta0, phi0 = fitGreatCircle(x, y, z)
	time_diff = time[-1] - time[0]
	angle_diff = np.arccos(np.dot([x[0],y[0],z[0]],[x[-1],y[-1],z[-1]])/(np.linalg.norm([x[0],y[0],z[0]])*np.linalg.norm([x[-1],y[-1],z[-1]])))
	angular_velocity = angle_diff/time_diff
	print(angle_diff,time_diff,angular_velocity)
	
	params = greatCircle(0, theta0, phi0)
	firstpoint = [x[-1],y[-1],z[-1]]
	start_angle = getAngle(params,firstpoint)
	t_array = np.arange(start_angle, start_angle+angle_diff, angular_velocity/sps)
	x_fit, y_fit, z_fit = greatCircle(t_array, theta0, phi0)
	

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.scatter(x,y,z)
	ax.scatter(x_fit, y_fit, z_fit, c='b', s=3)
	ax.set_aspect('equal')

	theta_fit, phi_fit = cartesian2sphere(x_fit, y_fit, z_fit)
	theta_fit_deg = np.rad2deg(theta_fit)
	phi_fit_deg = np.rad2deg(phi_fit)

	if(np.abs(given_elevation[0]-theta_fit_deg[0])>np.abs(given_elevation[0]-theta_fit_deg[-1])):
		theta_fit_deg = np.flip(theta_fit_deg,0)
	if(np.abs(given_azimuth[0]-phi_fit_deg[0])>np.abs(given_azimuth[0]-phi_fit_deg[-1])):
		phi_fit_deg = np.flip(phi_fit_deg,0)

	#sys.exit()
	plt.show()