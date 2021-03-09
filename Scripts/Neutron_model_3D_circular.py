import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt 
import random
import math
import time
import argparse

def randomangle_Cone(alpha):
	theta = math.pi/2 -np.arccos(random.random()*np.cos(math.pi/2-alpha/2))
	phi   = random.random()*2*math.pi
	A = np.sqrt(1+np.tan(phi)**2)
	if math.pi/2 < phi <= 3*math.pi/2:
		ZX = - np.arctan(np.tan(theta)/A)
		ZY = - np.arctan(np.tan(theta)*np.tan(phi)/A)
	else:
		ZX = np.arctan(np.tan(theta)/A)
		ZY = np.arctan(np.tan(theta)*np.tan(phi)/A)
	return(ZX,ZY,theta,phi)

def position_uniform(lower_generator,upper_generator):
	return (lower_generator + random.random()*(upper_generator-lower_generator),lower_generator + random.random()*(upper_generator-lower_generator))

def format_radians_label(float_in):
	# Converts a float value in radians into a
	# string representation of that float
	string_out = str(float_in / (np.pi))+"PI"
	
	return string_out

def convert_polar_xticks_to_radians(ax):
	# Converts x-tick labels from degrees to radians
	
	# Get the x-tick positions (returns in radians)
	label_positions = ax.get_xticks()
	
	# Convert to a list since we want to change the type of the elements
	labels = list(label_positions)
	
	# Format each label (edit this function however you'd like)
	labels = [format_radians_label(label) for label in labels]
	
	ax.set_xticklabels(labels)

def convert_to_alphaXY(theta,phi):
	A = np.sqrt(1+np.tan(phi)**2)
	if math.pi/2 < phi <= 3*math.pi/2:
		ZX = - np.arctan(np.tan(theta)/A)
		ZY = - np.arctan(np.tan(theta)*np.tan(phi)/A)
	else:
		ZX = np.arctan(np.tan(theta)/A)
		ZY = np.arctan(np.tan(theta)*np.tan(phi)/A)
	return(ZX,ZY)

def calculate_impact_to_object(dist_to_object,bursts,PosXX,PosYY,Theta,Phi):
	AngXX,AngYY = convert_to_alphaXY(Theta,Phi)
	#calculate intermediate point
	a = PosXX - dist_to_object * np.tan(AngXX)
	b = PosYY - dist_to_object * np.tan(AngYY)
	#shape condition, in this case, a circle centered in (10,10) with radius 5.
	centreX = 20
	centreY = 20
	if 100 > ((a - centreX)**2 + (b - centreY)**2) and math.floor(np.arctan((b-centreY)/(a-centreX))*bursts/math.pi)%2 == True:
		impact = True
	else:
		impact = False

	return impact

def neutron_transport(LD,upper_generator,lower_generator,number_of_channels,col_D, col_L, dist_to_col,dist_to_det,iterations,foldercount):

	Folder = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/PaviaC"+foldercount
	
	alpha = np.arctan(1/LD)
	dist_to_object = 10
	bursts = 128

	file_metadata = open(Folder+"/Metadata.txt","a")
	timestamp_start = time.time()
	Message = ""
	file_metadata.write("Description: "+ Message+"\n")
	file_metadata.write("Characteristics: \nLD: "+str(LD))
	file_metadata.write("\nExpected limit angle: "+str(alpha))
	file_metadata.write("\nupper generator: "+str(upper_generator))
	file_metadata.write("\nlower generator: "+str(lower_generator))
	file_metadata.write("\nnumber of channels: "+str(number_of_channels))
	file_metadata.write("\ncollimator L: "+str(col_L))
	file_metadata.write("\ncollimator D: "+str(col_D))
	file_metadata.write("\nDistance to collimator: "+str(dist_to_col))
	file_metadata.write("\nDistance to detector: "+str(dist_to_det))
	file_metadata.write("\nNumber of iterations: "+str(iterations))

	fileindex = 0
	file_positions = open(Folder+"/Testrun-of-full-sim-positionXY-index"+str(fileindex)+".txt","a")
	file_AnglesTheta = open(Folder+"/Testrun-of-full-sim-AngleTheta-index"+str(fileindex)+".txt","a")
	file_AnglesPhi = open(Folder+"/Testrun-of-full-sim-AnglePhi-index"+str(fileindex)+".txt","a")
	hit = 0
	
	print(alpha)
	upper_collimator = col_D * (number_of_channels * 2)
	upper_collimator_wide = col_D * ((number_of_channels * 2) + 1)

	rangemin = 0
	rangemax = 50

	for it in range(iterations):
		if (it+1)%5000000 == 0: 
			print("Checkpoint working on file number "+str(fileindex)+", iterations so far: "+str(it+1))
		if hit >= 1000000:
			file_positions.close()
			file_AnglesTheta.close()
			file_AnglesPhi.close()
			fileindex+=1
			file_positions = open(Folder+"/Testrun-of-full-sim-positionXY-index"+str(fileindex)+".txt","a")
			file_AnglesTheta = open(Folder+"/Testrun-of-full-sim-AngleTheta-index"+str(fileindex)+".txt","a")
			file_AnglesPhi = open(Folder+"/Testrun-of-full-sim-AnglePhi-index"+str(fileindex)+".txt","a")
			hit = 0
		
		ZX,ZY,theta,phi = randomangle_Cone(alpha)
		neutron_position_X,neutron_position_Y = position_uniform(lower_generator,upper_generator)
		LD_neutron_X = np.tan(ZX)
		LD_neutron_Y = np.tan(ZY)
		posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y

		hli_X = dist_to_col * LD_neutron_X + neutron_position_X
		hlf_X = (dist_to_col + col_L) * LD_neutron_X + neutron_position_X
		hli_Y = dist_to_col * LD_neutron_Y + neutron_position_Y
		hlf_Y = (dist_to_col + col_L) * LD_neutron_Y + neutron_position_Y
		
		if (hli_X-20)**2 + (hli_Y-20)**2 < 225 and (hlf_X-20)**2 + (hlf_Y-20)**2 < 225:
			posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
			impact = calculate_impact_to_object(dist_to_object,bursts,posX,posY,theta,phi)
			if rangemin < posX < rangemax and rangemin < posY < rangemax and impact == False:
				file_positions.write(str([posX,posY])+",")
				file_AnglesTheta.write(str(theta)+",")
				file_AnglesPhi.write(str(phi)+",")
				hit += 1

	file_positions.close()
	file_AnglesTheta.close()
	file_AnglesPhi.close()

	file_metadata.write("\nTime required for execution: "+str(time.time()-timestamp_start)+"s")
	file_metadata.close()

	message = "Everything's fine"
	return message

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-x","--LD",type=str,default="27.2",help="Divergence of the initial beam, in LD")
	parser.add_argument("-i","--iterations",type=str,default="200000000",help="Number of iterations")
	parser.add_argument("-c","--channels",type=str,default="20",help="Number of channels")
	parser.add_argument("-D","--col_D",type=str,default="2.5",help="Collimator's channel width")
	parser.add_argument("-L","--col_L",type=str,default="5000",help="Collimator's channel length")
	parser.add_argument("-t","--dist_to_col",type=str,default="500",help="Distance from source to collimator")
	parser.add_argument("-f","--dist_to_det",type=str,default="50",help="Distance collimator to detector")
	parser.add_argument("-n","--filename",type=str,default="Collimator-3D-model_",help="Basic title for the saved files")
	parser.add_argument("-g","--foldercount",type=str,default="0",help="Adjoint number for the folder count")
	args=parser.parse_args()
	
	#LD of the initial beam
	#LD = 27.2
	LD = float(args.LD)
	#Monkey control
	if LD == 0:
		LD = 5
	#upper height of the generator
	upper = 70
	lower = -20
	
	number_of_channels = int(args.channels)
	
	if number_of_channels <= 2:
		number_of_channels = 2
	
	#col_D = 2.5
	#col_L = 110
	col_D = float(args.col_D)
	col_L = float(args.col_L)

	dist_to_col = int(args.dist_to_col)
	dist_to_det = int(args.dist_to_det)
	iterations = int(args.iterations)

	maintitle = args.filename
	foldercount = args.foldercount
	print("--> Starting run ",0,"<--")
	start_time = time.time()

	#histogram = neutron_transport(LD,upper,lower,number_of_channels,col_D, col_L, dist_to_col,dist_to_det,iterations)
	message = neutron_transport(LD,upper,lower,number_of_channels,col_D, col_L, dist_to_col,dist_to_det,iterations,foldercount)
	print(message)
	print("--- %s seconds ---" % (time.time() - start_time))
	'''
	with open("Test3/Testrun-of-full-sim-positionXY.txt", 'r') as file:
		position = np.array(eval(file.read())) # read list string and convert to array
	file.close()
	
	plot3D(position[:,0],position[:,1],run,True)
	
	X = []
	Y = []	
	for ble in range(len(data["histogram_white"])):
		X.append(data["histogram_white"][ble][0])
		Y.append(data["histogram_white"][ble][1])
	plot3D(X,Y,run,False)
	'''
	print("--> Completed run ",0,"<--")
	print("\n","============================","\n")

