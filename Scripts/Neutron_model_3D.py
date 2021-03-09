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

def randomangle_Square(alpha):
	X = random.random() * alpha - (alpha / 2)
	Y = random.random() * alpha - (alpha / 2)
	return (X,Y)

def position_uniform(lower_generator,upper_generator):
	return (lower_generator + random.random()*(upper_generator-lower_generator),lower_generator + random.random()*(upper_generator-lower_generator))

def position_gaussian(lower_generator,upper_generator,upper_collimator_wide):
	neutron_position_X,neutron_position_Y = 10000,10000
	while(neutron_position_X > upper_generator or neutron_position_Y > upper_generator or neutron_position_X < lower_generator or neutron_position_Y < lower_generator):
		rng1 = random.random()
		rng2 = random.random()
		neutron_position_X = np.sqrt(-2*np.log(rng1))*np.cos(2*math.pi*rng2) * (upper_generator-lower_generator)/5 + upper_collimator_wide/2 
		neutron_position_Y = np.sqrt(-2*np.log(rng1))*np.sin(2*math.pi*rng2) * (upper_generator-lower_generator)/5 + upper_collimator_wide/2
	return (neutron_position_X,neutron_position_Y)

def neutron_transport(LD,upper_generator,lower_generator,number_of_channels,col_D, col_L, dist_to_col,dist_to_det,iterations,foldercount):
	#histogram = []
	#histogram_white = []
	#X_distribution = []
	#Y_distribution = []
	'''
	data = {}
	data["histogram"] = []
	data["histogram_white"] = []
	data["X_distribution"] = []
	data["Y_distribution"] = []
	data["X_distribution_white"] = []
	data["Y_distribution_white"] = []
	'''
	Folder = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/OPavia"+foldercount
	
	alpha = np.arctan(1/LD)

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
	#if upper_collimator > upper_generator:
	#	upper_generator = upper_collimator - lower_generator
	#print(LD,upper_generator,lower_generator,upper_collimator,col_D,col_L)
	
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
		
		#Define the position and angle of the neutron generated in 2D
		neutron_position_X,neutron_position_Y = position_uniform(lower_generator,upper_generator)
		#neutron_position_X,neutron_position_Y = position_gaussian(lower_generator,upper_generator,upper_collimator_wide)
		
		#neutron_angle_X,neutron_angle_Y = randomangle_Square(alpha)
		neutron_angle_X,neutron_angle_Y,theta,phi = randomangle_Cone(alpha)
		
		LD_neutron_X = np.tan(neutron_angle_X)
		LD_neutron_Y = np.tan(neutron_angle_Y)
		
		#Calculate entering an exiting positions in the collimator
		hli_X = dist_to_col * LD_neutron_X + neutron_position_X
		hlf_X = (dist_to_col + col_L) * LD_neutron_X + neutron_position_X
		hli_Y = dist_to_col * LD_neutron_Y + neutron_position_Y
		hlf_Y = (dist_to_col + col_L) * LD_neutron_Y + neutron_position_Y
		hlmid_X,hlmid_Y = 0,0
		bordercase_X,bordercase_Y = False,False
		
		
		if hli_X < upper_collimator_wide < hlf_X and col_D < hli_Y < upper_collimator:
			hlmid_X = upper_collimator
			bordercase_X = True
		elif hli_X < 0 < hlf_X and col_D < hli_Y < upper_collimator:
			hlmid_X = 0
			bordercase_X = True
		elif hlf_X < upper_collimator_wide < hli_X and col_D < hli_Y < upper_collimator:
			hlmid_X = upper_collimator
			bordercase_X = True
		elif hlf_X < 0 < hli_X and col_D < hli_Y < upper_collimator:
			hlmid_X = 0
			bordercase_X = True
		
		if hli_Y < upper_collimator_wide < hlf_Y and col_D < hli_X < upper_collimator:
			hlmid_Y = upper_collimator
			bordercase_Y = True
		elif hli_Y < 0 < hlf_Y and col_D < hli_X < upper_collimator:
			hlmid_Y = 0
			bordercase_Y = True
		elif hlf_Y < upper_collimator_wide < hli_Y and col_D < hli_X < upper_collimator:
			hlmid_Y = upper_collimator
			bordercase_Y = True
		elif hlf_Y < 0 < hli_Y and col_D < hli_X < upper_collimator:
			hlmid_Y = 0
			bordercase_Y = True
			
		#if 10 < neutron_position_X < 30 and 10 < neutron_position_Y < 30:
			#data["histogram_white"].append([(dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y])					
			#data["X_distribution_white"].append(neutron_angle_X)
			#data["Y_distribution_white"].append(neutron_angle_Y)
		#print("border X: ",bordercase_X," border Y: ",bordercase_Y)
		#Start the magic: general case: if the neutron is inside the collimator's area, and the start and end cell
		#are the same, the neutron is evaluated.
		if (0 <= hli_X <= upper_collimator_wide) and (0 <= hli_Y <= upper_collimator_wide) and (bordercase_X == False) and (bordercase_Y == False):
			if (math.floor(hli_X/col_D) % 2) != 0:
				if (math.floor(hli_Y/col_D) % 2) == 0:
					if math.floor(hlf_X/col_D) == math.floor(hli_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hli_Y/col_D):
						posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
						if rangemin < posX < rangemax and rangemin < posY < rangemax:
							file_positions.write(str([posX,posY])+",")
							file_AnglesTheta.write(str(theta)+",")
							file_AnglesPhi.write(str(phi)+",")
							hit+=1
			else:
				if (math.floor(hli_Y/col_D) % 2) != 0: 
					if math.floor(hlf_X/col_D) == math.floor(hli_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hli_Y/col_D):
						posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
						if rangemin < posX < rangemax and rangemin < posY < rangemax:
							file_positions.write(str([posX,posY])+",")
							file_AnglesTheta.write(str(theta)+",")
							file_AnglesPhi.write(str(phi)+",")
							hit+=1
			#This elif determines the case in which the neutron enters either the left or right border
			#there are 2 casualities: either enters within and exits in the middle, or enters in the middle and exits properly
		elif bordercase_X == True:
			if 0 < hli_X < upper_collimator_wide:
				#print("Border X from inside")
				if (math.floor(hli_X/col_D) % 2) != 0:
					if (math.floor(hli_Y/col_D) % 2) == 0:
						if math.floor(hlmid_X/col_D) == math.floor(hlmid_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hli_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
				else:
					if (math.floor(hli_Y/col_D) % 2) != 0: 
						if math.floor(hlmid_X/col_D) == math.floor(hli_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hli_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
			else:
				#print("Border X from outside")
				if (math.floor(hlmid_X/col_D) % 2) != 0:
					if (math.floor(hli_Y/col_D) % 2) == 0:
						if math.floor(hlf_X/col_D) == math.floor(hlmid_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hli_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
				else:
					if (math.floor(hli_Y/col_D) % 2) != 0: 
						if math.floor(hlf_X/col_D) == math.floor(hlmid_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hli_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
				
			#The same process is followed for upper and lower borders:
		elif bordercase_Y == True:
			if 0 < hli_Y < upper_collimator_wide:
				#print("Border Y from inside ","hli_X: ",hli_X,"hli_Y: ",hli_Y,"hlmid_Y: ",hlmid_Y,"hlf_Y: ",hlf_Y)
				if (math.floor(hli_X/col_D) % 2) != 0:
					if (math.floor(hli_Y/col_D) % 2) == 0:
						if math.floor(hlf_X/col_D) == math.floor(hli_X/col_D) and math.floor(hlmid_Y/col_D) == math.floor(hli_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
				else:
					if (math.floor(hli_Y/col_D) % 2) != 0: 
						if math.floor(hlf_X/col_D) == math.floor(hli_X/col_D) and math.floor(hlmid_Y/col_D) == math.floor(hli_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
				
			else: 
				#print("Border Y from outside")
				if (math.floor(hli_X/col_D) % 2) != 0:
					if (math.floor(hlmid_Y/col_D) % 2) == 0:
						if math.floor(hlf_X/col_D) == math.floor(hli_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hlmid_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
				else:
					if (math.floor(hlmid_Y/col_D) % 2) != 0: 
						if math.floor(hlf_X/col_D) == math.floor(hli_X/col_D) and math.floor(hlf_Y/col_D) == math.floor(hlmid_Y/col_D):
							posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
							if rangemin < posX < rangemax and rangemin < posY < rangemax:
								file_positions.write(str([posX,posY])+",")
								file_AnglesTheta.write(str(theta)+",")
								file_AnglesPhi.write(str(phi)+",")
								hit+=1
				
			#of course, if the neutron is outside of the collimator, nothing can stop it (unlimited power!)
		else:
			posX,posY = (dist_to_col+col_L+dist_to_det)*LD_neutron_X + neutron_position_X,(dist_to_col+col_L+dist_to_det)*LD_neutron_Y + neutron_position_Y
			if rangemin < posX < rangemax and rangemin < posY < rangemax:
				file_positions.write(str([posX,posY])+",")
				file_AnglesTheta.write(str(theta)+",")
				file_AnglesPhi.write(str(phi)+",")
				hit+=1
	#print("Completed")
	file_positions.close()
	file_AnglesTheta.close()
	file_AnglesPhi.close()

	file_metadata.write("\nTime required for execution: "+str(time.time()-timestamp_start)+"s")
	file_metadata.close()

	message = "Everything's fine"
	return message
#Plotting function
def plot3D(x,y,argument,filename):
	title1 = filename+"Hist_"+str(argument)+".txt"
	title2 = filename+"Xax_"+str(argument)+".txt"
	title3 = filename+"Yax"+str(argument)+".txt"
	
	f1 = open(title1,"a")
	f2 = open(title2,"a")
	f3 = open(title3,"a")
	
	hist,xax,yax,image = plt.hist2d(x,y,bins=200,range=[[-10,30],[-10,30]])
	print("Hist dimensions pre: ",len(hist))
	for item in (hist):
		for elem in item:
			f1.write(str(elem))
			f1.write(",")
	f1.close()
	for item in (xax):
		f2.write(str(item))
		f2.write(",")
	f2.close()
	for item in (yax):
		f3.write(str(item))
		f3.write(",")
	f3.close()
	print("file(s) saved")

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument("-x","--LD",type=str,default="27.2",help="Divergence of the initial beam, in LD")
	parser.add_argument("-i","--iterations",type=str,default="200000000",help="Number of iterations")
	parser.add_argument("-c","--channels",type=str,default="20",help="Number of channels")
	parser.add_argument("-D","--col_D",type=str,default="2.5",help="Collimator's channel width")
	parser.add_argument("-L","--col_L",type=str,default="400",help="Collimator's channel length")
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

