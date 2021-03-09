import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt 
import random
import math
import time
import argparse
import os

def plot3D(x,y):
	hist,xax,yax,image = plt.hist2d(x,y,bins=200,range=[[0,20],[0,20]])
	return hist

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

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-x","--Xbias",type=str,default="0",help="X movement")
	parser.add_argument("-y","--Ybias",type=str,default="0",help="Y movement")
	parser.add_argument("-d","--distance",type=str,default="10",help="distance object-detector")
	parser.add_argument("-f","--foldernumber",type=str,default="0",help="Folder number")
	args=parser.parse_args()
	Xbias = float(args.Xbias)
	Ybias = float(args.Ybias)
	distance_to_object = float(args.distance)
	runup = float(args.foldernumber)
	finalposition = []
	for j in range(9):
		Foldername = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/Pavia"+str(j)
		number_of_files = int(len(os.listdir(Foldername)[1:])/3)
		for elem in range(number_of_files):
			countpos,countneg = 0,0
			titleTheta = Foldername+"/Testrun-of-full-sim-AngleTheta-index"+str(elem)+".txt"
			with open(titleTheta, 'r') as file:
				AnTheta = np.array(eval(file.read())) # read list string and convert to array
			file.close()

			titlePhi = Foldername+"/Testrun-of-full-sim-AnglePhi-index"+str(elem)+".txt"
			with open(titlePhi, 'r') as file:
				AnPhi = np.array(eval(file.read())) # read list string and convert to array
			file.close()

			titlePos =  Foldername+"/Testrun-of-full-sim-positionXY-index"+str(elem)+".txt"
			with open(titlePos, 'r') as file:
				Position = np.array(eval(file.read())) # read list string and convert to array
			file.close()

			for i in range(len(AnPhi)):
				impact = calculate_impact_to_object(distance_to_object,128,Position[i,0]-Xbias,Position[i,1]-Ybias,AnTheta[i],AnPhi[i])
				if impact == False:
					finalposition.append(Position[i,:])
					countpos +=1
				else: 
					countneg +=1
			print("Iteration "+str(elem)+" Number of points: "+str(len(finalposition))+" Positive: "+str(countpos)+" Negative "+str(countneg))
	PosXX_after = np.zeros(len(finalposition))
	PosYY_after = np.zeros(len(finalposition))
	for i in range(len(finalposition)):
		PosXX_after[i] = finalposition[i][0] - Xbias
		PosYY_after[i] = finalposition[i][1] - Ybias
	#hist,xax,yax,image=plt.hist2d(PosXX_after,PosYY_after,400,range=[[-5, 15], [0, 20]],cmap="hot")
	hist,xax,yax,image=plt.hist2d(PosXX_after,PosYY_after,1600,range=[[0,40],[0,40]])
	title = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/ResultsPavia"+str(int(runup))+"/X"+str(int(Xbias*10))+"Y"+str(int(Ybias*10))+".txt"
	file2 = open(title,"a")
	for item in hist:
		for point in item:
			file2.write(str(point))
			file2.write(",")
	file2.close()