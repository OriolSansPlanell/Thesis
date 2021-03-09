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
	hist,xax,yax,image = plt.hist2d(x,y,bins=200,range=[[10,30],[10,30]])
	return hist, xax, yax

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-f","--foldernumber",type=str,default="0",help="Folder number")
	parser.add_argument("-n","--filesinfolder",type=str,default="0",help="File number inside folder")
	parser.add_argument("-x","--indexX",type=str,default="100",help="X coordinate")
	parser.add_argument("-y","--indexY",type=str,default="75",help="Y coordinate")
	args=parser.parse_args()
	foldernumber = int(args.foldernumber)
	filesinfolder = int(args.filesinfolder)
	IndexX = int(args.indexX)
	IndexY = int(args.indexY)
	Foldername = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/Cosine"+str(foldernumber)
	folderrange = len(os.listdir(Foldername))
	bins = np.zeros(20)
	for item in range(filesinfolder):
		title = Foldername+"/Testrun-of-full-sim-positionXY-index"+str(item)+".txt"
		print("I'm here")
		with open(title, 'r') as file:
			position = np.array(eval(file.read())) # read list string and convert to array
			file.close()
		hist,xax,yax = plot3D(position[:,0],position[:,1])
		print("Position loaded")
		titletheta = Foldername + "/Testrun-of-full-sim-AngleTheta-index"+str(item)+".txt"
		titlephi = Foldername + "/Testrun-of-full-sim-AnglePhi-index"+str(item)+".txt"

		with open(titletheta,"r") as file:
			angletheta = np.array(eval(file.read()))
			file.close()
		with open(titlephi,"r") as file:
			anglephi = np.array(eval(file.read()))
			file.close()
		print("Angles loaded")
		indices_X = [idx for idx,val in enumerate(position[:,0]) if xax[IndexX]< val <= xax[IndexX+25]]
		indices_Y = [idx for idx,val in enumerate(position[:,1]) if yax[IndexY]< val <= yax[IndexY+25]]
		set_IndicesX = set(indices_X)
		set_IndicesY = set(indices_Y)
		common = []
		if set_IndicesX & set_IndicesY:
			common = set_IndicesX & set_IndicesY
		Distribution_theta = np.zeros(len(common))
		Distribution_phi = np.zeros(len(common))
		i = 0
		for elem in common:
			Distribution_theta[i] = angletheta[elem]
			Distribution_phi[i] = anglephi[elem]
			i+=1

		titletheta = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/ResultsCosine/Theta"+str(foldernumber)+"-"+str(item)+".txt"
		titlephi = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/ResultsCosine/Phi"+str(foldernumber)+"-"+str(item)+".txt"

		ftheta = open(titletheta,"a")
		for item in (Distribution_theta):
			ftheta.write(str(item))
			ftheta.write(",")
		ftheta.close()

		fphi = open(titlephi,"a")
		for item in (Distribution_phi):
			fphi.write(str(item))
			fphi.write(",")
		fphi.close()

