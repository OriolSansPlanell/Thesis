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

if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-x","--foldernumber",type=str,default="0",help="Folder number")
	args=parser.parse_args()
	foldernumber = int(args.foldernumber)
	Foldername = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/Pavia"+str(foldernumber)
	folderrange = len(os.listdir(Foldername))
	finalhist = np.zeros([200,200])
	bins = np.zeros(20)
	for fileindex in range(int(folderrange/3)):
		title = Foldername+"/Testrun-of-full-sim-positionXY-index"+str(fileindex)+".txt"
		i,run = 0,0
		with open(title, 'r') as file:
			position = np.array(eval(file.read())) # read list string and convert to array
			file.close()
		hist = plot3D(position[:,0],position[:,1])
		finalhist += hist
	title = "/eos/home-o/osanspla/SWAN_projects/Neutron Collimator Geometry/ResultsPavia/"+str(foldernumber)+".txt"
	file2 = open(title,"a")
	for item in finalhist:
		for elem in item:
			file2.write(str(elem))
			file2.write(",")
	file2.close()