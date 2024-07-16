#this is the code for Nathaniel Wolff's assignment 4 in BME163

#imports
import numpy as np 
import matplotlib.pyplot as plt
import argparse
import matplotlib.patches as mplpatches  

#stylesheet
plt.style.use('BME163') 

#setting up flags for input and output files 
parser = argparse.ArgumentParser()
parser.add_argument('--qualityFile', '-q', default = 'BME163_Input_Data_4.quals', type=str, action = 'store', help ='cell quality file goes here')
parser.add_argument('--coverageFile', '-c', default = 'BME163_Input_Data_4.cov', type = str, action = 'store', help = 'cell coverage file goes here')
parser.add_argument('--outputFile', '-o', default = 'Wolff_Nathaniel_BME163_Assignment_Week4.png', type=str, action = 'store', help ='output file goes here')
parser.add_argument('--identityFile', '-i', default = 'BME163_Input_Data_4.ident', type=str, action = 'store', help ='cell identity file goes here')

args = parser.parse_args()

qualityFile = args.qualityFile
coverageFile = args.coverageFile
identityFile = args.identityFile
outputFile = args.outputFile

#setting up the figure

figureWidth=5
figureHeight=5

plt.style.use('BME163')

plt.figure(figsize=(figureWidth,figureHeight))

panelWidth=3
panelHeight=4

panel1 = plt.axes([.5/figureWidth,.75/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])

#color list, ordered blue, green, yellow, orange
listOfColors = [(44/255,86/255,134/255), (32/255,100/255,113/255), (248/255,174/255,51/255), (230/255,87/255,43/255)]

#swarmplot, written as a function (now horizontal value, vertical shift)
#contains default values from lecture that won't be used unless a flag is set

def swarmplot_function(values = np.random.normal(80,1,100), yPos = 80, ptColor = 'black', linePosCenter = 1):
	#xPos=1


	xmin=0
	xmax=100

	ymin=75
	ymax=100

	width=1

	xrange=xmax-xmin
	yrange=ymax-ymin

	markersize=1

	placedPoints=[]
	minDist=markersize/72
	shift=((minDist/20)*xrange)/panelWidth

	medianOne = np.median(values)
	panel1.vlines(x = medianOne, ymin = yPos-width, ymax = yPos +width, colors = 'red')
	#panel1.axvline(x=medianOne, ymin = 0, ymax = .5, color = 'r') 
	
#plotting loop 

	for x1 in values:
		placed=False
		numberMoves = 0
		if len(placedPoints)==0:
			placedPoints.append((x1,yPos))
		else:
			for move in np.arange(0,width,shift):
				numberMoves +=1  
				newMove = move*(-1)**(numberMoves)
				y1=yPos+newMove
				distList=[]
				for coords2 in placedPoints:
					x2,y2 = coords2[0],coords2[1]
					xdist=(np.abs(x1-x2)/xrange)*panelWidth 
					ydist=(np.abs(y1-y2)/yrange)*panelHeight
					distance=(xdist**2+ydist**2)**0.5
					distList.append(distance)

				if min(distList)>minDist:
					placed=True
					break
			if placed:
				placedPoints.append((x1,y1))
			else:
				print("no more space is available to plot points")
				break
				


	
	for coords in placedPoints:
		
		x,y = coords[0],coords[1]
		
		panel1.plot(x,y,marker='o',mew=0, mfc=ptColor,ms=markersize)


panel1.tick_params(bottom=True, labelbottom=True,\
                   left=True, labelleft=True, \
                   right=False, labelright=False,\
                   top=False, labeltop=False)

panel1.set_xlim(75,100)
panel1.set_ylim(0.5, 12.5)

panel1.set_xticks(np.arange(75,101,5))
panel1.set_yticks([2,5,8,11], ['1-3', '4-6', '7-9', '>10'])

plt.xlabel('Identity (%)')
plt.ylabel('Subread Coverage')


#reading in the data from the two files

identityList = []
subreadCoverageList = []

readName_IdentityDict = {}
readIdentityFile = open(identityFile, 'r')
for line in readIdentityFile:
	line.strip()
	splitLine = line.split("\t")
	readName = splitLine[0]
	fixedIdentity = float(splitLine[1].rstrip())
	identityList.append(fixedIdentity)
	readName_IdentityDict.update({readName:fixedIdentity})
	


readName_SubreadCoverageDict = {}
readCoverageFile = open(coverageFile, 'r')
for line in readCoverageFile:
	line.strip()
	splitLine = (line.split("\t"))
	readName = splitLine[0]
	fixedSC = float(splitLine[1].rstrip())
	subreadCoverageList.append(fixedSC)
	readName_SubreadCoverageDict.update({readName:fixedSC})




#use for XC
#read name 1st column, read quality 2nd column 
readName_ReadQualityDict = {}
readQualityFile = open(qualityFile, 'r')
for line in readQualityFile:
	line.strip()
	splitLine = line.split("\t")
	readName = splitLine[0]
	fixedRQ = splitLine[1].rstrip()
	
	readName_ReadQualityDict.update({ readName: fixedRQ }) 

#consolidating the identities and subread coverages into a single dictionary and sorting it by subread coverage (values)
identity_SubreadCoverageDict = dict(zip(identityList,subreadCoverageList))
orderedIdentity_SubreadCoverageDict = dict(sorted(identity_SubreadCoverageDict.items(), key = lambda kv:kv[1]))

workingID_SBCDict = orderedIdentity_SubreadCoverageDict


#binning the data by subread coverage

allBins = []

counter = 0 
for binUL in range(1,10,3):
	#counter+=1
	upperBin = binUL+2
	lowerBin = binUL
	
	
	binPoints = []
	lastBinPoints = []

	for item in workingID_SBCDict.items():
		counter +=1
		if item[1] in range(lowerBin, upperBin+1, 1):
			binPoints.append(item)
	allBins.append(binPoints)

#doing top bin collection
lastBinPoints = []
for item in workingID_SBCDict.items():
	if item[1] >= 10: 
		lastBinPoints.append(item)
			
allBins.append(lastBinPoints)


counter = 0
for bin in allBins:
	counter +=1
	if counter <= 3:
		pointList = []
		srcList = []
		for individualPoint in bin:
			coord,src = individualPoint
			pointList.append(coord)
			srcList.append(src)
		swarmplot_function(pointList,  min(srcList)+0.5*(max(srcList) - min(srcList)), listOfColors[counter-1], min(srcList)+0.5*(max(srcList) - min(srcList)) )
	elif counter == 4:
		pointList = []
		srcList = []
		for individualPoint in bin:
			coord,src = individualPoint
			pointList.append(coord)
			srcList.append(src)
		swarmplot_function(pointList,  11, listOfColors[counter-1], 11 )




plt.savefig(outputFile, dpi=600)
