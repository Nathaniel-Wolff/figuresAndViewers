#this is the code for Nathaniel Wolff's assignment 5 in BME163
#only use tabs for indented blocks!!

#imports
import numpy as np 
import matplotlib.pyplot as plt
import argparse
import matplotlib.patches as mplpatches  
import matplotlib.image as mplimg
import time

#can't install mappy yet
#import mappy as mp 


#stylesheet
plt.style.use('BME163') 

#setting up the figure and panel
figureWidth = 5
figureHeight = 2

panelHeight = 0.5
panelWidth = 1.5


plt.figure(figsize = (figureWidth, figureHeight))
panel1 = plt.axes([.5/figureWidth, .6/figureHeight, panelWidth/figureWidth, panelHeight/figureHeight])
panel2 = plt.axes([2.5/figureWidth, .6/figureHeight, panelWidth/figureWidth, panelHeight/figureHeight])

panel1.set_xlim(-10, 10)
panel1.set_ylim(0, 2)

panel2.set_xlim(-10, 10)
panel2.set_ylim(0, 2)





listOfPanels = [panel1, panel2]

#panel1.title("Distance to Splice Site")
plt.ylabel("Bits")
plt.title("5'SS")


for panel in listOfPanels:
	panel.set_xticks([-10, -5, 0, 5, 10])
	panel.set_yticks([0, 1, 2])
	

#5' splice site
panel1.tick_params(bottom=True, labelbottom=True,\
                   left=True, labelleft=True, \
                   right=False, labelright=False,\
                   top=False, labeltop=False)
#3' splice site
panel2.tick_params(bottom=True, labelbottom=True,\
                   left=True, labelleft=True, \
                   right=False, labelright=False,\
                   top=False, labeltop=False)


#parsing in the input files from the user
parser = argparse.ArgumentParser()
parser.add_argument('--inputFasta', '-i', default = "Splice_Sequences.fasta", type = str, action = 'store', help = 'input FASTA file goes here')
parser.add_argument('--outputFile', '-o', default = "Wolff_Nathaniel_Assignment_Week5.png", type = str, action = 'store', help = 'output file name and extension goes here')

args = parser.parse_args()

inputFasta = args.inputFasta
outputFile = args.outputFile

#reading and parsing fasta file (standin for mappy), storing each header and sequence in a dictionary 
readFastaFile = open(inputFasta, 'r')


fastaDict5Prime = {}
fastaDict3Prime = {}
counter = 0 
oneSequence = []
headerList = []


for line in readFastaFile:
	if line.startswith('>5'):
		header = line[1:].strip()
		newHeader = str(header)
		#print(type(header))
		if oneSequence: 
			combinedSequence = str(''.join(oneSequence))
			fastaDict5Prime[str(newHeader)] = combinedSequence
			oneSequence = []
	elif line.startswith('>3'):
		header = line[1:].strip()
		newHeader = str(header)
		#print(type(header))
		if oneSequence: 
			combinedSequence = str(''.join(oneSequence))
			fastaDict3Prime[str(newHeader)] = combinedSequence
			oneSequence = []
	else: 
		oneSequence.append(line.strip())
	

if fastaDict5Prime == fastaDict3Prime:
	print("same")

#getting the relative frequencies of each of the letters in the respective positions (Convert to -10,10)
#later I'll make this dynamic?
listOfPositions = []
for number in range(0,20,1):
	listOfPositions.append(number)


#5' SS frequencies 

bases = ["A","C","T","G"]
listOfPositionDicts = [] 
positionDict = {}
 

for position in listOfPositions:
	for base in bases:
		positionDict[base] = 0 
	listOfPositionDicts.append(positionDict)

countPositionDictList = []

#IMPORTANT! base order: A is 0, C is 1, T is 2, G is 3
#NOTE: the current positions are given 0-19. These will have to be converted to -10-10 later.

positionCountsDict = {0:[0,0,0,0],1:[0,0,0,0],2:[0,0,0,0],3:[0,0,0,0],4:[0,0,0,0],5:[0,0,0,0],6:[0,0,0,0],7:[0,0,0,0],8:[0,0,0,0],9:[0,0,0,0],10:[0,0,0,0],11:[0,0,0,0],12:[0,0,0,0],13:[0,0,0,0],14:[0,0,0,0],15:[0,0,0,0],16:[0,0,0,0],17:[0,0,0,0],18:[0,0,0,0],19:[0,0,0,0]}

frequencyBaseDict = {0:[0,0,0,0],1:[0,0,0,0],2:[0,0,0,0],3:[0,0,0,0],4:[0,0,0,0],5:[0,0,0,0],6:[0,0,0,0],7:[0,0,0,0],8:[0,0,0,0],9:[0,0,0,0],10:[0,0,0,0],11:[0,0,0,0],12:[0,0,0,0],13:[0,0,0,0],14:[0,0,0,0],15:[0,0,0,0],16:[0,0,0,0],17:[0,0,0,0],18:[0,0,0,0],19:[0,0,0,0]}



for position in listOfPositions: 	
	for sequence in fastaDict5Prime.values():
		whichBase = sequence[position]
		if whichBase == "A":
			positionCountsDict[position][0]+=1
		elif whichBase == "C":
			positionCountsDict[position][1]+=1
		elif whichBase == "T":
			positionCountsDict[position][2]+=1
		elif whichBase  == "G":
			positionCountsDict[position][3]+=1

	sumOfCounts = sum(positionCountsDict[position])

#calculating relative frequencies of each base in each position (same order) and saving those results to the frequencyBaseDict
	
	frequencyBaseDict = positionCountsDict
	
	for listElement in range(0,4,1):
		frequencyBaseDict[position][listElement] = (frequencyBaseDict[position][listElement]/sumOfCounts)

#calculating uncertainty, H	
#translating the Wikipedia page's parameters: i is position (0-19), b is which base (0-3), and log2 is np.log2	
uncertaintyDict = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0, 18:0, 19:0}	
for position in listOfPositions: 
	uncertainty = 0
	positionFrequencies = frequencyBaseDict[position]
	for baseFrequency in positionFrequencies:
		individualSubH = -baseFrequency * np.log2(baseFrequency)
		uncertainty += individualSubH
	uncertaintyDict[position] = uncertainty

#calculating small sample correction, e
e = (1/.69314718056) * 3/len(frequencyBaseDict)

#calculating Ri and the heights for positions 0-19


heightDict = {}
for position in listOfPositions:
	RPos =  np.log2(20) - (uncertaintyDict[position]-e)	
	orderedHeightList = []
	for individualFrequency in frequencyBaseDict[position]:
		individualHeight = individualFrequency * RPos
		orderedHeightList.append(individualHeight)
	
	heightDict[position] = orderedHeightList

#reading in the base images

imgDict = {'A': mplimg.imread('A.png'), "C": mplimg.imread('C.png'), "T": mplimg.imread('T.png'), "G" : mplimg.imread('G.png')}



#defining a plotting loop for the images
for position in listOfPositions:
	onePlotDict = {'A':heightDict[position][0], 'C':heightDict[position][1], 'T':heightDict[position][2], 'G':heightDict[position][3]}
	#print("olddict", onePlotDict)
	sortedOnePlotDict = dict(sorted(onePlotDict.items(), key = lambda kv:kv[-1]))
	#print("sort dict", sortedOnePlotDict)
	plottedHeightShift = 0 
	for baseHeight in sortedOnePlotDict.items():
		whichBase,height = baseHeight
		#print(height)
		panel1.imshow(imgDict[whichBase], extent = [position-10, position-9, (plottedHeightShift)/2, (plottedHeightShift + height)/2], aspect = 'auto', origin = 'upper')
		plottedHeightShift+=height


#3' SS		
	
bases = ["A","C","T","G"]
listOfPositionDicts = [] 
positionDict = {}
 

for position in listOfPositions:
	for base in bases:
		positionDict[base] = 0 
	listOfPositionDicts.append(positionDict)

countPositionDictList = []

#IMPORTANT! base order: A is 0, C is 1, T is 2, G is 3
#NOTE: the current positions are given 0-19. These will have to be converted to -10-10 later.

positionCountsDict = {0:[0,0,0,0],1:[0,0,0,0],2:[0,0,0,0],3:[0,0,0,0],4:[0,0,0,0],5:[0,0,0,0],6:[0,0,0,0],7:[0,0,0,0],8:[0,0,0,0],9:[0,0,0,0],10:[0,0,0,0],11:[0,0,0,0],12:[0,0,0,0],13:[0,0,0,0],14:[0,0,0,0],15:[0,0,0,0],16:[0,0,0,0],17:[0,0,0,0],18:[0,0,0,0],19:[0,0,0,0]}

frequencyBaseDict = {0:[0,0,0,0],1:[0,0,0,0],2:[0,0,0,0],3:[0,0,0,0],4:[0,0,0,0],5:[0,0,0,0],6:[0,0,0,0],7:[0,0,0,0],8:[0,0,0,0],9:[0,0,0,0],10:[0,0,0,0],11:[0,0,0,0],12:[0,0,0,0],13:[0,0,0,0],14:[0,0,0,0],15:[0,0,0,0],16:[0,0,0,0],17:[0,0,0,0],18:[0,0,0,0],19:[0,0,0,0]}



for position in listOfPositions: 	
	for sequence in fastaDict3Prime.values():
		whichBase = sequence[position]
		if whichBase == "A":
			positionCountsDict[position][0]+=1
		elif whichBase == "C":
			positionCountsDict[position][1]+=1
		elif whichBase == "T":
			positionCountsDict[position][2]+=1
		elif whichBase  == "G":
			positionCountsDict[position][3]+=1

	sumOfCounts = sum(positionCountsDict[position])

#calculating relative frequencies of each base in each position (same order) and saving those results to the frequencyBaseDict
	
	frequencyBaseDict = positionCountsDict
	
	for listElement in range(0,4,1):
		frequencyBaseDict[position][listElement] = (frequencyBaseDict[position][listElement]/sumOfCounts)

#calculating uncertainty, H	
#translating the Wikipedia page's parameters: i is position (0-19), b is which base (0-3), and log2 is np.log2	
uncertaintyDict = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0, 18:0, 19:0}	
for position in listOfPositions: 
	uncertainty = 0
	positionFrequencies = frequencyBaseDict[position]
	for baseFrequency in positionFrequencies:
		individualSubH = -baseFrequency * np.log2(baseFrequency)
		uncertainty += individualSubH
	uncertaintyDict[position] = uncertainty

#calculating small sample correction, e
e = (1/.69314718056) * 3/len(frequencyBaseDict)

#calculating Ri and the heights for positions 0-19


heightDict = {}
for position in listOfPositions:
	RPos =  np.log2(20) - (uncertaintyDict[position]-e)	
	orderedHeightList = []
	for individualFrequency in frequencyBaseDict[position]:
		individualHeight = individualFrequency * RPos
		orderedHeightList.append(individualHeight)
	
	heightDict[position] = orderedHeightList

#reading in the base images

imgDict = {'A': mplimg.imread('A.png'), "C": mplimg.imread('C.png'), "T": mplimg.imread('T.png'), "G" : mplimg.imread('G.png')}



#defining a plotting loop for the images
for position in listOfPositions:
	onePlotDict = {'A':heightDict[position][0], 'C':heightDict[position][1], 'T':heightDict[position][2], 'G':heightDict[position][3]}
	#print("olddict", onePlotDict)
	sortedOnePlotDict = dict(sorted(onePlotDict.items(), key = lambda kv:kv[-1]))
	#print("sort dict", sortedOnePlotDict)
	plottedHeightShift = 0 
	for baseHeight in sortedOnePlotDict.items():
		whichBase,height = baseHeight
		#print(height)
		panel2.imshow(imgDict[whichBase], extent = [position-10, position-9, (plottedHeightShift)/2, (plottedHeightShift + height)/2], aspect = 'auto', origin = 'upper')
		plottedHeightShift+=height




#exec('positionDict{:s}[{:s}] = 0'.format(position, base))
plt.savefig(outputFile, dpi=600)

print("program complete")