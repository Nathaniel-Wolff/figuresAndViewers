 #this is the code for Nathaniel Wolff's final assignment in BME 163

#imports
import numpy as np 
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict
import matplotlib.patches as mplpatches  
import matplotlib.image as mplimg
import time

#stylesheet
plt.style.use('BME163') 


#parsing in the input files from the user
parser = argparse.ArgumentParser()
parser.add_argument('--inputPSL', '-i', default = "BME163_Input_Data_6.psl", type = str, action = 'store', help = 'input PSL file goes here')
parser.add_argument('--inputGTF', '-g', default = "gencode.vM12.annotation.gtf", type = str, action = 'store', help = 'input GTF goes here')
parser.add_argument('--coords', '-c', default = "chr7:45232000-45241000", type = str, action = 'store', help = 'coordinates for chosen gene locus goes here') 

parser.add_argument('--outputFile', '-o', default = "FinalAssignment.png", type = str, action = 'store', help = 'output file name and extension goes here')

args = parser.parse_args()

inputPSL = args.inputPSL
inputGTF = args.inputGTF
coords = args.coords
outputFile = args.outputFile


#separating the coords into three variables, specifying the chromosome and span within it

coords = coords.split("-")
coords.extend(coords[0].split(":"))
coords.pop(0)



		
startSpan = coords[2]
chr = coords[1]
endSpan = coords[0]
	

#setting up the panels - complete 
figureWidth = 5
figureHeight = 5

panelHeight = 1.5
panelWidth = 4


plt.figure(figsize = (figureWidth, figureHeight))
panelBottom = plt.axes([.1/figureWidth, .1/figureHeight, panelWidth/figureWidth, panelHeight/figureHeight])
panelMiddle = plt.axes([.1/figureWidth, 1.7/figureHeight, panelWidth/figureWidth, panelHeight/figureHeight])
panelTop = plt.axes([.1/figureWidth, 3.3/figureHeight, panelWidth/figureWidth, panelHeight/figureHeight])

listOfPanels = [panelTop, panelMiddle, panelBottom]

#setting the x limits dynamically as a function of user input
#fix the y limits later!!


for panel in listOfPanels:
	panel.set_xlim(int(startSpan)-100, int(endSpan)+100)
	panel.tick_params(bottom=False, labelbottom=False,
                   left=False, labelleft=False, 
                   right=False, labelright=False,
                   top=False, labeltop=False)

#parsing in input data to prepare it for plotting 


#parsing in gtf file 
gtfPairDict = defaultdict(list)
with open(inputGTF) as inputGTF:
	for line in inputGTF:
		line = line.split('\t')
		if line[0] == chr:
			
			if line[2] in ['CDS', 'exon']:
				start = line[3]
				end = line[4]
				exonOrCDS = line[2]
				readList = [start, end, exonOrCDS]
				transcriptID = line[8].split(' transcript_id "')[1].split('"')[0]
				gtfPairDict[transcriptID].append(readList)

	

	#computing some bullshit


	listOfReads = []
	#list of reads will contain a read for each transcript, represented by a list of the following ordered elements: [minStart,maxEnd, false,  blockStarts, blockEnds, blockTypes]
	#CAUTION: blockEnds is NOT blockwidths, as is given for the psl file. you will need to get the actual width of the rectangle from the difference of start and end
	for transcript,reads in gtfPairDict.items():
		minStart = 0
		maxEnd = 0
		
		blockStarts = []
		blockEnds = []
		blockWidths = []
		blockTypes = []

		
		for read in reads:
			blockStarts.append(int(read[0]))
			blockEnds.append(int(read[1]))
			blockWidths.append(int(read[1]) - int(read[0]))
			
			blockTypes.append(read[2])
		minStart = int(min(blockStarts))
		maxEnd = int(max(blockEnds))
		
		if int(startSpan)<int(minStart)<int(endSpan) or int(startSpan)<int(maxEnd)<int(endSpan):
			readList = [minStart, maxEnd, False, blockStarts, blockWidths, blockTypes]
			listOfReads.append(readList)

	
listOfReads.sort(key = lambda kv:kv[0])

counter = 0


yPos = 1
baseList = []



for yPos in range(0,len(listOfReads), 1):
	plottedReadEnd = 0 
	for read in listOfReads:
		
		startAlign = read[0]
		endAlign = read[1]
		wasPlotted = read[2]
		blockStarts = read[3]
		blockWidths = read[4]
		whichType = read[5]

		if not wasPlotted:

			if plottedReadEnd<startAlign:
			
				sameLineRectangle = mplpatches.Rectangle([startAlign,yPos-.05], endAlign-startAlign, .1, facecolor = 'grey', edgecolor = 'black', linewidth = 0.25)
				panelTop.add_patch(sameLineRectangle)
				plottedReadEnd = endAlign
				finalY = yPos
				heightSameLineRect = .5
				for index in range(0,len(blockStarts),1):
					typeBlock = whichType[index]
					if typeBlock == 'CDS':
						blockRectangle = mplpatches.Rectangle([blockStarts[index],yPos-.25], blockWidths[index], .5, facecolor = 'grey', edgecolor = 'black', linewidth = 0.25)
						panelTop.add_patch(blockRectangle)
					elif typeBlock == 'exon':
						blockRectangle = mplpatches.Rectangle([blockStarts[index],yPos-.125], blockWidths[index], .25, facecolor = 'grey', edgecolor = 'black', linewidth = 0.25)
						panelTop.add_patch(blockRectangle)		
				read[2]=True
			
	
				


	
panelTop.set_ylim(-1.25,finalY*1.1)
panelTop.tick_params(bottom=False, labelbottom=False,
                   left=False, labelleft=False, 
                   right=False, labelright=False,
                   top=False, labeltop=False)









#psl file parsed below
#Note that here, a block = exon.
#each line refers to a single read already aligned to a certain chromosome, positions specified


numberOfChrReads = 0 
listReads = []
histList = []



#dictionary links the beginning of a read to the end of an read as a key value pair
listOfReads = []
listOfBlockStarts = []
listOfBlockWidths = []
with open(inputPSL, 'r') as f:
	for line in f:
		splitLine = line.split('\t')
		if splitLine[13] == chr:
			numberOfChrReads += 1  
			line = line.split('\t')
			yStart = 0 
			startAlign = int(line[15])
			endAlign = int(line[16]) 
			numberExons = int(line[17])
			blockWidths = np.array(line[18].split(',')[:-1], dtype = int)
			blockStarts = np.array(line[20].split(',')[:-1], dtype = int)
			
			height = .1

			#picking the proper reads for the given chr/coordinates and storing them in a dictionary 
			#make sure to use the same range for the endAlign



			#if startAlign in range(int(startSpan)+1, int(endSpan)+1) and endAlign<int(endSpan):
				#readStrEndDict[startAlign] = endAlign


			if int(startSpan)<int(startAlign)<int(endSpan) or int(startSpan)<int(endAlign)<int(endSpan):
				read = [startAlign,endAlign,False,blockStarts,blockWidths]
				listOfBlockStarts.append(blockStarts)
				listOfBlockWidths.append(blockWidths)
				listOfReads.append(read)
	
	#starting the stacking loop: comparing the new read to be plotted (reminder: that read is defined by start and end align) to the read that has previously been plotted and the last x limit

listOfReads.sort(key = lambda kv:kv[0])

counter = 0

iBlue = (88/255,85/255,120/255)
	

yPos = 1

for yPos in range(0,len(listOfReads), 1):
	plottedReadEnd = 0 
	for read in listOfReads:
		
		startAlign = read[0]
		endAlign = read[1]
		wasPlotted = read[2]
		blockStarts = read[3]
		blockWidths = read[4]


		if not wasPlotted:

			if plottedReadEnd<startAlign:
			
				sameLineRectangle = mplpatches.Rectangle([startAlign,yPos-.05], endAlign-startAlign, .1, facecolor = iBlue, edgecolor = 'black', linewidth = 0.0)
				panelMiddle.add_patch(sameLineRectangle)
				plottedReadEnd = endAlign
				finalY = yPos
				heightSameLineRect = .5
				for index in range(0,len(blockStarts),1):
					blockRectangle = mplpatches.Rectangle([blockStarts[index],yPos-.25], blockWidths[index], heightSameLineRect, facecolor = iBlue, edgecolor = 'black', linewidth = 0.0)
					for base in range(blockStarts[index], blockStarts[index]+blockWidths[index], 1):
						histList.append(base)
							
					panelMiddle.add_patch(blockRectangle)	
				read[2]=True

bins = range(int(startSpan),int(endSpan), 1)
xHistogram,Bins = np.histogram(histList, bins)


for i in range(0,len(xHistogram), 1):
	left=Bins[i]
	bottom = 0 
	width  = Bins[i+1]-left
	height = xHistogram[i]
	
	rectangle = mplpatches.Rectangle([left,bottom], width, height,
                                         facecolor = iBlue, 
                                         edgecolor=iBlue, 
                                         linewidth=0.25
                                         )

	panelBottom.add_patch(rectangle)
	

		


	#panel.set_xlim(int(startSpan)-100, int(endSpan)+100)
panelMiddle.set_ylim(-.5,finalY*1.1)
panelMiddle.tick_params(bottom=False, labelbottom=False,
                   left=False, labelleft=False, 
                   right=False, labelright=False,
                   top=False, labeltop=False)
panelBottom.set_ylim(-.5,finalY*1.1)
panelBottom.tick_params(bottom=False, labelbottom=False,
                   left=False, labelleft=False, 
                   right=False, labelright=False,
                   top=False, labeltop=False)









#selection of reads to display from argparse input 





plt.savefig(outputFile, dpi=1200)
print("program complete")