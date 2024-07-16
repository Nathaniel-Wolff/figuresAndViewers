#This is the code for Nathaniel Wolff's assignment 2 in BME163



import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches 
import numpy as np
import argparse 


#stylesheet
plt.style.use('BME163') 

#setting up iflag for input
parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', '-i', default = 'BME163_Input_Data_1.txt', type=str, action = 'store', help ='input file goes here')
parser.add_argument('--outputFile', '-o', default = 'Assignment2.png', type=str, action = 'store', help ='output file goes here')

#args refers to 
args = parser.parse_args()

inputFile = args.inputFile
outputFile = args.outputFile

figureHeight = 2
figureWidth = 5

mainPanelHeight1 = 1
mainPanelWidth1 = 1

topPanelHeight = .25
topPanelWidth = 1

sidePanelHeight = 1
sidePanelWidth = .25

panel2Height = 1
panel2Width = 1

plt.figure(figsize = (figureWidth,figureHeight))


#setting up the panels and placing the proper positions
firstMainPanel = plt.axes([.7/figureWidth,.375/figureHeight,mainPanelWidth1/figureWidth,mainPanelHeight1/
figureHeight])
firstSidePanel = plt.axes([.375/figureWidth, .375/figureHeight, sidePanelWidth/figureWidth, sidePanelHeight/figureHeight])
firstTopPanel = plt.axes([.7/figureWidth, 1.45/figureHeight, topPanelWidth/figureWidth, topPanelHeight/figureHeight])

listOfFirstPanels = [firstMainPanel,firstSidePanel,firstTopPanel]

#setting ticks of first main panel
firstMainPanel.tick_params(bottom=True, labelbottom=True,
left=False, labelleft=False, 
right=False, labelright=False,
top=False, labeltop=False)
firstMainPanel.set_xticks(np.arange(0,16,5))

#seting ticks of first top panel
firstTopPanel.tick_params(bottom=False, labelbottom=False,
left=True, labelleft=True, 
right=False, labelright=False,
top=False, labeltop=False)

#setting ticks of the first side panel
firstSidePanel.tick_params(bottom=True, labelbottom=True,
left=True, labelleft=True, 
right=False, labelright=False,
top=False, labeltop=False)

#setting the x and y limits of the first top, main, and side panels

firstMainPanel.set_xlim(0,15)
firstMainPanel.set_ylim(0,15)

firstTopPanel.set_xlim(0,15)
firstTopPanel.set_ylim(0,20)

firstSidePanel.set_xlim(20,0)
firstSidePanel.set_ylim(0,15)

#extra credit portion, same as above
secondMainPanel = plt.axes([2.5/figureWidth,.375/figureHeight,mainPanelWidth1/figureWidth,mainPanelHeight1/
figureHeight])
secondSidePanel = plt.axes([2.175/figureWidth, .375/figureHeight, sidePanelWidth/figureWidth, sidePanelHeight/figureHeight])
secondTopPanel = plt.axes([2.5/figureWidth, 1.45/figureHeight, topPanelWidth/figureWidth, topPanelHeight/figureHeight])


#setting ticks of second main panel
secondMainPanel.tick_params(bottom=True, labelbottom=True,
left=False, labelleft=False, 
right=False, labelright=False,
top=False, labeltop=False)
firstMainPanel.set_xticks(np.arange(0,16,5))

#seting ticks of first top panel
secondTopPanel.tick_params(bottom=False, labelbottom=False,
left=True, labelleft=True, 
right=False, labelright=False,
top=False, labeltop=False)

#setting ticks of the first side panel
secondSidePanel.tick_params(bottom=True, labelbottom=True,
left=True, labelleft=True, 
right=False, labelright=False,
top=False, labeltop=False)

secondMainPanel.set_xlim(0,15)
secondMainPanel.set_ylim(0,15)

secondTopPanel.set_xlim(0,15)
secondTopPanel.set_ylim(0,20)

secondSidePanel.set_xlim(20,0)
secondSidePanel.set_ylim(0,15)



#parsing the input file

inputFileHandle=open(inputFile,'r')

xValues = []
yValues = []

#constructing x and y values from the parsed data
for line in inputFileHandle:
	xValues.append(np.log2(float(line.rstrip().split('\t')[1])+1))
	yValues.append(np.log2(float(line.rstrip().split('\t')[2])+1))
inputFileHandle.close()


#plotting the points properly

iBlue=(88/255,85/255,120/255)
iYellow=(248/255,174/255,51/255)
iGreen=(120/255,172/255,145/255)

firstMainPanel.plot(xValues, yValues, marker ='o', markerfacecolor = iYellow, markersize = 1.45, markeredgewidth = 0, linewidth = 0, alpha = 0.1)


#forming histograms for the x and y 

bins = []
for number in range(0,len(xValues)):
	bins.append(number*.5)
	
	
xHistogram,Bins = np.histogram(xValues, bins)
yHistogram,Bins = np.histogram(yValues, bins)

for i in range(0,len(xHistogram), 1):
	left=bins[i]
	bottom = 0 
	width  = bins[i+1]-left
	height = np.log2(xHistogram[i]+1)
	
	rectangle = mplpatches.Rectangle([left,bottom], width, height,
                                         facecolor=(120/255,172/255,145/255), 
                                         edgecolor='black', 
                                         linewidth=0.25
                                         )

	firstTopPanel.add_patch(rectangle)
	

for i in range(0,len(xHistogram), 1):
	left=bins[i]
	bottom = 0 
	width  = bins[i+1]-left
	height = np.log2(xHistogram[i]+1)
	
	rectangle = mplpatches.Rectangle([left,bottom], width, height,
                                         facecolor=(120/255,172/255,145/255), 
                                         edgecolor='black', 
                                         linewidth=0.25
                                         )

	secondTopPanel.add_patch(rectangle)




for i in range(0,len(yHistogram),1):
	left=0
	bottom = bins[i] 
	width  = np.log2(yHistogram[i]+1) 
	height = bins[i+1]-bottom
	rectangle = mplpatches.Rectangle([left,bottom], width, height,
                                         facecolor=(88/255,85/255,120/255), 
                                         edgecolor='black', 
                                         linewidth=0.25
                                         )

	firstSidePanel.add_patch(rectangle)	

for i in range(0,len(yHistogram),1):
	left=0
	bottom = bins[i] 
	width  = np.log2(yHistogram[i]+1) 
	height = bins[i+1]-bottom
	rectangle = mplpatches.Rectangle([left,bottom], width, height,
                                         facecolor=(88/255,85/255,120/255), 
                                         edgecolor='black', 
                                         linewidth=0.25
                                         )

	secondSidePanel.add_patch(rectangle)

#setting up the legend panel
legendRight = plt.axes([4/figureWidth, .375/figureHeight, .1/figureWidth, sidePanelHeight/figureHeight])
#seting ticks of legend
legendRight.tick_params(bottom=False, labelbottom=False,
left=True, labelleft=True, 
right=False, labelright=False,
top=False, labeltop=False)

legendRight.set_xlim(0,.1)
legendRight.set_ylim(0,20)
legendRight.set_yticks([0,20],['0','>20'])

#color map tuple pairs, viridis values

viridisValues1 = [(68/255, 1/255, 84/255), (59/255, 82/255, 139/255)]
viridisValues2 = [(59/255, 82/255, 139/255), (33/255, 145/255, 140/255)]
viridisValues3 = [(33/255, 145/255, 140/255),(94/255, 201/255, 98/255)]
viridisValues4 = [(94/255, 201/255, 98/255), (253/255, 231/255, 37/255)]

#color map tuple pair linspaces, viridis values
vvLin1Red = np.linspace(68/255, 59/255, 5)
vvLin2Red = np.linspace(59/255, 33/255, 6)
vvLin3Red = np.linspace(33/255, 94/255, 6)
vvLin4Red = np.linspace(94/255, 253/255, 6)


vvLin1Green = np.linspace(1/255, 82/255, 5)
vvLin2Green = np.linspace(82/255, 145/255, 6)
vvLin3Green = np.linspace(145/255, 201/255, 6)
vvLin4Green = np.linspace(201/255, 231/255, 6)

vvLin1Blue = np.linspace(84/255, 139/255, 5)
vvLin2Blue = np.linspace(139/255, 140/255, 6)
vvLin3Blue = np.linspace(140/255, 98/255, 6)
vvLin4Blue = np.linspace(98/255, 37/255, 6)


plLin4Red = np.linspace(245/255, 237/255, 5)
plLin3Red = np.linspace(190/255, 245/255, 6)
plLin2Red = np.linspace(87/255, 190/255, 6)
plLin1Red = np.linspace(15/255, 87/255, 6)

plLin4Green = np.linspace(135/255,252/255, 5)
plLin3Green = np.linspace(48/255, 135/255, 6)
plLin2Green = np.linspace(0/255, 48/255, 6) 
plLin1Green = np.linspace(0/255, 0/255, 6)

plLin4Blue = np.linspace(48/255, 27/255, 5)
plLin3Blue = np.linspace(101/255, 48/255,  6)
plLin2Blue = np.linspace(151/255, 101/255, 6)
plLin1Blue = np.linspace(118/255, 151/255, 6)


#total linspaces for all tuple pairs, viridis values
vvListOfRedLins = list(vvLin1Red)+list(vvLin2Red)+list(vvLin3Red)+list(vvLin4Red)
vvListOfGreenLins = list(vvLin1Green)+list(vvLin2Green)+list(vvLin3Green)+list(vvLin4Green)
vvListOfBlueLins = list(vvLin1Blue)+list(vvLin2Blue)+list(vvLin3Blue)+list(vvLin4Blue)

orderedVVRed = list(dict.fromkeys(vvListOfRedLins))
orderedVVGreen = list(dict.fromkeys(vvListOfGreenLins))
orderedVVBlue = list(dict.fromkeys(vvListOfBlueLins))


#viridis heatmaps into the legend panel

for index in range(0,20,1):
	colorPaletteVV = (orderedVVRed[index], orderedVVGreen[index], orderedVVBlue[index])
	vvGradeRect = mplpatches.Rectangle([0,index], .1, 8, facecolor=colorPaletteVV,edgecolor = 'black',linewidth = 0)
		
	legendRight.add_patch(vvGradeRect)

secondMainPanel.plot(xValues, yValues, marker ='o', markerfacecolor = iYellow, markersize = 1.45, markeredgewidth = 0, linewidth=0,alpha=0.1)

#listOfBins = []
#for i in np.arange(0,15,.25): 
	#listOfBins.append(i)
	

#establishing a histogram of bin size .25 for the second main panel	
#xOtherHistogram,Bins = np.histogram(xValues, listOfBins, density = False)
#yOtherHistogram,Bins = np.histogram(yValues, listOfBins, density = False)



#for xIndividualBlock in np.clip(xOtherHistogram, 0, 20): 
	#for yIndividualBlock in yOtherHistogram:
		#clippedIndex=np.clip(yIndividualBlock, 0, 20)
		#print(clippedIndex)
		#secondMainPanel.plot(xValues, yValues, marker ='o', markerfacecolor = (orderedVVRed[clippedIndex], #orderedVVGreen[clippedIndex], orderedVVBlue[clippedIndex]), markersize = 1.45, markeredgewidth = 0, linewidth = 0, alpha = 0.1)



plt.savefig(outputFile, dpi=600)
print("complete")