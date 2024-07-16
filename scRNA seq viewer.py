#This is the code for Nathaniel Wolff's assignment 3 in BME163



import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.patches as mplpatches 
import numpy as np
import argparse 


#stylesheet
plt.style.use('BME163') 

#setting up flags for cell type, cell position, and output files 

parser = argparse.ArgumentParser()

parser.add_argument('--cellTypeFile', '-c', default = 'BME163_Input_Data_Week3.celltype.tsv', type=str, action = 'store', help ='cell type csv goes here')

parser.add_argument('--positionFile', '-p', default = 'BME163_Input_Data_Week3.position.tsv', type = str, action = 'store', help = 'cell position csv goes here')

parser.add_argument('--outputFile', '-o', default = 'Wolff_Nathaniel_BME163_Assignment_Week3.png', type=str, action = 'store', help ='output file goes here')

parser.add_argument('--expressionFile', '-e', default = 'BME163_Input_Data_Week3.expression.csv', type=str, action = 'store', help ='expression file goes here')

parser.add_argument('--gene', '-g', default = 'CD7', type=str, action = 'store', help ='gene of interest goes here')


args = parser.parse_args()

cellTypeFile = args.cellTypeFile
positionFile = args.positionFile
outputFile = args.outputFile
expressionFile = args.expressionFile
gene = args.gene

#setting up and plotting figure/panels
figureHeight = 4
figureWidth = 8

mainPanelHeight1 = 2
mainPanelWidth1 = 2

panelHeight = .25
panelWidth = 1

plt.figure(figsize = (figureWidth,figureHeight))

legendPanelWidth = .2
legendPanelHeight = .35



legendPanel = plt.axes([5.1/figureWidth, 1.5/figureHeight, legendPanelWidth/figureWidth, legendPanelHeight/figureHeight])



legendPanel.tick_params(bottom=False, labelbottom=False,
left=False, labelleft=False, 
right=False, labelright=True,
top=False, labeltop=False)



panel2 = plt.axes([3/figureWidth, .5/figureHeight,mainPanelWidth1/figureWidth,mainPanelHeight1/
figureHeight])

panel2.set_xlim(left=-30, right=30)
panel2.set_ylim(bottom=-40, top=30)

panel1 = plt.axes([.5/figureWidth,.5/figureHeight,mainPanelWidth1/figureWidth,mainPanelHeight1/
figureHeight])

panel1.set_xlim(left=-30, right=30)
panel1.set_ylim(bottom = -40, top = 30)

#setting the ticks of the first panel
panel1.tick_params(bottom=True, labelbottom=True,
left=True, labelleft=True, 
right=False, labelright=False,
top=False, labeltop=False)

panel1.set_xticks(np.arange(-30,32,10))
panel1.set_yticks(np.arange(-40,31,10))
panel1.set_xlabel("tSNE 2")
panel1.set_ylabel("tSNE 1")


panel2.set_xticks(range(-30,31,10))
panel2.set_yticks(range(-40,31,10))
panel2.set_xlabel("tSNE 2")
panel2.set_ylabel("tSNE 1")

panel2.tick_params(bottom=True, labelbottom=True,
left=True, labelleft=True, 
right=False, labelright=False,
top=False, labeltop=False)



#opening the position file and parsing its input into a dictionary called barcode_XYDict

readPositionFile=open(positionFile,'r')

barcode_XYDict = {}
for line in readPositionFile:
	line.strip()
	splitLine = line.split()
	barcode_XYDict.update({ splitLine[0]: (splitLine[1],splitLine[2]) }) 

#opening the cell type file and parsing its input into another dictionary called barcode_cellTypeDict

readCellTypeFile = open(cellTypeFile, 'r')

barcode_cellTypeDict = {}
for line in readCellTypeFile:
	line.strip()
	splitLine = line.split()
	barcode_cellTypeDict.update({splitLine[2]: splitLine[1] })


#setting up the colors for each cell type
colorDict = {"tCell":(0,128/255,128/255), "monocyte":(128/255,128/255,128/255), "bCell":(160/255, 32/255, 240/255)}

#plotting all of the points with their respective colors 

howmanypoints = 0
for barcode in barcode_XYDict.keys():
	x,y = barcode_XYDict[barcode]
	theColor = colorDict[barcode_cellTypeDict[barcode]]
	
	plt.plot(float(x),float(y), marker ='o', markerfacecolor = theColor, markersize = 4.1, markeredgewidth = 0.1, markeredgecolor = 'black', linewidth = 0, alpha = 1)
	howmanypoints +=1


#consolidating each x,y into its respective cluster
#each barcode corresponds to an x,y and each barcode corresponds to a cell type

xValuesmonocyte = []
xValuestCell = []
xValuesbCell = []
yValuesmonocyte = []
yValuestCell = []
yValuesbCell = []


for barcode in barcode_XYDict.keys():
	oneCoord = barcode_XYDict[barcode]
	whichCellType = barcode_cellTypeDict[barcode]
	exec('xValues{:s}.append({:s})'.format(whichCellType, barcode_XYDict[barcode][0]))
	exec('yValues{:s}.append({:s})'.format(whichCellType, barcode_XYDict[barcode][1]))

#establishing the median centers of each cluster

centermonocyte = (np.median(xValuesmonocyte), np.median(yValuesmonocyte))
centertCell = (np.median(xValuestCell), np.median(yValuestCell))
centerbCell = (np.median(xValuesbCell), np.median(yValuesbCell))	

centersDict = {'centermonocyte':centermonocyte, 'centertCell': centertCell, 'centerbCell': centerbCell}


#plotting the center text

for name in ['centermonocyte', 'centertCell', 'centerbCell']:
	panel1.annotate(name[6:], centersDict[name], xycoords = 'data', path_effects = [pe.withStroke(linewidth=1, foreground = "white")], horizontalalignment='center', verticalalignment='center')



#This is the code for the extra credit portion of Nathaniel Wolff's assignment 3 in BME163

#parsing the csv file 

expressionCSV = open(expressionFile, "r")


listOfBarcodes = []
listOfExpressionLevels = []


iterator = 0 
listLine = []
listOfExpressionLevels = []

for line in expressionCSV:
	if iterator == 0:
		listOfBarcodes = line.strip().replace('"','').split(",")
		iterator+=1
		continue
	rowList = line.strip().replace('"','').split(",")
	iterator+=1
	if iterator != 0:
		if rowList[0] == gene:
			listOfExpressionLevels = line.strip().replace('"','').split(",")[0:]
		
			
			

#making a dictionary linking barcodes and expression levels, given the user's gene input
barcode_ExpressionLevelDict = dict(zip(listOfBarcodes, listOfExpressionLevels)) 


#making a max and min expression level variable


maxExpression = float(max(listOfExpressionLevels[2:]))
minExpression = float(min(listOfExpressionLevels[2:]))
expressionRange = round(maxExpression-minExpression, 2)


legendPanel.set_xlim(0,.1)
legendPanel.set_ylim(0,20)
legendPanel.set_yticks([0,20],['0',str(round(maxExpression, 2))])	
		
#plasmaheatmap stuff from Chris		
plasma5 = (237/255, 252/255, 27/255)
plasma4 = (245/255, 135/255, 48/255)
plasma3 = (190/255, 48/255, 101/255)
plasma2 = (87/255, 0/255, 151/255)
plasma1 = (15/255, 0/255, 118/255)

R1=np.linspace(plasma1[0],plasma2[0],25)
G1=np.linspace(plasma1[1],plasma2[1],25)
B1=np.linspace(plasma1[2],plasma2[2],25)

R2=np.linspace(plasma2[0],plasma3[0],25)
G2=np.linspace(plasma2[1],plasma3[1],25)
B2=np.linspace(plasma2[2],plasma3[2],25)

R3=np.linspace(plasma3[0],plasma4[0],25)
G3=np.linspace(plasma3[1],plasma4[1],25)
B3=np.linspace(plasma3[2],plasma4[2],25)

R4=np.linspace(plasma4[0],plasma5[0],26)
G4=np.linspace(plasma4[1],plasma5[1],26)
B4=np.linspace(plasma4[2],plasma5[2],26)

R=np.concatenate((R1,R2,R3,R4),axis=None)
G=np.concatenate((G1,G2,G3,G4),axis=None)
B=np.concatenate((B1,B2,B3,B4),axis=None)

#plotting the cells with the proper plasma facecolors

howmanypointsright = 0

color_XYDict = {}
 

for barcode in barcode_XYDict.keys():
	a = barcode_XYDict[barcode]
	expressionLevelBarcodeNM = 100*(float(barcode_ExpressionLevelDict[barcode])/expressionRange)
	roundedELB = round(expressionLevelBarcodeNM)
	plasmaColorTuple = (R[roundedELB], G[roundedELB], B[roundedELB])
	
	color_XYDict.update({ a: plasmaColorTuple}) 





orderedColor_XYDict = dict(sorted(color_XYDict.items(), key = lambda kv:kv[1]))


for pair in orderedColor_XYDict.items():
	coord = 0 

	
	plasmaColorTuple = pair[1]
	xyCoord = pair[0]
	x = float(xyCoord[0])
	y = float(xyCoord[1])

	
	
	panel2.plot(x,y, marker ='o', markerfacecolor = plasmaColorTuple, markersize = 4.1, markeredgewidth = 0.1, markeredgecolor = 'black', linewidth = 0, alpha = 1)
	howmanypointsright+=1
	

panel2.annotate(gene, (-23.5, -35), xycoords ='data', path_effects = [pe.withStroke(linewidth=0, foreground = "black")], horizontalalignment='center', verticalalignment='center')

for index in range(0,101,1):
	colorPalettePlas = (R[index], G[index], B[index])
	vvGradeRect = mplpatches.Rectangle([0,index], .1, .8, facecolor=colorPalettePlas,edgecolor = 'black',linewidth = 0)


#couldn't get the legend panel fully in place :(

		
	legendPanel.add_patch(vvGradeRect)


plt.savefig(outputFile, dpi=600)


