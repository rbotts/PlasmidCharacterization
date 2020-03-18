#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 19:42:20 2019

@author: ghovis
"""

import sys
import re
import csv
import copy
import os
import linecache
from collections import defaultdict
from Bio import SeqIO

csv.field_size_limit(100000000)

def rowsFromProkka(inFile, dbDelim = ".faa"):
	### Note that infile is a .csv file of output from Prokka
    	# Input 
	#	inFile = csv file for reading Prokka output
	# 	dbDelim = the key indicating user defined database in the Prokka output
	# Output
	#	outlist = a list of each CDS's data formatted as
	# 	(start, stop, arg, arg type, bg, inc croup, plasmidName, keep (Boolean indicating if it was from a passed in database)
	with open(inFile, "r") as csvfile:

		csvreader = csv.reader(csvfile, delimiter= '\t')

		#create an empty dictionary to temporarily hold the information for one gene	
		tempout = {"start" : "NA", "stop" : "NA", "geneName" : "NA", "plasmidName" : "NA", 'keep' : False, "plasmidFileId" : "NA", "plasmidNumContigs" : "NA"}
		outlist = []

		for line in csvreader: 
			# identify lines specifying the start of a new plasmid
			if len(line) == 1 and line[0].find(">")==0:
				# identify when the plasmids are switching and write out the previous record
				outlist.append([tempout[I] for I in tempout])
				# clear tempout
				tempout = {"start" : "NA", "stop" : "NA", "resName" : "NA", "resType" : "NA", "geneName" : "NA","incGroup" : "NA",\
				 "plasmidName" : line[0][line[0].index(" "):], 'keep' : False, "plasmidNumContigs" : "NA"}

			# identify rows starting a new CDS
			elif len(line) == 3 and line[2] == 'CDS':
				outlist.append([tempout[I] for I in tempout])
				# replace values in tempout, note that we only need to 
				# replace the plasmid name when we get to a new plasmid
				tempout["start"] = line[0]
				tempout["stop"] = line[1]
				tempout["geneName"] = "NA" 
				tempout["resName"] = "NA" 
				tempout["resType"] = "NA"
				tempout["incGroup"] = "NA"
				tempout["plasmidFileId"] = "NA"
				tempout["plasmidNumContigs"] = "NA"
				tempout['keep'] = False
		
			# identify if the CDS is from database and tag for keeping
			elif len(line)>4 and line[3] == 'inference' and ((line[4].find(dbDelim)>-1) or (line[4].find('ISfinder'))):
				tempout['keep'] = True
		
			# identify the gene name if it has one
			elif len(line)>4 and line[3] == 'product':
				if "group" not in line[4] and "RES" not in line[4]:
					tempout['geneName'] = line[4]
					tempout['incGroup'] = "NA"
					tempout['resName'] = "NA"
					tempout["resType"] = "NA"                    
				if "group" in line[4]:
					tempout['incGroup'] = line[4]
					tempout['geneName'] = "NA"
					tempout['resName'] = "NA"
					tempout["resType"] = "NA"                    
				if "RES" in line[4]:
					tempout['resName'] = line[4]
					tempout['geneName'] = "NA"
					tempout['incGroup'] = "NA"
					tempout["resType"] = "NA"
				if "_:_" in line[4]:
					resToList = line[4].split("_:_")
					tempout['resName'] = resToList[0]
					tempout['geneName'] = "NA"
					tempout['incGroup'] = "NA"
					tempout['resType'] = resToList[1]					                    
	return outlist

def cleanList(CDSlist):
	### function removes rows that were not from the current database
	# Input
	#	CDSlist = a list of the information contained for each CDS in one row. 
	#		Last element determines if it was from user defined database
	# Out
	#	returns a list of only the rows from CDS's from our database, 
	#	removes column identifying keeper rows
	return [row[:-1] for row in CDSlist if row[-1]]

#########################findNeighbors.py#######################
#takes gene table and finds the nearest backbone genes up and downstream of every ARG
#input file:file created by prokkaReformat.py
#output file is .csv of ARG (with ARG type), upstream BG, downstream BG, inc group, and plasmid name 
def findNeighbors(filePath):
    with open(filePath,'r') as  csvreadfile:
        csvreader = csv.reader(csvreadfile)
        tempN1 = ""
        N1 = ""
        N2 = ""
        outList = []
        resList = []
        afterRES = False 
        for line in csvreader:
        	#set gene in BG column to upstream neighbor temporarily
            if not line[4] == "NA" and afterRES == False:
                tempN1 = line[4]
            #if theres an ARG, set temp upstream neighbor to neighbor1 and assign ARG
            if not line[2] == "NA":
                N1 = tempN1
                resName = line[2]
                resType = line[3]#inserted by Leslie
                afterRES = True 
                #resLIst accounts for the possibility of multiple ARGS in row
                resList.append([resName, resType])#resType inserted by Leslie 
            #first BG that occurs after the res gene is neighbor 2 (downstream)
            if not line[4] == "NA" and afterRES == True:
                N2 = line[4]
                afterRES = False
                #Account for absent data.
                if len(line) > 9:
                    Plasmid= line[6]
                    Inc= line[14]
                else:
                    Plasmid= "NA"
                    Inc= "NA"
                if not line[14] == None:
                    Inc= line[14]
                else:
                    Inc= "NA"
                #every ARG in a cluster will share the same up/downstream BG neighbors
                for gene in resList:
                    #Careful to run this function only after obtaining the MOB table
                    outList.append([gene,N1,N2,Inc,Plasmid])#resType inserted by Leslie, Gabrielle changed line[5] to line[14] to incorporate reliable inc groups from MOB-suite
                tempN1 = line[4]
                resList = []
    return outList

######################geneTableEditor.py#############################
#function takes in a parameter, the parameter will be an absolute file pathway
# to the gene table file (should be .csv file from prokkaReformat.py)
#(ex:/Users/marielelensink/Documents/H2/GeneTable1Edit.csv)
#output is  another gene table that is compatible with R
#Gabrielle altered to add incompatibility group 
#leslie just changes all the numbers up one basically
def geneTableEditor(filePath):
    with open(filePath, 'r') as csvreadfile:
        csvreader = csv.reader(csvreadfile)
        incList= []
        #Create default dictionary for incompatibility (inc) groups 
        dDict= defaultdict(list)
        outList = []
        plasmid = ""
        for line in csvreader:
            if not (line[0] == "NA" and line[1] == "NA"):
                plasmid = line[6]
                plasmid = plasmid.strip()
                if not line[5] == 'NA':
                    incOld= line[5]
                    #Splice inc group name to make it easier to read
                    incSave= incOld[5:10]
                    #Add plasmid and inc group pairs to list
                    incList.append([plasmid, incSave])
                #Add all the inc groups to each plasmid key in the default dictionary
                for plasmid, incSave in incList:
                    dDict[plasmid].append(incSave)
    #Write the data for each plasmid to a new file
    with open(filePath, 'r') as csvreadfile:
        csvreader = csv.reader(csvreadfile)  
        resName1 = ""
        for line in csvreader:
            #nullrow is a count variable to indicate if values were assigned to a row for writing to the new spreadsheet.
            nullrow= 1
            #Avoid entering a gene without values (specifically the start and end positions)
            if not (line[0] == "NA" and line[1] == "NA"):
                start= line[0]
                stop = line[1]
                if "RES" in line[2]:
                    resName1 = line[2]
                    resName = resName1[3:]
                    resType = line[3]
                else:
                    resName = line[2]
                    resType = line[3]
                    geneName = line[4]
                    parseinc = line[5]
                    plasmid = line[6]
                    plasmidNumContigs = line[7]
                    parseinc = parseinc[5:10]
                    nullrow= 0
                    #Ensure that all the inc groups will be listed for each plasmid 
                    countplasmid = 0
                    for i in dDict.keys():
                        key= str(plasmid)
                        key= key.strip()
                        if (key == str(i)):
                            countplasmid += 1
                            curincs= dDict[key]
                            lencurincs= len(curincs)
                            count= 0
                            if (lencurincs > 1):
                                for j in range(0,lencurincs):
                                    oneinc= str(curincs[j])
                                    if (oneinc == parseinc):
                                        count += 1
                                if (count == 0):
                                    incgroup= ",".join(dDict[key])
                                    break
                            else:
                                if (str(curincs) != str(parseinc)):
                                    incgroup= ",".join(dDict[key])
                                    break
                    if (countplasmid == 0):
                        incgroup= "other"
                    if (nullrow == 0):
                        outList.append([start,stop,resName,resType,geneName,incgroup,plasmid,plasmidNumContigs])
    return outList    

		
############################ResGeneEdit.py##########################
#takes output.csv file from GeneTableEdit.py and edits to make more compatible with R
#changes instances of "NA" to "-", cleans ARG names 
#output .csv file of the same format 
#Leslie added ResType and moved most index numbers up 1 
def resGeneEdit(inputFile):
    with open(inputFile,'r') as csvreadfile:
        csvreader = csv.reader(csvreadfile)
        start = ""
        stop = ""
        ResGene1 = ""
        ResGene = ""
        ResType = ""
        BackboneGene = ""
        IncGroup = ""
        PlasmidName = ""
        plasmidNumContigs= ""
        outlist = []
        for line in csvreader:
            start = line[0]
            stop = line[1]
            #cleans/standardizes ARGS by removing subgroups 
            ResGene1 = line[2]
            ResType = line[3]
            ResGene = ResGene1[:4]
            #change "NA's" to hyphen for easier viewing
            if ResGene == "NA":
                ResGene = "-"
            if ResType == "NA":
                ResType = "-"
            BackboneGene = line[4]
            if BackboneGene == "NA":
                BackboneGene = "-"
            IncGroup = line[5]
            PlasmidName = line[6]
            plasmidNumContigs = line[7]
            outlist.append([start,stop,ResGene,ResType,BackboneGene,IncGroup,PlasmidName,plasmidNumContigs])
    return outlist

############################isolatePlasmid.py##########################
#Created by Gabrielle for mob_suite data
#Creates separate files for each plasmid sequence (mob_suite requires separate files for each plasmid) 
#Input: one fasta file with data from multiple plasmids
#Output: list of plasmid sequences, list of plasmid names associated with each sequence
def isolatePlasmid (plasmidFile):
    plasmidList= []
    plasmidNames= []
    organismDict= {}
    organism= ""
    #Open file containing all plasmids
    with open(plasmidFile) as handle:
        #Isolate sequence of one single plasmid
        for plasmidrec in SeqIO.parse(handle,'fasta'):
            sequence= ">"+str(plasmidrec.id)+"\n"+str(plasmidrec.seq)
            plasmidList.append(sequence)
            #Separate ID of plasmid from file denotation and save name of plasmid for future use
            plasmidID= str(plasmidrec.id)
            if (plasmidID[-6:] == ".fasta"):
                plasmidNames.append(plasmidID[:-6])
            elif (plasmidID[-4:] == ".faa"):
                plasmidNames.append(plasmidID[:-4])
            else:
                plasmidNames.append(plasmidID)
            #Make a list of all the host organisms
            
            #descript= [x.strip() for x in plasmidrec.description.split(", ")]
            #print(descript)
            #organism= str(descript[0])
            #organism= organism[14:]
            #print(organism)
            organism1= plasmidrec.description.split()[1]
            organism2= plasmidrec.description.split()[2]
            organism= organism1 + " " + organism2
            #organism= plasmidrec.description.split()[1]
            organismDict.update({plasmidID: organism})
        #print(organismDict)
    return plasmidList,plasmidNames,organismDict

############################mobTyper.py##########################
#Created by Gabrielle for mob_suite data
#Runs mob_suite on the plasmid sequences then inputs the data from each of the mob_suite output files into a multidimensional list
#Input: lists from isolatePlasmid function output
#Output: multidimensional list of mob_suite data (each new entry contains data for a different plasmid)
def mobTyper (plasmidList, plasmidNames):
    #Desired output folder (in mobsuite directory)
    outFolder= "MobOut"
    #Path for location of the output folder
    mobpath= "/Users/ghovis/mob-suite/mob_suite/"+str(outFolder)
    #Path for plasmid file
    filepath= "/Users/ghovis/Documents/Research2019-20/Databases/"
    #Create list for mob outputs
    mobOutputList= []
    
    for i in range(0,len(plasmidList)):
        #Create a new file to write sequence for one plasmid parsed from the file containing all plasmids
        outUnoPlasmid= open(str(filepath)+"unoPlasmid.fasta","w")
        outUnoPlasmid.write(plasmidList[i])
        os.chdir(filepath)
        os.rename("unoPlasmid.fasta","outPlasmid"+str(i)+".fasta")
        #Run mob_typer on the file containing the sequence of one plasmid (from isolatePlasmid output)
        #Output goes into output directory specified above
        os.system("mob_typer --infile %s --outdir %s" % (str(filepath)+"outPlasmid"+str(i)+".fasta", str(mobpath)))
        #Open mob output file
        newfilepath= str(mobpath)+"/mobtyper_"+"outPlasmid"+str(i)+".fasta_report.txt"
        lineNum= 2
        #Pull mob-suite output data from the file
        line= linecache.getline(newfilepath, lineNum)
        mobOutputList.append(line)
        
    #Edit lines for later conversion to CSV format
    mobOutputList= [entry.replace(",", ";") for entry in mobOutputList]
    mobOutputList= [entry.replace("\t", ",") for entry in mobOutputList]
    return mobOutputList

############################cleanMobTyper.py##########################
#Created by Gabrielle for mob_suite data
#A default mob_suite run creates output files using the file name as the plasmid identifier. The file name must be replaced with the name of the plasmid when this data is utilized.
#Replaces file name in plasmid identifier column with the name of the plasmid
#Input: csv file from mobTyper function, list of plasmid names created in isolatePlasmid function
#Output: new list with edited mob_suite data (file name replaced with plasmid name)
def cleanMobTyper (mobFile, plasmidNames):
    outFilePath= "/Users/ghovis/Documents/Research2019-20/Databases/"
    outFileName= "outMob.csv"
    #Create new list for edited output from mobTyper function
    newMobOutput= []
    with open(outFilePath+outFileName,"r") as openmob:
        mobreader= csv.reader(openmob,skipinitialspace= "True")
        i= 0
        for line in mobreader:
            #Clear empty entries
            if (len(line) == 1):
                line= ""
            else:
                #Replace name of file with name of associated plasmid
                newPName= str(plasmidNames[i])
                line[0]= newPName
                newMobOutput.append(line)
                if (i < (len(plasmidNames)-1)):
                    i += 1
    for i in range(0,len(plasmidNames)):
        os.unlink(outFilePath+"outPlasmid"+str(i)+".fasta")
    return newMobOutput

############################addMobData.py##########################
#Created by Gabrielle for mob-suite data
#Adds data from the csv file containing mob-suite output to new columns for each plasmid containing the previous data
#Input: csv file from resGeneEdit function, csv file from cleanMobTyper function
#Output: list with previous data and the new mob-suite data for each plasmid
def addMobData(inputFile, mobFile, organismDict):
    mobDict= {}
    species= ""
    handle= open("newFile.csv", "w")
    newFile= csv.writer(handle)
    #Open file with mob-suite data and read rows
    with open(mobFile,"r") as readcsvfile:
        readMobFile= csv.reader(readcsvfile)
        #For each row, add a new dictionary entry; {(plasmid name) : (the row of mob-suite data in a list form)}
        for row in readMobFile:
            key= row[0].strip()
            value= list(row)
            mobDict.update({key: value})
    #Open output csv file from resGeneEdit function and read rows
    with open(inputFile, "r") as readFile:
        for line in csv.reader(readFile):
            #Check to ensure that there are enough entries to implement the code
            if not (len(line) < 9):
                if (line[6].strip() in organismDict):
                    species= organismDict[str(line[6].strip())]
                    line.append(species)
                #If the plasmid name from inputFile is in the mob-suite dictionary
                if (line[6].strip() in mobDict):
                    #Add the mob-suite data to the inputFile row for the corresponding plasmid and write to the new csv file
                    newFile.writerow(line+mobDict[str(line[6]).strip()])
                #When the plasmid name is not referenced in the mob-suite data, add only the row from inputFile to the new csv file
                else:
                    newFile.writerow(line)
    handle.close()
    return 


###################step1(ProkkaReformat.py)#########################
#INPUT: Prokka will output 10 files, this function takes in the file from Prokka ending in .tbl.
#OUTPUT: Function will output a .csv file with the start index, stop index, ARG name, Backbone Gene Name, 
#fName begins each file name to identify output files
fName = "NewRun"
with open(fName+'TableOutput.csv','w') as csvfile: #name of output file here
    csvwriter = csv.writer(csvfile)
    for i in range(1,19):
        cleanedUp = (cleanList(rowsFromProkka("PROKKA_06282019 ("+str(i)+").tbl")))
        csvwriter.writerows(cleanedUp)

#################step2(geneTableEditor.py)#########################
editedTable = geneTableEditor("/Users/ghovis/Documents/Research2019-20/Databases/"+fName+"TableOutput.csv")#name of input file from prokkareformat
with open(fName+'EditedGeneTable.csv','w') as csvwritefile:#need to change the input here
    csvwriter = csv.writer(csvwritefile)
    csvwriter.writerows(editedTable)

######################step3(ResGeneEdit.py)########################
geneOutputList = resGeneEdit(fName+'EditedGeneTable.csv')
with open(fName+'RCompatibleTable.csv','w') as csvwritefile: #write in output file name 
    csvwriter = csv.writer(csvwritefile)
    csvwriter.writerows(geneOutputList)

######################step4(isolatePlasmid.py)########################
#Path for plasmid file
filepath= "/Users/ghovis/Documents/Research2019-20/Databases/"
#Name of file with multiple plasmids
filename= "1794plasmidSeq_9-6.fasta"

plasmidList, plasmidNames, organismDict= isolatePlasmid(str(filepath)+str(filename))

######################step5(mobTyper.py)########################
mobRun= mobTyper(plasmidList, plasmidNames)
with open("outMob.csv","a",newline= '') as csvfile:
    #Make new csv file for mob_suite data of all the plasmids
    writecsv= csv.writer(csvfile, delimiter= "\t", lineterminator= "'", skipinitialspace= "True")
    writecsv.writerow(mobRun)
    
######################step6(cleanMobTyper.py)########################
#Name of csv file from mobRun output
mobfile= "outMob.csv"
cleanMob= cleanMobTyper(filepath+mobfile,plasmidNames)
with open("replaceFile.csv","a",newline= '') as replaceFile:
    writenew= csv.writer(replaceFile, skipinitialspace= "True")
    for i in range(0,len(cleanMob)):
        writenew.writerow(cleanMob[i])
    #Replace the mobRun output file with the edited spreadsheet from the cleanMob function
    os.rename(filepath+"replaceFile.csv","outMob.csv")
    
######################step7(addMobData.py)################################
#Name of csv file from last run.
inputFilePath= "/Users/ghovis/Documents/Research2019-20/Databases/"
inputFile= fName+"TableOutput.csv"
#Name of csv file with mob_suite data from each plasmid
outMobPath= "/Users/ghovis/Documents/Research2019-20/Databases/"
outMobName= "outMob.csv"
addMobData(inputFilePath+inputFile,outMobPath+outMobName,organismDict)
os.rename(inputFilePath+"newFile.csv",fName+"mobTable.csv")

#################step8(findNeighbors.py)##########################
newNeighbors = findNeighbors('/Users/ghovis/Documents/Research2019-20/Databases/'+fName+'mobTable.csv')
with open(fName+'Neighbors.csv','w') as csvwritefile:
    csvwriter = csv.writer(csvwritefile)
    csvwriter.writerows(newNeighbors)