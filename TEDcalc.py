#!/usr/bin/env python
import argparse
import os
import sys
from Bio import SeqIO
import pybedtools

def mainArgs():
	parser = argparse.ArgumentParser(description='Calculate TE dispersion distance.',
									 prog='TEDcalc')
	# Input options
	parser.add_argument('--genome',type=str,required=True,default=None,help='Fasta file containing sequences referenced by gff file.')
	parser.add_argument('--gff',type=str,required=True,help='GFF/GTF formatted repeat annotations.')
	args = parser.parse_args()
	return args

def getGFF(infile=None):
	a = pybedtools.BedTool(infile)
	b = a.merge()
	featureDict = dict()
	IDcounter = 1
	for x in b:
		if x[0] not in featureDict.keys():
			featureDict[x[0]] = list()
			featureDict[x[0]].append((IDcounter,int(x[1]),int(x[2])))
			IDcounter += 1
		else:
			featureDict[x[0]].append((IDcounter,int(x[1]),int(x[2])))
			IDcounter += 1
	return featureDict 


def getLens(infile=None):
	lenDict = dict()
	for rec in SeqIO.parse(infile, "fasta"):
		lenDict[rec.id] = int(len(rec.seq))
	return lenDict

def featureSpace(intervals):
	# From list of intervals get total length covered by all features
	T = dict()
	for chrm in intervals.keys():
		TEtotal = 0
		for x in intervals[chrm]:
			TEtotal += abs(int(x[1]) - int(x[2]))
		T[chrm] = TEtotal
	return T

def getdist(r1, r2):
	# sort the two ranges such that the range with smaller first element
	# is assigned to x and the bigger one is assigned to y
	x, y = sorted((r1, r2))
	#now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
	#then the ranges are not overlapping and return the differnce of y[0] and x[1]
	#otherwise return 0 
	if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
		return y[0] - x[1]
	return 0

def getStatbyChrm(intervals,TElens,chrlens):
	chrmScores = list()
	for chrm in intervals.keys():
		print(chrm)
		T = TElens[chrm]
		S = chrlens[chrm]
		Tprop = T / S
		print("T = %s S = %s Tprop = %s" % (str(T),str(S),str(Tprop)) )
		rowMeans = list()
		if len(intervals[chrm]) >1:
			# For each feature i in chrom
			for i in intervals[chrm]:
				iSum = 0
				# Look at all other features j in chrom
				for j in intervals[chrm]:
					# If not self
					if i[0] != j[0]:
						print("compare %s to %s" % (str(i[0]),str(j[0])))
						# Get distance from self, scaled to chrom len
						d = getdist((i[1],i[2]),(j[1],j[2])) / S
						print("d = %s" % str(d))
						# Multiply the distance by the weight of the target
						# Add to cumulative weighted score for feature i
						jw = (abs(j[1] - j[2]) / T)
						print("Target weight = %s " % str(jw))
						iSum += d * jw
				# Append mean weighted distance for feature i to 
				rowMeans.append(iSum/(len(intervals[chrm])-1))
			# Get mean dist across all features (rows) on chrom
			globalMean = sum(rowMeans) / len(rowMeans)
			# Update chrm stats with mean weighted dist and prop. of chrm that is covered by feature class.
			chrmScores.append((chrm,globalMean,Tprop))
		else:
			# To avoid div by 0 errors chromsomes with only 1 TE (or none, but we'll fix that later) get a dispersion score of 0
			chrmScores.append((chrm,0,Tprop))
	return chrmScores

def main():
	# Get cmd line args
	args = mainArgs()
	# Import and flatten features
	if args.gff:
		if os.path.isfile(args.gff):
			intervals = getGFF(infile=args.gff)
		else:
			print("Input file not found: %s" % args.gff)
			sys.exit(1)
	# Get chrom lens
	chrlens = getLens(infile=args.genome)
	# Get Total TE space within chroms
	TElens = featureSpace(intervals)
	# Get (chrmName, dispersion stat, prop TEs) per chromosome
	chrmScores = getStatbyChrm(intervals,TElens,chrlens)
	# Print results
	print('\t'.join(['Chrm','D-stat','TE-prop','ChrmLen']))
	for name,score,prop in chrmScores:
		print('\t'.join([name,str(score),str(prop),str(chrlens[name])]))
