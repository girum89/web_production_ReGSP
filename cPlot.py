#!/usr/bin/python3
#-*- coding: utf-8 -*-
#
#
# ---------------------------------------------------------------------------------------------------------------------------
#  cPlot (comparison Plot)
#  Original Author		: Mingeun Ji (Dept. of Multimedia Engineering at Dongguk University, Seoul, Korea)
#  E-mail		: mingeun@mme.dongguk.edu
#  Version		: 0.1
#  Rev. Date	: Oct. 06, 2019.

#  cPlot Adjusted to handle ,ultiple references for use in ReGSP
#  Nodified by: Girum Fitihamlak Ejigu
#  E-mail: girum89@mju.ac.kr
#  Required software packages to use cPlot
#   - Python3
#
# ---------------------------------------------------------------------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import os
import sys
import time
import itertools
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import gc
import psutil
from collections import OrderedDict
GLOBAL_CONF = {}
"""
	Web Conf
"""

GLOBAL_CONF["DIR"] = ""								#output directory
#GLOBAL_CONF["THRESHOLD"] = 13
GLOBAL_CONF["FILE_NAME"] = {}				#file names dictionary
GLOBAL_CONF["FILE_NAME"]["REFER"] = []		#reference file names
GLOBAL_CONF["FILE_NAME"]["QUERY"] = []		#query file names

GLOBAL_CONF["FILE_DATA"] = {}
GLOBAL_CONF["FILE_DATA"]["REFER"] = {}		#GLOBAL_CONF["FILE_DATA"]["REFER"]
GLOBAL_CONF["FILE_DATA"]["QUERY"] = {}
GLOBAL_CONF["FILE_DATA"]["REVER"] = {}

GLOBAL_CONF["LENGTH"] = {}
GLOBAL_CONF["LENGTH"]["REFER"] = {}			# length of each reference
GLOBAL_CONF["LENGTH"]["QUERY"] = {}

GLOBAL_CONF["COUNT_REFER"] = 0
GLOBAL_CONF["COUNT_QUERY"] = 0

GLOBAL_CONF["SLICING_DATA"] = {}			# sliced data dictonary
GLOBAL_CONF["SLICING_DATA"]["REFER"] = {}	#
GLOBAL_CONF["SLICING_DATA"]["QUERY"] = {}	#
GLOBAL_CONF["SLICING_DATA"]["REVER"] = {}	#

GLOBAL_CONF["INTERSECTION_DATA"] = {}			#슬라이싱된 레퍼런스 데이터를 기준으로 QUERY와 RC 쿼리의 교집합 read에 대한 딕셔너리
GLOBAL_CONF["INTERSECTION_DATA"]["FORWARD"] = {}  #쿼리 read와 RC 쿼리 read와 교집합인 레퍼런스의 read 리스트 { REF_IDX: {QUERY_IDX: ['ACG', 'CGT', ... ] }, ...}
GLOBAL_CONF["INTERSECTION_DATA"]["BACKWARD"] = {}	#레퍼런스 read와 쿼리 read간의 교집합인 read 리스트 {REF_IDX: {RCQUERY_IDX: ['ACG', CGT', ...]}, ...}

GLOBAL_CONF["INDEX"] = {}						#각 데이터의 read마다 인덱스 부여한 딕셔너리
GLOBAL_CONF["INDEX"]["FORWARD"] = {}				#GLOBAL_CONF["INDEX"]["REVER"][REF_IDX][QUERY_IDX]["ABC"] = [ 0 , 2 , 3 , ... ]
GLOBAL_CONF["INDEX"]["BACKWARD"] = {}				#GLOBAL_CONF["INDEX"]["REVER"][REF_IDX][RECQUERY_IDX]["ABC"] = [ 0, 2, 3, ...]

GLOBAL_CONF["MAX_INDEX_READ"] = {}			#각 Refer, Query, RC-Query의 중복되지 않는 각 read에 대한 개수를 값으로 가지는 딕셔너리
GLOBAL_CONF["MAX_INDEX_READ"]["FORWARD"] = {}	#레퍼런스내 read의 개수 GLOBAL_CONF["MAX_INDEX_READ"]["REFER"][REF_IDX][QUERY_IDX]["ABC"] = 1 ...
GLOBAL_CONF["MAX_INDEX_READ"]["BACKWARD"] = {} #쿼리내 read의 개수 GLOBAL_CONF["MAX_INDEX_READ"]["QUERY"][REF_IDX][QUERY_IDX]["ABC"] = 1 ....

GLOBAL_CONF["BEFORE_MERGE"] = {}
GLOBAL_CONF["BEFORE_MERGE"]["FORWARD"] = {}		#GLOBAL_CONF["BEFORE_MERGE"]["FORWARD"][REF_IDX][QUERY_IDX] = {'ACG_1':(start_x_pos, start_y_pos, length),... }
GLOBAL_CONF["BEFORE_MERGE"]["BACKWARD"] = {}
			
GLOBAL_CONF["SLOPE"] = {}
GLOBAL_CONF["SLOPE"]["FORWARD"] = {}			#GLOBAL_CONF["SLOPE"]["FORWARD"][REF_IDX][QUERY_IDX] = {기울기:['ACG_1', ... ], ... }
GLOBAL_CONF["SLOPE"]["BACKWARD"] = {}

GLOBAL_CONF["MERGE"] = {}
GLOBAL_CONF["MERGE"]["FORWARD"] = {}			#GLOBAL_CONF["MERGE"]["FORWARD"][REF_IDX][QUERY_IDX] = {'ACG_1' : (start_x_pos, start_y_pos, length), ... }
GLOBAL_CONF["MERGE"]["BACKWARD"] = {}

GLOBAL_CONF["BEFORE_OPTIMAL"] = {}				#MERGE-FORWARD, MERGE-BACKWARD 통합된 딕셔너리, GLOBAL_CONF["BEFORE_OPTIMAL"][REF_IDX][QUERY_IDX][START_FIRST_X] = { 'ACG_1', : (start_x_pos, start_y_pos, length), ...}

GLOBAL_CONF["X_POS_READ"] = {}					#GLOBAL_CONF["X_POS_READ"][REF_IDX][QUERY_IDX][X_POS_DIGIT] = ['ACG_1' , ... ] #X_POS_DIGIT = X_POS_0, X_POS_1, ...
GLOBAL_CONF["X_POS_READ_VALUE"] = {}			#GLOBAL_CONF["X_POS_READ_VALUE"][REF_IDX][QUERY_IDX][X_POS_DIGIT] = {'ACG_1':value, ...}
GLOBAL_CONF["LIST_OPTIMAL"] = []
GLOBAL_CONF["OPTIMAL"] = {}						#GLOBAL_CONF["OPTIMAL"][REF_IDX][QUERY_IDX] = { 'ACG_1':(start_x_pos, start_y_pos, length), ... }

GLOBAL_CONF["SHIFT"] = {}						#GLOBAL_CONF["SHIFT"][REF_IDX][QUERY_IDX] = { 파일번호:long_read_가중치}

GLOBAL_CONF["SHIFT_DATA"] = {}


def READ_FILE(fpFile, listFileName, listData, length):
	strContent = ""
	idx = 0
	fileName = ""
	for lines in fpFile.readlines():
		if lines[0] == '>':
			if strContent != "":
				listData[idx] = strContent
				length[idx] = len(strContent)
				idx += 1
				strContent = ""
				
			fileName = lines.split(' ')[0][1:]
			listFileName.append(fileName)
				
		else:
			strContent += lines.strip()

	listData[idx] = strContent
	length[idx] = len(strContent)

	
def REVERSE_QUERY_DATA(listQueryData, listRCQueryData):
	#print(">>REVERS", end="\t")
	start = time.time()
	for idx in range(0, len(listQueryData) ):
		strRevData = listQueryData[idx][::-1]
		strRCData = ""
		for char in strRevData:
			if char == 'A':
				strRCData += 'T'
			elif char == 'C':
				strRCData += 'G'
			elif char == 'G':
				strRCData += 'C'
			elif char == 'T':
				strRCData += 'A'
		listRCQueryData[idx] = strRCData
	end = time.time() - start
	#print(end)

def SLICING_BY_SEED(listData, listSlicingData,kmer=13):
	#print(">>SLICING", end="\t")
	start = time.time()
	seed = kmer
	for idx, data in listData.items():
		listSlicingData[idx] = {}
		for pos in range(0, len(data)-seed+1):
			strRead = data[pos:pos+seed]
			if 'N' in strRead:
				continue
			if not strRead in listSlicingData[idx].keys():
				listSlicingData[idx][strRead] = []

			listSlicingData[idx][strRead].append(pos)
	end = time.time() - start
	#print(end)

def REVERSE_SLICING_BY_SEED(listData, listSlicingData,kmer=13):
	#print(">>REVERSE SLICING", end="\t")
	start = time.time()
	seed = kmer
	for idx, data in listData.items():
		listSlicingData[idx] = {}
		max_idx = len(data) - seed
		for pos in range(0, max_idx):
			strRead = data[pos:pos+seed]
			if 'N' in strRead:
				continue
			#print(strRead, max_idx-pos+seed-1)
			if not strRead in listSlicingData[idx].keys():
				listSlicingData[idx][strRead] = []
			listSlicingData[idx][strRead].append(max_idx-pos+seed-1)
	end = time.time() - start
	#print(end)

def INTERSECTION(listSrcData, listDestData, listIntersection):
	#print(">>Intersection", end="\t")
	start = time.time()
	for idxSrc, dataSrc in listSrcData.items():
		#print('>'+str(idxSrc), dataSrc)
		listIntersection[idxSrc] = {}
		for idxDest, dataDest in listDestData.items():
			listSrcKey = list(dataSrc.keys())
			listDestKey = list(dataDest.keys())

			listIntersection[idxSrc][idxDest] = list( set(listSrcKey).intersection(listDestKey))
	end = time.time() - start
	#print(end)

def PREPROCESSING_BEFORE_MERGE(listSrcData, listDestData, listIntersectionData, listIndex, listIndexCount, listBeforeMerge, listSlope,kmer=13):
	#print(">>PREPROCESSING",end='\t')
	start = time.time()
	seed = kmer
	length_refer = 0
	for idxRef, idxValue in listIntersectionData.items():
		length_query = 0
		for idxQuery, data in idxValue.items():
			for intersectionRead in data:
				if intersectionRead in listSrcData[idxRef].keys() and intersectionRead in listDestData[idxQuery].keys():
					product_list = list(itertools.product(listSrcData[idxRef][intersectionRead], listDestData[idxQuery][intersectionRead]))
					if not idxRef in listIndex.keys():
						listIndex[idxRef] = {}
					if not idxQuery in listIndex[idxRef].keys():
						listIndex[idxRef][idxQuery] = {}

					listIndex[idxRef][idxQuery][intersectionRead] = product_list

					for item in product_list:
						#print('>>>item', intersectionRead, item)
						if not idxRef in listBeforeMerge.keys():
							listBeforeMerge[idxRef] = {}
						if not idxQuery in listBeforeMerge[idxRef].keys():
							listBeforeMerge[idxRef][idxQuery] = {}

						#	MAX_INDEX_READ[REF_IDX][QUERY_IDX]에서 카운팅된 READ가 있는지 확인하여 포지션에서 라벨링함
						if not idxRef in listIndexCount.keys():
							listIndexCount[idxRef] = {}
						if not idxQuery in listIndexCount[idxRef].keys():
							listIndexCount[idxRef][idxQuery] = {}

						if not intersectionRead in listIndexCount[idxRef][idxQuery].keys():
							listIndexCount[idxRef][idxQuery][intersectionRead] = 1
						else:
							listIndexCount[idxRef][idxQuery][intersectionRead] += 1

						#print(idxRef, idxQuery, intersectionRead)
						label = str(intersectionRead) + '_' + str(listIndexCount[idxRef][idxQuery][intersectionRead])
						length = 0
						slope = 0

						if listBeforeMerge is GLOBAL_CONF["BEFORE_MERGE"]["FORWARD"]:
							length = seed
							#slope = item[1]-item[0]#y=x+b b=y-x --> b=item[1]-item[0]
						else:
							length = -seed
							#slope = item[1]+item[0]#y=-x+b b=y+x --> b=item[1]+item[0]
						#ref_x = item[0] + length_refer
						#query_y = item[1] + length_query
						slope = GET_SLOPE(item[0], item[1], length)
						listBeforeMerge[idxRef][idxQuery][label] = (item[0], item[1], length)
						
						if not idxRef in listSlope.keys():
							listSlope[idxRef] = {}
						if not idxQuery in listSlope[idxRef].keys():
							listSlope[idxRef][idxQuery] = {}
						if not slope in listSlope[idxRef][idxQuery].keys():
							listSlope[idxRef][idxQuery][slope] = []
						listSlope[idxRef][idxQuery][slope].append(label)
			length_query += GLOBAL_CONF["LENGTH"]["QUERY"][idxQuery]

		length_refer += GLOBAL_CONF["LENGTH"]["REFER"][idxRef]
								
	del(listIntersectionData)
	del(listIndex)
	gc.collect()
	end = time.time() - start
	#print(end)
    

def GET_SLOPE(x, y, length):
	return (((y+length)-y)/((x+abs(length)-x)))

def GET_POINT(first):
	x_pos = first[0]
	y_pos = first[1]
	length = first[2]

	if length > 0:
		return (x_pos, y_pos), (x_pos+length-1, y_pos+length-1)
	else:
		return (x_pos, y_pos), (x_pos-length-1, y_pos+length+1)

def CHECK_COND_MERGE(first, second):
	flag = False
	if first[0] > second[0]:
		#swap
		temp = first
		first = second
		second = temp
		flag = True

	a_1, a_2 = GET_POINT(first)
	b_1, b_2 = GET_POINT(second)
	if a_2[0] == b_1[0]:
		return (a_1[0], a_1[1], 1, flag)
	elif a_2[0] > b_1[0]:
		if b_2[0] <= a_2[0]:
			return (a_1[0], a_1[1], -1, flag)
		else:
			if a_1[0] == b_1[0]:
				return (b_1[0], b_1[1], -2, flag)
			else:
				return (a_1[0], a_1[1], a_2[0]-b_1[0]+1, flag)

	else:
		return None
	
"""
	MERGE sequnces that have overlaping sequences
"""
def MERGE(fileData, listBeforeMerge, listSlope, listIndexCount, listMerge, listBeforeOptimal,threshold=13):
	#print(">>MERGE",end='\t')
	start = time.time()
	listCopyBeforeMerge = copy.deepcopy(listBeforeMerge)
	listBeforeOrderedMerge = {}
	listMergeInprocess = {}
	listAvailableMerge = []
	listOrderedMerge = OrderedDict()
	for idxRef, idxValue in listSlope.items():
		if not idxRef in listBeforeOptimal.keys():
			listBeforeOptimal[idxRef] = {}
		for idxQuery,slopeData in idxValue.items():
			if not idxQuery in listBeforeOptimal[idxRef].keys():
				listBeforeOptimal[idxRef][idxQuery] = {}
			for slope, data  in slopeData.items():
				listBeforeOrderedMerge = {}
				listOrderedMerge = {}

				for read in data:
					listBeforeOrderedMerge[read] = listCopyBeforeMerge[idxRef][idxQuery][read]
					#print(read, listBeforeOrderedMerge[read])
				
				ordered = sorted(listBeforeOrderedMerge.items(), key=lambda k: k[1][0])

				for read in ordered:
					listOrderedMerge[read[0]] = read[1]
				firstRead = ""
				startPosX = -1
				startPosY = -1
				endPosX = -1
				length = -1
				secondRead = ""
				newStartPosX = -1
				newStartPosY = -1
				newLength = -1
				first = None
				second = None

				idxMax = len(listOrderedMerge) - 1
				idx = 0
				if not idxRef in listMerge.keys():
					listMerge[idxRef] = {}
				if not idxQuery in listMerge[idxRef].keys():
					listMerge[idxRef][idxQuery] = {}

				keys = list(listOrderedMerge.keys())
				flag = False	#합쳐졌는지 아닌지 구분
				mergeList = []
				while True:
					if idxMax == 0 :
						read, value = list(listOrderedMerge.items())[0]
						lengthRead = abs(value[2])
						read = read[:-2]
						if lengthRead < threshold:
							break
						if not read in listIndexCount[idxRef][idxQuery].keys():
							listIndexCount[idxRef][idxQuery][read] = 1
							read += '_1'
						else:
							listIndexCount[idxRef][idxQuery][read] += 1
							read += '_' + str(listIndexCount[idxRef][idxQuery][read])
						listMerge[idxRef][idxQuery][read] = value
						if not value[0] in listBeforeOptimal[idxRef][idxQuery].keys():
							listBeforeOptimal[idxRef][idxQuery][value[0]] = {}
						listBeforeOptimal[idxRef][idxQuery][value[0]][read] = value
							
						break
						
					if idx > idxMax:
						newEndPosX = newStartPosX + abs(newLength) - 1
						read = fileData[idxRef][startPosX:newEndPosX]
						newLength = len(read)
						if abs(newLength) < threshold:
							break

						if slope < 0 :
							newLength = -newLength
						
						if not read in listIndexCount[idxRef][idxQuery].keys():
							listIndexCount[idxRef][idxQuery][read] = 1
							read += '_1'
						else:
							listIndexCount[idxRef][idxQuery][read] += 1
							read += '_' + str(listIndexCount[idxRef][idxQuery][read])
						listMerge[idxRef][idxQuery][read] = (startPosX, startPosY, newLength)
						if not startPosX in listBeforeOptimal[idxRef][idxQuery].keys():
							listBeforeOptimal[idxRef][idxQuery][startPosX] = {}
						listBeforeOptimal[idxRef][idxQuery][startPosX][read] = (startPosX, startPosY, newLength)
						#print('read2', idx, idxMax, read, listMerge[idxRef][idxQuery][read], slope)
						break
					
					if startPosX == -1:
						firstRead = keys[idx]
						#처음 비교값이 없을 경우
						startPosX = listOrderedMerge[firstRead][0]
						startPosY = listOrderedMerge[firstRead][1]
						length = abs(listOrderedMerge[firstRead][2])
						endPosX = startPosX + length - 1
						#print('read', firstRead, listOrderedMerge[firstRead], slope)
						idx += 1
					else:	
						#비교값이 있을 경우
						secondRead = keys[idx]
						newStartPosX = listOrderedMerge[secondRead][0]
						newStartPosY = listOrderedMerge[secondRead][1]
						newLength = abs(listOrderedMerge[secondRead][2])
						newEndPosX = newStartPosX + newLength - 1
						if startPosX <= newStartPosX and newStartPosX <= endPosX:
							endPosX = newEndPosX
							flag = True
							idx += 1
						else:
							"""

							"""
							if flag is False:
								first = (startPosX, startPosY, length)
								second = (newStartPosX, newStartPosY, newLength)
								if not first in mergeList:
									firstRead = fileData[idxRef][startPosX:endPosX]
									length = len(firstRead)
									if not abs(length) < threshold:
										if slope < 0 :
											length = -length
										mergeList.append(first)
										if not firstRead in listIndexCount[idxRef][idxQuery].keys():
											listIndexCount[idxRef][idxQuery][firstRead] = 1
											firstRead += '_1'
										else:
											listIndexCount[idxRef][idxQuery][firstRead] += 1
											firstRead += '_' + str(listIndexCount[idxRef][idxQuery][firstRead])
										first = (first[0], first[1], first[2])
										listMerge[idxRef][idxQuery][firstRead] = first

										if not first[0] in listBeforeOptimal[idxRef][idxQuery].keys():
											listBeforeOptimal[idxRef][idxQuery][first[0]] = {}
										listBeforeOptimal[idxRef][idxQuery][first[0]][firstRead] = first

								if not second in mergeList:
									secondRead = fileData[idxRef][newStartPosX:newEndPosX]
									newLength = len(secondRead)
									if not abs(newLength) < threshold:
										if slope < 0 :
											newLength = -newLength
										mergeList.append(second)
										if not secondRead in listIndexCount[idxRef][idxQuery].keys():
											listIndexCount[idxRef][idxQuery][secondRead] = 1
											secondRead += '_1'
										else:
											listIndexCount[idxRef][idxQuery][secondRead] += 1
											secondRead += '_' + str(listIndexCount[idxRef][idxQuery][secondRead])
										second = (second[0], second[1], second[2])
										listMerge[idxRef][idxQuery][secondRead] = second

										if not second[0] in listBeforeOptimal[idxRef][idxQuery].keys():
											listBeforeOptimal[idxRef][idxQuery][second[0]] = {}
										listBeforeOptimal[idxRef][idxQuery][second[0]][secondRead] = second
								idx += 1
							else:
								secondRead = keys[idx-1]
								newStartPosX = listOrderedMerge[secondRead][0]
								newStartPosY = listOrderedMerge[secondRead][1]
								newLength = abs(listOrderedMerge[secondRead][2])
								newEndPosX = newStartPosX + newLength - 1
								mergeRead = fileData[idxRef][startPosX:newEndPosX]
								mergeLength = len(mergeRead)
								#print('merge', mergeRead, len(mergeRead))
								if not abs(mergeLength) < threshold:
									if slope < 0 :
										mergeLength = -mergeLength
									merge = (startPosX, startPosY, mergeLength)
									
									if not merge in mergeList:
										mergeList.append(merge)
										if not mergeRead in listIndexCount[idxRef][idxQuery].keys():
											listIndexCount[idxRef][idxQuery][mergeRead] = 1
											mergeRead += '_1'
										else:
											listIndexCount[idxRef][idxQuery][mergeRead] += 1
											mergeRead += '_' + str(listIndexCount[idxRef][idxQuery][mergeRead])
										listMerge[idxRef][idxQuery][mergeRead] = merge
										if not merge[0] in listBeforeOptimal[idxRef][idxQuery].keys():
											listBeforeOptimal[idxRef][idxQuery][merge[0]] = {}
										listBeforeOptimal[idxRef][idxQuery][merge[0]][mergeRead] = merge
								

							flag = False
							startPosX = -1
							startPosY = -1
							endPosX = -1
							length = -1						
				listOrderedMerge.clear()
			if listBeforeOptimal[idxRef][idxQuery] == {}:
				del(listBeforeOptimal[idxRef][idxQuery])
	del(fileData)
	del(listBeforeMerge)
	del(listSlope)
	del(listIndexCount)
	gc.collect()
	end = time.time() - start
	#print(end)

def GET_VALUE_FOR_OPTIMIZATION(item):
	weight_length = GLOBAL_CONF["WEIGHT_LENGTH"]
	weight_pos = GLOBAL_CONF["WEIGHT_POSITION"]
	weight_overlap = GLOBAL_CONF["WEIGHT_OVERLAP"]	#오버랩 가중치: Optimal 계산할 때 각 리드당 오버랩 안된 점수에 대한 가중치

	a_1, b_1 = GET_POINT(item)
	length = abs(item[2])

def GET_SCORE_FOR_OPTIMIZATION(value):
	alpha = GLOBAL_CONF["WEIGHT_LENGTH"]
	beta = GLOBAL_CONF["WEIGHT_POSITION"]
	gamma = value[2] <= 0 and GLOBAL_CONF["WEIGHT_DIRECTION"] or 1-GLOBAL_CONF["WEIGHT_DIRECTION"]
	length_value = abs(value[2])
	first_point, second_point = GET_POINT(value)
	a = -((second_point[1]-first_point[1])/(second_point[0]-first_point[0]))
	b = 1
	c = 0
	central_point = ( (second_point[0]+first_point[0])/2, (second_point[1]+first_point[1])/2)
	pos_value = abs(a*central_point[0]+b*central_point[1]+c)/(((a**2)+(b**2))**0.5)

	total = (length_value*alpha-pos_value*beta)*gamma
	
	return (length_value*alpha-pos_value*beta)*gamma

def CHECK_COND_OVERLAP(first, second):
	a_1, a_2 = GET_POINT(first)
	b_1, b_2 = GET_POINT(second)
		
	flag = None
	if a_2[0] == b_1[0]:
		flag = True
	elif a_2[0] > b_1[0]:
		flag = False

	return flag

def CHECK_COND_OPTIMAL(first, second):
	a_1, a_2 = GET_POINT(first)
	b_1, b_2 = GET_POINT(second)


	flag = None
	if a_2[0] == b_1[0]:
		flag = True
	elif a_2[0] > b_1[0]:
		if b_2[0] <= a_2[0]:
			flag = False
		else:
			if a_1[0] == b_1[0]:
				flag = False
			else:
				flag = True
	else:
		flag = None
	
	return flag

def GET_POWERSET(listSrc):
	size = len(listSrc)
	powerSet = []
	for idx in range(2**size):
		flag = bin(idx)[2:].zfill(size)
		subset = [listSrc[nested_idx] for nested_idx in range(size) if flag[nested_idx] == '1']
		if len(subset) > 0 :
			powerSet.append(subset)
	return powerSet

def CHECK_OVERLAP(listSrc, listCopyBeforeOptimal):
	sum = 0
	for firstItem in listSrc:
		first_value = list(listCopyBeforeOptimal[firstItem].values())[0]
		for secondItem in listSrc:
			if firstItem == secondItem:
				break
					
			second_value = list(listCopyBeforeOptimal[secondItem].values())[0]
					
			flag_overlap = CHECK_COND_OVERLAP(first_value, second_value)
			if flag_overlap == True:
				pass
			else:
				return -1
		sum += GET_SCORE_FOR_OPTIMIZATION(first_value)
	return sum

def CALCULATE_OPTIMAL(listSrc, listCopyBeforeOptimal):
	size = len(listSrc)
	powerSet = []
	best_score = None
	best_subset = None
	for idx in range(2**size):
		flag = bin(idx)[2:].zfill(size)
		subset = [listSrc[nested_idx] for nested_idx in range(size) if flag[nested_idx] == '1']
		if len(subset) > 0 :
			flag_overlap = CHECK_OVERLAP(subset, listCopyBeforeOptimal)

			if flag_overlap != -1:
				if best_score == None or best_score < flag_overlap:
					best_score = flag_overlap
					best_subset = subset
				
	return best_subset, best_score

def OPTIMIZATION(listBeforeOptimal, listOptimal, listOptimalScore):
	#print(">>OPTIMIZATION",end='\t')
	start = time.time()
	listCopyBeforeOptimal = copy.deepcopy(listBeforeOptimal)
	for idxRef, idxValue in listCopyBeforeOptimal.items():
		if not idxRef in listOptimal.keys():
			listOptimal[idxRef] = {}
		#print(len( listCopyBeforeOptimal[idxRef]))
		for idxQuery,optimalData in idxValue.items():
			if not idxQuery in listOptimal[idxRef].keys():
				listOptimal[idxRef][idxQuery] = {}
			#"SCOREING"
			optimalDataList = list(optimalData.keys())
			optimalDataList.sort()
			for startIdx  in optimalDataList:
				best_score = None
				best_read = ""
				keys = list(optimalData[startIdx].keys())
				for idx in range(0, len(keys)):
					read = keys[idx]
					value = optimalData[startIdx][read]

					score = GET_SCORE_FOR_OPTIMIZATION(value)

					if best_score == None or best_score < score:

						if best_score != None:
							del listCopyBeforeOptimal[idxRef][idxQuery][startIdx][best_read]
						best_score = score
						best_read = read
						
					else:
						del listCopyBeforeOptimal[idxRef][idxQuery][startIdx][read]

			#"OPTIMIZATION"
			srcRead = ""
			srcValue = None
			listCandidateOptimalRead = []
			keys = list( listCopyBeforeOptimal[idxRef][idxQuery].keys())
			keys.sort()
			idx = 0
			total_score = 0
			while True:
				if idx == len(keys):
					if len(listCandidateOptimalRead) == 1:
						start_point = listCandidateOptimalRead[0]
						read, value = list(listCopyBeforeOptimal[idxRef][idxQuery][start_point].items())[0]
						listOptimal[idxRef][idxQuery][read] = value
						
					break
				start_point = keys[idx]
				read, value = list(listCopyBeforeOptimal[idxRef][idxQuery][start_point].items())[0]

				if srcRead == "":
					srcRead = read
					srcValue = value
					listCandidateOptimalRead.append(srcValue[0])
					idx+=1
					continue

				flag = CHECK_COND_OPTIMAL(srcValue, value)

				if flag == False:
					pass
				elif flag == True:
					listCandidateOptimalRead.append(srcValue[0])
					listCandidateOptimalRead.append(value[0])
				else:
					#flag == None
					"""-> listCandidateOptimalRead Optimal"""
					best_set, best_score = CALCULATE_OPTIMAL(listCandidateOptimalRead, listCopyBeforeOptimal[idxRef][idxQuery])
					if not best_set is None:
						for start_pos in best_set:
							read, value = list(listCopyBeforeOptimal[idxRef][idxQuery][start_pos].items())[0]
							listOptimal[idxRef][idxQuery][read] = value
							#print('best', len(read), value)
						total_score += best_score
					srcRead = ""
					srcValue = None
					listCandidateOptimalRead.clear()
					continue
				idx += 1
			
	slope_list = {}
	slope_list_reads = []
	f = open(GLOBAL_CONF["DIR"] + "/optimal.output", "w")
	for idxRef, idxValue in listCopyBeforeOptimal.items():
		for idxQuery,optimalData in idxValue.items():
			ordered = sorted( listOptimal[idxRef][idxQuery].items(), key=lambda k:k[1][0] )
			for read, value in ordered:
				if not idxQuery in GLOBAL_CONF["SHIFT_DATA"].keys():
					GLOBAL_CONF["SHIFT_DATA"][idxQuery] = ""
				if value[2] > 0: #
					y_threshold = int(GLOBAL_CONF["LENGTH"]["QUERY"][idxQuery]/GLOBAL_CONF["LENGTH"]["REFER"][idxRef]*value[0])
					temp = listOptimal[idxRef][idxQuery][read]
					listOptimal[idxRef][idxQuery][read] = (value[0], y_threshold, value[2])
					f.write( ">" + GLOBAL_CONF["FILE_NAME"]["QUERY"][idxQuery] + "_" + str(temp[0]) + "_" + str(temp[1]) + "_" + str(y_threshold) + "_" + str(value[2]) + "\n")
					f.write(read[:-2] +"\n")

					GLOBAL_CONF["SHIFT_DATA"][idxQuery] += ">" + GLOBAL_CONF["FILE_NAME"]["QUERY"][idxQuery] + "_"+ str(temp[0]) + "_" + str(temp[1]) + "_" + str(y_threshold) + "_" + str(value[2]) + "\n" + read[:-2] + "\n";

					if len(slope_list) > 0 :
						limit = listOptimal[idxRef][idxQuery][read][1]-1
						slope_list_len = len(slope_list)
						for idx in range(0, int(slope_list_len)):
							temp = slope_list[slope_list_reads[idx]]
							gap = int(slope_list[slope_list_reads[0]][1]) - int(temp[1])
							adjust=limit-gap
							if(adjust >= GLOBAL_CONF["LENGTH"]["QUERY"][idxQuery] or adjust<0):
								if ((limit+gap)>= GLOBAL_CONF["LENGTH"]["QUERY"][idxQuery] or (limit+gap)<0 ):
									adjust=abs(gap)
								else:
									adjust=limit+gap
							listOptimal[idxRef][idxQuery][slope_list_reads[idx]] = (temp[0], adjust, temp[2]) 
							f.write( ">" + GLOBAL_CONF["FILE_NAME"]["QUERY"][idxQuery] + "_" + str(temp[0]) + "_" + str(temp[1]) + "_" + str(adjust) + "_" + str(temp[2]) + "\n")
							f.write(slope_list_reads[idx][:-2] +"\n")
							GLOBAL_CONF["SHIFT_DATA"][idxQuery] += ">" + GLOBAL_CONF["FILE_NAME"]["QUERY"][idxQuery] + "_"+ str(temp[0]) + "_" + str(temp[1]) + "_" + str(adjust) + "_" + str(temp[2]) + "\n" + slope_list_reads[idx][:-2] + "\n";
						slope_list = {}
						slope_list_reads = []
							
				elif value[2] < 0:
					y_threshold = int(GLOBAL_CONF["LENGTH"]["QUERY"][idxQuery]/GLOBAL_CONF["LENGTH"]["REFER"][idxRef]*value[0])
					slope_list[read] = listOptimal[idxRef][idxQuery][read]
					slope_list_reads.append(read)

			if len(slope_list) > 0 :
				limit = listOptimal[idxRef][idxQuery][read][1] - 1
				slope_list_len = len(slope_list)
				for idx in range(0, int(slope_list_len)):
					temp = slope_list[slope_list_reads[idx]]
					gap = int(slope_list[slope_list_reads[0]][1]) - int(temp[1])
					adjust=limit-gap
					if(adjust >= GLOBAL_CONF["LENGTH"]["QUERY"][idxQuery] or adjust<0):
						if ((limit+gap)>= GLOBAL_CONF["LENGTH"]["QUERY"][idxQuery] or (limit+gap)<0 ):
							adjust=abs(gap)
						else:
							adjust=limit+gap
					listOptimal[idxRef][idxQuery][slope_list_reads[idx]] = (temp[0], adjust, temp[2]) 
					f.write( ">" + GLOBAL_CONF["FILE_NAME"]["QUERY"][idxQuery] + "_" + str(temp[0]) + "_" + str(temp[1]) + "_" + str(adjust) + "_" + str(temp[2]) + "\n")
					f.write(slope_list_reads[idx][:-2] + "\n")
					GLOBAL_CONF["SHIFT_DATA"][idxQuery] += ">" + GLOBAL_CONF["FILE_NAME"]["QUERY"][idxQuery] + "_"+ str(temp[0]) + "_" + str(temp[1]) + "_" + str(adjust) + "_" + str(temp[2]) + "\n" + slope_list_reads[idx][:-2] + "\n";
				slope_list = {}
				slope_list_reads = []
				
	end = time.time() - start
	#print(end)

def SORT_QUERY_IDX(listOptimal, listSortedQueryIdx):
	#print(">>SORTING",end='\t')
	start = time.time()
	listCopyOptimal = copy.deepcopy(listOptimal)
	listBeforeSortedQueryIdx = {}
	#print(listOptimal)
	for idxRef, idxValue in listOptimal.items():
		if not idxRef in listSortedQueryIdx.keys():
			listBeforeSortedQueryIdx[idxRef] = {}
		if not idxRef in listSortedQueryIdx.keys():
			listSortedQueryIdx[idxRef] = []
		for idxQuery, data in idxValue.items():
			if not idxQuery in listBeforeSortedQueryIdx[idxRef].keys():
				listBeforeSortedQueryIdx[idxRef][idxQuery] = {}
			long_point = None
			long_read = ""
			for read, pos in data.items():
				value = pos[2]
				if long_point == None or abs(long_point[2]) < value:
					long_point = pos
					long_read = read
			
			listBeforeSortedQueryIdx[idxRef][idxQuery][long_read] = long_point
		ordered = sorted( listBeforeSortedQueryIdx[idxRef].items(), key=lambda k:list(k[1].values())[0][0] )
		for idx in ordered:
			listSortedQueryIdx[idxRef].append(idx[0])

		
	end = time.time() - start
	#print(end)

def DRAW_MERGE_READ(t,outfilename,kmert):
	global GLOBAL_CONF

	color_ = [ 'red', 'royalblue', 'deepskyblue', 'chartreuse', 'black', 'm', 'tan', 'rosybrown', 'pink', 'gold', 'mediumpurple']
	line_style = {'Solid': ":", 'Dashed': "-.", 'Dash-dot': "--", "Dotted": "-"}
	marker_style = {"point":".", "pixel":",", "circle":"o", "octagon":"8", "square":"s", "pentagon":"p", "hexagon":"h", "plus":"+", "plus_filled":"P", "x":"x", "x_filled":"X", "star":"*"}


	listMergeForward = GLOBAL_CONF["MERGE"]["FORWARD"]
	listMergeBackward = GLOBAL_CONF["MERGE"]["BACKWARD"]

	listOptimal = GLOBAL_CONF["OPTIMAL"]
	listShift = GLOBAL_CONF["SHIFT"]
						
	seq_count_refer = GLOBAL_CONF["COUNT_REFER"]
	seq_count_query = GLOBAL_CONF["COUNT_QUERY"]
	title="cPlot Result"

	color_idx = 0
	COLOR_LINE_SET = []
	COLOR_MARKER_SET = []

	for _idx in range(0, seq_count_query):
		if color_idx == 10:
			color_idx = 0
		COLOR_LINE_SET.append(color_[color_idx])
		COLOR_MARKER_SET.append(color_[color_idx])
		color_idx += 1

	fig,subPlot_shift = plt.subplots(1,seq_count_refer,sharey='row',gridspec_kw = {'wspace':0, 'hspace':0})

	#reads_thickness = GLOBAL_CONF["READS_THICKNESS"]/10
	reads_thickness =2/10

	#marker_thickness = GLOBAL_CONF["MARKER_THICKNESS"]/10
	marker_thickness = 2/10
	marker_ms = marker_style["circle"]
	axis_thickness = 2/10
	axis_color = "black"
	axis_ls = line_style["Solid"]

	list_length_shift = []
	optimal_query_ticks = {}
	optimal_query_ticks_label = {}
	optimal_query_ticks_major = {}
	optimal_query_ticks_major_labels = {}
	length_query_=0
	lines_optimal = {}
	#print(GLOBAL_CONF["SHIFT_DATA"].keys())

	f = open(GLOBAL_CONF["DIR"] + "/shift.output", "w")
	for idx in range(0, seq_count_refer):
		if len(listShift) >0:
			ticks_val = 0
			optimal_value = 0
			query_major_lables=[]
			query_major_ticks=[]
			query_minticks=[]
			query_minticks_label=[]
			line_legend=[]
			for shift_idx in listOptimal[idx]: #listShift[idx]
				for k, v in listOptimal[idx][shift_idx].items():
					a_1, a_2 = GET_POINT(v)
					x_1 = a_1[0]
					x_2 = a_2[0]
					y_1 = a_1[1] + optimal_value
					y_2 = a_2[1] + optimal_value 
					if seq_count_refer==1:
						plot_lines, =subPlot_shift.plot( [x_1, x_2], [y_1, y_2], c=COLOR_LINE_SET[shift_idx], markerfacecolor=COLOR_MARKER_SET[shift_idx], markersize=marker_thickness, marker=marker_ms, lw=reads_thickness, zorder=2, label=GLOBAL_CONF["FILE_NAME"]["QUERY"][shift_idx])
					else:
						plot_lines, =subPlot_shift[idx].plot( [x_1, x_2], [y_1, y_2], c=COLOR_LINE_SET[shift_idx], markerfacecolor=COLOR_MARKER_SET[shift_idx], markersize=marker_thickness, marker=marker_ms, lw=reads_thickness, zorder=2, label=GLOBAL_CONF["FILE_NAME"]["QUERY"][shift_idx])
				
				query_major_lables.append(GLOBAL_CONF["FILE_NAME"]["QUERY"][shift_idx])
				list_length_shift.append(optimal_value)
				len_ = GLOBAL_CONF["LENGTH"]["QUERY"][shift_idx]
				optimal_value += GLOBAL_CONF["LENGTH"]["QUERY"][shift_idx]
				add_val = int(GLOBAL_CONF["LENGTH"]["QUERY"][shift_idx]/4)
				ticks_label_val = 0 
				line_legend.append(plot_lines)
				while True:
					ticks_label_val += add_val
					query_minticks.append(ticks_val+ticks_label_val)
					query_minticks_label.append(str(ticks_label_val))
					if ticks_label_val+add_val > len_ - add_val :
						break
				query_major_ticks.append(ticks_val)
				ticks_val += GLOBAL_CONF["LENGTH"]["QUERY"][shift_idx]
				if seq_count_refer==1:
					subPlot_shift.axhline(y=ticks_val, color=axis_color, linewidth=axis_thickness, ls=axis_ls, zorder=1)
				else:
					subPlot_shift[idx].axhline(y=ticks_val, color=axis_color, linewidth=axis_thickness, ls=axis_ls, zorder=1)
				#print(GLOBAL_CONF["SHIFT_DATA"][shift_idx])
				f.write(GLOBAL_CONF["SHIFT_DATA"][shift_idx])
			optimal_query_ticks_major[idx] = query_major_ticks
			optimal_query_ticks_major_labels[idx]=query_major_lables
			optimal_query_ticks[idx] = query_minticks
			optimal_query_ticks_label[idx] = query_minticks_label
			length_query_=optimal_value
			lines_optimal[idx]=line_legend
	f.close()
	optimal_name = []
	for idx in range(0, seq_count_refer):
		for idx_ in listOptimal[idx]: #listShift[idx]
			optimal_name.append( GLOBAL_CONF["FILE_NAME"]["QUERY"][idx_] )

	ticks_val = 0

	for idx in GLOBAL_CONF["LENGTH"]["REFER"]:
		if seq_count_refer==1:
			subPlot_shift.set_xlabel(GLOBAL_CONF["FILE_NAME"]["REFER"][idx])
		else:
			subPlot_shift[idx].set_xlabel(GLOBAL_CONF["FILE_NAME"]["REFER"][idx])
	

	for idx in GLOBAL_CONF["LENGTH"]["REFER"]:
		if seq_count_refer==1:
			#pass
			subPlot_shift.tick_params(axis='x', which='major', direction='inout', colors='gray')
			subPlot_shift.legend(handles=lines_optimal[idx], labels=optimal_name,bbox_to_anchor=(-0.4,-0.13), loc="lower left", ncol=len(optimal_name), fontsize='x-small')
		else:
			subPlot_shift[idx].tick_params(axis='x', which='major', direction='inout', colors='gray', rotation=315)
			subPlot_shift[0].legend(handles=lines_optimal[idx], labels=optimal_name,bbox_to_anchor=(-0.4,-0.125), loc="lower left", ncol=len(optimal_name), fontsize='x-small')

	
	for idx in GLOBAL_CONF["LENGTH"]["REFER"]:
		if seq_count_refer==1:
			subPlot_shift.tick_params(axis='x', which='minor', direction='inout', colors='gray')
			subPlot_shift.tick_params(axis='x', which='major', labelsize=6)
		else:
			subPlot_shift[idx].tick_params(axis='x', which='minor', direction='inout', colors='gray')
			subPlot_shift[idx].tick_params(axis='x', which='major', labelsize=6)

	if seq_count_refer==1:
		subPlot_shift.set_ylabel('Query')
	else:
		subPlot_shift[0].set_ylabel('Query')
	
	for idx in GLOBAL_CONF["LENGTH"]["REFER"]:
		if seq_count_refer==1:
			subPlot_shift.yaxis.set_major_locator(ticker.FixedLocator(optimal_query_ticks_major[idx]))
			subPlot_shift.yaxis.set_major_formatter(ticker.FixedFormatter(optimal_query_ticks_major_labels[idx]))
			subPlot_shift.tick_params(axis='y', which='major', direction='inout', width=2, length=10)
		else:
			subPlot_shift[idx].yaxis.set_major_locator(ticker.FixedLocator(optimal_query_ticks_major[idx]))
			subPlot_shift[idx].yaxis.set_major_formatter(ticker.FixedFormatter(optimal_query_ticks_major_labels[idx]))
			subPlot_shift[idx].tick_params(axis='y', which='major', direction='inout', width=2, length=10)
	
	for idx in GLOBAL_CONF["LENGTH"]["REFER"]:
		if seq_count_refer==1:
			subPlot_shift.yaxis.set_minor_locator(ticker.FixedLocator(optimal_query_ticks[idx]))
			subPlot_shift.yaxis.set_minor_formatter(ticker.FixedFormatter(optimal_query_ticks_label[idx]))
			subPlot_shift.tick_params(axis='y', which='minor', direction='inout', colors='gray')
			subPlot_shift.tick_params(axis='y', which='minor', labelsize=6)
			subPlot_shift.tick_params(axis='y', which='major', labelsize=8)

		else:
			subPlot_shift[idx].yaxis.set_minor_locator(ticker.FixedLocator(optimal_query_ticks[idx]))
			subPlot_shift[idx].yaxis.set_minor_formatter(ticker.FixedFormatter(optimal_query_ticks_label[idx]))
			subPlot_shift[idx].tick_params(axis='y', which='minor', direction='inout', colors='gray')
			subPlot_shift[idx].tick_params(axis='y', which='minor', labelsize=6)
			subPlot_shift[idx].tick_params(axis='y', which='major', labelsize=8)
	

	for idx in GLOBAL_CONF["LENGTH"]["REFER"]:
		if seq_count_refer==1:
			subPlot_shift.set_xlim(0, GLOBAL_CONF["LENGTH"]["REFER"][idx])
			subPlot_shift.set_ylim(0, length_query_)

		else:
			subPlot_shift[idx].set_xlim(0, GLOBAL_CONF["LENGTH"]["REFER"][idx])
			subPlot_shift[idx].set_ylim(0, length_query_)

	
	plt.suptitle('cPlot for $kmer='+str(kmert)+'$')
	
	#plt.show()
	fig = plt.gcf()
	fig.set_size_inches(12.1, 9.1, forward=True)
	
	fig.tight_layout(rect=[0, 0.01, 1, 0.98])
	fig.savefig(GLOBAL_CONF["DIR"] + '/'+outfilename+'.cPlot.pdf', format='pdf')
	fig.savefig(GLOBAL_CONF["DIR"] + '/'+outfilename+'.cPlot.png', format='png',dpi=250)
	#GLOBAL_CONF["FILE_NAME"]["QUERY"]=[]
	#GLOBAL_CONF["FILE_NAME"]["REFER"]=[]


if __name__ == '__main__':
	process = psutil.Process(os.getpid())
	cpu_before = process.cpu_percent(interval=1)
	mem_berfore = process.memory_info().rss / 1024 / 1024
	GLOBAL_CONF["FILE_REFER"] = open(sys.argv[1], 'r')
	GLOBAL_CONF["FILE_QUERY"] = open(sys.argv[2], 'r')
	GLOBAL_CONF["DIR"] = sys.argv[3]
	outfilename=sys.argv[4]
	kmert=int(sys.argv[5])
	start = time.time()


	READ_FILE(GLOBAL_CONF["FILE_REFER"], GLOBAL_CONF["FILE_NAME"]["REFER"], GLOBAL_CONF["FILE_DATA"]["REFER"], GLOBAL_CONF["LENGTH"]["REFER"])
	READ_FILE(GLOBAL_CONF["FILE_QUERY"], GLOBAL_CONF["FILE_NAME"]["QUERY"], GLOBAL_CONF["FILE_DATA"]["QUERY"], GLOBAL_CONF["LENGTH"]["QUERY"])
	GLOBAL_CONF["COUNT_REFER"] = len(GLOBAL_CONF["LENGTH"]["REFER"])
	GLOBAL_CONF["COUNT_QUERY"] = len(GLOBAL_CONF["LENGTH"]["QUERY"])

	REVERSE_QUERY_DATA(GLOBAL_CONF["FILE_DATA"]["QUERY"], GLOBAL_CONF["FILE_DATA"]["REVER"])
	
	SLICING_BY_SEED(GLOBAL_CONF["FILE_DATA"]["REFER"], GLOBAL_CONF["SLICING_DATA"]["REFER"],kmert)
	SLICING_BY_SEED(GLOBAL_CONF["FILE_DATA"]["QUERY"], GLOBAL_CONF["SLICING_DATA"]["QUERY"],kmert)
	REVERSE_SLICING_BY_SEED(GLOBAL_CONF["FILE_DATA"]["REVER"], GLOBAL_CONF["SLICING_DATA"]["REVER"],kmert)

	INTERSECTION(GLOBAL_CONF["SLICING_DATA"]["REFER"], GLOBAL_CONF["SLICING_DATA"]["QUERY"], GLOBAL_CONF["INTERSECTION_DATA"]["FORWARD"])
	INTERSECTION(GLOBAL_CONF["SLICING_DATA"]["REFER"], GLOBAL_CONF["SLICING_DATA"]["REVER"], GLOBAL_CONF["INTERSECTION_DATA"]["BACKWARD"])

	PREPROCESSING_BEFORE_MERGE(GLOBAL_CONF["SLICING_DATA"]["REFER"], GLOBAL_CONF["SLICING_DATA"]["QUERY"], GLOBAL_CONF["INTERSECTION_DATA"]["FORWARD"], GLOBAL_CONF["INDEX"]["FORWARD"], GLOBAL_CONF["MAX_INDEX_READ"]["FORWARD"], GLOBAL_CONF["BEFORE_MERGE"]["FORWARD"], GLOBAL_CONF["SLOPE"]["FORWARD"],kmert)
	PREPROCESSING_BEFORE_MERGE(GLOBAL_CONF["SLICING_DATA"]["REFER"], GLOBAL_CONF["SLICING_DATA"]["REVER"], GLOBAL_CONF["INTERSECTION_DATA"]["BACKWARD"], GLOBAL_CONF["INDEX"]["BACKWARD"], GLOBAL_CONF["MAX_INDEX_READ"]["BACKWARD"], GLOBAL_CONF["BEFORE_MERGE"]["BACKWARD"], GLOBAL_CONF["SLOPE"]["BACKWARD"],kmert)

	
	MERGE(GLOBAL_CONF["FILE_DATA"]["REFER"], GLOBAL_CONF["BEFORE_MERGE"]["FORWARD"], GLOBAL_CONF["SLOPE"]["FORWARD"], GLOBAL_CONF["MAX_INDEX_READ"]["FORWARD"], GLOBAL_CONF["MERGE"]["FORWARD"], GLOBAL_CONF["BEFORE_OPTIMAL"],kmert)
	
	MERGE(GLOBAL_CONF["FILE_DATA"]["REFER"], GLOBAL_CONF["BEFORE_MERGE"]["BACKWARD"], GLOBAL_CONF["SLOPE"]["BACKWARD"], GLOBAL_CONF["MAX_INDEX_READ"]["BACKWARD"], GLOBAL_CONF["MERGE"]["BACKWARD"], GLOBAL_CONF["BEFORE_OPTIMAL"],kmert)

	OPTIMIZATION(GLOBAL_CONF["BEFORE_OPTIMAL"],GLOBAL_CONF["OPTIMAL"], None)

	SORT_QUERY_IDX(GLOBAL_CONF["OPTIMAL"], GLOBAL_CONF["SHIFT"])

	a = time.time()-start
	DRAW_MERGE_READ(a,outfilename,kmert)