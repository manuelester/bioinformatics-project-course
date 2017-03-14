def windowautomation(windowsize,crossvalnum,cvalue):
	print('running SVM model for window size of: ' + str(windowsize) + ' and ' + str(crossvalnum) + '-fold cross-validation')
	import matplotlib
	matplotlib.use('Agg')
	# testingfile feb 26 new PSSMs from offline PSI-BLAST run
	# aa order: A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V

	import numpy as np

	# define sigmoid function to normalize substitution matrix scores
	import math
	def sigmoid(x):
		return 1 / (1 + math.exp(-x))

	############################ open test file and and create empty lists #########################################
	opentest = open('dssp_8_state.3line.txt', 'r+')
	opentest = opentest.readlines()
	protname = []
	proteinseq = []
	windows = []
	secondstruc = []
	pssmlist = []
	lst = []




	############################################# user chooses window size #########################################
	win = windowsize
	zero = win/2 - 0.5
	zero = int(zero)




	################################################################################################################ append protein names, sequence segments (according to window size), and secondary structure elements into empty lists from above
	# open PSSM file corresponding to protein ID so that order of windows in window list matches order of secondary structures in secondary structure list
	################################################################################################################
	for i in range(0, len(opentest), 3):
		protid = opentest[i].rstrip('\n')
		if protid in protname:
			protid = opentest[i].rstrip('\n') + '_2'
			if protid in protname:
				protid = opentest[i].rstrip('\n') + '_3'
				if protid in protname:
					protid = opentest[i].rstrip('\n') + '_4'
		protname.append(protid)
		# print(protid)
		protseq = opentest[i+1].rstrip('\n')
		proteinseq.append(protseq)
		second = opentest[i+2].rstrip('\n')
		secondstruc.append(second)




	################################## randomize proteins before splitting into windows ############################
	from sklearn.utils import shuffle
	non_random = list(zip(protname, secondstruc))
	randomized = shuffle(non_random)
	randomized_protname = []
	randomized_secondstruc = []
	for i in range(0,len(randomized)):
		randomized_protname.append(randomized[i][0])
		randomized_secondstruc.append(randomized[i][1])




	################################ split proteins into cross-validation sets #####################################
	cvnum = crossvalnum
	avg = len(randomized_protname)/float(cvnum)
	protname_test = []
	protname_train = []
	ss_test = []
	ss_train = []
	ss_test_new = []
	ss_train_new = []
	windows_test = []
	windows_train = []
	last = 0.0

	while last < len(randomized_protname):
		protname_test.append(randomized_protname[int(last):int(last + avg)])
		protname_train.append(randomized_protname[0:int(last)] + randomized_protname[int(last + avg):])
		ss_test.append(randomized_secondstruc[int(last):int(last + avg)])
		ss_train.append(randomized_secondstruc[0:int(last)] + randomized_secondstruc[int(last + avg):])
		last += avg

	#print(len(protname_test))
	#print(len(protname_train))
	#print(len(protname_test[0]))



	################################# make windows list for test sets ############################################
	pad = list(np.zeros(20))
	for i in range(0,len(protname_test)):
		tmplist = []
		pssmlist = []
		for j in range (0,len(protname_test[i])):
			pssmlist = []
			pssmfile = open('../input/PSSM/' + protname_test[i][j] + '.fasta.pssm', 'r+')
			pssmfile = pssmfile.readlines()

			for k in range(3, len(pssmfile)-6):
				seq = pssmfile[k].rstrip('\n')
				seq = seq.split()
				# seq = seq[2:22]
				seq = seq[22:42]

				# make sure that all the numbers are converted to float by iterating through the list
				for l in range (0,len(seq)):
					seq[l] = (float(seq[l]))/100
				pssmlist.append(seq)
			for k in range(0,len(pssmlist)):
				if k < zero:
					seqwin = pssmlist[0:win-(zero-k)]
					seqwin = [j for k in seqwin for j in k]
					window = (pad*(zero-k)) + seqwin
					#print(len(window))
					tmplist.append(window)
				elif k >= (len(pssmlist)-zero):
					addprot = pssmlist[(k-zero):win-(zero-k)+1]
					addprot = [j for k in addprot for j in k]
					#print(len(addprot))
					if len(addprot) != win*20:
						addprot.extend(pad*(int(win-(len(addprot)/20))))
					#print(len(addprot))
					tmplist.append(addprot)
				else:
					window = pssmlist[(k-zero):(k+zero+1)]
					window = [j for k in window for j in k]
					#print(len(window))
					tmplist.append(window)
		windows_test.append(tmplist)

		# use map to convert features to numbers
		map = {'H':1, 'G':2, 'I':3, 'E':4, 'B':5, 'T':6, 'S':7, 'C':8}
		ss2 = ''.join(ss_test[i])
		ss2 = [char for char in ss2]
		ss2 = [map[i] for i in ss2]
		ss_test_new.append(ss2)

		#print(len(ss_test_new[i]))





	#################################### make windows list for train sets ##########################################
	for i in range(0,len(protname_train)):
		tmplist = []
		pssmlist = []
		for j in range (0,len(protname_train[i])):
			pssmlist = []
			pssmfile = open('../input/PSSM/' + protname_train[i][j] + '.fasta.pssm', 'r+')
			pssmfile = pssmfile.readlines()
		
			# convert each line of PSSMfile into a list by splitting on white space, only take range [2:] because column 1 is the position and column 2 the original amino acid which are not needed
			for k in range(3, len(pssmfile)-6):
				seq = pssmfile[k].rstrip('\n')
				seq = seq.split()
				# seq = seq[2:22]
				seq = seq[22:42]

				# make sure that all the numbers are converted to float by iterating through the list
				for l in range (0,len(seq)):
					seq[l] = (float(seq[l]))/100
				pssmlist.append(seq)
			for k in range(0,len(pssmlist)):
				if k < zero:
					seqwin = pssmlist[0:win-(zero-k)]
					seqwin = [j for k in seqwin for j in k]
					window = (pad*(zero-k)) + seqwin
					#print(len(window))
					tmplist.append(window)
				elif k >= (len(pssmlist)-zero):
					addprot = pssmlist[(k-zero):win-(zero-k)+1]
					addprot = [j for k in addprot for j in k]
					#print(len(addprot))
					if len(addprot) != win*20:
						addprot.extend(pad*(int(win-(len(addprot)/20))))
					#print(len(addprot))
					tmplist.append(addprot)
				else:
					window = pssmlist[(k-zero):(k+zero+1)]
					window = [j for k in window for j in k]
					#print(len(window))
					tmplist.append(window)
		
		windows_train.append(tmplist)
		#print(len(tmplist))
	
		# use map to convert features to numbers
		map = {'H':1, 'G':2, 'I':3, 'E':4, 'B':5, 'T':6, 'S':7, 'C':8}
		ss2 = ''.join(ss_train[i])
		ss2 = [char for char in ss2]
		ss2 = [map[i] for i in ss2]
		ss_train_new.append(ss2)
	
		#print(len(ss_train_new[i]))





	########################################## LINEAR SVM NEW CODE WITH AUTO ITERATION ############################
	from sklearn import svm
	from sklearn.model_selection import cross_val_score
	cvalue = int(cvalue)
	predlist = []
	for i in range(0,len(windows_train)):
		lin_clf = svm.LinearSVC(C = cvalue)
		lin_clf.fit(windows_train[i], ss_train_new[i])
		prediction = lin_clf.predict(windows_test[i])
		#predictionprobs = lin_clf.predict_proba(X_test[0])
		predlist.append(list(prediction))

	#clf = svm.SVC(decision_function_shape = 'ovo')
	#clf.fit(X_train1, Y_train1)
	#B = clf.predict(X_test1)

	# compare output to truth to get accuracy and confusion matrix

	from sklearn.metrics import confusion_matrix
	from sklearn.metrics import accuracy_score
	from sklearn.metrics import f1_score
	from sklearn.metrics import roc_auc_score
	from sklearn.metrics import precision_score
	from sklearn.metrics import recall_score
	from sklearn.metrics import average_precision_score
	from sklearn.metrics import classification_report
	# i am now confused
	stat_matrix = []
	acc = []
	stat_matrix_1 = confusion_matrix(ss_test_new[0], predlist[0], labels=[1,2,3,4,5,6,7,8])
	stat_matrix_normalized = stat_matrix_1.astype('float') / stat_matrix_1.sum(axis=1)[:, np.newaxis]
	F1 = f1_score(ss_test_new[0], predlist[0], average = 'weighted')
	precision = precision_score(ss_test_new[0], predlist[0], average = 'weighted')
	recall = recall_score(ss_test_new[0], predlist[0], average = 'weighted')
	target_names=['H:1','G:2','I:3','E:4','B:5','T:6','S:7','C:8']
	class_report = classification_report(ss_test_new[0], predlist[0], target_names=target_names)
	for i in range(0,len(predlist)):
		# stat_matrix.append(confusion_matrix(ss_test_new[i],predlist[i],labels=[1,2,3,4,5,6,7,8]))
		acc.append(accuracy_score(ss_test_new[i], predlist[i]))	

	#'H':1, 'G':2, 'I':3, 'E':4, 'B':5, 'T':6, 'S':7, 'C':8

	# Do accuracy average and sensitivity/specificity calculations
	averageacc = sum(acc)/len(acc)

	#TP = stat_matrix_1[0,0]
	#FP = sum(stat_matrix_1[0,1:])
	#TN = float(sum(stat_matrix_1[1,1:])) + float(sum(stat_matrix_1[2,1:])) + float(sum(stat_matrix_1[3,1:])) + float(sum(stat_matrix_1[4,1:])) + float(sum(stat_matrix_1[5,1:])) + float(sum(stat_matrix_1[6,1:])) + float(sum(stat_matrix_1[7,1:]))
	#FN = sum(stat_matrix_1[1:,0])

	#MCC_H = (TP*FN - FP*TN) / (((TP+FP)(TP+FN)(TN+FP)(TN+FN))**-1)

	# Write statistics to a file
	printfile = open('SM_windowsize_' + str(win) + '_fold_validation_' + str(cvnum) + '_cvalue_' + str(cvalue) + '.txt', 'w+')
	printfile.write('Average accuracy is: ' + str(averageacc) + '\n')
	printfile.write('Accuracies are: ' + str(acc) + '\n')

	printfile.write('F1 score is: ' + str(F1) + '\n')
	printfile.write('Precision score is: ' + str(precision) + '\n')
	printfile.write('Recall score is: ' + str(recall) + '\n')
	printfile.write(class_report + '\n')





	######################################### plot beautiful confusion matrix ######################################
	import itertools
	import numpy as np
	import matplotlib.pyplot as plt

	from sklearn import svm, datasets
	from sklearn.model_selection import train_test_split
	from sklearn.metrics import confusion_matrix

	classes = [1,2,3,4,5,6,7,8]
	plt.figure()
	plt.imshow(stat_matrix_normalized, interpolation='nearest', cmap='PuBuGn')
	plt.title('SM_ws: ' + str(win) + ' fold-validation: ' + str(cvnum) + ' F1: ' + str(F1) + ' cvalue: ' + str(cvalue))
	plt.colorbar()
	tick_marks = np.arange(len(classes))
	plt.xticks(tick_marks, classes, rotation=45)
	plt.yticks(tick_marks, classes)

	thresh = stat_matrix_normalized.max() / 2.
	for i, j in itertools.product(range(stat_matrix_normalized.shape[0]), range(stat_matrix_normalized.shape[1])):
		plt.text(j, i, round(stat_matrix_normalized[i, j],3), horizontalalignment="center", color="white" if stat_matrix_normalized[i, j] > thresh else "black")

	plt.tight_layout()
	plt.ylabel('True label')
	plt.xlabel('Predicted label')
	plt.savefig('SM_windowsize_' + str(win) + 'fold_validation_' + str(cvnum) + '_cvalue_' + str(cvalue) + '.png')
	#plt.show()


windowautomation(windowsize = 7, crossvalnum = 5, cvalue = 1)
windowautomation(windowsize = 9, crossvalnum = 5, cvalue = 1)
windowautomation(windowsize = 11, crossvalnum = 5, cvalue = 1)
windowautomation(windowsize = 13, crossvalnum = 5, cvalue = 1)
windowautomation(windowsize = 15, crossvalnum = 5, cvalue = 1)
windowautomation(windowsize = 17, crossvalnum = 5, cvalue = 1)
