print('WARNING: In order to run this code, PSSM files for all the corresponding proteins in the dataset must first be available in a folder named PSSM')

def randomforest(windowsize):
	print('building random forest predictor for window size of: ' + str(windowsize) )
	import matplotlib
	matplotlib.use('Agg')
	# testingfile feb 26 new PSSMs from offline PSI-BLAST run
	# aa order: A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V

	import numpy as np

	# define sigmoid function to normalize substitution matrix scores
	import math
	def sigmoid(x):
		return 1 / (1 + math.exp(-x))

	############################ open train data and and create empty lists #########################################
	opentest = open('dssp_8_state.3line.txt', 'r+')
	opentest = opentest.readlines()
	protname = []
	proteinseq = []
	windows = []
	secondstruc = []
	pssmlist = []
	lst = []




	############################################# calculate padding width according to window size #########################################
	win = windowsize
	zero = win/2 - 0.5
	zero = int(zero)




	######################################################################### make lists for train dataset ########################################################################
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





	################################# make windows list for train set ############################################
	pad = list(np.zeros(20))
	
	windows_train = []
	tmplist = []
	pssmlist = []

	for j in range (0,len(protname)):
		pssmlist = []
		pssmfile = open('./PSSM/' + protname[j] + '.fasta.pssm', 'r+')
		pssmfile = pssmfile.readlines()

		for k in range(3, len(pssmfile)-6):
			seq = pssmfile[k].rstrip('\n')
			seq = seq.split()
			seq = seq[2:22]
			# seq = seq[22:42]

			# make sure that all the numbers are converted to float by iterating through the list
			for l in range (0,len(seq)):
				seq[l] = float(sigmoid(float(seq[l])))
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
#	windows_train.append(tmplist)
	windows_train = tmplist

	# use map to convert features to numbers

	map = {'H':1, 'G':2, 'I':3, 'E':4, 'B':5, 'T':6, 'S':7, 'C':8}
	ss2 = ''.join(secondstruc)
	ss2 = [char for char in ss2]
	ss2 = [map[i] for i in ss2]
	ss_new = ss2







	########################################## Run RandomForest Classifier ############################

	from sklearn.ensemble import RandomForestClassifier

	predlist = []

	clf = RandomForestClassifier()
	clf.fit(windows_train, ss_new)

	from sklearn.externals import joblib
	joblib.dump(clf, 'randomforest_SM_model.pkl')




randomforest(windowsize = 11)

