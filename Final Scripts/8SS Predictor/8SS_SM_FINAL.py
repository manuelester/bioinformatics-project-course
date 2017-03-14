def randomforest(windowsize, filename):
	print('running random forest predictor for window size of: ' + str(windowsize) )
	import matplotlib
	matplotlib.use('Agg')
	# testingfile feb 26 new PSSMs from offline PSI-BLAST run
	# aa order: A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V

	import numpy as np

	# define sigmoid function to normalize substitution matrix scores
	import math
	def sigmoid(x):
		return 1 / (1 + math.exp(-x))





	############################ open test data and and create empty lists #########################################
	filename = str(filename)	
	test = open(filename, 'r+')
	test = test.readlines()
	test_protname = []
	test_proteinseq = []
	test_windows = []
	test_secondstruc = []
	test_pssmlist = []
	test_lst = []



	############################################# calculate padding width according to window size #########################################
	win = windowsize
	zero = win/2 - 0.5
	zero = int(zero)
	pad = list(np.zeros(20))




	######################################################################### make lists for test dataset ########################################################################
	for i in range(0, len(test), 3):
		protid = test[i].rstrip('\n')
		if protid in test_protname:
			protid = test[i].rstrip('\n') + '_2'
			if protid in test_protname:
				protid = test[i].rstrip('\n') + '_3'
				if protid in test_protname:
					protid = test[i].rstrip('\n') + '_4'
		test_protname.append(protid)
		# print(protid)
		protseq = test[i+1].rstrip('\n')
		test_proteinseq = protseq





	#################################### make windows list for test data ##########################################

	windows_test = []
	tmplist = []
	pssmlist = []

	for j in range (0,len(test_protname)):
		pssmlist = []
		pssmfile = open(test_protname[j] + '.fasta.pssm', 'r+')
		pssmfile = pssmfile.readlines()
		
		# convert each line of PSSMfile into a list by splitting on white space, only take range [2:] because column 1 is the position and column 2 the original amino acid which are not needed
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
	
	windows_test = tmplist






	########################################## Run RandomForest Classifier ############################

	from sklearn.ensemble import RandomForestClassifier

	from sklearn.externals import joblib
	clf = joblib.load('randomforest_SM_model.pkl')
	prediction = clf.predict(windows_test)
	prediction_list = list(prediction)

	map = {1:'H', 2:'G', 3:'I', 4:'E', 5:'B', 6:'T', 7:'S', 8:'C'}
	prediction_string = [char for char in prediction_list]
	prediction_string = [map[i] for i in prediction_list]
	prediction_string = ''.join(prediction_string)


	# Write prediction and sequence to a file
	predictionfile = open(str(test_protname[0]) + '_prediction.txt', 'w+')
	predictionfile.write('protein ID:' + '\n')
	predictionfile.write(str(test_protname[0]) + '\n')

	predictionfile.write('protein sequence:' + '\n')
	predictionfile.write(str(test_proteinseq) + '\n')

	predictionfile.write('predicted secondary structure:' + '\n')
	predictionfile.write(str(prediction_string) + '\n')





import sys
file_name = sys.argv[1]
#file_name = '>d1m4tA.fasta'
randomforest(windowsize = 11, filename = file_name)

