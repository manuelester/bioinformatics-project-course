# testingfile

# open test file and create empty lists
opentest = open('dssp_8_state.3line.txt', 'r+')
opentest = opentest.readlines()
protname = []
proteinseq = []
windows = []
secondstruc = []

# user chooses window size
win = int(input('Please enter an odd number integer: '))
zero = win/2 - 0.5
zero = int(zero)

# append protein names, sequence segments (according to window size), and secondary structure elements into empty lists from above
for i in range(0, len(opentest), 3):
    protname.append(opentest[i].lstrip('>').rstrip('\n'))
for i in range(1, len(opentest), 3):
	seq = ('0'*(zero)) + opentest[i].rstrip('\n') + ('0'*(zero))

	for j in range(zero,len(seq)-zero):
		lst = seq[(j-(zero)):(j+1+(zero))]
		lst = lst.upper()
		windows.append(lst)	

for i in range(2, len(opentest), 3):
    secondstruc.append(opentest[i].rstrip('\n'))

# double check to confirm that all windows are of the specified length

count = 0
for i in range(0,len(windows)):
	if len(windows[i]) < win:
		count += 1

# print(count)

# use map to convert features to numbers
map = {'H':1, 'G':2, 'I':3, 'E':4, 'B':5, 'T':6, 'S':7, 'C':8}
ss2 = ''.join(secondstruc)
ss2 = [char for char in ss2]
ss2 = [map[i] for i in ss2]

#make list of feature-window pairs
protSStruc = list(zip(windows, ss2))
#print(len(protSStruc))

# convert amino acids from string to integer number
map = {'A':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7,
    'I':8, 'K':9, 'L':10, 'M':11, 'N':12, 'P':13, 'Q':14,
    'R':15, 'S':16, 'T':17, 'V':18, 'W':19, 'Y':20, 'X':21, '0':0}

windowsnew = []
for element in windows:
	tmp = []
	for character in element:
		tmp2 = map[character]
		tmp.append(tmp2)
	windowsnew.append(tmp)

# print(windowsnew)

# convert amino acid integers to sparse encoding

from sklearn import preprocessing
enc = preprocessing.OneHotEncoder()
X_master = windowsnew
enc.fit(X_master)
X_master = enc.transform(X_master).toarray()

# X input is X_master and Y input is ss2
