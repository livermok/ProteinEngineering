'''Calclates a score for predicted orthagonality of PhoB in our Bacterial System'''
#open files and save as a list
f = open('alignments') '''File that contains clustal alignments converted to one line per sequence'''
lines = f.readlines()
PhoBFile = open('PhoB')
PhoB = PhoBFile.readlines() '''PhoB extracted from alignment'''
#initilize variables
length = len(PhoB[0])
roi = list(range(4, length))
samples = range(len(lines))
score=list(range(len(lines)))
#Set residues to score
surface = [4,5,7]
core = [6,8]
#loop to produce scores
for j in samples:
	score[j] = 0
	for i in surface:
		if unicode(PhoB[0][i]) != unicode(lines[j][i]): ''' unicode to recognize period as a period'''
			score[j] = score[j] + 1                 ''' Rewards divergenc from PhoB specificity Residues'''
	for i in core:
		if unicode(PhoB[0][i]) == unicode(lines[j][i]): ''' Rewards core structural similarities to PhoB'''
			score[j] = score[j] + 1
#output scores
print score
