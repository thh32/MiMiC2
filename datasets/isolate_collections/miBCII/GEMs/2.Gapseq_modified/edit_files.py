import sys
import shutil
import os

print ("File of interest; ", sys.argv[1])

deleting = 0
turn_of_deletion = 0
group = 0

outputting = open('tmp.txt','w')

for line in open(sys.argv[1],'r'):
	if 'subsys_PWY_6151' in line:
		if '<groups:group sboTerm="SBO:0000411" groups:id="subsys_PWY_6151" groups:id="subsys_PWY_6151" groups:name="PWY-6151" groups:kind="partonomy">' in line:
			group+=1
		if group >1:
			deleting = 1
	

	if deleting == 0:       
		if 'groups:id="subsys_PWY_6151" groups:id="subsys_PWY_6151"' in line:
			outputting.write(line.replace('groups:id="subsys_PWY_6151" groups:id="subsys_PWY_6151"', 'groups:id="subsys_PWY_6151"'))
		else:
			outputting.write(line)
	

	if line == '      </groups:group>\n':
		deleting = 0



outputting.close()


shutil.move(os.path.join('./' , 'tmp.txt'), os.path.join('./', sys.argv[1]))
