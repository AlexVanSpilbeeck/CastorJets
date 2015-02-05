import os
import re
import sys
import time

# define our method
def replace_all(text, dic):
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text

nruns = 4
nevents = 250000
date = '20150205_no_gen_cut'

os.system("make clean")
#os.system('cd /scratch/; rm -rf CMSSW_4_2_10_patch2')


gen_rad = [5, 5]
det_rad = [5, 7]

for gen in gen_rad:
    for det in det_rad:

        newRadii = "ak" + str(gen) + "_ak" + str(det) + "_" + date + "_" + str(nevents)
#	reps = {'radii': newRadii, 'MainAnalyzer':'MainAnalyzer_2'}

	# .cc
#        ccfile = open("src/JetAnalyzer_radii.cc")
#        outcc = open("src/JetAnalyzer_" + newRadii + ".cc", "w")
#        for line in ccfile:
#	    outcc.write( line.replace("JetAnalyzer_radii","JetAnalyzer_ak" + str(gen) + "_ak" + str(det)))
#           outcc.write( replace_all(line, reps) )
#        outcc.close()

	# .h
#        hfile = open("src/JetAnalyzer_radii.h")
#        outh = open("src/JetAnalyzer_" + newRadii + ".h", "w")
#        for line in hfile:
#
#            outh.write( line.replace("JetAnalyzer_radii","JetAnalyzer_ak" + str(gen) + "_ak" + str(det)))
#           outh.write( replace_all(line, reps) )
#        outh.close()

	# MainAnalyzer.cc
#	MAfile = open("src/MainAnalyzer.cc")
#        outMA = open("src/MainAnalyzer_2.cc","w")
#	for line in MAfile:
#	    if 'radii' in line:
#   	        outMA.write( line.replace("radii","ak" + str(gen) + "_ak" + str(det)))
#	    elif 'MainAnalyzer' in line:
# 	        outMA.write( line.replace("MainAnalyzer", "MainAnalyzer_2"))
#	    else:
#		outMA.write(line)	
#            if  'radii' in line:
#	        print "Found lin\t" + line
#	    outMA.write( replace_all(line, reps) )
#	outMA.close()
#os.system('mv src/MainAnalyzer_2.cc src/MainAnalyzer.cc')

        # MainAnalyzer.h
#        MAhfile = open("src/MainAnalyzer.h")
#       outMAh = open("src/MainAnalyzer_2.h","w")
#        for line in MAhfile:
#            if 'radii' in line:
#                outMAh.write( line.replace("radii","ak" + str(gen) + "_ak" + str(det)))
#            elif 'MainAnalyzer' in line:
#                outMAh.write( line.replace("MainAnalyzer", "MainAnalyzer_2"))
#	    else: 
#		outMAh.write(line)
#	    outMAh.write( replace_all(line, reps) )
#        outMAh.close()
#	os.system('src/MainAnalyzer_2.h src/MainAnalyzer.h')

	#Run.cc
#	Runfile = open("Run_backup.cc")
#	outRun = open("Run_2.cc","w")
#        for line in Runfile:
#            if 'radii' in line:
#	        outRun.write( line.replace("radii","ak" + str(gen) + "_ak" + str(det)))
#            elif "MainAnalyzer.h" in line:
#	        outRun.write('#include "./src/MainAnalyzer_2.h"\n')
#	    else:
#	        outRun.write(line)
#	    outRun.write( replace_all(line, reps) )
#        outRun.close()
#	os.system('mv Run_2.cc Run.cc')
#	os.system('cp Run.cc Run' + str(gen) + str(det) )

	print "before make"

#        print str(gen) + str(det)
        os.system('cmsenv')
        
	print "done cmsenv"

        os.system('make')

	print "done make"

        shfile = open("script_3.sh")
        out_sh = open("script_" + newRadii +".sh", "w")

        for line in shfile:        
            if 'date=' in line:
                out_sh.write("date='" + date + "'\n")
            elif 'gen_radius=' in line:
                out_sh.write("gen_radius='ak" + str(gen) + "'\n")
            elif 'det_radius=' in line:
                out_sh.write("det_radius='ak" + str(det) + "'\n")
            elif 'totalEvents=' in line:
                out_sh.write("totalEvents=" + str(nevents) + "\n")
            else:
	        out_sh.write(line)
#                out_sh.write( line.replace("radii", newRadii ))
        out_sh.close()
        print "script_ak" + str(gen) + "_ak" + str(det) +".sh" 
#        os.system('mv ' + 'script_ak' + str(gen) + '_ak' + str(gen) + '.sh  script.sh')

        jobname = str(gen) + str(det) + str (nevents)

        os.system('qsub -q localgrid@cream02 -o ' + jobname + '.stdout -e ' + jobname + '.stderr -N ' + jobname + ' script_' + newRadii + '.sh')
        print '\t\tqsub -q localgrid@cream02 -o ' + jobname + '.stdout -e ' + jobname + '.stderr -N ' + jobname + ' script_ak' + str(gen) + '_ak' + str(det) + '.sh'
	print './Run_2 Pythia6Z2star_new 7000 JetAnalyzer_radii ' + str(gen) + ' ' + str(det) + ' ' + str(nevents)

	time.sleep(10)
#os.system('rm -rf src/MainAnalyzer_2*')
