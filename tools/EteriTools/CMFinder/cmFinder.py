import glob, os, sys
import shutil
from shutil import copyfile
import subprocess
from os import system

def sh(script):
    system("bash -c '%s'" % script)


#print sh("pwd")


#input_dir = sys.argv[1]
model_tree_stk = sys.argv[1]
cmfinder_fa = sys.argv[2]
gapCmd = ""
gapVal = ""
if len(sys.argv) > 3:
    gapCmd = sys.argv[3]
    gapVal = sys.argv[4]


#sh("mkdir  CLUSTER")
#sh("cd CLUSTER/")
cmd = " cp -f " + model_tree_stk  +  " model.cmfinder.stk"
     #print cmd
sh(cmd)

alifoldCmd = "perl alifold.pl -file " + model_tree_stk
sh(alifoldCmd)

cmd_stk = "perl mloc2stockholm.pl -file model.cmfinder.stk  -split_input yes  --con_struct " + model_tree_stk + ".alifold"

sh(cmd_stk)

# sh("rm -f model.cmfinder.stk");

model_tree_stk = "model.cmfinder.stk.sth"


os.popen("cmfinder " +  gapCmd + " " + gapVal + " -a " + model_tree_stk +  " " + cmfinder_fa + " " +  " output > model.cmfinder.stk")# + directory +"/model.cmfinder.stk && rm " + directory + "/output")
if os.path.isfile('output') :
    ##copyfile("output", "model.cmfinder.stk")
    sh("rm output")
    sh("rm model.cmfinder.stk.sth")
else :
    copyfile("model.cmfinder.stk.sth", "model.cmfinder.stk")
    sh("rm model.cmfinder.stk.sth")
