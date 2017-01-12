import os
import sys
from shutil import copyfile
from os import system


def sh(script):
    system("bash -c '%s'" % script)


model_tree_stk = sys.argv[1]
cmfinder_fa = sys.argv[2]
path = sys.argv[3]

gapCmd = ""
gapVal = ""
if len(sys.argv) > 4:
    gapCmd = sys.argv[4]
    gapVal = sys.argv[5]


cmd = " cp -f %s model.cmfinder.stk" % (model_tree_stk)
sh(cmd)

alifoldCmd = "%salifold.pl -file  %s" % (path, model_tree_stk)
# alifoldCmd = "perl " + path + "/alifold.pl -file " + model_tree_stk
sh(alifoldCmd)

cmd_stk = "%smloc2stockholm.pl -file model.cmfinder.stk  -split_input yes --con_struct %s.alifold" % (path, model_tree_stk)
# cmd_stk = "perl " + path + "/mloc2stockholm.pl -file model.cmfinder.stk  -split_input yes --con_struct " + model_tree_stk + ".alifold"
sh(cmd_stk)

model_tree_stk_sth = "model.cmfinder.stk.sth"
x = "cat " + model_tree_stk_sth
sh("mv model.cmfinder.stk.sth model.tree.stk")

sh("cmfinder %s %s -a model.tree.stk %s output > model.cmfinder.stk" % (gapCmd, gapVal, cmfinder_fa))
# sh("cmfinder " + gapCmd + " " + gapVal + " -a model.tree.stk" + " " + cmfinder_fa + " " + " output > model.cmfinder.stk")

if os.path.isfile('output'):
    sh("rm output")
else:
    copyfile("model.tree.stk", "model.cmfinder.stk")
