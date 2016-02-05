import os
import random

random.seed()

if os.path.exists("HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_ALIGN.txt"):
    command = "rm HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_ALIGN.txt"
    os.system(command)

if os.path.exists("HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_VALIDATE.txt"):
    command = "rm HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_VALIDATE.txt"
    os.system(command)

file = open("temp.txt")

command = "echo \"CastorPool=cmscaf\" >> HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_ALIGN.txt"
os.system(command)
command = "echo \"CastorPool=cmscaf\" >> HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_VALIDATE.txt"
os.system(command)

for line in file.readlines():
    line=line.replace("\n","")

    i = random.randrange(1,100,1)

    if i <= 7:
        command = "echo \"" + str(line) + "\" >> HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_ALIGN.txt"
    else:
        command = "echo \"" + str(line) + "\" >> HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_VALIDATE.txt"

    os.system(command)
