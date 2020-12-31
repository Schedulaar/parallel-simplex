import json

parTimesFile = open("par_times.json")
parTimesJson = "[" + parTimesFile.read() + "]"
parTimesFile.close()

seqTimesFile = open("seq_times.json")
seqTimesJson = "[" + seqTimesFile.read() + "]"
seqTimesFile.close()




