import json

parTimesFile = open("par_times.json", encoding='utf-8')
parJson = json.loads(parTimesFile.read())
parTimesFile.close()

seqTimesFile = open("seq_times.json", encoding='utf-8')
seqJson = json.loads(seqTimesFile.read())
seqTimesFile.close()

output = ""
speedup_output = ""

coordinates = ""
speedup_coordinates = ""
last_n = 0
for i in range(len(seqJson)):
    if seqJson[i]["n"] > last_n:
        coordinates += "(%i, %f)\n" % (seqJson[i]["n"], seqJson[i]["t"])
        speedup_coordinates += "(%i, %f)\n" % (seqJson[i]["n"], 1)
        last_n = seqJson[i]["n"]


output += r'''
\addplot coordinates { %% seq
    %s
};
''' % coordinates
speedup_output += r'''
\addplot coordinates { %% seq
    %s
};
''' % speedup_coordinates


def find(pred, iterable):
    for element in iterable:
        if pred(element):
            return element
    return None


for N in range(1, 5):
    last_n = 0
    coordinates = ""
    speedup_coordinates = ""
    for i in range(len(parJson)):
        if parJson[i]["N"] == N and parJson[i]["n"] > last_n:
            coordinates += "(%i, %f)\n" % (parJson[i]["n"], parJson[i]["t"])
            seq_t = find(lambda l: l["n"] == parJson[i]["n"], seqJson)["t"]
            speedup_coordinates += "(%i, %f)\n" % (parJson[i]["n"], seq_t / parJson[i]["t"])
            last_n = parJson[i]["n"]
    output += r'''
        \addplot coordinates { %% %ix%i
            %s
        };
    ''' % (N, N, coordinates)
    speedup_output += r'''
        \addplot coordinates { %% %ix%i
            %s
        };
    ''' % (N, N, speedup_coordinates)

outputFile = open("times.tex", "w")
outputFile.write(output)
outputFile.close()

outputFile = open("times_speedup.tex", "w")
outputFile.write(speedup_output)
outputFile.close()
