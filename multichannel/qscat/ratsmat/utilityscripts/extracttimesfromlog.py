SRCHSTR = "nroots"
OUTFILE = "nroots"

with open("cal.log", "r") as fin:
  with open(OUTFILE+".txt", "w+") as fout:
    for line in fin:
      ind = line.find(SRCHSTR + " Time Taken: ")
      if ind != -1:
        fout.write(line[ind+len(SRCHSTR + " Time Taken: "):])