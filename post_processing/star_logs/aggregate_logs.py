logs=open('logs.txt','r').read().strip().split('\n')
outf=open('aggregate_logs.txt','w')
for line in logs:
    tokens=line.split('\t')
    sample=tokens[0]
    logfile=open(tokens[1],'r').read().strip().split('\n')
    outf.write(sample)
    for line in logfile:
        outf.write('|'+line)
    outf.write('\n')
    
