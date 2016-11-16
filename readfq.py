def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

if __name__ == "__main__":
    import sys
    n, slen, qlen, GCs, Q20s, Q30s = 0, 0, 0, 0, 0, 0
    for name, seq, qual in readfq(sys.stdin):
        n += 1
        slen += len(seq)
        qlen += qual and len(qual) or 0
        for base in seq:
            if base in "GCgc":
                GCs += 1
        if qlen > 0:
            for quality in qual:
                if (ord(quality)-33) >= 20:
                    Q20s += 1
                if (ord(quality)-33) >= 30:
                    Q30s += 1
    gc = round(float(GCs)/float(slen) * 100, 2)
    if qlen > 0:
        q20 = round(float(Q20s) / float(qlen) * 100, 2)
        q30 = round(float(Q30s) / float(qlen) * 100, 2)
        print('\t'.join(('Reads Number'.center(10), 'Base Number'.center(10), 'Quality Number'.center(10), 'GC%'.center(10), 'Q20%'.center(10), 'Q30%'.center(10))))
        print('\t'.join((str(n).center(10), str(slen).center(10), str(qlen).center(10), str(gc).center(10), str(q20).center(10), str(q30).center(10))))
    else:
        print('\t'.join(('Reads Number'.center(10), 'Base Number'.center(10), 'GC%'.center(10))))
        print('\t'.join((str(n).center(10), str(slen).center(10), str(gc).center(10))))
