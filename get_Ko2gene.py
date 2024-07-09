import sys


def get_hmm_info(hmmout):
    gene2kolis = {}
    flag = 0
    with open(hmmout) as f_in:
        for line in f_in:
            if "Query:" in line:
                flag = 1
            elif "Domain annotation" in line:
                flag = 0
            if flag == 1:
                line = line.strip()
                if line == '':
                    pass
                else:
                    if line[0:3] != '---':
                        if line[0:3] != 'Sco':
                            if line[0:3] != 'E-v':
                                if line[0:3] != '[No':
                                    line = line.strip().split()
                                    if line[0] == 'Query:':
                                        ids = line[1]
                                    else:
                                        if float(line[0]) < 1e-5:
                                            gid = line[8]
                                            score = float(line[1])
                                            if gid in gene2kolis:
                                                gene2kolis[gid][ids] = score
                                            else:
                                                gene2kolis[gid] = {ids:score}

    return gene2kolis



def get_besthit(gene2ko):
    for gene,vlis in gene2ko.items():
        koid = max(vlis,key=vlis.get)
        print(f"{gene}\t{koid}")


if __name__ == '__main__':
    gene2kolis = get_hmm_info(sys.argv[1])
    get_besthit(gene2kolis)

