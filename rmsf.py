from pylbtc.vmd import *
import numpy as np
import path_sys

def select_rmsf(NAME, INPUT, SELECTIONS, OUTPUT, EXTRA='', AV=True):

# Some verifications
    if len(SELECTIONS) != len(OUTPUT):
            raise IndexError("SELECTIONS and OUTPUTS must have the same length.")

    path = path_sys.path(NAME)
    psf = NAME+".0.psf"
    pdb = NAME+".0.pdb"

    mol = System(NAME)
    mol.load(path+psf)
    mol.load(path+pdb)

    for sel_str, out in zip(SELECTIONS, OUTPUT):
        if 'site' in sel_str:
            AV = False
        else:
            fd1 = open(out.format('-1'), 'w')
            fd2 = open(out.format('-2'), 'w')
        fd = open(out.format(''), 'w')

        sel = mol.selectAtoms(sel_str+EXTRA)
        ref = []
        for atoms in sel:
            tmp = mol.selectAtoms('index {}'.format(atoms))
            ref.append((tmp['segname'][0], tmp['resid'][0]))

        ref = sorted(ref, key=lambda x: (x[0], int(x[1])))
        rmsf = open(INPUT).readlines()
        seg1 = []
        seg2 = []
        indexAxisX1 = 0
        indexAxisX2 = 0
        indexAxisXav = 0
        for ref_line in ref:
            for line in rmsf:
                line_piece = line.split()
                if line_piece[0] == ref_line[0] and int(line_piece[1]) == ref_line[1]:
                    fd.write('   '.join(line_piece)+'\n')
                    if 'site' in sel_str:
                        continue
                    elif '1' in line_piece[0]:
                        indexAxisX1 += 1
                        line_piece.insert(0, str(indexAxisX1))
                        fd1.write('   '.join(line_piece)+'\n')
                        seg1.append(line_piece)
                    elif '2' in line_piece[0]:
                        indexAxisX2 += 1
                        line_piece.insert(0, str(indexAxisX2))
                        fd2.write('   '.join(line_piece)+'\n')
                        seg2.append(line_piece)
        if AV:
            fd_av = open(out.format('-av'), 'w')
            for line1,line2 in zip(seg1, seg2):
                tmp_l = []
                rmsf_line1 = float(line1[3])
                rmsf_line2 = float(line2[3])
                tmp_l.append(rmsf_line1)
                tmp_l.append(rmsf_line2)
                line1[3] = '{:.4f}'.format(np.mean(tmp_l))
                line1[4] = '{:.4f}'.format(np.std(tmp_l))
                line1[1] = line1[1][:3]
                #indexAxisXav += 1
                #line_piece.insert(0, str(indexAxisXav))
                fd_av.write('   '.join(line1)+'\n')

