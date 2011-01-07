from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pdz_exp_data

def calculate(pdz_seq, ligand_seq):

  binary_data_score = 0
  affinity_data_score = 0

  for relative_pos_pdz, pos_ligand in pdz_exp_data.position_pairs:
    #note the correctio nfor string indexing
    score1, score2 = pdz_exp_data.exp_data[(
      relative_pos_pdz, 
      pos_ligand, 
      pdz_seq[pdz_exp_data.a1synPDZ_key_positions[relative_pos_pdz] - 1], 
      ligand_seq[pos_ligand - 1])]
    
    binary_data_score += score1
    affinity_data_score += score2
    

  return binary_data_score,affinity_data_score
  

# a1synPDZ
pdz_prot = Seq("RRRVTVRKADAGGLGISIKGGRENKMPILISKIFKGLAADQTEALFVGDAILSVNGEDLSSATHDEAVQALKKTGKEVVLEVKYMK", IUPAC.protein)

ligand_terminal_residues = [
  ("kek1",Seq("CGTDV", IUPAC.protein)),
  ("kek2",Seq("ETSDI", IUPAC.protein)),
  ("kek3",Seq("DIFKS", IUPAC.protein)),
  ("kek4",Seq("VDISI", IUPAC.protein)),
  ("kek5",Seq("DGTEV", IUPAC.protein)),
  ("kek6",Seq("EFVSL", IUPAC.protein)),
  ("gliotactin",Seq("PQTSV", IUPAC.protein)),
  ("egfr",Seq("TETRV", IUPAC.protein)),
  ("example",Seq("GPDRDRESIV", IUPAC.protein))
]

print 'protein 1\tprotein 2\t    score    \t    score'
print '         \t         \t(binary data)\t(affinity data)'
print '----------------------------------------------------\n'
for name, seq in ligand_terminal_residues:
  binary_data_score, affinity_data_score = calculate(pdz_prot, seq)
  print "a1synPDZ\t%s\t\t%f\t%f" % (name, binary_data_score, affinity_data_score)


