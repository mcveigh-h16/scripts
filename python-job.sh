virtualenv myenv 
source myenv/bin/activate
pip install biopython
pip install pandas

python FindIdenticalSeqMixedOrg.py missing_bacteria_nocyanos.gbk missing_bacteria_nocyanos.unique
