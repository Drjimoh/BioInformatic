from Bio import ExPASy
from Bio import SwissProt


class uniPROT:
	"""
	author: @_drjimoh
	descr: This is a module for reverse trascription and reverse translation of
	protein sequence pulled from uniPROT
	The mission of UniProt is to provide the scientific community with a comprehensive, 
	high-quality and freely accessible resource of protein sequence and functional information. 
	"""
	rna_codon = [{ 'codon':'UUU', 'protein':'F'},      {'codon':'CUU', 'protein':'L'}, 
				    {'codon':'AUU', 'protein':'I'},
				     {'codon':'GUU', 'protein':'V'},
				   {'codon':'UUC', 'protein':'F'}, 
				   {'codon':'CUC', 'protein':'L'},  {'codon':'AUC', 'protein':'I'},
				   {'codon':'GUC', 'protein':'V'},
							{'codon':'UUA', 'protein':'L'},      {'codon':'CUA', 'protein':'L'}, 
							 {'codon':'AUA', 'protein':'I'},  {'codon':'GUA', 'protein':'V'},
							{'codon':'UUG', 'protein':'L'},  {'codon':'CUG', 'protein':'L'}, 
							{'codon':'AUG', 'protein':'M'},   
							{'codon':'GUG', 'protein':'V'}, 
							{'codon':'UCU' , 'protein':'S'}, {'codon':'CCU', 'protein':'P'},
							 {'codon':'ACU' , 'protein':'T'}, {'codon':'GCU', 'protein':'A'},
							{'codon':'UCC' , 'protein':'S'},     
							{'codon':'CCC', 'protein':'P'},      
							{'codon':'ACC','protein':'T'},
							 {'codon':'GCC', 'protein':'A'},
							{'codon':'UCA' , 'protein':'S'},     
							{'codon':'CCA', 'protein':'P'},      
							{'codon':'ACA','protein':'T'},
							 {'codon':'GCA', 'protein':'A'},
							{'codon':'UCG' , 'protein':'S'},     
							{'codon':'CCG', 'protein':'P'},      
							{'codon':'ACG','protein':'T'},
							 {'codon':'GCG', 'protein':'A'},
							{'codon':'UAU', 'protein':'Y'}, {'codon':'CAU', 'protein':'H'},     
							{'codon':'AAU', 'protein':'N'},      
							{'codon':'GAU', 'protein':'D'},
							{'codon':'UAC', 'protein':'Y'}, 
							{'codon':'CAC', 'protein':'H'},     
							{'codon':'AAC', 'protein':'N'},      
							{'codon':'GAC', 'protein':'D'},
							{'codon':'UAA' , 'protein':'Stop'},   {'codon':'CAA', 'protein': 'Q'},      
							{'codon':'AAA', 'protein': 'K'},      {'codon':'GAA', 'protein': 'E'},
							{'codon':'UAG' , 'protein':'Stop'},   {'codon':'CAG', 'protein': 'Q'},      
							{'codon':'AAG', 'protein': 'K'},     {'codon':'GAG', 'protein':'E'},
							{'codon':'UGU', 'protein':'C'},  {'codon':'CGU', 'protein':'R'}, 
							{'codon':'AGU' , 'protein':'S'},  {'codon':'GGU', 'protein':'G'},
							{'codon':'UGC', 'protein':'C'},  {'codon':'CGC', 'protein':'R'} , 
							{'codon':'AGC','protein':'S'},      {'codon':'GGC','protein':'G'},
							{'codon':'UGA' , 'protein':'Stop'},   
							{'codon':'CGA', 'protein': 'R'},    
							  {'codon':'AGA', 'protein': 'R'},      {'codon':'GGA', 'protein':'G'},
							{'codon':'UGG', 'protein':'W'},{'codon':'CGG' ,'protein':'R'},  
							    {'codon':'AGG' ,'protein':'R'},    
							   {'codon':'GGG', 'protein':'G'}]

	def __init__(self, prot_id):
		self.handle = ExPASy.get_sprot_raw(prot_id)
		self.record = SwissProt.read(self.handle)

	def r_translation(self):
		protein = self.record.sequence
		self.rna = ''
		for x in protein:
			for r in self.rna_codon:
				if x == r['protein']:
					self.rna += r['codon']
		return self.rna

	def r_trascription(self):
		self.r_translation()		
		self.dna = ''
		for x in self.rna:
			if x == 'A':
				self.dna += 'A'
			elif x == 'U':
				self.dna += 'T'
			elif x == 'G':
				self.dna += 'G'
			elif x == 'C':
				self.dna += 'C'
		return self.dna


if __name__ == '__main__':
	first_protein = uniPROT('Q5SLP9')
	print(first_protein.r_translation())
	print(first_protein.r_trascription())

