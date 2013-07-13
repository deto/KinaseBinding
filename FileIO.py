import DataManipulation;

def read_fasta(filename):
    fin = open(filename);
    fcontents = fin.read();
    fin.close();
    
    records = fcontents.split('>')[1:];
    records = [r.splitlines() for r in records];
    output_pairs = [[r[0], ''.join(r[1:])] for r in records];
	
    return output_pairs;
    
def read_kinase(filename):
	fin = open(filename);
	fin_contents = fin.read().splitlines();
	fin.close();
	
	#stip out the first few lines
	fin_contents = fin_contents[4:];
	fin_contents_split = [x.split('\t') for x in fin_contents];

	fin_extract = [[x[0], x[4], x[13]] for x in fin_contents_split];
	return fin_extract;

def load_specific_kinase(kinase, start_center):
    fin_extract = read_kinase('Kinase_Substrate_Dataset');
    
    #select specidfied kinase and organism
    fg_read = [x[2] for x in fin_extract if x[0] == kinase and x[1] == 'human'];
    data_A = DataManipulation.list_to_formA(fg_read); 
    data_A = DataManipulation.subset_formA(data_A,0,start_center);
    
    return dict([('seq',fg_read),
                 ('formA',data_A)]);
    

def read_lines(filename):
    fin = open(filename);
    fcontents = fin.read().splitlines();
    fin.close();
    
    return fcontents;
    
