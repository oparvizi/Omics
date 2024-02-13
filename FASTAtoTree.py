#//$menu window=main;popup=Scripts
# export Flu seqs from BN for V1 validation. 
# 1 file containing all 
# 1 file for each subtype containing selected sequences + refseqs, aligned & trimmed
# All outputs into user selected I:\\Sequencing-Respiratory Reference\\Flu\\NGS run transfers\\<subfolder>
# Created:	08/09/2020 Steve Platt
# REQUIRES:	BioEdit, MEGAX (standard installations for both available from IT software centre)

import bns
from Bio import SeqIO
import m_SeqHelpers as sh
import m_DialogBoxes as dbxs
import numpy as np
import time, subprocess

messagebox = bns.Util.Program.MessageBox

# QUIT if no entries are selected
if bns.Database.Db.Selection.__len__() == 0:
	messagebox("Information", "Please select some entries and rerun this script",'')
	__bnscontext__.Stop()

exper='HA_WGS'
header_fields=['Key','subtype (script)','Genetic group (script)']
temp_dir=bns.Database.Db.Info.TempDir
ref_dict={}
ref_dict['H1']={'present':False,'tmpfile':'H1_run_yymmdd.fas'}# outdir 'I:/Sequencing-Respiratory Reference/Flu/H1N1 (2009) WGS'}
ref_dict['H3']={'present':False,'tmpfile':'H3_run_yymmdd.fas'}# outdir 'I:/Sequencing-Respiratory Reference/Flu/H3N2 WGS'}
ref_dict['B']={'present':False,'tmpfile':'Bb_run_yymmdd.fas'}# outdir 'I:/Sequencing-Respiratory Reference/Flu/B_Influenza_WGS'}

# DIALOG for run yymmdd selection/entry [seqexp (always HA_WGS), header label fields (always header_fields above)]
dlg = dbxs.TextInputDlg("simple","Enter run date (yymmdd)")
if not dlg.Show():
	__bnscontext__.Stop()
run_date = dlg.inputSimple.GetValue()
# SET tmpfile in ref_dict
for k in ref_dict.keys():
	ref_dict[k]['tmpfile'] = ref_dict[k]['tmpfile'].replace('yymmdd',run_date)

# DIALOG for dest folder selection in I:\Sequencing-Respiratory Reference\Flu\NGS run transfers\
dlg = dbxs.SelectFileDir(style='directory', header='Select Run Folder')
if not dlg.Show():
	__bnscontext__.Stop()
dest_dir = dlg.fname+"\\"

refdata_path = "I:\\Sequencing-Respiratory Reference\\Flu\\NGS run transfers\\"
# DISPLAY settings
retval = messagebox("Information","SETTINGS FOR THIS ANALYSIS\n\nSequence experiment: "+exper+
"\n\nDatabase fields to label sequences: \n-"+'\n-'.join(header_fields)+
"\n\nRun date: "+run_date+
"\n\nRun folder location: \n"+dest_dir+
"\n\nReference data specified in: \n"+refdata_path+"\n-Ref_alignment_directories.txt\n-FASTA_seq_export_Flu_V1_validation_infer_NJ_nucleotide.mao"
"\n\n\nIf any of this information is incorrect please 'Cancel' and correct accordingly",
"question+okcancel")
if retval != 1:
	__bnscontext__.Stop()


# READ FILE of refseq & aligned output locations
with open(refdata_path+"Ref_alignment_directories.txt") as reffile:
			data=reffile.read().splitlines()
for line in data:
	if not line.startswith('#'):
		ref_dict[line.split('\t')[0]]['refset']=line.split('\t')[1]

#__bnscontext__.Stop()

# GATHER select seqs
all_dict={}
for entry in bns.Database.Db.Selection:
	seqstr = sh.GetSequence(entry.Key,exper)
	seq_header_list=[]
	subtype=''
	for header in header_fields:
		if header == 'Key':
			seq_header_list.append(entry.Key)
		else:
			seq_header_list.append(entry.Field(header).Content)
		if header == 'subtype (script)':
			subtype=entry.Field(header).Content
			if subtype=='H1':
				ref_dict[subtype]['present']=True
			if subtype=='H3':
				ref_dict[subtype]['present']=True
			if subtype.startswith('B'):
				ref_dict['B']['present']=True
	seqheader = '_'.join(seq_header_list)
	all_dict[seqheader]=[subtype,seqstr]

# WRITE all to fasta in dest folder
all_out = open(dest_dir+"Run_"+run_date+".fas",'w')
for k in all_dict.keys():
	all_out.write('>'+k+'\n'+all_dict[k][1]+'\n')
all_out.close()

# for EACH SUBTYPE specific (files in BN temp dir)
for stype in ref_dict.keys():
	if ref_dict[stype]['present']==True:
		stype_outfile_path = temp_dir+'\\'+ref_dict[stype]['tmpfile']
		stype_outfile = open(stype_outfile_path,'w')
		
		# EXTRACT, MERGE with ref seqs & write file
		__bnscontext__.SetBusy("Subtype: "+stype+"..merging refseqs")
		for k in all_dict.keys():
			if all_dict[k][0].startswith(stype):
				stype_outfile.write('>'+k+'\n'+all_dict[k][1]+'\n')
		ref_dict[stype]['num_seqs'] = sum(1 for line in open(ref_dict[stype]['refset']) if line.startswith('>')) # count number of refseqs
		with open(ref_dict[stype]['refset']) as refset:
			stype_outfile.write(refset.read())
		stype_outfile.close()

		# ALIGN with CLUSTALW in Bioedit apps folde)
		__bnscontext__.SetBusy("Subtype: "+stype+"..merged refseqs..aligning")
		program='C:/BioEdit/apps/clustalw.exe'
		ret = subprocess.call(program+' "'+stype_outfile_path+'"') #, shell=True)
		if ret == 1:	# clustalw.exe path error
			messagebox("Error","Cannot find "+program,'exclamation')
			__bnscontext__.Stop()
		if ret == 2:	# infile path error
			messagebox("Error","Infile path error for subtype: "+stype+"\n\n"+stype_outfile_path+"\n\nNot Found!",'exclamation')
			__bnscontext__.Stop()
		if ret != 0: # generic error
			messagebox("Error","An unspecified error prevented sequence alignment from running for subtype "+stype+"\n\nReturn code: "+str(ret),'exclamation')
			__bnscontext__.Stop()
		
		# CREATE a numpy matrix from the alignment file
		__bnscontext__.SetBusy("Subtype: "+stype+"..merged refseqs..aligned..trimming")
		alignfile = stype_outfile_path.replace('.fas','.aln')
		seqidlist = []
		seqlist = []
		for seq_record in SeqIO.parse(open(alignfile,"rU"),"clustal"):
			seqidlist.append(seq_record.id)
			seqlist.append(seq_record.seq)
		seqMatrix = np.array(seqlist)
		#(rows,cols) = seqMatrix.shape
		
		# DETERMINE how many starting columns where all refseqs are -
		if stype != 'B':
			maxtrim=0
			for x in range(100):
				if list(seqMatrix[:,x]).count('-') >= ref_dict[stype]['num_seqs']:
					maxtrim = x
			# REMOVE those columns from numpy array
			seqMatrix = np.delete(seqMatrix,np.s_[0:maxtrim+1],axis=1)
		# WRITE output to FASTA file
		# trimmed_file = ref_dict[stype]['outdir']+'/'+stype+'_Run_'+run_date+'.fas'					# LIVE ... use for final script version
		trimmed_file = dest_dir+'/'+stype+'_Run_'+run_date+'.fas'					# LIVE ... use for final script version
		with open(trimmed_file,'w') as outfile:
			(rows2,cols2)= seqMatrix.shape
			for r in range(rows2):
				outfile.write('>'+seqidlist[r]+'\n')
				outfile.write(('').join(seqMatrix[r])+'\n')
		
		# AUTOMATE TREE CREATION WITH MEGA - requires MEGAX <<<<<<<<<< ####### <<<<<<<< ####### <<<<<<<<<< ####### <<<<<<<<<< ####### <<<<<<<<<<<<
		__bnscontext__.SetBusy("Subtype: "+stype+"..merged refseqs..aligned..trimmed..tree creation")
		ret2 = subprocess.call('megacc -a "'+refdata_path+'FASTA_seq_export_Flu_V1_validation_infer_NJ_nucleotide.mao" -d "'+trimmed_file+'" -o "'+temp_dir+'\\'+stype+'"')
		if ret2 != 0:
			messagebox("Error","Tree generation with MEGAX failed for subtype "+stype+"\n\nReturn code: "+str(ret2),'exclamation')
			__bnscontext__.Stop()

		__bnscontext__.SetBusy("")
		
		# DISPLAY nwk tree in MEGA for manual manipulation/saving ... 
		# call DOESN'T WORK AS HOPED, script waits for MEGAX to finish before progressing	#######################################
		# USERS: reroot & save PDF trees as per SOP, close MEGAX & wait for next tree to open	#######################################
		program2 = "C:/Program Files/MEGA-X/MEGAX64.exe"
		ret3 = subprocess.call(program2+' "'+temp_dir+'\\'+stype+'.nwk"')
		if ret3 != 0:
			messagebox("Error","Displaying tree with MEGAX failed for subtype "+stype+"\n\nReturn code: "+str(ret3.returncode),'exclamation')
			__bnscontext__.Stop()

messagebox("Closing","Alignment, trimming and tree generation complete",'')












