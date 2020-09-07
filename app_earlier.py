from flask import Flask, render_template,request,g,flash,redirect,url_for
from flask_wtf import FlaskForm
from wtforms import StringField, IntegerField,  PasswordField, TextAreaField, RadioField, SelectField
from wtforms.validators import InputRequired,Length,AnyOf,ValidationError    
from subprocess	import	call
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

import pandas as pd
import sqlite3
import string
import random
import primer3
import uuid

#import pyodbc



app = Flask(__name__)
app.config['SECRET_KEY'] = 'Feludaisjatayu!'

#@app.teardown_appcontext



class jatayu(FlaskForm):

    variation = SelectField('Type of Mutation', choices=[('A', 'A'), ('C', 'C'),('G', 'G'), ('T', 'T') ]   )
    organism = SelectField('Organism', choices=[('Homo Sapiens', 'Homo Sapiens') ]   )

def base_check(form, field):
    
    pattern = "[^ACTG]+"
    if re.findall(pattern, field.data):
               #flash('Provide Valid DNA sequence','error')

               raise ValidationError('Provide Valid DNA sequence')
     
class MyForm(jatayu):
    sequence = TextAreaField('Sequence', validators=[InputRequired('Please provide the genomic sequence'),base_check,Length(min=20,max=30,message='Must be between 20 and 30') ] ) 
    position = IntegerField('Position', validators=[InputRequired('Position of the mutation')])
    
@app.route('/', methods=['GET', 'POST'])
def form():
    form = MyForm()
    if request.method == 'POST' :
       s =list(form.sequence.data)
       print(form.sequence.data)
       mut_pam = ("".join((s[form.position.data-1:form.position.data+4 ] )) )

       if (form.variation.data == s[form.position.data-1]) :
         flash('Mutation is same as wild type','error')
         return redirect(url_for('form'))
       if not mut_pam.endswith ('GG'): 
           
         flash ('Enter the mutation followed by NNGG','error')
         return redirect(url_for('form'))
                
       else:
        mutation_type=  s[form.position.data-1] + ">" + form.variation.data   
        s[form.position.data-1]=form.variation.data
        mutation_sequence="".join(s)
        sq_len=len(form.sequence.data)
        if form.validate_on_submit():
         if  align(form.sequence.data,sq_len,mutation_sequence,form.position.data,form.variation.data) == "Not a Homo Sapiens Sequence":
          human= "Not A Homo Sapiens Sequence, Try A Valid Sequence"
          return render_template('results_failure.html',sequence=human)
         else:    
          mapped,primer3_desg =align(form.sequence.data,sq_len,mutation_sequence,form.position.data,form.variation.data) 
          print (mapped)
          print (primer3_desg)
          return render_template('results.html', sequence=form.sequence.data, position=form.position.data,variation_base=form.variation.data, variation=mutation_type,mutation=mutation_sequence,organism=form.organism.data,
          tables=[mapped.to_html(classes='sgrna',index=None),primer3_desg.to_html(classes='primer',index=None)   ],
          titles = [ 'na' ,'sgRNA', 'Primers']  )
        
    return render_template('form.html', form=form)

def sgrna_des(muta):
    sg=list(muta)
    #print (sg[2])
    #return sg
    puri=['A','G']
    pyri=['C','T']
    if sg[14]  == 'T' or sg[14] == 'C':
       #return sg[5]   
       sg[14] =random.choice(puri)
       sgrna="".join(sg)
       return sgrna
    else:
       sg[14] =random.choice(pyri)
       sgrna="".join(sg)
       return sgrna

tab = str.maketrans("ACTG", "TGAC")

def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]

def desg_primer_wd_pam(oligo):
    primer =primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'random',
        'SEQUENCE_TEMPLATE': oligo,
        'SEQUENCE_INCLUDED_REGION': [0,490]
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 20,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 50.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 75.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_LEFT_NUM_RETURNED':2,
        'PRIMER_RIGHT_NUM_RETURNED':2,
        'PRIMER_PRODUCT_SIZE_RANGE': [[489,500]],
    })
    data = {'Left Primer': [primer['PRIMER_LEFT_0_SEQUENCE'], primer['PRIMER_LEFT_1_SEQUENCE'] ],
       'Left Primer Tm': [primer['PRIMER_LEFT_0_TM'], primer['PRIMER_LEFT_1_TM'] ],     
            
            'Right Primer': [primer['PRIMER_RIGHT_0_SEQUENCE'], primer['PRIMER_RIGHT_1_SEQUENCE']],
     'Right Primer Tm': [primer['PRIMER_RIGHT_0_TM'], primer['PRIMER_RIGHT_1_TM']],       
           
           
           }
    primer_df=pd.DataFrame.from_dict(data)

    return primer_df

def sgrna_sequence(genome,out,out1,position,variation):

 file1=open(out,'r') 
 file2=open(out1,'w')

 records = SeqIO.to_dict(SeqIO.parse(open(genome), 'fasta'))
 for i  in file1:
         i=i.strip()
         line=re.split('\t',i)
         print (line[0])
         long_seq_record = records[line[0]]
         long_seq = long_seq_record.seq
         #alphabet = long_seq.alphabet
         #short_seq = str(long_seq)[int(line[1])-1:int(line[2])]

         if line[3] == '+':
             short_seq20 = str(long_seq)[int(line[1]) -1 + int(position)  -19:int(line[1])+int(position)]
             short_seq23 = str(long_seq)[int(line[1]) -1 + int(position)  -19:int(line[1])+int(position) + 3 ]

             line1="\t".join(line)
             end=int(line[1])+int(position)
             start=int(line[1]) -1 + int(position) -19
             ss=list(short_seq20)
             print (len(ss))

             ss[18]=variation
             sg_sequence="".join(ss)
             print (sg_sequence)
             
             sgrna_des1=sgrna_des(sg_sequence)
             if short_seq23.endswith('GG'):
                 print (short_seq23)  
                 oligo = str(long_seq) [int(line[1]) -252 + int(position)  :int(line[1])+int(position) +250   ]
                 file2.write(line1+"\t"+ str(start)  + "\t" +   str(end) + "\t" +  sg_sequence + "\t" +  sgrna_des1 +"\t" +  oligo   + "\n")	
             
             """
             else:
                 oligo1 = str(long_seq) [int(line[1]) -200 :int(line[1])+ int(position) +1 ]
                 oligo2 = str(long_seq) [int(line[1])+ int(position) +1 :int(line[1])+int(position) +17   ]
                 pam='GG'
                 oligo=oligo1+pam +oligo2  
                 file2.write(line1+"\t"+ str(start)  + "\t" +   str(end) + "\t" +  sg_sequence + "\t" +  sgrna_des1 +"\t" +  oligo   + "\n")    
             """          
         elif line[3] == '-' : 
             short_seq20_rev = str(long_seq)[int(line[2]) - int(position) -2:int(line[2])  - int(position) +19-1 ]
             short_seq20_rev1 =  reverse_complement_table(short_seq20_rev) 
             short_seq23 = str(long_seq)[int(line[1]) -3 + int(position)  -19:int(line[1])+int(position) + 3 ]

             line1="\t".join(line)
             end=int(line[2])  - int(position) +19-1
             start=int(line[2]) - int(position) -2
             ss=list(short_seq20_rev1)
             print (len(ss))
             ss[18]=variation
             sg_sequence="".join(ss)
             print (sg_sequence)
             sgrna_des1=sgrna_des(sg_sequence)
             if short_seq23.endswith('GG'):
                 print (short_seq23)  
                 oligo = str(long_seq) [int(line[1]) -252 + int(position)  :int(line[1])+int(position) +250   ]
                 file2.write(line1+"\t"+ str(start)  + "\t" +   str(end) + "\t" + "\t" +  sgrna_des1 + sg_sequence +"\n")	



         
        

	 
def align(wtype,len1,mutation_sequence,position,variation):
    my_seqs = SeqRecord(Seq(wtype,generic_dna), id = "randomsequence")
    filename = uuid.uuid4().hex + ".fasta"
    filename_sam = filename + ".sam"
    filename_bed = filename_sam + ".bed"
    filename_bed_gene = filename_bed + "_gene"
    filename_bed_gene1 = filename_bed_gene + "_1"
    filename_bed_gene1_sgrna = filename_bed_gene1 + "_sgrna"
    SeqIO.write(my_seqs, filename, "fasta")
    #cmdArgs=['bwa', 'aln' ,'-N','-l', '40', '-k', '0'  ,'-t', '4'  ,'GRCh38.p13.genome.fa' , 'example.fasta',  '-f', 'example.sai'  ] 
    cmdArgs=['bwa' ,'mem'  ,'-t' ,'8', '-k', str(len1) , '-T' , '20', '-a' ,'GRCh38.p13.genome.fa',filename, '-o' ,filename_sam] # 'bwa', 'mem' ,'a','-t', '4'  ,'GRCh38.p13.genome.fa' , 'example.fasta',  '-o', 'example.sam'  
    call(cmdArgs)
    #cmdArgs=['bwa' ,'mem'  ,'-t' ,'8', '-k', str(len1) , '-T' , '20', '-a' ,'GRCh38.p13.genome.fa','example.fasta', '-o' ,'example.sam'] # 'bwa', 'mem' ,'a','-t', '4'  ,'GRCh38.p13.genome.fa' , 'example.fasta',  '-o', 'example.sam'  
    #call(cmdArgs)
    #cmdArgs1=        [ 'bwa' ,'samse', '-f', 'example.sam'  ,'GRCh38.p13.genome.fa' ,'example.sai' , 'example.fasta'  ]
    #call(cmdArgs1) 
    #cmdArgs2=['grep', 'random', 'example.sam' ]
    #print (cmdArgs2)
    #call(cmdArgs2,stdout=open('actual_seq','w'))
      
    df=pd.read_csv(filename_sam,sep='\t',header=None,comment="@")
    df=df.iloc[:,0:12]
    df.columns=['QNAME','Strand','Chr','Start','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','NM']
    if df['Chr'][0] == '*':
        print ("NOT HUMAN SEQ")
        return "Not a Homo Sapiens Sequence"
    else:    
    #df['CIGAR']=df['CIGAR'].replace('M','')
     df['CIGAR'] = df['CIGAR'].str.replace(r'\D', '').astype(int)    
     df['End']=df['Start'] + df['CIGAR']   
     df=df.loc[df['NM']=='NM:i:0', ['Chr','Start', 'End','Strand'] ]
    #df['SEQ'] = wt_seq
     df['mutation_seq'] = mutation_sequence
     df.loc[df['Strand'] == 0, 'Strand'] = '+'
     df.loc[df['Strand'] == 256, 'Strand'] = '+'
     df.loc[df['Strand'] == 16, 'Strand'] = '-'
     df.loc[df['Strand'] == 272, 'Strand'] = '-'
     df.to_csv(filename_bed,sep='\t',index=None,header=None)
     cmdArgs3=['bedtools', 'intersect', '-a', filename_bed,'-b', 'gencode.v32.annotation_genes_filtered.gtf' , '-wb'   ]
     call(cmdArgs3,stdout=open(filename_bed_gene,'w'))
     cmdArgs4=[ 'cut',  '-f' , '1-6,9', filename_bed_gene]
     call(cmdArgs4,stdout=open(filename_bed_gene1,'w'))
     sgrna_sequence('GRCh38.p13.genome.fa',  filename_bed_gene1,filename_bed_gene1_sgrna ,position,variation)
     df3=pd.read_csv(filename_bed_gene1_sgrna,sep='\t',header=None)
     df3.columns=['Chr','Start','End', 'Strand'  ,'WT-Seq', 'Disease_Seq',  'Gene','Start_sgrna','End_sgrna','sgrna','sgrna_seq','oligo']
     print (df3['oligo'][0])
     primer3_desg = desg_primer_wd_pam(df3['oligo'][0])
    #primer3_desg = df3.apply(lambda row : desg_primer_wd_pam(row['oligo']), axis = 1)
    #df3['sgrna_seq']= df3.apply(lambda row : sgrna_des(row['sgrna']), axis = 1)
     df3=df3.drop(['Start','End', 'Strand'  ,'WT-Seq','Disease_Seq','sgrna','oligo'],axis=1   )   
     print (type(primer3_desg))
     return df3,primer3_desg
 

if __name__ == '__main__':
    app.run(debug=True)